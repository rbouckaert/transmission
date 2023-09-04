package transmission.operator;


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import transmission.distribution.TransmissionTreeLikelihood;
import transmission.distribution.Validator;
import transmission.logger.ColouredTreeLogger;

@Description("Operator that randomly picks an infection and moves it elsewhere")
public class InfectionMover2 extends Operator {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
	final public Input<TransmissionTreeLikelihood> likelihoodInput = new Input<>("likelihood", "transmission treelikelihood containing the colouring", Validate.REQUIRED);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private TransmissionTreeLikelihood likelihood;
    private TreeInterface tree;
    private int [] colourAtBase;

    @Override
	public void initAndValidate() {
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	likelihood = likelihoodInput.get();
    	tree = likelihood.treeInput.get();
	}

    
    final static boolean debug = false;
    int [][] stats = new int[3][3];
    
	@Override
	public double proposal() {
		
		
		int pre = blockCount.getValue(0);
		double logHR = doproposal();
		logHR = 0;
		int post = blockCount.getValue(0);
		updateStats(pre, post);

		if (debug) {
			colourAtBase = likelihood.getFreshColouring();		
			Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
			if (!validator.isValid(colourAtBase)) {
				validator.isValid(colourAtBase);
				System.err.println("Invalid state proposed");
			}
		}
		
		// candidate set of infections to move may change after moving the infection, 
		// so the HR must include that information		
		// Therefore, we determine number of eligible infections after proposal
		int afterEligbleInfectionCount = 0;
		for (int i : blockCount.getValues()) {
			if (i == 0) {
				afterEligbleInfectionCount += 1;
			} else if (i > 0) {
				afterEligbleInfectionCount += 2;
			}
		}
		logHR += Math.log(afterEligbleInfectionCount) - Math.log(eligbleInfectionCount);
		
		
//		ColouredTreeLogger logger = new ColouredTreeLogger();
//		logger.initByName("tree", tree, "blockcount", blockCount, "blockstart", blockStartFraction, "blockend", blockEndFraction);
//		String newick = logger.toString();
//		System.out.println(newick);
		
		return logHR;
	}
	
	private void updateStats(int pre, int post) {
		stats[1+pre][1+post]++;
		if (stats[0][0] > 0 && (stats[0][0] + stats[0][1] + stats[1][0] + stats[1][1] + stats[1][2] + stats[2][1] + stats[2][2]) % 1000000 == 0) {
			StringBuilder b = new StringBuilder();
			b.append(stats[0][0] + " " + stats[0][1] + ";");
			b.append(stats[1][0] + " " + stats[1][1] + " " + stats[1][2] + ";");
			b.append(stats[2][1] + " " + stats[2][2] + ";");

			b.append(stats[0][0] + stats[0][1] + " " );
			b.append(stats[1][0] + stats[1][1] + stats[1][2]+ " " );
			b.append(stats[2][1] + stats[2][2]);
			b.append("\n");
			
			System.out.println(b.toString());
		}
	}

	private double doproposal() {
		colourAtBase = likelihood.getColouring();
		int leafCount = tree.getLeafNodeCount();
		boolean [] topOfBlock = new boolean[1];
		int i = randomlySelectInfection(topOfBlock);

		int bc = blockCount.getValue(Math.abs(i));
		if (bc == 0) {
			Node parent = tree.getNode(i).getParent();
			if (colourAtBase[i] < leafCount && parent != null && colourAtBase[parent.getNr()] < leafCount) {
				// infection is on border of two sampled colours
				return moveInfectionOnBorderOfTwoSampledColours(i, colourAtBase[parent.getNr()]);
			} else {
				// infection is on border of at least one unsampled colour
				return moveInfectionAtSegment(i);
			}
		} else {
			if (topOfBlock[0]) {
				return moveInfectionAtTopOfSegment(i);
			} else {
				return moveInfectionAtBottomOfSegment(i);
			}
		}
	}

	
	enum segmentType {top, bottom, all};
	DecimalFormat format = new DecimalFormat("#.##");
	class Segment {
		int nodeNr;
		double length;
		segmentType type;
		Segment(int nodeNr, double length, segmentType type) {
			this.nodeNr = nodeNr;
			this.length = length;
			this.type = type;
		}
		
		@Override
		public String toString() {
			return  nodeNr + ":" + type + format.format(length);
		}
	}
	
	private double moveInfectionAtTopOfSegment(int nodeNr) {
		// 1. remove infection & remove uniformly part of block at end
		blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) - 1);
		double deltaLength = 0;
		if (blockCount.getValue(nodeNr) == 0) {
			deltaLength = (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr)) * tree.getNode(nodeNr).getLength();
			blockEndFraction.setValue(nodeNr, blockStartFraction.getValue(nodeNr));
			
		} else {
			double newBlockEndFraction = blockEndFraction.getValue(nodeNr) - Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr));
			deltaLength = (blockEndFraction.getValue(nodeNr) - newBlockEndFraction) * tree.getNode(nodeNr).getLength();
			blockEndFraction.setValue(nodeNr, newBlockEndFraction);
		}

		// 2. determine segment
		int parent = tree.getNode(nodeNr).getParent().getNr();
		int parentColour = colourAtBase[parent];
		List<Segment> segments = getSegments(parentColour);
		segments.addAll(getSegments(colourAtBase[nodeNr]));

		// 3. insert randomly in segment
		double length = insertInfectionIntoSegments(segments);
		return Math.log((length - deltaLength)/length);
	}

	private double moveInfectionAtBottomOfSegment(int nodeNr) {
		// 1. remove infection & remove uniformly part of block at start
		blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) - 1);
		double deltaLength = 0;
		if (blockCount.getValue(nodeNr) == 0) {
			deltaLength = (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr)) * tree.getNode(nodeNr).getLength();
			blockStartFraction.setValue(nodeNr, blockEndFraction.getValue(nodeNr));
		} else {
			double newBlockStartFraction = blockStartFraction.getValue(nodeNr) + Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) -blockStartFraction.getValue(nodeNr));
			deltaLength = (newBlockStartFraction - blockStartFraction.getValue(nodeNr)) * tree.getNode(nodeNr).getLength();
			blockStartFraction.setValue(nodeNr, newBlockStartFraction);
		}

		// 2. determine segment
		List<Segment> segments = getSegments(colourAtBase[nodeNr]);
		int parent = tree.getNode(nodeNr).getParent().getNr();
		int parentColour = colourAtBase[parent];
		segments.addAll(getSegments(parentColour));
		
		// 3. insert randomly in segment
		double length = insertInfectionIntoSegments(segments);
		return Math.log((length - deltaLength)/length);
	}

	// infection is on border of at least one unsampled colour
	private double moveInfectionAtSegment(int nodeNr) {
		// 2. determine segment
		List<Segment> segments = getSegments(colourAtBase[nodeNr]);
		int parent = tree.getNode(nodeNr).getParent().getNr();
		int parentColour = colourAtBase[parent];
		segments.addAll(getSegments(parentColour));
		
//System.out.print(Arrays.toString(segments.toArray()));		

		// 1. remove infection
		blockCount.setValue(nodeNr, -1);

		// 3. insert randomly in segment
		insertInfectionIntoSegments(segments);
		
		return 0;
	}

	private List<Segment> getSegments(int colour) {
		List<Segment> list = new ArrayList<>();
		int n = tree.getNodeCount();
		for (int j = 0; j < n; j++) {
			if (colourAtBase[j] == colour) {
				if (blockCount.getValue(j) == -1) {
					list.add(new Segment(j, tree.getNode(j).getLength(), segmentType.all));
				} else {
					list.add(new Segment(j, tree.getNode(j).getLength() * blockStartFraction.getValue(j), segmentType.bottom));
				}
				if (!tree.getNode(j).isLeaf()) {
					Node left = tree.getNode(j).getLeft();
					int leftIndex = left.getNr();
					if (blockCount.getValue(leftIndex) != -1) {
						list.add(new Segment(left.getNr(), left.getLength() * (1.0-blockEndFraction.getValue(leftIndex)), segmentType.top));					
					}
					Node right = tree.getNode(j).getRight();
					int rightIndex = right.getNr();
					if (blockCount.getValue(rightIndex) != -1) {
						list.add(new Segment(right.getNr(), right.getLength() * (1.0-blockEndFraction.getValue(rightIndex)), segmentType.top));					
					}
				}
			}
		}
		return list;
	}

	private double insertInfectionIntoSegments0(List<Segment> segments) {
		// randomly select non-zero length segment
		int i = Randomizer.nextInt(segments.size());
		Segment segment = segments.get(i);
		while (segment.length <= 0) {
			i = Randomizer.nextInt(segments.size());
			segment = segments.get(i);
		}
		double r = Randomizer.nextDouble() * segment.length;
		blockCount.setValue(segment.nodeNr, blockCount.getValue(segment.nodeNr) + 1);
		Node node = tree.getNode(segment.nodeNr);
		switch (segment.type) {
		case top:
			blockEndFraction.setValue(segment.nodeNr, 1.0 - r / node.getLength());
			if (blockCount.getValue(segment.nodeNr) == 0) {
				blockStartFraction.setValue(segment.nodeNr, 1.0 - r / node.getLength());
			}
			break;
		case bottom:
			blockStartFraction.setValue(segment.nodeNr, r / node.getLength());
			if (blockCount.getValue(segment.nodeNr) == 0) {
				blockEndFraction.setValue(segment.nodeNr, r / node.getLength());
			}
			break;
		case all:
			if (blockCount.getValue(segment.nodeNr) != 0) {
				throw new RuntimeException("Expected block count to be 0");
			}
			blockStartFraction.setValue(segment.nodeNr, r / node.getLength());
			blockEndFraction.setValue(segment.nodeNr, r / node.getLength());
			break;
		}
		return segment.length;
	}

	private double insertInfectionIntoSegments(List<Segment> segments) {
		double length = 0;
		for (Segment segment : segments) {
			length += segment.length;
		}
		double r = Randomizer.nextDouble() * length;
		for (Segment segment : segments) {
			if (segment.length > r) {
				blockCount.setValue(segment.nodeNr, blockCount.getValue(segment.nodeNr) + 1);
				Node node = tree.getNode(segment.nodeNr);
				switch (segment.type) {
				case top:
					blockEndFraction.setValue(segment.nodeNr, 1.0 - r / node.getLength());
					if (blockCount.getValue(segment.nodeNr) == 0) {
						blockStartFraction.setValue(segment.nodeNr, 1.0 - r / node.getLength());
					}
					break;
				case bottom:
					blockStartFraction.setValue(segment.nodeNr, r / node.getLength());
					if (blockCount.getValue(segment.nodeNr) == 0) {
						blockEndFraction.setValue(segment.nodeNr, r / node.getLength());
					}
					break;
				case all:
					if (blockCount.getValue(segment.nodeNr) != 0) {
						throw new RuntimeException("Expected block count to be 0");
					}
					blockStartFraction.setValue(segment.nodeNr, r / node.getLength());
					blockEndFraction.setValue(segment.nodeNr, r / node.getLength());
					break;
				}
				return length;
			}
			r = r - segment.length;
		}
		throw new RuntimeException("Progammer error 1: should never get here");
	}

	private double moveInfectionOnBorderOfTwoSampledColours(int nodeNr, int otherColour) {
		// 1. remove infection
		blockCount.setValue(nodeNr, -1);

		// 2. determine path between node i and node colourAtBase[parent.getNr()]
		int colour = colourAtBase[nodeNr];
		List<Node> path = getPathExcludingMRCA(colour, otherColour);

		// 3. insert infection uniform randomly at path between node i and node colourAtBase[parent.getNr()]
		double length = 0;
		for (Node n : path) {
			length += n.getLength();
		}
		double r = Randomizer.nextDouble() * length;
		for (Node n : path) {
			if (n.getLength() > r) {
				if (blockCount.getValue(n.getNr()) != -1) {
					throw new RuntimeException("Expected block count to be -1");
				}
				blockCount.setValue(n.getNr(), 0);
				blockStartFraction.setValue(n.getNr(), r / n.getLength());
				blockEndFraction.setValue(n.getNr(), r / n.getLength());
				return 0;
			}
			r = r - n.getLength();
		}
		throw new RuntimeException("Progammer error 2: should never get here");
	}

	// randomly pick an infection that can be moved (i.e. at top or bottom of a block, but not insider a block)
	// returns node number with infection
	// use topOfBlock to indicate using top of block (true), or for bottom of block (false)
	private int eligbleInfectionCount = 0;
	private int randomlySelectInfection(boolean [] topOfBlock) {
		// 1. determine number of eligible infections 
		eligbleInfectionCount = 0;
		for (int i : blockCount.getValues()) {
			if (i == 0) {
				eligbleInfectionCount += 1;
			} else if (i > 0) {
				eligbleInfectionCount += 2;
			}
		}

		// 2. pick one uniformly at random;
		int k = Randomizer.nextInt(eligbleInfectionCount);
		for (int i = 0; i < blockCount.getDimension(); i++) {
		//for (int i = blockCount.getDimension()-1; i >= 0; i--) {
			int bc = blockCount.getValue(i);
			if (bc == 0) {
				k--;
				if (k < 0) {
					topOfBlock[0] = true;
					return i;
				}
			} else if (bc > 0) {
				k -= 2;
				if (k < 0) {
					topOfBlock[0] = true;
					return i;
				}
				if (k < -1) {
					topOfBlock[0] = false;
					return i;
				}
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}

	
    protected List<Node> getPathExcludingMRCA(int nodeIndex1, int nodeIndex2) {
    	Node n1 = tree.getNode(nodeIndex1); 
    	Node n2 = tree.getNode(nodeIndex2);
    	List<Node> path = new ArrayList<>();
    	path.add(n1);
    	path.add(n2);
    	
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            path.add(n1);
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            path.add(n2);
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();
                }
	            path.add(n);
	        }
        }
        // remove MRCA itself
        while (path.remove(n1)) {}
        return path;
    }

}
