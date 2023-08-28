package transmission.operator;


import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import transmission.distribution.TransmissionTreeLikelihood;
import transmission.distribution.Validator;

@Description("Operator that randomly picks an infection and moves it elsewhere")
public class InfectionMover extends Operator {
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

	@Override
	public double proposal() {
		double logHR = doproposal();
		colourAtBase = likelihood.getFreshColouring();
		Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
		if (!validator.isValid(colourAtBase)) {
			validator.isValid(colourAtBase);
			System.err.println("Invalid state proposed");
		}
		return logHR;
	}
	
	private double doproposal() {
		colourAtBase = likelihood.getColouring();
		int leafCount = tree.getLeafNodeCount();
		
		// TODO: candidate set may change after moving the infection, so the HR must include that information		
		boolean [] topOfBlock = new boolean[1];
		int i = randomlySelectInfection(topOfBlock);
		int bc = blockCount.getValue(Math.abs(i));
		if (bc == 0) {
			Node parent = tree.getNode(i).getParent();
			if (colourAtBase[i] < leafCount && parent != null && colourAtBase[parent.getNr()] < leafCount) {
				// infection is on border of two sampled colours
				return moveInfectionOnBorderOfTwoSampledColours(i, parent);
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
			return  nodeNr + ":" + type;
		}
	}
	
	private double moveInfectionAtTopOfSegment(int i) {
		// 2. determine segment
		int parent = tree.getNode(i).getParent().getNr();
		int parentColour = colourAtBase[parent];
		List<Segment> segments = getSegments(parentColour);
		
		// 1. remove infection & remove uniformly part of block at end
		blockCount.setValue(i, blockCount.getValue(i) - 1);
		double deltaLength = 0;
		if (blockCount.getValue(i) == 0) {
			deltaLength = blockEndFraction.getValue(i) - blockStartFraction.getValue(i);
			blockEndFraction.setValue(i, blockStartFraction.getValue(i));
			
		} else {
			double newBlockStartFraction = blockEndFraction.getValue(i) - Randomizer.nextDouble() * (blockEndFraction.getValue(i) - blockStartFraction.getValue(i));
			deltaLength = blockEndFraction.getValue(i) - newBlockStartFraction;
			blockEndFraction.setValue(i, newBlockStartFraction);
		}

		// 3. insert randomly in segment
		double length = insertInfectionIntoSegments(segments);
		return (length + deltaLength)/length;
	}

	private double moveInfectionAtBottomOfSegment(int i) {
		// 2. determine segment
		List<Segment> segments = getSegments(colourAtBase[i]);
		
		// 1. remove infection & remove uniformly part of block at start
		blockCount.setValue(i, blockCount.getValue(i) - 1);
		double deltaLength = 0;
		if (blockCount.getValue(i) == 0) {
			deltaLength = blockEndFraction.getValue(i) - blockStartFraction.getValue(i);
			blockStartFraction.setValue(i, blockEndFraction.getValue(i));
		} else {
			double newBlockEndFraction = blockStartFraction.getValue(i) + Randomizer.nextDouble() * (blockEndFraction.getValue(i) -blockStartFraction.getValue(i));
			deltaLength = newBlockEndFraction - blockStartFraction.getValue(i);
			blockStartFraction.setValue(i, newBlockEndFraction);
		}
		
		// 3. insert randomly in segment
		double length = insertInfectionIntoSegments(segments);
		return (length + deltaLength)/length;
	}

	// infection is on border of at least one unsampled colour
	private double moveInfectionAtSegment(int i) {
		// 2. determine segment
		List<Segment> segments = getSegments(colourAtBase[i]);
		int parent = tree.getNode(i).getParent().getNr();
		int parentColour = colourAtBase[parent];
		segments.addAll(getSegments(parentColour));

		// 1. remove infection
		blockCount.setValue(i, -1);

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
					Node right = tree.getNode(j).getRight();
					int rightIndex = right.getNr();
					if (blockCount.getValue(rightIndex) != -1) {
						list.add(new Segment(right.getNr(), right.getLength() * (1.0-blockEndFraction.getValue(rightIndex)), segmentType.top));					
					}
					Node left = tree.getNode(j).getLeft();
					int leftIndex = left.getNr();
					if (blockCount.getValue(leftIndex) != -1) {
						list.add(new Segment(left.getNr(), left.getLength() * (1.0-blockEndFraction.getValue(leftIndex)), segmentType.top));					
					}
				}
			}
		}
		return list;
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
		throw new RuntimeException("Progammer error 2: should never get here");
	}

	private double moveInfectionOnBorderOfTwoSampledColours(int i, Node parent) {
		// 1. remove infection
		blockCount.setValue(i, -1);

		// 2. determine path between node i and node colourAtBase[parent.getNr()]
		List<Node> path = getPathExcludingMRCA(i, colourAtBase[parent.getNr()]);

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
	private int randomlySelectInfection(boolean [] topOfBlock) {
		// 1. determine number of eligible infections 
		int eligbleInfectionCount = 0;
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
			int j = blockCount.getValue(i);
			if (j == 0) {
				k--;
				if (k <= 0) {
					topOfBlock[0] = true;
					return i;
				}
			} else if (j > 0) {
				k -= 2;
				if (k == 0) {
					topOfBlock[0] = true;
					return i;
				}
				if (k < 0) {
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
