package transmission.operator;


import java.util.ArrayList;
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

    
    final static boolean debug = false;
        
	@Override
	public double proposal() {
		//System.out.print(blockCount.getValue());
		
		// randomly pick two distinct leafs
		int n = tree.getLeafNodeCount();
		int i = Randomizer.nextInt(n);
		int j = Randomizer.nextInt(n);
		while (j == i) {
			j = Randomizer.nextInt(n);
		}
		
		// get path between two leafs
		List<Node> path = getPathExcludingMRCA(i, j);
		
		
		removeInfectionFromPath(path);
		insertInfectionToPath(path);
		
		
		// make sure the colouring is valid
		colourAtBase = likelihood.getFreshColouring();		
		Validator validator = new Validator((Tree)tree, colourAtBase, blockCount, blockStartFraction, blockEndFraction);
		if (!validator.isValid(colourAtBase)) {
			return Double.NEGATIVE_INFINITY;
		}

		//System.out.print("=>"+blockCount.getValue());

		return 0;
	}
	
	private void removeInfectionFromPath(List<Node> path) {
		// 1. determine number of eligible infections 
		int eligbleInfectionCount = 0;
		for (Node node : path) {
			int i = blockCount.getValue(node.getNr());
			if (i == 0) {
				eligbleInfectionCount += 1;
			} else if (i > 0) {
				eligbleInfectionCount += 2;
			}
		}

		// 2. pick one uniformly at random;
		int k = Randomizer.nextInt(eligbleInfectionCount);
		for (Node node : path) {
			int nodeNr = node.getNr();
			int bc = blockCount.getValue(nodeNr);
			if (bc == 0) {
				k--;
				if (k < 0) {
					// remove from top of block
					blockCount.setValue(nodeNr, -1);
					return;
				}
			} else if (bc > 0) {
				k -= 2;
				if (k < 0) {
					// remove from bottom of block
					blockCount.setValue(node.getNr(), blockCount.getValue(node.getNr()) - 1);
					if (blockCount.getValue(nodeNr) == 0) {
						blockStartFraction.setValue(nodeNr, blockEndFraction.getValue(nodeNr));						
					} else if (blockCount.getValue(nodeNr) > 0) {
						double newBlockStartFraction = blockStartFraction.getValue(nodeNr) + Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) -blockStartFraction.getValue(nodeNr));
						blockStartFraction.setValue(nodeNr, newBlockStartFraction);
					}
					return;
				}
				if (k < -1) {
					// remove from top of block
					blockCount.setValue(node.getNr(), blockCount.getValue(node.getNr()) - 1);
					if (blockCount.getValue(nodeNr) == 0) {
						blockEndFraction.setValue(nodeNr, blockStartFraction.getValue(nodeNr));						
					} else if (blockCount.getValue(nodeNr) > 0) {
						double newBlockEndFraction = blockEndFraction.getValue(nodeNr) - Randomizer.nextDouble() * (blockEndFraction.getValue(nodeNr) - blockStartFraction.getValue(nodeNr));
						blockEndFraction.setValue(nodeNr, newBlockEndFraction);
					}
					return;
				}
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}	

	private void insertInfectionToPath(List<Node> path) {
		// insert infection uniform randomly on path 
		double length = 0;
		for (Node n : path) {
			length += n.getLength();
		}
		
		double r = Randomizer.nextDouble() * length;
		for (Node node : path) {
			if (node.getLength() > r) {
				int nodeNr = node.getNr();
				blockCount.setValue(nodeNr, blockCount.getValue(nodeNr) + 1);
				if (blockCount.getValue(nodeNr) == 0) {
					blockStartFraction.setValue(nodeNr, r / node.getLength());
					blockEndFraction.setValue(nodeNr, r / node.getLength());
					return;
				} else { // blockCount > 0
					if (blockStartFraction.getValue(nodeNr) > r) {
						blockStartFraction.setValue(nodeNr, r / node.getLength());
						return;
					}
					if (blockEndFraction.getValue(nodeNr) < r) {
						blockEndFraction.setValue(nodeNr, r / node.getLength());						
					}
					// (blockStartFraction.getValue(nodeNr) < r && r < blockEndFraction.getValue(nodeNr))
					return;
				}
			}
			r = r - node.getLength();
		}
		throw new RuntimeException("Progammer error 2: should never get here");
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
