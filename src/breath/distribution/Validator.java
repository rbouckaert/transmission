package breath.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

/** make sure a valid colouring/blockcount assignment is done **/
public class Validator {

	public Validator(Tree tree, int [] colourAtBase, IntegerParameter blockCount, RealParameter blockStartFraction, RealParameter blockEndFraction) {
		this.tree = tree;
		this.colourAtBase = colourAtBase;
		this.blockCount = blockCount;
		this.blockStartFraction = blockStartFraction;
	    this.blockEndFraction = blockEndFraction; 
	}

	private boolean [] nodesTraversed;
	private boolean allowInfectionsAfterSampling = false;
    private int nseen;
    private Tree tree;
    private int [] colourAtBase;
    private IntegerParameter blockCount;
    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    
	/** check whether the colouring is valid, that is
	 * o each leaf i has colour i
	 * o each branch with blockcount > 0 has different colour at base than at parent
	 * o each coloured segment is connected
	 */
	public boolean isValid(int [] colourAtBase) {
		// each leaf i has colour i
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			if (colourAtBase[i] != i) {
				return false;
			}
		}
		
		// each branch with blockcount > 0 has different colour at base than at parent
		for (int i = 0; i < tree.getLeafNodeCount() - 1; i++) {
			if (blockCount.getValue(i) >= 0) {
				Node node = tree.getNode(i);
				int baseColour = colourAtBase[node.getNr()];
				int parentColour = colourAtBase[node.getParent().getNr()];
				if (baseColour == parentColour) {
					return false;
				}
			}
		}
		
		// each block of size 0 has no length and of size > 0 has positive length
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			double start = blockStartFraction.getValue(i);
			double end = blockEndFraction.getValue(i);
			if (blockCount.getValue(i) == 0) {
				if (Math.abs(end-start) > 1e-6) {
					return false;
				}
			}
			if (blockCount.getValue(i) > 0) {
				if (end-start < 0) {
					return false;
				}
			}
		}		
		
		// each coloured segment is connected
		List<Node> [] segments = new List[tree.getNodeCount()];
		for (int i = 0; i < tree.getNodeCount(); i++) {
			int colour = colourAtBase[i];
			if (segments[colour] == null) {
				segments[colour] = new ArrayList<>();
			}
			segments[colour].add(tree.getNode(i));			
		}
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (segments[i] != null) {
				if (segments[i].size() > 1) {
					nseen = 0;
					nodesTraversed = new boolean[tree.getNodeCount()];
					List<Node> nodes = segments[i];
			        Node cur = nodes.get(0);

			        for (int k = 1; k < nodes.size(); ++k) {
			            cur = getCommonAncestor(cur, nodes.get(k));
			        }
			        if (nseen != nodes.size()) {
			        	return false;
			        }
				}
			}
		}

		// do not allow infections after sampling
		// so coloured tree does not exceed sampling time
//		if (!allowInfectionsAfterSampling) {
//			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//				xx
//			}
//		}
	
		return true;
	}
    
    private Node getCommonAncestor(Node n1, Node n2) {
        // assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
            nodesTraversed[n1.getNr()] = true;
            nseen += 1;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
            nseen += 1;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	                nseen += 1;
	            }
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	                nseen += 1;
	            }
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
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	                nseen += 1;
	            } 
	        }
        }
        return n1;
    }

	protected Node getMRCA(Tree tree, Set<String> taxa) {
		List<Node> leafs = new ArrayList<>();
		for (Node node : tree.getExternalNodes()) {
			if (taxa.contains(node.getID())) {
				leafs.add(node);
			}
		}

        nodesTraversed = new boolean[tree.getNodeCount()];
        nseen = 0;
        Node cur = leafs.get(0);

        for (int k = 1; k < leafs.size(); ++k) {
            cur = getCommonAncestor(cur, leafs.get(k));
        }
		return cur;
	}

}
