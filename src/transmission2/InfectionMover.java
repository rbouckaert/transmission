package transmission2;



import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;

@Description("Operator that randomly picks an infection and moves it elsewhere")
public class InfectionMover extends Operator {
    final public Input<TransmissionSet> transmissionsInput = new Input<>("transmissions", "represents set of transmission on the tree", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

    private TreeInterface tree;
    private TransmissionSet transmissions;

    @Override
	public void initAndValidate() {
    	tree = treeInput.get();
    	transmissions = transmissionsInput.get();
	}

    
    final static boolean debug = false;
        
	@Override
	public double proposal() {
		// randomly pick transmission
		int transIndex = Randomizer.nextInt(transmissions.nodeNrInput.get().getDimension());
		
		// pick a leaf below the branch
		int branchNr = transmissions.nodeNrInput.get().getValue(transIndex);
		List<Node> leafs = new ArrayList<>();
		tree.getNode(branchNr).getAllLeafNodes(leafs);
		int k = Randomizer.nextInt(leafs.size());
		int i = leafs.get(k).getNr();
		

		// pick a leaf from above the branch
		boolean [] forbiddenLeafs = new boolean[tree.getLeafNodeCount()];
		for (Node node : leafs) {
			forbiddenLeafs[node.getNr()] = true;
		}
		k = Randomizer.nextInt(tree.getLeafNodeCount() - leafs.size());
		int j = 0;
		while (forbiddenLeafs[j]) {
			j++;
		}
		while (k > 0) {
			j++;
			while (forbiddenLeafs[j]) {
				j++;
			}
			k--;
		}
		
		// get path between two leafs
		List<Node> path = getPathExcludingMRCA(i, j);
		
		moveInfectionAlongPath(path, transIndex);

		return 0;
	}
	
	


	private void moveInfectionAlongPath(List<Node> path, int transIndex) {
		// insert infection uniform randomly on path 
		double length = 0;
		for (Node n : path) {
			length += n.getLength();
		}
		
		double r = Randomizer.nextDouble() * length;
		for (Node node : path) {
			if (node.getLength() > r) {
				int nodeNr = node.getNr();
				transmissions.nodeNrInput.get().setValue(transIndex, nodeNr);
				transmissions.branchFractionInput.get().setValue(transIndex, r / node.getLength());
				return;
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

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(transmissionsInput.get().branchFractionInput.get());
        list.add(transmissionsInput.get().nodeNrInput.get());
        return list;
    }
}
