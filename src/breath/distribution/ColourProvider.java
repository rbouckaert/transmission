package breath.distribution;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;

public class ColourProvider {

	
	// initialise colourAtBase
	// return true if a valid colouring can be found, 
	// return false if there is a path between leafs without a transmission
	static public boolean getColour(
		     Node root,
		     IntegerParameter blockCount,
		     int leafCount,
		     int [] colourAtBase
			) {
		
		colourAtBase[root.getNr()] = root.getNr();
		calcColourAtBase(root, colourAtBase, blockCount);
		// normalise colours so leaf i has colour i
		// but unsampled nodes remain at their colour number
		int n = colourAtBase.length;
		int [] permutation = new int[n];
		for (int i = 0; i < n; i++) {
			permutation[i] = i;
		}
		for (int i = 0; i < leafCount; i++) {
			int j = colourAtBase[i];
			if (j >= leafCount && permutation[j] < leafCount) {
				// we already assigned permutation[j] to another leaf
				// so there must be a path without transmission between that leaf
				// and leaf i, i.e. this is not a valid colouring
				return false;
			}
			permutation[j] = i;
		}
//		int j = leafCount;
//		for (int i = 0; i < n; i++) {
//			if (permutation[i] == 2*n+i) {
//				permutation[i] = j++;
//			}
//		}
		for (int i = 0; i < n; i++) {
			colourAtBase[i] = permutation[colourAtBase[i]];
		}
		return true;
	}

	static private void calcColourAtBase(Node node, int [] colourAtBase, IntegerParameter blockCount) {
		if (!node.isRoot()) {
			int k = node.getNr();
			if (blockCount.getArrayValue(node.getNr()) < 0) {
				colourAtBase[k] = colourAtBase[node.getParent().getNr()];
			} else {
				colourAtBase[k] = node.getNr();
			}
		}
		for (Node child : node.getChildren()) {
			calcColourAtBase(child, colourAtBase, blockCount);
		}
	}
}
