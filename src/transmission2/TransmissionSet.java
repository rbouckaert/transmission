package transmission2;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import breath.distribution.ColourProvider;

@Description("Helper class to convert set of transmissions encoded as node numbes and branch fractions to per node lists of transmissions")
public class TransmissionSet extends CalculationNode {
    final public Input<IntegerParameter> nodeNrInput = new Input<>("nodeNr", "number of node having tranmission in the branch above", Validate.REQUIRED);
    final public Input<RealParameter> branchFractionInput = new Input<>("branchFraction", "fraction of branch length at which the transmission happens", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

    
    private IntegerParameter nodeNr;
    private RealParameter branchFraction;
    private TreeInterface tree;
    
	private double [][] transmissionsForNode;
	private boolean needsUpdate = true;

	private DecimalFormat f = new DecimalFormat("#.####");

	@Override
	public void initAndValidate() {
		nodeNr = nodeNrInput.get();
		branchFraction = branchFractionInput.get();
		tree = treeInput.get();
		transmissionsForNode = new double[tree.getNodeCount()][0];
	}
	
	
	
	public double [] getTransmissionForNode(int nodeNr) {
		if (needsUpdate) {
			update();
		}
		
		return transmissionsForNode[nodeNr];
	}
	
	
	private void update() {
		Integer [] counts = new Integer[transmissionsForNode.length];
		// allocate memory
		for (int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		for (Integer i : this.nodeNr.getValues()) {
			counts[i]++;
		}
		for (int i = 0; i < counts.length; i++) {
			if (transmissionsForNode[i].length != counts[i]) {
				transmissionsForNode[i] = new double[counts[i]];
			}
		}
		for (int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		for (int k = 0; k < this.nodeNr.getDimension(); k++) { 
			Integer i = this.nodeNr.getValue(k);
			transmissionsForNode[i][counts[i]] = branchFraction.getValue(k);
			counts[i]++;
		}
		for (int k = 0; k < transmissionsForNode.length; k++) {
			sort(transmissionsForNode[k]);
		}
		needsUpdate = false;
	}



	private void sort(double[] x) {
		switch (x.length) {
		case 0:
			break;
		case 1:
			break;
		case 2:
			double tmp;
			if (x[0]>x[1]) {tmp = x[0]; x[0] = x[1]; x[1] = tmp;}
			break;
		case 3:
			if (x[0] > x[1]) {tmp = x[0]; x[0] = x[1]; x[1] = tmp;}
			if (x[0] > x[2]) {tmp = x[0]; x[0] = x[2]; x[2] = tmp;}
			if (x[1] > x[2]) {tmp = x[1]; x[1] = x[2]; x[2] = tmp;}
			break;
		case 4:
			if (x[0]>x[1]) {tmp = x[0]; x[0] = x[1]; x[1] = tmp;}
			if (x[2]>x[3]) {tmp = x[2]; x[2] = x[3]; x[3] = tmp;}
			if (x[0]>x[2]) {tmp = x[0]; x[0] = x[2]; x[2] = tmp;}
			if (x[1]>x[3]) {tmp = x[1]; x[1] = x[3]; x[3] = tmp;}
			if (x[1]>x[2]) {tmp = x[1]; x[1] = x[2]; x[2] = tmp;}
			break;
		default:
			Arrays.sort(x);
		}
	}



	// initialise colourAtBase
	// return true if a valid colouring can be found, 
	// return false if there is a path between leafs without a transmission
	public boolean getColour(
		     int [] colourAtBase
			) {
		int leafCount = tree.getLeafNodeCount();
		Integer [] counts = new Integer[leafCount * 2 - 1];
		for (int i = 0; i < counts.length; i++) {
			counts[i] = -1;
		}
		for (Integer i : nodeNr.getValues()) {
			counts[i]++;
		}
	    IntegerParameter blockCount = new IntegerParameter(counts);
	    return ColourProvider.getColour(tree.getRoot(), blockCount, leafCount, colourAtBase);
	}
	
	

	@Override
	protected void store() {
		super.store();
	}
	
	@Override
	protected void restore() {
		needsUpdate = true;
		super.restore();
	}
	
	@Override
	protected boolean requiresRecalculation() {
		needsUpdate = true;
		return super.requiresRecalculation();
	}



	public double getMostRecentTransmission(int nr) {
		if (needsUpdate) {
			update();
		}
		return transmissionsForNode[nr][0];
	}



	public double getLeastRecentTransmission(int nr) {
		if (needsUpdate) {
			update();
		}
		return transmissionsForNode[nr][transmissionsForNode[nr].length-1];
	}



	public int infectionCount(int nr) {
		if (needsUpdate) {
			update();
		}
		return transmissionsForNode[nr].length;
	}
	

	@Override
	public String toString() {
		if (needsUpdate) {
			update();
		}
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < transmissionsForNode.length; i++) {
			if (transmissionsForNode[i].length > 0) {
				b.append(i +":");
				for (double d : transmissionsForNode[i]) {
					b.append(f.format(d) + " ");
				}
				b.append("\n");
			}
		}
		return b.toString();
	}


}
