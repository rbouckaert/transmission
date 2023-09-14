package transmission.util;

import java.io.PrintStream;
import java.util.Arrays;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.IntegerParameter;
import transmission.distribution.TransmissionTreeLikelihood;

@Description("Reports infector of leaf nodes, or -1 if infected from unsampled case")
public class InfectorOfLogger extends BEASTObject implements Loggable {
	final public Input<TransmissionTreeLikelihood> likelihoodInput = new Input<>("likelihood", "transmission treelikelihood containing the colouring", Validate.REQUIRED);
	final public Input<Boolean> directOnlyInput = new Input<>("directOnly", "consider direct infections only, if false block counts are ignored", true);

	private TransmissionTreeLikelihood likelihood;
	private TreeInterface tree;
	private boolean directOnly;
	
	@Override
	public void initAndValidate() {
		likelihood = likelihoodInput.get();
		tree = likelihood.treeInput.get();
		directOnly = directOnlyInput.get();
	}

	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			out.print("infectorOf." + (i+1) + "\t");
		}
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		int [] colourAtBase = likelihood.getColouring();
		
    	// determine who infected who
    	int n = tree.getLeafNodeCount();
    	int [] infectedBy = new int[n];
    	Arrays.fill(infectedBy, -1);
    	IntegerParameter blockCount = likelihood.blockCountInput.get();
    	for (int i = 0; i < 2 * n - 2; i++) {
    		Node node = tree.getNode(i);
    		Node parent = node.getParent();
    		if (colourAtBase[node.getNr()] < n && colourAtBase[parent.getNr()] < n && 
    				colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
    			if (!directOnly || blockCount.getValue(node.getNr()) == 0) {
    				infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
    			}
    		}
    	}
    	
    	for (int i = 0; i < n; i++) {
    		out.printf(infectedBy[i] + "\t");
    	}
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

}
