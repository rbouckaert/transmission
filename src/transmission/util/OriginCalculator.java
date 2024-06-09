package transmission.util;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;

@Description("calculates origin time from tree and ")
public class OriginCalculator extends CalculationNode implements Function, Loggable {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");
    final public Input<Function> deltaStartTimeInput = new Input<>("deltaStartTime", "time at which the study start till root of tree", new Constant("0.0"));

    private TreeInterface tree;
    private Function deltaStartTime;
    
	@Override
	public void initAndValidate() {
		deltaStartTime = deltaStartTimeInput.get();
		tree = treeInput.get();
	}

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		return tree.getRoot().getHeight() + deltaStartTime.getArrayValue();
	}

	@Override
	public void init(PrintStream out) {
		out.print(getID()+ "\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(getArrayValue() + "\t");
	}

	@Override
	public void close(PrintStream out) {
	}

}
