package transmission2;

import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

@Description("Provide valid initialisation for transmission set")
public class InfectionInitialiser extends BEASTObject implements StateNodeInitialiser {
    final public Input<TransmissionSet> transmissionsInput = new Input<>("transmissions", "represents set of transmission on the tree", Validate.REQUIRED);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");

	@Override
	public void initAndValidate() {
		initStateNodes();
	}

	@Override
	public void initStateNodes() {
		int n = treeInput.get().getNodeCount();
		IntegerParameter nodeNrs = transmissionsInput.get().nodeNrInput.get();
		RealParameter branchFraction = transmissionsInput.get().branchFractionInput.get();
		nodeNrs.setDimension(n-1);
		branchFraction.setDimension(n-1);
		for (int i = 0; i < n-1; i++) {
			nodeNrs.setValue(i, i);
			branchFraction.setValue(i, 0.5);
		}
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(transmissionsInput.get().branchFractionInput.get());
		stateNodes.add(transmissionsInput.get().nodeNrInput.get());
	}

}
