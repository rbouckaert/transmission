package transmission2;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;

@Description("Log block count in TreeWithMetaDataLogger as in original representation")
public class BlockCountLogger extends BEASTObject implements Function {
    final public Input<IntegerParameter> nodeNrInput = new Input<>("nodeNr", "number of node having tranmission in the branch above", Validate.REQUIRED);
	
	private IntegerParameter nodeNr;

	@Override
	public void initAndValidate() {
		nodeNr = nodeNrInput.get();
	}

	@Override
	public int getDimension() {
		return Integer.MAX_VALUE;
	}

	@Override
	public double getArrayValue(int dim) {
		int sum = -1;
		for (int i =0; i < nodeNr.getDimension(); i++) {
			if (nodeNr.getArrayValue(i) == dim) {
				sum++;
			}
		}
		return sum;
	}
	
	@Override
	public String getID() {
		String id = super.getID();
		if (id == null) {
			return "blockcount";
		} else {
			return id;
		}
	}

}
