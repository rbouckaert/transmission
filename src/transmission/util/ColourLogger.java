package transmission.util;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import transmission.distribution.TransmissionTreeLikelihood;

@Description("Reports colouring of a transmission tree")
public class ColourLogger extends BEASTObject implements Function {
	final public Input<TransmissionTreeLikelihood> likelihoodInput = new Input<>("likelihood", "transmission treelikelihood containing the colouring", Validate.REQUIRED);

	private TransmissionTreeLikelihood likelihood;
	
	@Override
	public void initAndValidate() {
		likelihood = likelihoodInput.get();		
	}

	@Override
	public int getDimension() {
		return likelihood.treeInput.get().getNodeCount();
	}

	@Override
	public double getArrayValue(int dim) {
		return likelihood.getColour(dim);
	}

}
