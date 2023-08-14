package transmission.distribution;

import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.GammaDistribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Function.Constant;
import beast.base.core.Input.Validate;

@Description("Defines a hazard function based on a gamma distribution")
public class GammaHazardFunction extends HazardFunction {
    final public Input<Function> A_trInput = new Input<>("A", "scale of gamma hazard function", new Constant("1.0"));
    final public Input<Function> p_trInput = new Input<>("p", "shape of gamma hazard function", Validate.REQUIRED);
    final public Input<Function> constantInput = new Input<>("C", "constant of tranmission process", new Constant("1.0"));

    private GammaDistribution samplingDist = new GammaDistributionImpl(1, 1);
    
    private Function scale;
	private Function shape;
	private Function constant;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		scale = A_trInput.get();
		shape = p_trInput.get();
		constant = constantInput.get();
	}
	
	@Override
	public double logS(double t, double d) {
		samplingDist.setAlpha(scale.getArrayValue());
		samplingDist.setBeta(shape.getArrayValue());
		final double logS = -constant.getArrayValue() * samplingDist.density(t - d);
		return logS;
	}

	@Override
	public double logTr(double t, double d) {
		samplingDist.setAlpha(scale.getArrayValue());
		samplingDist.setBeta(shape.getArrayValue());
		final double logTr = Math.log(constant.getArrayValue()) + samplingDist.logDensity(t - d);
		return logTr;
	}

}
