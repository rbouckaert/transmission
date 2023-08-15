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
    final public Input<Function> scaleInput = new Input<>("scale", "scale of gamma hazard function");
    final public Input<Function> shapeInput = new Input<>("shape", "shape of gamma hazard function");
    final public Input<Function> rate_trInput = new Input<>("rate", "rate of gamma hazard function", Validate.XOR, scaleInput);
    final public Input<Function> constantInput = new Input<>("C", "constant of tranmission process", new Constant("1.0"));

    private GammaDistribution samplingDist = new GammaDistributionImpl(1, 1);
    
    private Function scale;
	private Function shape;
	private Function rate;
	private Function constant;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		scale = scaleInput.get();
		shape = shapeInput.get();
		rate = rate_trInput.get();
		constant = constantInput.get();
	}
	
	@Override
	public double logS(double t, double d) {
		update();
		final double logS = -constant.getArrayValue() * samplingDist.density(t - d);
		return logS;
	}

	private void update() {
		samplingDist.setBeta(shape.getArrayValue());
		if (rate != null) {
			samplingDist.setAlpha(1.0/rate.getArrayValue());
		} else {
			samplingDist.setAlpha(scale.getArrayValue());
		}
	}

	@Override
	public double logTr(double t, double d) {
		update();
		final double logTr = Math.log(constant.getArrayValue()) + samplingDist.logDensity(t - d);
		return logTr;
	}

}
