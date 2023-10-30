package transmission.distribution;

import org.apache.commons.math.distribution.GammaDistributionImpl;

import java.text.DecimalFormat;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Function.Constant;
import beast.base.core.Input.Validate;
import beast.base.util.Randomizer;

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
	
	DecimalFormat f = new DecimalFormat("#.####");
	DecimalFormat f4 = new DecimalFormat("#.######");
	
	@Override
	public double logS(double t, double d) {
		update();
		double logS = 0;
		try {
			logS = -constant.getArrayValue() * samplingDist.cumulativeProbability(t - d);
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
System.err.println("logS" + getID() + "(" + f.format(t) + "-" + f.format(d) + "=" + f.format(t-d) + ") \t= " + f4.format(logS));		
//System.err.println("logS" + getID() + "(" + f.format(t-d) + ") # = " + f4.format(logS));		
		return logS;
	}

	private void update() {
		samplingDist.setAlpha(shape.getArrayValue());
		if (rate != null) {
			samplingDist.setBeta(1.0/rate.getArrayValue());
		} else {
			samplingDist.setBeta(scale.getArrayValue());
		}
	}

	@Override
	public double logH(double t, double d) {
		update();
		final double logH = Math.log(constant.getArrayValue()) + samplingDist.logDensity(t - d);
System.err.println("logh" + getID() + "(" + f.format(t) + "-" + f.format(d) + "=" + f.format(t-d) + ") = " + f4.format(logH));		
//System.err.println("logh" + getID() + "(" + f.format(t-d) + ") # = " + f4.format(logH));		
		return logH;
	}

	
	@Override
	public double simulate() throws MathException {
		double p = Randomizer.nextDouble();
		update();
		double t = samplingDist.inverseCumulativeProbability(p);
		return t;
	}
	
	public double getRate() {
		if (rate != null) {
			return rate.getArrayValue();
		} else {
			return 1.0/scale.getArrayValue();
		}
	}
}
