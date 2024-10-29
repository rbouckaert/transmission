package breath.distribution;


import org.apache.commons.math.MathException;

import beast.base.inference.CalculationNode;

abstract public  class HazardFunction extends CalculationNode {

	@Override
	public void initAndValidate() {
	}
	
	public abstract double logS(double t, double d);
	public abstract double logH(double t, double d);
	public abstract double simulate() throws MathException;

}
