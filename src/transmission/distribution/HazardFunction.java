package transmission.distribution;

import org.apache.commons.math.MathException;

import beast.base.core.BEASTObject;

abstract public  class HazardFunction extends BEASTObject {

	@Override
	public void initAndValidate() {
	}
	
	public abstract double logS(double t, double d);
	public abstract double logH(double t, double d);
	public abstract double simulate() throws MathException;

}
