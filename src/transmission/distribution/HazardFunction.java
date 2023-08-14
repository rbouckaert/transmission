package transmission.distribution;

import beast.base.core.BEASTObject;

abstract public  class HazardFunction extends BEASTObject {

	@Override
	public void initAndValidate() {
	}
	
	public abstract double logS(double t, double d);
	public abstract double logTr(double t, double d);

}
