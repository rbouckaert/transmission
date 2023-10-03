package transmission2;



import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.IntegerDistribution;
import org.apache.commons.math3.util.FastMath;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;


@Description("Distribution of the number X of Bernoulli trials needed to get one success.")
public class GeometricDistribution extends ParametricDistribution {
    final public Input<Function> pInput = new Input<>("p", "probability of successfull trial, defaults to 0.5");
    
    private double probabilityOfSuccess, log1mProbabilityOfSuccess;

    org.apache.commons.math.distribution.Distribution dist = new IntegerDistribution() {
		@Override
		public double probability(int x) {
	        if (x < 0) {
	            return 0.0;
	        } else {
	            return FastMath.exp(log1mProbabilityOfSuccess * x) * probabilityOfSuccess;
	        }
		}
		
		@Override
		public int inverseCumulativeProbability(double p) throws MathException {
			throw new RuntimeException("Not implemented yet");
		}
		
		@Override
		public double cumulativeProbability(int x0, int x1) throws MathException {
			return cumulativeProbability(x1) - cumulativeProbability(x0);
		}
		
		@Override
		public double cumulativeProbability(int x) throws MathException {
			return 1.0 - Math.pow(1-probabilityOfSuccess, x);
		}
		
		@Override
		public double probability(double x) {
	        if (x < 0) {
	            return 0.0;
	        } else {
	            return FastMath.exp(log1mProbabilityOfSuccess * x) * probabilityOfSuccess;
	        }
		}

		@Override
		public double cumulativeProbability(double x) throws MathException {
			return 1.0 - Math.pow(1-probabilityOfSuccess, x);
		}

		@Override
		public double cumulativeProbability(double x0, double x1) throws MathException {
			return cumulativeProbability(x1) - cumulativeProbability(x0);
		}

	};
    	
    
		
    		

    // Must provide empty constructor for construction by XML. Note that this constructor DOES NOT call initAndValidate();
    public GeometricDistribution() {
    }

    public GeometricDistribution(RealParameter p) {
        try {
            initByName("p", p);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initByName p parameter when constructing GeometricDistribution instance.");
        }
    }

    @Override
    public void initAndValidate() {
        refresh();
    }

    
    public double calcLogP(final Function fun) {
        double logP = 0;
        for (int i = 0; i < fun.getDimension(); i++) {
            final double x = fun.getArrayValue(i);
            logP += logDensity(x);
        }
        
//System.err.println(fun.getArrayValue() + " " + logP);       
        return logP;
    }

    @Override
    public double logDensity(double x) {
        final double offset = offsetInput.get();
        x -= offset;
    	return log1mProbabilityOfSuccess * x + FastMath.log(probabilityOfSuccess);
    }
    
    @Override
    public double density(double x) {
        final double offset = offsetInput.get();
        x -= offset;
        return FastMath.exp(log1mProbabilityOfSuccess * x) * probabilityOfSuccess;
    }
    
    /**
     * make sure internal state is up to date *
     */
	void refresh() {
        if (pInput.get() == null) {
        	probabilityOfSuccess = 0.5;
        } else {
        	probabilityOfSuccess = pInput.get().getArrayValue();
            if (probabilityOfSuccess < 0 || probabilityOfSuccess > 1) {
            	probabilityOfSuccess = 0.5;
            }
        }
    	log1mProbabilityOfSuccess = FastMath.log(1-probabilityOfSuccess);
    }

    @Override
    public org.apache.commons.math.distribution.Distribution getDistribution() {
        refresh();
        return dist;
    }
    
    @Override
    public double getMeanWithoutOffset() {
    	refresh();
    	return 1.0/probabilityOfSuccess;
    }
    
} // class GeometricDistribution
