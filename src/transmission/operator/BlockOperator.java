package transmission.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Operator that moves block parameters of a transmission tree")
public class BlockOperator extends Operator {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;

    @Override
	public void initAndValidate() {
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
	}

	@Override
	public double proposal() {
		int i = Randomizer.nextInt(blockStartFraction.getDimension());
		
		if (Randomizer.nextBoolean()) {
			switch (blockCount.getValue()) {
			case -1:
				// nothing to do since start and end fractions are ignored
				break;
			case 0:
				// make sure start == end fraction after proposal
				double f = Randomizer.nextDouble();
				blockStartFraction.setValue(i, f);
				blockEndFraction.setValue(i, f);
				break;
			default:
				// boundary move only
				f = Randomizer.nextDouble();
				double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
				if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
					return Double.NEGATIVE_INFINITY;
				}
				blockStartFraction.setValue(i, mid - f / 2.0);
				blockEndFraction.setValue(i, mid + f / 2.0);
			}
			return 0;
		}
		
		
		
		switch (blockCount.getValue(i)) {
		case -1:
			if (Randomizer.nextBoolean()) {
				blockCount.setValue(i, 0);
				double f = Randomizer.nextDouble();
				blockStartFraction.setValue(i, f);
				blockEndFraction.setValue(i, f);
			} else {
				// do nothing
			}
			break; 
			
		case 0:
			if (Randomizer.nextBoolean()) {
				// remove infection
				blockCount.setValue(i, -1);
				return Math.log(2.0);
			} else {
				// add infection
				blockCount.setValue(i, 1);

				double f = Randomizer.nextDouble();
				double mid = blockStartFraction.getValue(i);
				if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
					return Double.NEGATIVE_INFINITY;
				}
				blockStartFraction.setValue(i, mid - f / 2.0);
				blockEndFraction.setValue(i, mid + f / 2.0);
			}
			break;
			
		case 1:
			if (Randomizer.nextBoolean()) {
				// remove infection
				blockCount.setValue(i, 0);
				double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
				blockStartFraction.setValue(i, mid);
				blockEndFraction.setValue(i, mid);
				return -Math.log(2.0);
			} else {
				// add infection
				blockCount.setValue(i, 2);

				double f = Randomizer.nextDouble();
				double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
				if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
					return Double.NEGATIVE_INFINITY;
				}
				blockStartFraction.setValue(i, mid - f / 2.0);
				blockEndFraction.setValue(i, mid + f / 2.0);				
			}
			break;
		
		default:
			if (Randomizer.nextBoolean()) {
				// remove infection
				blockCount.setValue(i, blockCount.getValue()-1);
			} else {
				// add infection
				blockCount.setValue(i, blockCount.getValue()+1);
			}
			double f = Randomizer.nextDouble();
			double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
				return Double.NEGATIVE_INFINITY;
			}
		}
		
		return 0;
	}

}
