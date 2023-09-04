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
    final public Input<Boolean> keepConstantCountInput = new Input<>("keepconstantcount", "if true, for every deleting there is an insertion to keep total sum of block counts constant", false);

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
		
		if (Randomizer.nextBoolean()) {
			int i = Randomizer.nextInt(blockStartFraction.getDimension());
			// only move start and end fraction but not block count
			switch (blockCount.getValue(i)) {
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
		
		if (keepConstantCountInput.get()) {
			int i = Randomizer.nextInt(blockCount.getDimension());
			if (blockCount.getValue(i) == -1) {
				return -1;
			}
			double logHR = removeInfection(i);
			i = Randomizer.nextInt(blockCount.getDimension());
			logHR += insertInfection(i);
if (blockCount.getValue(0) + blockCount.getValue(1) != 0) {
	int h = 3;
	h++;
}
			return logHR;
		} else	if (Randomizer.nextBoolean()) {
			int i = Randomizer.nextInt(blockCount.getDimension());
			return removeInfection(i);
		} else {
			int i = Randomizer.nextInt(blockCount.getDimension());
			return insertInfection(i);
		}
		
	}

	private double insertInfection(int i) {
		switch (blockCount.getValue(i)) {
		case -1:
			blockCount.setValue(i, 0);
			double f = Randomizer.nextDouble();
			blockStartFraction.setValue(i, f);
			blockEndFraction.setValue(i, f);
			return 0;
			
		case 0:
			// add infection
			blockCount.setValue(i, 1);

			f = Randomizer.nextDouble();
			double mid = blockStartFraction.getValue(i);
			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
				return Double.NEGATIVE_INFINITY;
			}
			blockStartFraction.setValue(i, mid - f / 2.0);
			blockEndFraction.setValue(i, mid + f / 2.0);
			return 0;
			
		case 1:
			// add infection
			blockCount.setValue(i, 2);

			f = Randomizer.nextDouble();
			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
				return Double.NEGATIVE_INFINITY;
			}
			blockStartFraction.setValue(i, mid - f / 2.0);
			blockEndFraction.setValue(i, mid + f / 2.0);				
			return Math.log(2.0);
		
		default:
			// add infection
			blockCount.setValue(i, blockCount.getValue(i)+1);
			f = Randomizer.nextDouble();
			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
				return Double.NEGATIVE_INFINITY;
			}
		}
		return 0;
		
	}

	private double removeInfection(int i) {
		switch (blockCount.getValue(i)) {
		case -1:
			// do nothing
			return 0; 
			
		case 0:
			// remove infection
			blockCount.setValue(i, -1);
			return Math.log(2.0);
			
		case 1:
			// remove infection
			blockCount.setValue(i, 0);
			double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
			blockStartFraction.setValue(i, mid);
			blockEndFraction.setValue(i, mid);
			return -Math.log(2.0);
		
		default:
			// remove infection
			blockCount.setValue(i, blockCount.getValue(i)-1);
			double f = Randomizer.nextDouble();
			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
				return Double.NEGATIVE_INFINITY;
			}
		}
		
		return 0;
	}

}
