package transmission.util;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;

@Description("Counts the number of infections")
public class InfectionCount extends CalculationNode implements Function, Loggable{
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);

    private IntegerParameter blockCount;
    private int n;

	
	@Override
	public void initAndValidate() {
		blockCount = blockCountInput.get();
		n = blockCount.getDimension();
		if (n % 2 == 1) {
			n--;
		}
	}


	@Override
	public void init(PrintStream out) {
		out.print("infectionCount\n");
	}


	@Override
	public void log(long sample, PrintStream out) {
		int infectionCount = 0;
		for (int i = 0; i < n; i++) {
			infectionCount += blockCount.getValue(i) + 1; 
		}
		out.print(infectionCount + "\t");
	}


	@Override
	public void close(PrintStream out) {
	}


	@Override
	public int getDimension() {
		return 1;
	}


	@Override
	public double getArrayValue(int dim) {
		int infectionCount = 0;
		for (int i = 0; i < n; i++) {
			infectionCount += blockCount.getValue(i) + 1; 
		}
		return infectionCount;
	}
	
	

}
