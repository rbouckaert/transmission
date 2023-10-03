package transmission2;

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
    final public Input<IntegerParameter> nodeNrInput = new Input<>("nodeNr", "number of node having tranmission in the branch above", Validate.REQUIRED);


	
	@Override
	public void initAndValidate() {
	}


	@Override
	public void init(PrintStream out) {
		out.print("infectionCount\n");
	}


	@Override
	public void log(long sample, PrintStream out) {
		out.print(nodeNrInput.get().getDimension() + "\t");
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
		return nodeNrInput.get().getDimension();
	}
}
