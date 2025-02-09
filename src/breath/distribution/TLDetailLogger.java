package breath.distribution;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;

@Description("Logger of component details of the transmission tree likelihood")
public class TLDetailLogger extends CalculationNode implements Loggable {
	public Input<TransmissionTreeLikelihood> tlInput = new Input<>("breathLikelihood", "BREATH transmission likelihood to log details from", Validate.REQUIRED);

	private TransmissionTreeLikelihood tl;
	
	@Override
	public void initAndValidate() {
		tl = tlInput.get();
	}

	@Override
	public void init(PrintStream out) {
		out.print("Coalescent\t");
		out.print("Sampled_hosts\t");
		out.print("Unsampled_hosts\t");
		out.print("Blocks\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(tl.calculateCoalescent() + "\t");
		out.print(tl.calculateSampledHostContribution() + "\t");
		out.print(tl.calculateUnsampledHostContribution() + "\t");
		out.print(tl.calculateBlockContribution() + "\t");
	}

	@Override
	public void close(PrintStream out) {
	}

}
