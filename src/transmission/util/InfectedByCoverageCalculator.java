package transmission.util;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.inference.Runnable;
import beast.base.util.HeapSort;
import beastfx.app.tools.Application;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.util.LogFile;
import beastfx.app.util.OutFile;

@Description("Calculate coverage of who infected who")
public class InfectedByCoverageCalculator extends Runnable {
	public Input<LogFile> truthInput = new Input<>("truth", "trace file with true infection information", Validate.REQUIRED);
	public Input<File> logFilePrefixInput = new Input<>("logFilePrefix", "log file name without the number and '.log' missing. It is assumed there are as many log files as there are entries in the truth file", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified", new OutFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);	
	final public Input<Double> coverageInput = new Input<>("coverage", "percentage of coverage to be tested against (between 0 and 100, default 95)", 95.0);	
	final public Input<String> tagInput = new Input<>("tag", "name of the entry in log files containing infector of information. "
			+ "If not specified, assume all entries must be used");	
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		LogAnalyser trueTrace = new LogAnalyser(truthInput.get().getPath(), 0, true, false);

		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}
		String tag = tagInput.get();
	
		int n = trueTrace.getLabels().size();
		if (tag != null) {
			n = 0;
			for (String str : trueTrace.getLabels()) {
				if (str.startsWith(tag + ".")) {
					n++;
				}
			}
		}
		
		out.print("Sample\t");
		for (int i = 0; i < n; i++) {
			out.print("Covered"+i+"\t");
		}
		out.println();
		
				
		
		for (int i = 0; i < trueTrace.getTrace(0).length; i++) {
			String filename = logFilePrefixInput.get().getPath() + i + ".log";
			LogAnalyser trace = new LogAnalyser(filename, burnInPercentageInput.get(), true, false);
			if (i % 10 == 0) {
				Log.warning.print("|");
			} else {
				Log.warning.print(".");
			}
			
			out.print(i + "\t");
			int offset = 0;
			if (tag != null) {
				offset = trace.indexof(tag + ".1") - 1;
			}
			for (int j = 1; j <= n; j++) {
				// get true source value
				int trueSource = 0;
				trueSource = (n + 1 + (int)(double)trueTrace.getTrace(j)[i]) % n;
				
				// collect info from trace
	        	double [] infectedBy = new double[n+1];
	        	Arrays.fill(infectedBy, 0);
	        	Double [] currenttrace = trace.getTrace(offset + j);
	        	for (Double d : currenttrace) {
	        		if (d < -1 || d >= n) {
	        			infectedBy[n]++;
	        		} else {
	        			infectedBy[(n+1+(int)(double)d) % n]++;
	        		}
	        	}
	        	
	        	// determine 95% coverage
	        	int[] index = new int[n+1];
	        	for (int k = 0; k <= n; k++) {
	        		index[k] = k;
	        	}
	        	HeapSort.sort(infectedBy, index);

	        	{
//System.out.print(i*10+j + "\t");
//int sum = 0;
//for (int k = 0; k <= n; k++) {
//	sum += infectedBy[index[k]];
//	System.out.print(infectedBy[index[k]] + "\t");
//}
//System.out.println(sum);
			}

	        	// is true value in the 95% coverage set?
	        	int threshold = (int)(currenttrace.length * coverageInput.get())/ 100;
	        	int sum = 0;
	        	int k = n;
	        	boolean covered = false;
	        	while (sum < threshold) {
	        		if (index[k] == trueSource) {
	        			covered = true;
	        			break;
	        		}
	        		sum += infectedBy[index[k]];
	        		k--;
	        	}
				out.print((covered?1:0) + "\t");
			}
			out.println();
		}
				
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
		Log.warning("\nDone");
	}
	
	
	public static void main(String[] args) throws Exception {
		new Application(new InfectedByCoverageCalculator(), "InfectedByCoverageCalculator", args);
	}

}
