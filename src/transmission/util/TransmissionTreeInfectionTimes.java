package transmission.util;


import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Function.Constant;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;

@Description("Convert transmission tree log into bar plot with transmission times")
public class TransmissionTreeInfectionTimes extends Runnable {	
	final public Input<TreeFile> srcInput = new Input<>("in", "source tree (set) file with transmission tree annotations");
	final public Input<String> partitionInput = new Input<>("partition", "name of the partition appended to `blockcount, blockend and blockstart`");
	final public Input<Function> endTimeInput = new Input<>("endTime", "end time of the study", new Constant("1.0"));
	final public Input<OutFile> pngFileInput = new Input<>("png", "name of file to write bar-chart plot", new OutFile("[[none]]"));	

	
	@Override
	public void initAndValidate() {
	}

	double endTime;
	@Override
	public void run() throws Exception {
		endTime = endTimeInput.get().getArrayValue();
		double [] bins = new double[100];
		

        // read trees one by one, adjust tip heights and write out relabeled tree in newick format
        MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(srcInput.get().getAbsolutePath(), 0);
        trees.reset();
    	Tree tree = trees.next();
    	int n = tree.getLeafNodeCount();
    	
    	
    	
        trees.reset();
        while (trees.hasNext()) {
        	tree = trees.next();
        	
        	// extract meta data from tree
        	Double [] start = new Double[n*2-1];
        	Double [] end = new Double[n*2-1];
        	Integer [] count = new Integer[n*2-1];
        	for (int i = 0; i < tree.getNodeCount(); i++) {
        		Node node = tree.getNode(i);
        		Object o = node.getMetaData("start");
        		if (o == null) {
        			o = node.getMetaData("blockstart");
        		}
        		if (o == null) {
        			o = node.getMetaData("blockstart.t:" + partitionInput.get());
        		}
        		start[i] = (Double) o;

        		o = node.getMetaData("end");
        		if (o == null) {
        			o = node.getMetaData("blockend");
        		}
        		if (o == null) {
        			o = node.getMetaData("blockend.t:" + partitionInput.get());
        		}
        		end[i] = (Double) o;
        		        		
        		o = node.getMetaData("blockcount");
        		if (o == null) {
        			o = node.getMetaData("blockcount.t:" + partitionInput.get());
        		}
        		count[i] = o == null ? 0 : (int)(double) o;
        	}
        	
        	for (int i = 0; i < n*2-1; i++) {
        		int c = count[i];
        		Node node = tree.getNode(i);
        		if (c == 0) {
        			add(node.getHeight() + node.getLength() * start[i], bins);
        		} else if (c > 0) {
        			for (int j = 0; j <= count[i]; j++) {
            			add(node.getHeight() + node.getLength() * (start[i] + j * (end[i]-start[i])/count[i]), bins);
        			}	
        		}
        	}
        }
        
        double sum = 0;
        for (double d : bins) {
        	sum += d;
        }
        for (int i = 0; i < bins.length; i++) {
        	bins[i] /= sum;
        }
        
        System.out.println(Arrays.toString(bins));
        
        
        Log.warning("Done");
	}


    private void add(double d, double[] bins) {
		int i = (int)(bins.length * d / endTime);
		bins[i]++;
	}

	public static void main(String[] args) throws Exception {
		new Application(new TransmissionTreeInfectionTimes(), "Transmission Tree Infection Times", args);
	}
}
