package transmission.util;

import java.io.PrintStream;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import transmission.distribution.ColourProvider;

@Description("Convert transmission tree log into trace log with who infected who. "
		+ "If leaf i is infected by leaf j, the log contains j at position i. "
		+ "If leaf is infected by an unsampled node, the log contains -1.")
public class TranmissionTree2InfectedByLog extends Runnable {	
	final public Input<TreeFile> srcInput = new Input<>("in", "source tree (set) file with transmission tree annotations");
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified", new OutFile("[[none]]"));

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}

        // read trees one by one, adjust tip heights and write out relabeled tree in newick format
        MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(srcInput.get().getAbsolutePath(), 0);
        trees.reset();
    	Tree tree = trees.next();
    	
    	
    	out.print("Sample\t");
    	int n = tree.getLeafNodeCount();
    	for (int i = 0; i < n; i++) {
    		out.print(tree.getNode(i).getID() + "\t");
    	}
    	out.println();
    	
        trees.reset();
        int k = 0;
        while (trees.hasNext()) {
        	tree = trees.next();
        	
        	// extract meta data from tree
        	Double [] start = new Double[n*2-1];
        	Double [] end = new Double[n*2-1];
        	Integer [] count = new Integer[n*2-1];
        	for (int i = 0; i < tree.getNodeCount(); i++) {
        		Node node = tree.getNode(i);
        		start[i] = (Double) node.getMetaData("start");
        		end[i] = (Double) node.getMetaData("end");
        		Object o = node.getMetaData("blockcount");
        		count[i] = o == null ? 0 : (int)(double) o;
        	}
            RealParameter blockStartFraction = new RealParameter(start);
            RealParameter blockEndFraction = new RealParameter(end);
            IntegerParameter blockCount = new IntegerParameter(count);

            // calculate colouring
        	int [] colourAtBase = new int[n*2-1];
        	ColourProvider.getColour(tree, blockStartFraction, blockEndFraction, blockCount, colourAtBase);
        	
        	// determine who infected who
        	int [] infectedBy = new int[tree.getLeafNodeCount()];
        	Arrays.fill(infectedBy, -1);
        	for (int i = 0; i < 2 * n - 2; i++) {
        		Node node = tree.getNode(i);
        		Node parent = node.getParent();
        		if (colourAtBase[node.getNr()] < n && colourAtBase[parent.getNr()] < n && 
        				colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
        			infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
        		}
        	}
        	
        	// log the result
        	out.print(k + "\t");
        	for (int i = 0; i < n; i++) {
        		out.print(infectedBy[i] + "\t");
        	}
        	out.println();

        	k++;
        }
		
		
		
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
        Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new TranmissionTree2InfectedByLog(), "TranmissionTree2InfectedByLog", args);
	}
}
