package transmission.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;

import org.graphper.api.Graphviz;
import org.graphper.api.Graphviz.GraphvizBuilder;
import org.graphper.api.Line;
import org.graphper.api.Node;

import beastfx.app.tools.Application;
import beastfx.app.util.LogFile;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import transmission.distribution.ColourProvider;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;

@Description("Create SVG files to visualise who infected who")
public class WIWVisualiser extends beast.base.inference.Runnable {
	final public Input<TreeFile> treeFile = new Input<>("trees", "tree file file with transmission trees.", new TreeFile("[[none]]"));
	final public Input<LogFile> inFile = new Input<>("log", "trace file containing infectorOf log. Ignored if tree file is specified", new LogFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("/tmp/wiw.svg"));
	final public Input<String> prefixInput = new Input<>("prefix", "prefix of infectorOf entry, e.g., infectorOf", "infectorOf");
	final public Input<Double> thresholdInput = new Input<>("threshold", "probability threshold below which edges will be ignored.", 0.1);
	final public Input<String> partitionInput = new Input<>("partition", "name of the partition appended to `blockcount, blockend and blockstart`");

	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		String [] nodeLabels;
		int n;
		
		double [][] transitions;
		if (treeFile.get() != null && !treeFile.get().getName().equals("[[none]]")) {
	        MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.get().getAbsolutePath(), burnInPercentageInput.get());
	        trees.reset();
	    	Tree tree = trees.next();
	    	n = tree.getLeafNodeCount();
			nodeLabels = new String[n];
			for (int i = 0; i < n; i++) {
				nodeLabels[i] = tree.getNode(i).getID();
			}
			
	        trees.reset();
			transitions = new double[n][n+1];
			int sampleCount = 0;
	        while (trees.hasNext()) {
	        	tree = trees.next();
	        	
	        	// extract meta data from tree
	        	Integer [] count = new Integer[n*2-1];
	        	for (int i = 0; i < tree.getNodeCount(); i++) {
	        		beast.base.evolution.tree.Node node = tree.getNode(i);
	        		Object o = node.getMetaData("blockcount");
	        		if (o == null) {
	        			o = node.getMetaData("blockcount.t:" + partitionInput.get());
	        		}
	        		count[i] = o == null ? 0 : (int)(double) o;
	        	}
	            IntegerParameter blockCount = new IntegerParameter(count);

	            // calculate colouring
	        	int [] colourAtBase = new int[n*2-1];
	        	ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);
	        	
	        	// determine who infected who
	        	int [] infectedBy = new int[tree.getLeafNodeCount()];
	        	Arrays.fill(infectedBy, -1);
	        	for (int i = 0; i < 2 * n - 2; i++) {
	        		beast.base.evolution.tree.Node node = tree.getNode(i);
	        		beast.base.evolution.tree.Node parent = node.getParent();
	        		if (colourAtBase[node.getNr()] < n && colourAtBase[parent.getNr()] < n && 
	        				colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
	        			if (blockCount.getValue(node.getNr()) == 0) {
	        				infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
	        			}
	        		}
	        	}
	        	
	        	
	        	// log the result
	        	for (int i = 0; i < n; i++) {
		        	transitions[i][(infectedBy[i] + n + 1) % (n+1)]++;
	        	}
	        	sampleCount++;
	        }			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n+1; j++) {
					transitions[i][j] /= sampleCount;
				}
			}
		} else {
			LogAnalyser trace = new LogAnalyser(inFile.get().getPath(), burnInPercentageInput.get(), true, false);
			
			// count number of taxa from trace file
			List<String> labels = trace.getLabels();
			n = 0;
			String prefix = prefixInput.get() + ".1";
			int offset = 0;
			for (int i = 0; i < labels.size(); i++) {
				String label = labels.get(i);
				if (label.equals(prefix)) {
					offset = i;
					prefix = prefixInput.get();
					while (labels.get(i+n).equals(prefix + "." + (n+1))) {
						n++;
					}
					break;
				}
			}
			
			// count transitions
			transitions = new double[n][n+1];
			int sampleCount = 0;
			for (int i = 0; i < n; i++) {
				Double [] infectorOf = trace.getTrace(offset + i + 1);
				sampleCount = infectorOf.length;
				for (Double d : infectorOf) {
					transitions[i][ (int)(d + n + 1) % (n+1)]++;
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n+1; j++) {
					transitions[i][j] /= sampleCount;
				}
			}
			nodeLabels = new String[n];
			for (int i = 0; i < n; i++) {
				nodeLabels[i] = (i+1) + "";
			}
		}		
		

		
		// build transition graph
		GraphvizBuilder dotty = Graphviz.digraph();
		//dotty = dotty.tempNode(Node.builder().shape(NodeShapeEnum.RECT).build());
		
		// add nodes
		DecimalFormat f = new DecimalFormat("#.##");
		Node [] nodes = new Node[n];

		boolean [] nodesInUse = new boolean[n];
		double threshold = thresholdInput.get();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i!=j && transitions[i][j] >= threshold) {
					nodesInUse[i] = true;
					nodesInUse[j] = true;
				}
			}
		}
		
		for (int i = 0; i < n; i++) {
			if (nodesInUse[i]) {
				Node node = Node.builder().label(nodeLabels[i] + " (" + f.format(1.0-transitions[i][n]-transitions[i][i]) + ")").build();
	//			System.err.println(transitions[i][i]);
				dotty = dotty.addNode(node);
				nodes[i] = node;
			}
		}
		
		// add edges
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i!=j && transitions[i][j] >= threshold) {
					dotty = dotty.addLine(Line.builder(nodes[j], nodes[i]).label(f.format(transitions[i][j])).build());
				}
			}
		}
		
		// create SVG
		String svg = dotty.build().toSvgStr();

		// output graph
		try {
			File svgFile = outputInput.get();
			Log.warning("Writing to file " + svgFile.getPath());
			FileWriter outfile = new FileWriter(svgFile);
			outfile.write(svg);
			outfile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("Done");	
	}
	
	
	public static void main(String[] args) throws Exception {
		new Application(new WIWVisualiser(), "Who-Infected-Who Visualiser", args);
	}
}
