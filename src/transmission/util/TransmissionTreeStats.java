package transmission.util;

import java.text.DecimalFormat;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beast.base.inference.parameter.IntegerParameter;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.TreeFile;
import transmission.distribution.ColourProvider;

@Description("Provide statistics of set of transmission trees")
public class TransmissionTreeStats extends Runnable {
	final public Input<TreeFile> treeFile = new Input<>("trees", "tree file file with transmission trees.",
			new TreeFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
			"percentage of trees to used as burn-in (and will be ignored). NB default 0", 0);
	final public Input<String> partitionInput = new Input<>("partition",
			"name of the partition appended to `blockcount, blockend and blockstart`");

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.get().getAbsolutePath(),
				burnInPercentageInput.get());
		trees.reset();
		Tree tree = trees.next();
		int n = tree.getLeafNodeCount();

		trees.reset();

		int sampleCount = 0;
		Integer[] count = new Integer[n * 2 - 1];
		double leafBranchLength = 0;
		double internalBranchLength = 0;
		double leafTransmissionCount = 0;
		double internalTransmissionCount = 0;
		DecimalFormat f = new DecimalFormat("#.##");
		while (trees.hasNext()) {
			tree = trees.next();

			// extract meta data from tree
			for (int i = 0; i < tree.getNodeCount(); i++) {
				beast.base.evolution.tree.Node node = tree.getNode(i);
				Object o = node.getMetaData("blockcount");
				if (o == null) {
					o = node.getMetaData("blockcount.t:" + partitionInput.get());
				}
				count[i] = o == null ? 0 : (int) (double) o;
			}
			// root count = -1
			count[count.length-1] = -1;
			IntegerParameter blockCount = new IntegerParameter(count);

			// calculate colouring
			int[] colourAtBase = new int[n * 2 - 1];
			ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);

			// determine who infected who
			int[] infectedBy = new int[tree.getLeafNodeCount()];
			Arrays.fill(infectedBy, -1);
			
			for (int i = 0; i < 2 * n - 2; i++) {
				beast.base.evolution.tree.Node node = tree.getNode(i);
				beast.base.evolution.tree.Node parent = node.getParent();
				if (colourAtBase[node.getNr()] < n && colourAtBase[parent.getNr()] < n
						&& colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
					if (blockCount.getValue(node.getNr()) == 0) {
						infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
					}
				}
			}

			// calc stats
			Node[] nodes = tree.getNodesAsArray();
			for (int i = 0; i < n; i++) {
				leafBranchLength += nodes[i].getLength();
				leafTransmissionCount += blockCount.getValue(i) + 1;
			}
			for (int i = n; i < n * 2 - 1; i++) {
				internalBranchLength += nodes[i].getLength();
				internalTransmissionCount += blockCount.getValue(i) + 1;
			}

			sampleCount++;
		}
		Log.info("average leaf Branch Length = " + f.format(leafBranchLength / (sampleCount * n)));
		Log.info("average internal Branch Length = " + f.format(internalBranchLength / (sampleCount * (n-2))));
		Log.info("average leaf transmission count = " + f.format(leafTransmissionCount / (sampleCount * n)));
		Log.info("average internal transmission count = " + f.format(internalTransmissionCount / (sampleCount * (n-2))));
		Log.info("average transmission count per tree = " + f.format((leafTransmissionCount+internalTransmissionCount) / sampleCount));

	}

	public static void main(String[] args) throws Exception {
		new Application(new TransmissionTreeStats(), "TransmissionTreeStats", args);
	}

}
