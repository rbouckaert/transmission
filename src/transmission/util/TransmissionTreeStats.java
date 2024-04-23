package transmission.util;

import java.io.PrintStream;
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

	final public Input<String> outputDirInput = new Input<>("out", "directory where to put files with tranmsision & sampling time stats", "/tmp");
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		MemoryFriendlyTreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.get().getAbsolutePath(),
				burnInPercentageInput.get());
		trees.reset();
		Tree tree = trees.next();
		int leafNodeCount = tree.getLeafNodeCount();
		
		PrintStream outTimeTillTransmission = new PrintStream(outputDirInput.get()+ "/timeTillTransmission.dat");
		outTimeTillTransmission.print("Sample\t");
		for (int i = 0; i < leafNodeCount; i++) {
			outTimeTillTransmission.print(tree.getNode(i).getID() + "\t");
		}
		outTimeTillTransmission.print("\n");
		PrintStream outTimeTillSampling = new PrintStream(outputDirInput.get()+ "/timeTillSampling.dat");
		outTimeTillSampling.print("Sample\t");
		for (int i = 0; i < leafNodeCount; i++) {
			outTimeTillSampling.print(tree.getNode(i).getID() + "\t");
		}
		outTimeTillSampling.print("\n");

		trees.reset();

		int sampleCount = 0;
		Integer[] count = new Integer[leafNodeCount * 2 - 1];
		Double[] blockStart = new Double[leafNodeCount * 2 - 1];
		Double[] blockEnd = new Double[leafNodeCount * 2 - 1];
		double[] infectionTimeHost = new double[leafNodeCount * 2 - 1];
		double[] firstInfecteeTimeByHost = new double[leafNodeCount * 2 - 1];
		
		double leafBranchLength = 0;
		double internalBranchLength = 0;
		double leafTransmissionCount = 0;
		double internalTransmissionCount = 0;
		double maxBlockCount = 0;
		
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
				
				o = node.getMetaData("blockstart");
				if (o == null) {
					o = node.getMetaData("blockstart.t:" + partitionInput.get());
				}
				blockStart[i] = o == null ? 1.0 : (double) o;

				o = node.getMetaData("blockend");
				if (o == null) {
					o = node.getMetaData("blockend.t:" + partitionInput.get());
				}
				blockEnd[i] = o == null ? 1.0 : (double) o;
			}
			// root count = -1
			count[count.length-1] = -1;
			IntegerParameter blockCount = new IntegerParameter(count);

			// calculate colouring
			int[] colourAtBase = new int[leafNodeCount * 2 - 1];
			ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);

			// determine who infected who
			int[] infectedBy = new int[tree.getLeafNodeCount()];
			Arrays.fill(infectedBy, -1);
			
			for (int i = 0; i < 2 * leafNodeCount - 2; i++) {
				beast.base.evolution.tree.Node node = tree.getNode(i);
				beast.base.evolution.tree.Node parent = node.getParent();
				if (colourAtBase[node.getNr()] < leafNodeCount && colourAtBase[parent.getNr()] < leafNodeCount
						&& colourAtBase[node.getNr()] != colourAtBase[parent.getNr()]) {
					if (blockCount.getValue(node.getNr()) == 0) {
						infectedBy[colourAtBase[node.getNr()]] = colourAtBase[parent.getNr()];
					}
				}
			}

			// calc stats
			Node[] nodes = tree.getNodesAsArray();
			outTimeTillSampling.print(sampleCount + "\t");
			for (int i = 0; i < leafNodeCount; i++) {
				leafBranchLength += nodes[i].getLength();
				leafTransmissionCount += blockCount.getValue(i) + 1;
				
				Node node = nodes[i];
				double timeToSampling = 0;
				while (!node.isRoot() && blockCount.getValue(node.getNr()) < 0) {
					timeToSampling += node.getLength();
					node = node.getParent();
				}
				timeToSampling += node.getLength() * blockStart[node.getNr()];
				infectionTimeHost[i] = node.getHeight() + node.getLength() * blockStart[node.getNr()];
				outTimeTillSampling.print(timeToSampling + "\t");
			}
			outTimeTillSampling.println();
			
			
			// determine firstInfecteeTimeByHost
			Arrays.fill(firstInfecteeTimeByHost, -1.0);
			for (int i = 0; i < 2 * leafNodeCount - 2; i++) {
				Node node = nodes[i];
				int host = colourAtBase[node.getParent().getNr()];
				if (blockCount.getValue(i) >= 0 && host < leafNodeCount) {
					double infectionTime = node.getHeight() + node.getLength() * blockEnd[i];
					if (infectionTime > firstInfecteeTimeByHost[host]) {
						firstInfecteeTimeByHost[host] = infectionTime;
					}
				}
			}
			
			outTimeTillTransmission.print(sampleCount + "\t");
			for (int i = 0; i < leafNodeCount; i++) {
				if (firstInfecteeTimeByHost[i] > 0) {
					outTimeTillTransmission.print((infectionTimeHost[i] - firstInfecteeTimeByHost[i]) + "\t");
				} else {
					outTimeTillTransmission.print("-1\t");
				}
			}
			outTimeTillTransmission.println();

			
			int maxBlockCountInTree = 0;
			for (int i = leafNodeCount; i < leafNodeCount * 2 - 1; i++) {
				internalBranchLength += nodes[i].getLength();
				internalTransmissionCount += blockCount.getValue(i) + 1;
				maxBlockCountInTree = Math.max(maxBlockCountInTree,  blockCount.getValue(i) + 1);
			}
			maxBlockCount += maxBlockCountInTree;

			sampleCount++;
		}
		Log.info("average leaf Branch Length = " + f.format(leafBranchLength / (sampleCount * leafNodeCount)));
		Log.info("average internal Branch Length = " + f.format(internalBranchLength / (sampleCount * (leafNodeCount-2))));
		Log.info("average leaf transmission count = " + f.format(leafTransmissionCount / (sampleCount * leafNodeCount)));
		Log.info("average internal transmission count = " + f.format(internalTransmissionCount / (sampleCount * (leafNodeCount-2))));
		Log.info("average transmission count per tree = " + f.format((leafTransmissionCount+internalTransmissionCount) / sampleCount));

		Log.info("average unsampled hosts per tree = " + f.format((leafTransmissionCount+internalTransmissionCount) / sampleCount - (leafNodeCount-1)));
		Log.info("average maximum block count = " + f.format(maxBlockCount / sampleCount));
		
		outTimeTillTransmission.close();
		outTimeTillSampling.close();
	}

	public static void main(String[] args) throws Exception {
		new Application(new TransmissionTreeStats(), "TransmissionTreeStats", args);
	}

}
