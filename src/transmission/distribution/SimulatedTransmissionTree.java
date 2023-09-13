package transmission.distribution;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Function.Constant;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.Runnable;
import beast.base.util.HeapSort;
import beast.base.util.Randomizer;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;

@Description("Simulates transmission tree with colouring and block counts")
public class SimulatedTransmissionTree extends Runnable {
	final public Input<Function> endTimeInput = new Input<>("endTime", "end time of the study", new Constant("1.0"));
	final public Input<Function> popSizeInput = new Input<>("popSize",
			"population size governing the coalescent process", new Constant("0.1"));

	final public Input<Function> sampleShapeInput = new Input<>("sampleShape",
			"shape parameter of the sampling intensity function", new Constant("2.0"));
	final public Input<Function> sampleRateInput = new Input<>("sampleRate",
			"rate parameter of the sampling intensity function", new Constant("5.0"));
	final public Input<Function> sampleConstantInput = new Input<>("sampleConstant",
			"constant multiplier of the sampling intensity function", new Constant("0.75"));

	final public Input<Function> transmissionShapeInput = new Input<>("transmissionShape",
			"shape parameter of the transmission intensity function", new Constant("2.5"));
	final public Input<Function> transmissionRateInput = new Input<>("transmissionRate",
			"rate parameter of the transmission intensity function", new Constant("10.0"));
	final public Input<Function> transmissionConstantInput = new Input<>("transmissionConstant",
			"constant multiplier of the transmission intensity function", new Constant("1.5"));

	public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified");
	public Input<Long> seedInput = new Input<>("seed","random number seed used to initialise the random number generator");
	final public Input<Integer> maxAttemptsInput = new Input<>("maxAttempts",
			"maximum number of attempts to generate coalescent sub-trees", 1000);
	public Input<Integer> taxonCountInput = new Input<>("taxonCount", "generate tree with taxonCount number of taxa. Ignored if negative", -1);
	public Input<Integer> treeCountInput = new Input<>("treeCount", "generate treeCount number of trees", 1);

	private Node root;
	private Map<Node, Integer> colourMap;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}
		
		int taxonCount = taxonCountInput.get();
		
		for (int i = 0; i < treeCountInput.get(); i++) {
			do {
				runOnce();
			} while (taxonCount > 0 && taxonCount != root.getAllLeafNodes().size());
			
			// convert to binary tree
			String newick = toNewick(root);
		
			// for debugging
			double h = root.getHeight();
			for (Node node : root.getAllLeafNodes()) {
				h = Math.min(h,  node.getHeight());
			}
			System.err.print( -h+ " ");
			System.err.println(toShortNewick(root, colourMap));
		
			out.println(newick);
		}
		
		if (outputInput.get() != null) {
			out.close();
		}		
		Log.warning("Done");
	}
		
	private void runOnce() throws MathException {	
		double endTime = endTimeInput.get().getArrayValue();
		double popSize = popSizeInput.get().getArrayValue();

		double sampleShape = sampleShapeInput.get().getArrayValue();
		double sampleRate = sampleRateInput.get().getArrayValue();
		double sampleConstant = sampleConstantInput.get().getArrayValue();

		double transmissionShape = transmissionShapeInput.get().getArrayValue();
		double transmissionRate = transmissionRateInput.get().getArrayValue();
		double transmissionConstant = transmissionConstantInput.get().getArrayValue();

		root = new Node();
		root.setHeight(endTime);
		List<Node> nodes = new ArrayList<>();
		nodes.add(root);

		PoissonDistribution poisson = new PoissonDistributionImpl(transmissionConstant);
		GammaDistribution sampleIntensity = new GammaDistributionImpl(sampleShape, 1.0 / sampleRate);
		GammaDistribution transmissionIntensity = new GammaDistributionImpl(transmissionShape, 1.0 / transmissionRate);
		ConstantPopulation popFun = new ConstantPopulation();
		popFun.initByName("popSize", popSize + "");

		List<Node> leafs = new ArrayList<>();
		int colour = 0;
		colourMap = new HashMap<>();
		colourMap.put(root, colour);
		
		while (nodes.size() > 0) {
			// continue with last node
			Node node = nodes.remove(nodes.size() - 1);

			// 1. draw number of events
			int n = poisson.inverseCumulativeProbability(Randomizer.nextDouble());

			// 2. draw whether colour will be sampled
			boolean sample = (Randomizer.nextDouble() < sampleConstant);

			// 3. simulate the time of sampling:
			// Note: do not need multiply by sampleConstant
			double sampletime = !sample ? 0
					: node.getHeight() - sampleIntensity.inverseCumulativeProbability(Randomizer.nextDouble());
			if (sampletime < 0) {
				sample = false;
			}

			// 4. Simulate the times when node infects the new infectees
			double[] times = new double[n];
			for (int i = 0; i < n; i++) {
				// Note: do not need multiply by transmissionConstant
				times[i] = node.getHeight()
						- transmissionIntensity.inverseCumulativeProbability(Randomizer.nextDouble());
			}
			Arrays.sort(times);
			// remove times that are invalid: after study time, or after sample time (if any)
			while (n > 0 && (times[n - 1] < 0 || times[n - 1] < sampletime)) {
				n--;
			}
			
			List<Node> current = new ArrayList<>();
			// create leaf node
			if (sample) {
				Node leaf = new Node();
				colourMap.put(leaf, colour);
				leaf.setHeight(sampletime);
				leafs.add(leaf);
				current.add(leaf);
				leaf.setID("t" + leafs.size());
			}
			// create internal (infection) nodes
			for (int i = 0; i < n; i++) {
				Node infectee = new Node();
				colourMap.put(infectee, colour);
				infectee.setHeight(times[i]);
				nodes.add(infectee);
				current.add(infectee);
			}

			
			// 5. Sampling a within-host phylogeny
			double currentHeight = sample ? sampletime : (n > 0 ? times[0] : 0);
		
			Node fragment = simulateCoalescent(current, popFun, currentHeight, node.getHeight(), maxAttemptsInput.get());
			// connect to node
			node.addChild(fragment);

			colourFragment(fragment, colour, colourMap);
			
			colour++;
		}
		
		
		// find all nodes to include in tree
		Set<Node> includedNodes = new HashSet<>();
		for (Node node : leafs) {
			while (node != root) {
				includedNodes.add(node);
				node = node.getParent();
			}
		}
		// remove nodes not included from tree
		traverse(root, includedNodes);
	}

	private Node simulateCoalescent(List<Node> current, ConstantPopulation popFun, double currentHeight,
			double height, int maxAttemptCount) {
		if (current.size() == 0) {
			Node node =  new Node();
			node.setHeight(currentHeight);
			return node;
		}
		
		List<Node> fragment;
		int attempt = 0;
		do {
			List<Node> currentCopy = new ArrayList<>(current);
			fragment = simulateCoalescent(currentCopy, popFun, currentHeight, height);
			if (fragment.size() == 1) {
				return fragment.get(0);
			}
			attempt++;
		} while (attempt < maxAttemptCount);
		throw new RuntimeException("Could not find a proper coalescent tree after " + maxAttemptCount + " attempts. "
				+ "Consider decreasing the population size or increasing maxAttempts.");
	}

	private void colourFragment(Node node, int colour, Map<Node, Integer> colourMap) {
		colourMap.put(node,  colour);
		for (Node child : node.getChildren()) {
			colourFragment(child, colour, colourMap);
		}
	}

	private String toShortNewick(Node node, Map<Node, Integer> colourMap) {
		if (node.isLeaf()) {
			return node.getID();
		} else {
			String newick = "(";
			for (int i = 0; i < node.getChildCount(); i++) {
				Node child = node.getChild(i);
				newick += toShortNewick(child, colourMap);
				newick += "[" + colourMap.get(child)+"]";
				if (i < node.getChildCount() - 1) {
					newick += ",";
				}
			}
			return newick + ")";
		}
	}
	
	private String toNewick(Node node) {
if (node.getID()!=null && node.getID().equals("t8")) {
	int h = 4;
	h--;
}
		switch(node.getChildCount()) {
		case 0: {// leaf
			double length = node.getLength();
			Node p = node.getParent();
			double blockStart = length;
			double blockEnd = length;
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				blockEnd = length;
				length += p.getLength();
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			return node.getID() + "[&blockcount=" + blockCount + (blockCount >= 0 ? ",blockstart=" + (blockStart/length) + ",blockend=" + (blockEnd/length): "") +",color=" + colourMap.get(node) + "]:" + length;
		}
		case 1:
			return toNewick(node.getChild(0));
		case 2:
			Node left = node.getLeft();
			String leftNewick = toNewick(left);
			Node right = node.getRight();
			String rightNewick = toNewick(right);

			double length = node.getLength();
			Node p = node.getParent();
			double blockStart = length;
			double blockEnd = length;
			int blockCount = -1;
			while (p != null && p.getChildCount() == 1) {
				blockEnd = length;
				length += p.getLength();
				p = p.getParent();
				if (colourMap.get(node) != colourMap.get(p)) {
					blockCount++;
				}
			}
			if (p == null) {
				return "(" + leftNewick + "," + rightNewick + ")";
			}
			return "(" + leftNewick + "," + rightNewick + ")" + "[&blockcount=" + blockCount + (blockCount >= 0 ? ",blockstart=" + (blockStart/length) + ",blockend=" + (blockEnd/length): "") +",color=" + colourMap.get(node) + "]:" + length; 
		}
		return null;
	}

	private void traverse(Node node, Set<Node> inclodeNodes) {
		List<Node> children = node.getChildren();
		for (int i = children.size()-1; i >= 0; i--) {
			Node child = children.get(i);
			if (!inclodeNodes.contains(child)) {
				node.removeChild(child);
			}
		}
		for (int i = node.getChildren().size()-1; i >= 0; i--) {
			Node child = children.get(i);
			traverse(child, inclodeNodes);
		}		
	}

	public List<Node> simulateCoalescent(final List<Node> nodes, final PopulationFunction demographic,
			double currentHeight, final double maxHeight) {
		// If only one node, return it
		// continuing results in an infinite loop
		if (nodes.size() == 1)
			return nodes;

		final double[] heights = new double[nodes.size()];
		for (int i = 0; i < nodes.size(); i++) {
			heights[i] = nodes.get(i).getHeight();
		}
		final int[] indices = new int[nodes.size()];
		HeapSort.sort(heights, indices);

		// node list
		List<Node> nodeList = new ArrayList<>();
		int activeNodeCount = 0;
		for (int i = 0; i < nodes.size(); i++) {
			nodeList.add(nodes.get(indices[i]));
		}
        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
            activeNodeCount += 1;
        }

		// get at least two tips
		while (activeNodeCount < 2) {
			currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
	        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
	            activeNodeCount += 1;
	        }
		}

		// simulate coalescent events
		double nextCoalescentHeight = currentHeight
				+ PopulationFunction.Utils.getSimulatedInterval(demographic, activeNodeCount, currentHeight);

		// while (nextCoalescentHeight < maxHeight && (getNodeCount() > 1)) {
		while (nextCoalescentHeight < maxHeight && (nodeList.size() > 1)) {

			if (nextCoalescentHeight >= getMinimumInactiveHeight(activeNodeCount, nodeList)) {
				currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
		        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
		            activeNodeCount += 1;
		        }
			} else {
				currentHeight = coalesceTwoActiveNodes(currentHeight, nextCoalescentHeight, nodeList, activeNodeCount);
				activeNodeCount--;
			}

			// if (getNodeCount() > 1) {
			if (nodeList.size() > 1) {
				// get at least two tips
				while (activeNodeCount < 2) {
					currentHeight = getMinimumInactiveHeight(activeNodeCount, nodeList);
			        while (getMinimumInactiveHeight(activeNodeCount, nodeList) <= currentHeight) {
			            activeNodeCount += 1;
			        }
				}

				// nextCoalescentHeight = currentHeight +
				// DemographicFunction.Utils.getMedianInterval(demographic,
				// getActiveNodeCount(), currentHeight);
				nextCoalescentHeight = currentHeight + PopulationFunction.Utils.getSimulatedInterval(demographic,
						activeNodeCount, currentHeight);
			}
		}

		return nodeList;
	}
	
    /**
     * Coalesce two nodes in the active list. This method removes the two
     * (randomly selected) active nodes and replaces them with the new node at
     * the top of the active list.
     * @param minHeight
     * @param height
     * @return
     */
    private double coalesceTwoActiveNodes(final double minHeight, double height, List<Node> nodeList, int activeNodeCount){
        final int node1 = Randomizer.nextInt(activeNodeCount);
        int node2 = node1;
        while (node2 == node1) {
            node2 = Randomizer.nextInt(activeNodeCount);
        }

        final Node left = nodeList.get(node1);
        final Node right = nodeList.get(node2);

        final Node newNode = new Node();
//		System.err.println(2 * m_taxa.get().getNrTaxa() - nodeList.size());
//        newNode.setNr(nextNodeNr++);   // multiple tries may generate an excess of nodes assert(nextNodeNr <= nrOfTaxa*2-1);
        newNode.setHeight(height);
        newNode.setLeft(left);
        left.setParent(newNode);
        newNode.setRight(right);
        right.setParent(newNode);

        nodeList.remove(left);
        nodeList.remove(right);

        activeNodeCount -= 2;

        nodeList.add(activeNodeCount, newNode);

        activeNodeCount += 1;

        
        if (getMinimumInactiveHeight(activeNodeCount, nodeList) < height) {
            throw new RuntimeException(
                    "This should never happen! Somehow the current active node is older than the next inactive node!\n"
            		+ "One possible solution you can try is to increase the population size of the population model.");
        }
        return height;
    }
    
    private double getMinimumInactiveHeight(int activeNodeCount, List<Node> nodeList) {
        if (activeNodeCount < nodeList.size()) {
            return (nodeList.get(activeNodeCount)).getHeight();
        } else
            return Double.POSITIVE_INFINITY;
    }


	public static void main(String[] args) throws Exception {
		new Application(new SimulatedTransmissionTree(), "SimulatedTransmissionTree", args);
	}

}
