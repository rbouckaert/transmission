package breath.operator;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import breath.distribution.ColourProvider;

@Description("Operator that moves block parameters of a transmission tree")
public class BlockOperator extends Operator {
	final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
    final public Input<Boolean> keepConstantCountInput = new Input<>("keepconstantcount", "if true, for every deleting there is an insertion to keep total sum of block counts constant", false);
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood", Validate.REQUIRED);

    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private TreeInterface tree;
    private double lowerStart, upperStart;
    private double lowerEnd, upperEnd;
    
    @Override
	public void initAndValidate() {
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	tree = treeInput.get();
    	
    	lowerStart = blockStartFraction.getLower();
    	if (lowerStart < 0) {
    		lowerStart = 0;
    	}
    	if (lowerStart > 1) {
    		throw new IllegalArgumentException("lower bound of block start should be less than 1");
    	}
    	upperStart = blockStartFraction.getUpper();
    	if (upperStart < 0) {
    		upperStart = 0;
    	}
    	if (upperStart > 1) {
    		throw new IllegalArgumentException("upper bound of block start should be less than 1");
    	}
    	if (upperStart < lowerStart) {
    		throw new IllegalArgumentException("upper bound of block start should be higher than lower bound");
    	}

    	lowerEnd = blockEndFraction.getLower();
    	if (lowerEnd < 0) {
    		lowerEnd = 0;
    	}
    	if (lowerEnd > 1) {
    		throw new IllegalArgumentException("lower bound of block end should be less than 1");
    	}
    	upperEnd = blockEndFraction.getUpper();
    	if (upperEnd < 0) {
    		upperEnd = 0;
    	}
    	if (upperEnd > 1) {
    		throw new IllegalArgumentException("upper bound of block end should be less than 1");
    	}
    	if (upperEnd < lowerEnd) {
    		throw new IllegalArgumentException("upper bound of block end should be higher than lower bound");
    	}
    	
    	if (lowerStart > lowerEnd) {
    		throw new IllegalArgumentException("lower bound of block start should be lower than lower bound of block end");
    	}
    	if (upperStart > upperEnd) {
    		throw new IllegalArgumentException("upper bound of block start should be lower than upper bound of block end");
    	}
    	
    }

	@Override
	public double proposal() {
//		if (true)
//			if (Randomizer.nextBoolean()) {
//				int [] i = chooseInfectionToRemove();
//				if (i == null) {
//					return Double.NEGATIVE_INFINITY;
//				}
//				return removeInfection0(i);
//			} else {
//				int k = chooseBlockToInsert();
//				return insertInfection0(k);
//			}

		
		
		if (Randomizer.nextBoolean()) {
			int i = Randomizer.nextInt(blockStartFraction.getDimension());

			int attempts = 0;
			while (blockCount.getValue(i) == -1 && attempts < 100) {
				i = Randomizer.nextInt(blockStartFraction.getDimension());
				attempts++;
			}
				
			// only move start and end fraction but not block count
			switch (blockCount.getValue(i)) {
			case -1:
				// nothing to do since start and end fractions are ignored
				break;
			case 0:
				// make sure start == end fraction after proposal
				double f = lowerStart + Randomizer.nextDouble() * (upperStart - upperEnd);
				blockStartFraction.setValue(i, f);
				blockEndFraction.setValue(i, f);
				break;
			default:
				double blockStart = Randomizer.nextDouble();
				double blockEnd = Randomizer.nextDouble();
				if (blockEnd < blockStart) {
					double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
				}
				blockStartFraction.setValue(i, blockStart);
				blockEndFraction.setValue(i, blockEnd);					
			}
			return 0;
		}
		
		if (keepConstantCountInput.get()) {
			int pre = blockCount.getValue(0);

			int [] i = chooseInfectionToRemove();
			double logHR = 0;
			if (i != null) {
				logHR += removeInfection(i);
			} else {
				return Double.NEGATIVE_INFINITY;
			}
			int k = chooseBlockToInsert();
			logHR += insertInfection(k);
			
			int post = blockCount.getValue(0);
			updateStats(pre, post);

			
			return 0*logHR;
		} else	if (Randomizer.nextBoolean()) {
			int [] i = chooseInfectionToRemove();
			if (i == null) {
				return Double.NEGATIVE_INFINITY;
			}
			return removeInfection(i);
		} else {
			int k = chooseBlockToInsert();
			return insertInfection(k);
		}
		
	}

    int [][] stats = new int[3][3];
	private void updateStats(int pre, int post) {
		if (true) return;
		stats[1+pre][1+post]++;
		int total = stats[0][0] + stats[0][1] + stats[1][0] + stats[1][1] + stats[1][2] + stats[2][1] + stats[2][2];
		if (stats[0][0] > 0 && (total) % 1000000 == 0) {
			StringBuilder b = new StringBuilder();
			double t = stats[0][0] + stats[0][1];
			b.append(stats[0][0]/t + " " + stats[0][1]/t + ";");
			t = stats[1][0] + stats[1][1] + stats[1][2];
			b.append(stats[1][0]/t + " " + stats[1][1]/t + " " + stats[1][2]/t + ";");
			t = stats[2][1] + stats[2][2];
			b.append(stats[2][1]/t + " " + stats[2][2]/t + ";");

			t = total;
			b.append((stats[0][0] + stats[0][1])/t + " " );
			b.append((stats[1][0] + stats[1][1] + stats[1][2])/t+ " " );
			b.append((stats[2][1] + stats[2][2])/t);
			b.append("\n");
			
			System.out.println(b.toString());
		}
	}

	
	private int eligbleInfectionCount = 0;
	
	
	private int [] calcEligbleInfectionCount() {
		int [] colourAtBase = new int[tree.getNodeCount()];
		int n = tree.getLeafNodeCount();
		ColourProvider.getColour(tree.getRoot(), blockCount, n, colourAtBase);
		
		eligbleInfectionCount = 0;
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				if (!(colourAtBase[i] < n && !tree.getNode(i).isRoot() && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					eligbleInfectionCount += 1;
				}
			} else if (blockCount.getValue(i) > 0) {
				eligbleInfectionCount += 2;
//			} else {
//				eligbleInfectionCount += blockCount.getValue(i) + 1;
			}
		}
		return colourAtBase;
	}
	
	private int[] chooseInfectionToRemove() {
		int [] colourAtBase = calcEligbleInfectionCount();
		int n = tree.getLeafNodeCount();
		if (eligbleInfectionCount == 0) {
			return null;
		}
		
		int k = Randomizer.nextInt(eligbleInfectionCount);
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				if (!(colourAtBase[i] < n && !tree.getNode(i).isRoot() && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					k--;
				}
			} else if (blockCount.getValue(i) > 0) {
				k -= 2;
			}
			if (k < 0) {
				return new int[] {i, k};
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}
	
	private int chooseBlockToInsert() {
//		{
//			int i = Randomizer.nextInt(tree.getNodeCount()-1);
//			if (true) return i;
//		}
		
		
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		double r = Randomizer.nextDouble() * length;
		int i = 0;
		while (r > 0) {
			Node node = tree.getNode(i);
			if (r < node.getLength()) {
				return i;
			}
			r = r - node.getLength();
			i++;
		}
		throw new RuntimeException("Programmer error: should not get here");
	}

	private double insertInfection(int i) {
		switch (blockCount.getValue(i)) {
		case -1:
			blockCount.setValue(i, 0);
			double f = Randomizer.nextDouble();
			blockStartFraction.setValue(i, f);
			blockEndFraction.setValue(i, f);
			break;
			
		case 0:
			// add infection
			blockCount.setValue(i, 1);
			
			double blockStart = Randomizer.nextDouble();
			double blockEnd = Randomizer.nextDouble();
			if (blockEnd < blockStart) {
				double tmp = blockEnd; blockEnd = blockStart; blockStart = tmp;
			}
			blockStartFraction.setValue(i, blockStart);
			blockEndFraction.setValue(i, blockEnd);					
			break;
		
		default:
			// add infection
			blockCount.setValue(i, blockCount.getValue(i)+1);
			
		}
	
		calcEligbleInfectionCount();
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		return Math.log(1.0/eligbleInfectionCount)
			   - Math.log(tree.getNode(i).getLength() / length)
				;
		
	} // insertInfection

	private double removeInfection(int [] infection) {
		int i = infection[0];
		switch (blockCount.getValue(i)) {
		case -1:
			// do nothing, should not get here
			return 0; 
			
		case 0:
			// remove infection
			blockCount.setValue(i, -1);
			break;
			
		case 1:
			// remove infection
			blockCount.setValue(i, 0);
			if (Randomizer.nextBoolean()) {
				blockStartFraction.setValue(i, blockEndFraction.getValue(i));
			} else {
				blockEndFraction.setValue(i, blockStartFraction.getValue(i));
			}
			break;
		default:
			// remove infection
			blockCount.setValue(i, blockCount.getValue(i)-1);
		}
		
		
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}
		return Math.log(tree.getNode(i).getLength() / length) 
				- Math.log(1.0/eligbleInfectionCount) 
				;
	} // removeInfection


	private double insertInfection0(int i) {
			// add infection
			blockCount.setValue(i, blockCount.getValue(i)+1);
			double k = blockCount.getValue(i) + 2;
			blockStartFraction.setValue(i, 1.0/k);
			blockEndFraction.setValue(i, (k-1.0)/k);
			return 0;
	} // insertInfection

	private double removeInfection0(int [] infection) {
		int i = infection[0];

		int bc = blockCount.getValue(i);
		if (bc == -1) {
			throw new RuntimeException("Programmer error");
		}
		// remove infection
		blockCount.setValue(i, blockCount.getValue(i)-1);

		double k = blockCount.getValue(i) + 2;
		blockStartFraction.setValue(i, 1.0/k);
		blockEndFraction.setValue(i, (k-1.0)/k);
		return 0;
	} // removeInfection

	@Override
	public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(blockCount);
        list.add(blockStartFraction);
        list.add(blockEndFraction);
		return list;
	}
}
