package transmission.operator;

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
import transmission.distribution.ColourProvider;

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
		if (Randomizer.nextBoolean()) {
			int i = Randomizer.nextInt(blockStartFraction.getDimension());
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
				// move one boundary only
				if (Randomizer.nextBoolean()) {
					f = lowerStart + Randomizer.nextDouble() * (Math.min(blockEndFraction.getValue(i), upperStart) - lowerStart);
					blockStartFraction.setValue(i, f);
				} else {
					f = upperEnd - Randomizer.nextDouble() * (upperEnd - Math.max(blockStartFraction.getValue(i), lowerEnd));
					blockEndFraction.setValue(i, 1-f);					
				}
			}
			return 0;
		}
		
		if (keepConstantCountInput.get()) {
			int pre = blockCount.getValue(0);

			int i = chooseInfectionToRemove();
			double logHR = 0;
			if (i != -1) {
				logHR += removeInfection(i);
			}
			i = chooseBlockToInsert();
			logHR += insertInfection(i);
			
			int post = blockCount.getValue(0);
			updateStats(pre, post);

			
			return 0*logHR;
		} else	if (Randomizer.nextBoolean()) {
			int i = chooseInfectionToRemove();
			if (i == -1) {
				return Double.NEGATIVE_INFINITY;
			}
			return removeInfection(i);
		} else {
			int i = chooseBlockToInsert();
			return insertInfection(i);
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
	
	private int chooseInfectionToRemove() {
		int [] colourAtBase = new int[tree.getNodeCount()];
		int n = tree.getLeafNodeCount();
		ColourProvider.getColour(tree.getRoot(), blockCount, n, colourAtBase);
		
		eligbleInfectionCount = 0;
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				if (!(colourAtBase[i] < n && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					eligbleInfectionCount += 1;
				}
			} else {
				eligbleInfectionCount += blockCount.getValue(i) + 1;
			}
		}
		if (eligbleInfectionCount == 0) {
			return -1;
		}
		
		
		int k = Randomizer.nextInt(eligbleInfectionCount);
		for (int i = 0; i < blockCount.getDimension(); i++) {
			if (blockCount.getValue(i) == 0) {
				if (!(colourAtBase[i] < n && colourAtBase[tree.getNode(i).getParent().getNr()] < n)) {
					k--;
				}
			} else {
				k -= blockCount.getValue(i) + 1;
			}
			if (k < 0) {
				return i;
			}
		}
		throw new RuntimeException("Programmer error: should not get here");
	}
	
	private int chooseBlockToInsert() {
//		int i = Randomizer.nextInt(tree.getNodeCount()-1);
//		return i;
		
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
			return 0;
			
		case 0:
			// add infection
			blockCount.setValue(i, 1);

			f = lowerStart + Randomizer.nextDouble() * (upperEnd - lowerStart);
			if (f < blockStartFraction.getValue(i)) {
				blockStartFraction.setValue(i, f);
			} else {
				blockEndFraction.setValue(i, f);				
			}
//			double mid = blockStartFraction.getValue(i);
//			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
//			blockStartFraction.setValue(i, mid - f / 2.0);
//			blockEndFraction.setValue(i, mid + f / 2.0);
			return 0;
			
//		case 1:
//			// add infection
//			blockCount.setValue(i, 2);
//
//			f = Randomizer.nextDouble();
//			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
//			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
//			blockStartFraction.setValue(i, mid - f / 2.0);
//			blockEndFraction.setValue(i, mid + f / 2.0);				
//			return Math.log(2.0);
		
		default:
			// add infection
			blockCount.setValue(i, blockCount.getValue(i)+1);
//			f = Randomizer.nextDouble();
//			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
//			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
		}
		return 0;
		
	}

	private double removeInfection(int i) {
		switch (blockCount.getValue(i)) {
		case -1:
			// do nothing, should not get here
			return 0; 
			
		case 0:
			// remove infection
			blockCount.setValue(i, -1);
			break;
			//return Math.log(2.0);
			
		case 1:
			// remove infection
			blockCount.setValue(i, 0);
			if (Randomizer.nextBoolean()) {
				blockStartFraction.setValue(i, blockEndFraction.getValue(i));
			} else {
				blockEndFraction.setValue(i, blockStartFraction.getValue(i));
				
			}
//			double mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
//			blockStartFraction.setValue(i, mid);
//			blockEndFraction.setValue(i, mid);
//			return -Math.log(2.0);
			break;
		default:
			// remove infection
			blockCount.setValue(i, blockCount.getValue(i)-1);
//			double f = Randomizer.nextDouble();
//			mid = (blockStartFraction.getValue(i) + blockEndFraction.getValue(i))/2.0;
//			if (mid- f / 2.0 < 0 || mid + f / 2.0 > 1.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
		}
		
		
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			length += node.getLength();
		}

		return Math.log(tree.getNode(i).getLength() / length) - Math.log(1.0/eligbleInfectionCount);
	}

	
	@Override
	public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(blockCount);
        list.add(blockStartFraction);
        list.add(blockEndFraction);
		return list;
	}
}
