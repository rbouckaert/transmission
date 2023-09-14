package transmission.logger;

import java.io.PrintStream;
import java.util.Arrays;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import transmission.distribution.ColourProvider;

@Description("Logs transmission tree with binary and single child nodes annotated with colour")
public class ColouredTreeLogger extends BEASTObject implements Loggable {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");
    final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);

    private Tree tree;
    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private int [] colourAtBase;

    @Override
	public void initAndValidate() {
    	tree = (Tree) treeInput.get();

    	int n = tree.getNodeCount();
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	colourAtBase = new int[n];
	}

	@Override
	public void init(PrintStream out) {
		tree.init(out);
	}

	@Override
	public void log(long sample, PrintStream out) {
		ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);
		
        out.print("tree STATE_" + sample + " = ");
        final String newick = toSortedNewick(tree.getRoot(), new int[1]);
        out.print(newick);
        out.print(";");
	}

	String toSortedNewick(Node node, int [] maxNodeInClade) {
		StringBuilder buf = new StringBuilder();
		if (!node.isLeaf()) {

            // General method for >2 children

            String[] childStrings = new String[node.getChildCount()];
            int[] maxNodeNrs = new int[node.getChildCount()];
            Integer[] indices = new Integer[node.getChildCount()];
            for (int i = 0; i < node.getChildCount(); i++) {
                childStrings[i] = toSortedNewick(node.getChild(i), maxNodeInClade);
                maxNodeNrs[i] = maxNodeInClade[0];
                indices[i] = i;
            }

            Arrays.sort(indices, (i1, i2) -> {
                if (maxNodeNrs[i1] < maxNodeNrs[i2])
                    return -1;

                if (maxNodeNrs[i1] > maxNodeNrs[i2])
                    return 1;

                return 0;
            });

            maxNodeInClade[0] = maxNodeNrs[maxNodeNrs.length - 1];

            buf.append("(");
            for (int i = 0; i < indices.length; i++) {
                if (i > 0) {
                    buf.append(",");
                }
                buf.append(childStrings[indices[i]]);
            }

            buf.append(")");

            if (getID() != null) {
                buf.append(node.getNr() + 1);
            }

    	} else {
    		maxNodeInClade[0] = node.getNr();
    		buf.append(node.getNr() + 1);
    	}

    	int i = node.getNr();
    	switch (blockCount.getValue(i)) {
    		case -1:
    	    	buf.append("[&colour=" + colourAtBase[i] + "]:");
    			buf.append(node.getLength());
    	    	return buf.toString();
    		case 0:    			
    	    	buf.append("[&colour=" + colourAtBase[i] + "]:");
    			buf.append(node.getLength() * blockStartFraction.getValue(i));
    			buf.append(")");
    	    	buf.append("[&colour=" + colourAtBase[node.getParent().getNr()] + "]:");
    			buf.append(node.getLength() * (1.0-blockEndFraction.getValue(i)));
    			return "(" + buf.toString();
    		default:
    	    	buf.append("[&colour=" + colourAtBase[i] + "]:");
    			buf.append(node.getLength() * blockStartFraction.getValue(i));
    			buf.append(")");
    	    	buf.append("[&colour=-1,count=" + blockCount.getValue(i)+ "]:");
    			buf.append(node.getLength() * (blockEndFraction.getValue(i) - blockStartFraction.getValue(i)));
    			buf.append(")");
    	    	buf.append("[&colour=" + colourAtBase[node.getParent().getNr()] + "]:");
    			buf.append(node.getLength() * (1.0-blockEndFraction.getValue(i)));
    			return "((" + buf.toString();
    	}
	}

	
	@Override
	public void close(PrintStream out) {
		tree.close(out);
	}
	
	@Override
	public String toString() {
		ColourProvider.getColour(tree.getRoot(), blockCount, tree.getLeafNodeCount(), colourAtBase);
        String newick = toSortedNewick(tree.getRoot(), new int[1]);
		return newick;
	}

}
