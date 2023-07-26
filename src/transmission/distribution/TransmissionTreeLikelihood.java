package transmission.distribution;



import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.IntervalList;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;

@Description("Likelihood of a transmission tree")
public class TransmissionTreeLikelihood extends TreeDistribution {
    final public Input<RealParameter> blockStartFractionInput = new Input<>("blockstart", "start of block in fraction of branch length", Validate.REQUIRED);
    final public Input<RealParameter> blockEndFractionInput = new Input<>("blockend", "end of block in fraction of branch length", Validate.REQUIRED);
    final public Input<IntegerParameter> blockCountInput = new Input<>("blockcount", "number of transitions inside a block", Validate.REQUIRED);
    final public Input<IntegerParameter> colourInput = new Input<>("colour", "colour of the base of the branch", Validate.REQUIRED);
    final public Input<PopulationFunction> popSizeInput = new Input<>("populationModel", "A population size model", Validate.REQUIRED);

    private Tree tree;
    private RealParameter blockStartFraction;
    private RealParameter blockEndFraction;
    private IntegerParameter blockCount;
    private IntegerParameter colourAtBase;
    private PopulationFunction popSizeFunction;
    
    @Override
    public void initAndValidate() {
    	tree = (Tree) treeInput.get();
    	if (tree == null) {
    		tree = treeIntervalsInput.get().treeInput.get();
    	}

    	int n = tree.getNodeCount();
    	blockStartFraction = blockStartFractionInput.get();
    	blockEndFraction = blockEndFractionInput.get();
    	blockCount = blockCountInput.get();
    	colourAtBase = colourInput.get();
    	
    	if (blockStartFraction.getDimension() != n) {
    		blockStartFraction.setDimension(n);
    		Log.warning("Setting dimension of parameter " + blockStartFraction.getID() + " to " + n);
    	}
    	if (blockEndFraction.getDimension() != n) {
    		blockEndFraction.setDimension(n);
    		Log.warning("Setting dimension of parameter " + blockEndFraction.getID() + " to " + n);
    	}
    	if (blockCount.getDimension() != n) {
    		blockCount.setDimension(n);
    		Log.warning("Setting dimension of parameter " + blockCount.getID() + " to " + n);
    	}
    	if (colourAtBase.getDimension() != n) {
    		colourAtBase.setDimension(n);
    		Log.warning("Setting dimension of parameter " + colourAtBase.getID() + " to " + n);
    	}
    	
    	popSizeFunction = popSizeInput.get();
    }
    
    @Override
    public double calculateLogP() {
    	logP = 0;
    	logP = calculateCoalescent();
    	return logP;
    }
    
    
    public double calculateCoalescent() {
		List<IntervalList> segments = collectSegments();
		double logP = 0;
		for (IntervalList intervals : segments) {
			if (intervals != null) {
				logP += calculateCoalescent(intervals, 0.0);
			}
		}
		return logP;
	}

    public List<Double> calculateCoalescents() {
		List<IntervalList> segments = collectSegments();
		List<Double> logP = new ArrayList<>();
		for (IntervalList intervals : segments) {
			if (intervals != null) {
				logP.add(calculateCoalescent(intervals, 0.0));
			}
		}
		return logP;
	}

    class SegmentIntervalList implements IntervalList  {

    	private List<Double> times = new ArrayList<>();
    	private List<IntervalType> events = new ArrayList<>();

    	private int intervalCount = 0;
        /**
         * The widths of the intervals.
         */
        private double[] intervals;
        /**
         * The number of uncoalesced lineages within a particular interval.
         */
        private int[] lineageCounts;


		@Override
		public int getIntervalCount() {
	        return intervalCount;
		}

		@Override
		public int getSampleCount() {
			throw new RuntimeException("Not implemented yet");
		}

		@Override
		public double getInterval(int i) {
	        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
	        return intervals[i];
		}

		@Override
		public int getLineageCount(int i) {
	        if (i >= intervalCount) throw new IllegalArgumentException();
	        return lineageCounts[i];
		}

		@Override
		public int getCoalescentEvents(int i) {
	        if (i >= intervalCount) throw new IllegalArgumentException();
	        if (i < intervalCount - 1) {
	            return lineageCounts[i] - lineageCounts[i + 1];
	        } else {
	            return lineageCounts[i] - 1;
	        }
		}

		@Override
		public IntervalType getIntervalType(int i) {
	        if (i >= intervalCount) throw new IllegalArgumentException();
	        int numEvents = getCoalescentEvents(i);

	        if (numEvents > 0) return IntervalType.COALESCENT;
	        else if (numEvents < 0) return IntervalType.SAMPLE;
	        else return IntervalType.NOTHING;
		}

		@Override
		public double getTotalDuration() {
			throw new RuntimeException("Not implemented yet");
		}

		@Override
		public boolean isBinaryCoalescent() {
			throw new RuntimeException("Not implemented yet");
		}

		@Override
		public boolean isCoalescentOnly() {
			throw new RuntimeException("Not implemented yet");
		}


		public void calculateIntervals() {
			double multifurcationLimit = 0.0;
			int nodeCount = events.size();
			
	        if (intervals == null || intervals.length != nodeCount) {
	            intervals = new double[nodeCount];
	            lineageCounts = new int[nodeCount];
	        }

	        // start is the time of the first tip
	        double start = times.get(0);
	        int numLines = 0;
	        int nodeNo = 0;
	        intervalCount = 0;
	        while (nodeNo < nodeCount) {

	            int lineagesRemoved = 0;
	            int lineagesAdded = 0;

	            double finish = times.get(nodeNo);
	            double next;

	            do {
	                final int childIndex = nodeNo;
	                final IntervalType type = events.get(childIndex);
	                // don't use nodeNo from here on in do loop
	                nodeNo += 1;
	                if (type == IntervalType.SAMPLE) {
	                    lineagesAdded++;
	                } else {
	                    lineagesRemoved++;

	                    // no mix of removed lineages when 0 th
	                    if (multifurcationLimit == 0.0) {
	                        break;
	                    }
	                }

	                if (nodeNo < nodeCount) {
	                    next = times.get(nodeNo);
	                } else break;
	            } while (Math.abs(next - finish) <= multifurcationLimit);

	            if (lineagesAdded > 0) {

	                if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
	                    intervals[intervalCount] = finish - start;
	                    lineageCounts[intervalCount] = numLines;
	                    intervalCount += 1;
	                }

	                start = finish;
	            }

	            // add sample event
	            numLines += lineagesAdded;

	            if (lineagesRemoved > 0) {

	                intervals[intervalCount] = finish - start;
	                lineageCounts[intervalCount] = numLines;
	                intervalCount += 1;
	                start = finish;
	            }
	            // coalescent event
	            numLines -= lineagesRemoved;
	        }
		}

		public void addEvent(double time, IntervalType type) {
			if (times.size() == 0 || times.get(times.size()-1) < time) {
				times.add(time);
				events.add(type);			
			} else {
				int index = Collections.binarySearch(times, time);
				times.add(index, time);
				events.add(index, type);
			}
		}
		
		@Override
		public String toString() {
			if (times == null) {
				return "empty SegmentIntervalList";
			}
			String str = "";
			for (int i = 0;i < times.size(); i++) {
				str += "(" + lineageCounts[i] + " " + (events.get(i) == IntervalType.SAMPLE ? "S": "C") + " " + times.get(i) + ") ";
			}
			return str;
		}
    	
    }
    
	private List<IntervalList> collectSegments() {
		List<IntervalList> segments = new ArrayList<>();
		int nodeCount = tree.getNodeCount();
		for (int i = 0; i < nodeCount; i++) {
			segments.add(null);
		}
		
		for (int i =  0; i < nodeCount; i++) {
			Node node = tree.getNode(i);
			int colour = colourAtBase.getValue(i);
			if (segments.get(colour) == null) {
				segments.set(colour, new SegmentIntervalList());
			}
			SegmentIntervalList intervals = (SegmentIntervalList) segments.get(colour);
			
			// add node event
			intervals.addEvent(node.getHeight(), node.isLeaf() ? IntervalType.SAMPLE : IntervalType.COALESCENT);
			if (!node.isRoot()) {
				int parentNr = node.getParent().getNr();
				int parentColour = colourAtBase.getValue(parentNr);
				if (colour != parentColour) {
					// add sampling event at top of block
					if (segments.get(parentColour) == null) {
						segments.set(parentColour, new SegmentIntervalList());
					}
					intervals = (SegmentIntervalList) segments.get(parentColour);
					double h = node.getHeight() + blockEndFraction.getValue(node.getNr()) * (node.getParent().getHeight() - node.getHeight());
					intervals.addEvent(h, IntervalType.SAMPLE);
				}
			}
		}
		
		for (IntervalList intervals : segments) {
			if (intervals != null) {
				((SegmentIntervalList)intervals).calculateIntervals();
			}
		}
		
		return segments;
	}

	private double calculateCoalescent(IntervalList intervals, double threshold) {
        double logL = 0.0;

        double startTime = 0.0;
        final int n = intervals.getIntervalCount();
        for (int i = 0; i < n; i++) {

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            final double intervalArea = popSizeFunction.getIntegral(startTime, finishTime);
            if (intervalArea == 0 && duration > 1e-10) {
            	/* the above test used to be duration != 0, but that leads to numerical issues on resume
            	 * (https://github.com/CompEvol/beast2/issues/329) */
                return Double.NEGATIVE_INFINITY;
            }
            final int lineageCount = intervals.getLineageCount(i);

            final double kChoose2 = Binomial.choose2(lineageCount);
            // common part
            logL += -kChoose2 * intervalArea;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                final double demographicAtCoalPoint = popSizeFunction.getPopSize(finishTime);

                // if value at end is many orders of magnitude different than mean over interval reject the interval
                // This is protection against cases where ridiculous infinitesimal
                // population size at the end of a linear interval drive coalescent values to infinity.

                if (duration == 0.0 || demographicAtCoalPoint * (intervalArea / duration) >= threshold) {
                    //                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
                    logL -= Math.log(demographicAtCoalPoint);
                } else {
                    // remove this at some stage
                    //  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
                    return Double.NEGATIVE_INFINITY;
                }
            }
            startTime = finishTime;
        }

        return logL;
  	}

	@Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(blockStartFractionInput.get().getID());
        conditions.add(blockEndFractionInput.get().getID());
        conditions.add(blockCountInput.get().getID());
        conditions.add(colourInput.get().getID());
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        arguments.add(treeInput.get().getID());
        return arguments;
    }

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
	}

}
