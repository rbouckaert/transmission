package breath.test;

import java.io.IOException;
import java.io.PrintStream;

import beast.base.util.Randomizer;

public class BlockOperatorTest {

	int blockCount = 0;
	double blockStartFraction = 0.5;
	double blockEndFraction = blockStartFraction;

	private void test() throws IOException {
		// sample single branch block 
		// on uniform distribution
		int N = 20000; // nr of steps
		
		PrintStream out = new PrintStream("/tmp/BlockOperatorTest.log");
		out.println("Sample\tBlockCount\tBlockStart\tBlockEnd\t");
		int acceptCount = 0;
		for (int i = 0; i < N; i++) {
			// store
			int oldBlockCount = blockCount;
			double oldStart = blockStartFraction;
			double oldEnd = blockEndFraction;
			
			double logAlpha = propose();
			double oldLogP = -(1+oldBlockCount) / 3.0; 
			double newLogP = -(1+blockCount) / 3.0; 
			logAlpha += newLogP - oldLogP;
			
            if (logAlpha >= 0 || (logAlpha != Double.NEGATIVE_INFINITY && Randomizer.nextDouble() < Math.exp(logAlpha))) {
            	// accept
            	acceptCount++;
            } else {
            	// reject
            	blockCount = oldBlockCount;
            	blockStartFraction = oldStart;
            	blockEndFraction = oldEnd;
            }

    		out.println(i + "\t"+ blockCount + "\t" + blockStartFraction + "\t" + 
    		blockEndFraction + "\t");
		}
		out.close();
		
		System.err.println("Accepts: " + acceptCount + " out of " + N);
		
	}
	
	
	
	
	private double propose() {
		// new proposal
		if (Randomizer.nextBoolean()) {
			// only move start and end fraction but not block count
			switch (blockCount) {
			case -1:
				// nothing to do since start and end fractions are ignored
				break;
			case 0:
				// make sure start == end fraction after proposal
				double f = Randomizer.nextDouble();
				blockStartFraction = f;
				blockEndFraction = f;
				break;
			default:
				// move one boundary only
				if (Randomizer.nextBoolean()) {
					f = Randomizer.nextDouble() * blockEndFraction;
					blockStartFraction = f;
				} else {
					f = 1-Randomizer.nextDouble() * (1 - blockStartFraction);
					blockEndFraction = f;					
				}
			}
			return 0;
		} else {
			if (Randomizer.nextBoolean()) {
				// reduce block count
				switch (blockCount) {
				case -1:
					// fail
					return Double.NEGATIVE_INFINITY;
				case 0:
					blockCount--;
					return 0;
				case 1:
					blockCount--;
//					double delta = Math.min(blockStartFraction, 1-blockEndFraction);
//					blockEndFraction = (blockEndFraction + blockStartFraction) / 2.0;
//					return -Math.log(delta);
					
					if (Randomizer.nextBoolean()) {
						double delta = 1 - blockStartFraction;
						blockEndFraction = blockStartFraction;
						return -Math.log(1.0/delta);
					} else {
						double delta = blockEndFraction;
						blockStartFraction = blockEndFraction;
						return -Math.log(1.0/delta);
					}
				default:
					blockCount--;
					return 0;
				}
			} else {
				// increase block count
				switch (blockCount) {
				case -1:
					blockCount++;
					double f = Randomizer.nextDouble();
					blockStartFraction = f;
					blockEndFraction = f;
					return 0;
				case 0:
					blockCount++;
					f = Randomizer.nextDouble();
//					double r = blockStartFraction;
//					double delta = Math.min(blockStartFraction, 1-blockEndFraction);
//					blockStartFraction = r - f * delta;
//					blockEndFraction = r + f * delta;
//					return Math.log(delta);

					
					if (Randomizer.nextBoolean()) {
						double delta = 1 - blockStartFraction;
						blockEndFraction = blockStartFraction + (1-blockStartFraction) * f;
						return Math.log(1.0/delta);
					} else {
						double delta = blockEndFraction;
						blockStartFraction = blockEndFraction * f;
						return Math.log(1.0/delta);
					}
				case 11114:
					return Double.NEGATIVE_INFINITY;
				default:
					blockCount++;
					return 0;
				}
			}
		}
	}




	public static void main(String[] args) throws IOException {
		new BlockOperatorTest().test();
	}

}
