package breath.test;



import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import breath.distribution.GammaHazardFunction;
import breath.distribution.HazardFunction;
import breath.distribution.TransmissionTreeLikelihood1;

public class TransmissionLikelihoodSpeedTest {

	@Test
	public void testSimpleCase1() {
		boolean approx = true;
		TreeParser tree = new TreeParser("(Bob:1.7,Eve:2.1);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        
        double bs1 = (1.7-0.5)/1.7;
        double bs2 = (2.1-0.65)/2.1;
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 2, "value",  bs1+" "+bs2);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 2, "value",  bs1+" "+bs2);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 2, "value", "0 0");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 3, "value", "0 1 2");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "1.0", "shape", "2.5", "rate", "10.0", "approx", approx);
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0", "approx", approx);
        
        TransmissionTreeLikelihood1 coal = new TransmissionTreeLikelihood1();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "1.0");
        
        coal.calcColourAtBase();
        
        long start = System.currentTimeMillis();
        for (int i = 0; i < 100; i++) {
        	double transmissionLikelihood = coal.calcTransmissionLikelihood();
        }
        long end = System.currentTimeMillis();
        System.err.println((end-start) + " ms");

        
        for (int b = 0; b <= 100; b++) {
        	blockcount.initByName("dimension", 2, "lower", -1, "value", b + " " + b + " " + b);
        	coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		//"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "1.0");

	        start = System.currentTimeMillis();
	        for (int i = 0; i < 10000; i++) {
	        	double transmissionLikelihood = coal.calcTransmissionLikelihood();
	        }
	        end = System.currentTimeMillis();
	        System.err.println(b + " " + (end-start) + " ms");
        }
	}

}
