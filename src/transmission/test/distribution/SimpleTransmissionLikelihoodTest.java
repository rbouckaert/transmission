package transmission.test.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import transmission.distribution.GammaHazardFunction;
import transmission.distribution.HazardFunction;
import transmission.distribution.TransmissionTreeLikelihood;

public class SimpleTransmissionLikelihoodTest {

	@Test
	public void testSimpleCase1() {
		TreeParser tree = new TreeParser("(Bob:1.7,Eve:2.1);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        
        double bs1 = (1.7-0.5)/1.7;
        double bs2 = (2.1-0.65)/2.1;
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 2, "value",  bs1+" "+bs2);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 2, "value",  bs1+" "+bs2);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 2, "value", 0);
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 3, "value", "0 1 2");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "1.0", "shape", "2.5", "rate", "10.0");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.5", "rate", "10.0");
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour,
        		"endTime", "0.1",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "1.0");
        
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        assertEquals(-24.13889, transmissionLikelihood, 1e-5);
	}
}
