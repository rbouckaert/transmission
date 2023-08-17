package transmission.test.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import javax.xml.crypto.NodeSetData;

import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.Node;
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
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "1.0");
        
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        assertEquals(-24.8279, transmissionLikelihood, 1e-5);
	}


	@Test
	public void testSimpleCase2() {
		TreeParser tree = new TreeParser("(Bob:1.7,Eve:2.1);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        
        double bs1 = (1.7-0.5)/1.7;
        double bs2 = (2.1-0.65)/2.1;
        double bs3 = (2.1-1.8)/2.1;
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 2, "value",  bs1+" "+bs3);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 2, "value",  bs1+" "+bs2);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 2, "value", "0 4");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 3, "value", "0 1 2");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "1.0", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "4.0");
        
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        assertEquals( -20.37201, transmissionLikelihood, 1e-5);
	}
	
	@Test
	public void testSimpleCase3() {
		TreeParser tree = new TreeParser("((t1:0.6587438122,t2:0.22):0.7863448577,(t3:0.3307722867,(t4:0.7084983373,t5:0.6330101104):0.6262222228):0.2);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        Node [] nodes = tree.getNodesAsArray();
        
        double h = tree.getRoot().getHeight();
        double e4 = (h-0.96 - nodes[4].getHeight()) / nodes[4].getLength(); 
        double e5 = (h-0.35 - nodes[5].getHeight()) / nodes[5].getLength(); 
        double e6 = (h-0.22 - nodes[6].getHeight()) / nodes[6].getLength(); 

        double s4 = (h-1.33 - nodes[4].getHeight()) / nodes[4].getLength(); 
        double s5 = (h-0.72 - nodes[5].getHeight()) / nodes[5].getLength(); 
        double s6 = (h-0.50 - nodes[6].getHeight()) / nodes[6].getLength(); 
        
        
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 8,       "value", "0.5 0.5 0.5 0.5 " + s4 + " " + s5 + " " + s6 + " 0.5" );
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 8,           "value", "0.5 0.5 0.5 0.5 " + e4 + " " + e5 + " " + e6 + " 0.5");
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 8, "value", "0 0 0 0 2 4 3 0");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 9,         "value", "0 1 2 3 4 0 3 2 8");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "1.0", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "4.0");
        
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        assertEquals( -26.10661, transmissionLikelihood, 1e-5);
	}	

	@Test
	public void testSimpleCase4() {
		TreeParser tree = new TreeParser("(t3:0.3307722867,(t4:0.7084983373,t5:0.6330101104):0.6262222228);");
		
        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(1.0));

        Node [] nodes = tree.getNodesAsArray();
        
        double h = tree.getRoot().getHeight();
        double e4 = (h-0.73 - nodes[2].getHeight()) / nodes[2].getLength(); 
        double e3 = (h-0.12 - nodes[0].getHeight()) / nodes[0].getLength(); 

        double s4 = (h-1.05 - nodes[2].getHeight()) / nodes[2].getLength(); 
        double s3 = e3; 
        
        
        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 4,       "value", s3 + " 0.5 " + s4 + " 0.5" );
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 4,           "value", e3 + " 0.5 " + e4 + " 0.5");
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 4, "value", "0 0 2 0");
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 5,         "value", "0 1 2 1 1");
        
        HazardFunction samplingHazard = new GammaHazardFunction();
        samplingHazard.initByName("C", "1.0", "shape", "2.5", "rate", "10.0");
        samplingHazard.setID("s");
        
        HazardFunction transmissionHazard = new GammaHazardFunction();
        transmissionHazard.initByName("C", "1.5", "shape", "2.0", "rate", "10.0");
        transmissionHazard.setID("tr");
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour,
        		"endTime", "0.0",
        		"samplingHazard", samplingHazard,
        		"transmissionHazard", transmissionHazard,
        		"lambda", "4.0");
        
        double transmissionLikelihood = coal.calcTransmissionLikelihood();
        
        assertEquals( -12.76696, transmissionLikelihood, 1e-5);
	}
}
