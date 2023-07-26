package transmission.test.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import test.beast.BEASTTestCase;
import transmission.distribution.TransmissionTreeLikelihood;


public class CoalescentTest extends BEASTTestCase {
    String[] trees = new String[]{"((A:1.0,B:1.0):1.0,C:2.0);", 
    		"(((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0);"};
    Alignment data;
    final double pop = 10000;

    @BeforeEach
    protected void setUp() throws Exception {
        data = getFourTaxaNoData();
    }

    @Test
    public void testConstantPopulation() throws Exception {
        // *********** 3 taxon **********
        Tree tree0 = getTree(data, trees[0]);

        ConstantPopulation cp = new ConstantPopulation();
        cp.initByName("popSize", Double.toString(pop));

        RealParameter blockStart = new RealParameter(); blockStart.initByName("dimension", 4, "value", 0.25);
        RealParameter blockEnd = new RealParameter(); blockEnd.initByName("dimension", 4, "value", 0.75);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 4, "value", 0);
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 4, "value", 0);
        
        
        TransmissionTreeLikelihood coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree0,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour);

        double logL = coal.calculateCoalescent();

        assertEquals(logL, -(4 / pop) - 2 * Math.log(pop), PRECISION);
        
        // *********** 4 taxon **********
        // contains 3 taxon tree with same colouring
        Tree tree1 = getTree(data, trees[1]);
        
        blockStart.initByName("dimension", 6, "value", 0.25);
        blockEnd.initByName("dimension", 6, "value", 0.75);
        blockcount.initByName("dimension", 6, "value", 0);
        colour.initByName("dimension", 6, "value", "0 0 0 1 0 0 2");
        
        
        coal = new TransmissionTreeLikelihood();
        coal.initByName(
        		"tree", tree1,
        		"populationModel", cp, 
        		"blockstart", blockStart, 
        		"blockend", blockEnd, 
        		"blockcount", blockcount, 
        		"colour", colour);

        List<Double> logLs = coal.calculateCoalescents();

        assertEquals(logLs.get(0), -(4 / pop) - 2 * Math.log(pop), PRECISION);
        assertEquals(logLs.get(1), 0.0, PRECISION);

    }

}
