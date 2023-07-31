package transmission.test.distribution;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import test.beast.BEASTTestCase;
import transmission.distribution.Validator;


public class ValidatorTest extends BEASTTestCase {
    String[] trees = new String[]{"((A:1.0,B:1.0):1.0,C:2.0);", 
    		"(((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0);"};
    Alignment data;

    @BeforeEach
    protected void setUp() throws Exception {
        data = getFourTaxaNoData();
    }

    @Test
    public void testThreeTaxonTree() throws Exception {
        Tree tree0 = getTree(data, trees[0]);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 4, "value", 0);
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 5, "value", "0 1 2 0 0");
        
        Validator validator = new Validator(tree0, colour, blockcount);
        assertEquals(validator.isValid(), true);
        blockcount.setValue(0, 1);
        assertEquals(validator.isValid(), false);
    }

    @Test
    public void testFourTaxonTree() throws Exception {
        Tree tree0 = getTree(data, trees[1]);
        IntegerParameter blockcount = new IntegerParameter(); blockcount.initByName("dimension", 6, "value", 0);
        IntegerParameter colour = new IntegerParameter(); colour.initByName("dimension", 7, "value", "0 1 2 3 0 0 0");
        
        Validator validator = new Validator(tree0, colour, blockcount);
        assertEquals(validator.isValid(), true);
        blockcount.setValue(0, 1);
        assertEquals(validator.isValid(), false);
    }
}
