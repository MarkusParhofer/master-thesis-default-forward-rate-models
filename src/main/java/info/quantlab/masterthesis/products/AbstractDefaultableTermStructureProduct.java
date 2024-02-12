package info.quantlab.masterthesis.products;

import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;

/**
* This class is a wrapper for Defaultable TermStructure Products. The valuation should always be performed on a 
* TermStructureMonteCarloSimulationModel that is associated with a DefaultableLIBORMarketModel or a MultiLIBORVectorModel.
*/
public abstract class AbstractDefaultableTermStructureProduct extends AbstractTermStructureMonteCarloProduct {

	public AbstractDefaultableTermStructureProduct() {
	}

	public AbstractDefaultableTermStructureProduct(String currency) {
		super(currency);
	}


}
