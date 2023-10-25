package info.quantlab.masterthesis.products;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.stochastic.RandomVariable;

/**
 * This class is a wrapper for Defaultable TermStructure Products. The valuation should always be performed on a 
 * TermStructureMonteCarloSimulationModel that is associated with a DefaultableLIBORMarketModel. 
 * We differ between different cases:
 * <ul>
 * <li> <b>AbstractBuyerDefaultableTermStructureProduct:</b> The model assumes the buyer is defaultable. 
 * 		This of course constitutes that the buyer has an obligation in the future, that they might not be able to meet. 
 * 		<b>Caution!:</b> As the buyer will sell the contract, if it is positive at maturity, even products, where the buyer 
 * 		has to invest money to make use of the option (e.g. a call option) are not considered buyer-defaultable.
 * 		Could be used for non transferable options though.</li>
 * <li> <b>AbstractSellerDefaultableTermStructureProduct (this class):</b> A seller that is defaultable. E.g. typical vanilla options, 
 * 		where the buyer has the right but not the obligation to buy/sell the underlying.</li>
 * <li> <b>AbstractDefaultableTermStructureProduct:</b> The buyer as well as the seller are assumed to be defaultable. E.g. Swaps, Swaptions 
 * 		and Forwards.</li>
 * </ul>
 * @author Markus Parhofer
 *
 */
public abstract class AbstractSellerDefaultableTermStructureProduct extends AbstractTermStructureMonteCarloProduct {

	public AbstractSellerDefaultableTermStructureProduct() {
		super();
	}

	public AbstractSellerDefaultableTermStructureProduct(String currency) {
		super(currency);
	}

	/**
	 * Gets the value of the product. The value is always F<sub>t</sub> measurable. Furthermore, it is always given that default has not 
	 */
	@Override
	public abstract RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException;

}
