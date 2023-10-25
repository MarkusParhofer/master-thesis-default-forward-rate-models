package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

/**
 * Class for the valuation of a defaultable caplet or floorlet, hence a product with the payoff 
 * <ul> 
 * <li><i>(L - K)<sup>+</sup> * (T - S) * N </i> for a caplet and</li>
 * <li><i>(K - L)<sup>+</sup> * (T - S) * N </i> for a floorlet,</li>
 * </ul>
 * where <i>L</i> is the forward rate, <i>K</i> is the strike-rate, <i>T</i> is the payment date, <i>S</i> the fixing date of the forward rate and <i>N</i> is the notional.<p/>
 * If the valuation Term Structure Model is non defaultable, this class behaves like a normal caplet/floorlet.
 * 
 * @author Markus Parhofer
 *
 */
public class DefaultableCaplet extends AbstractSellerDefaultableTermStructureProduct {

	final double _strikeRate;
	final double _fixingTime;
	final double _periodLength;
	final double _notional;
	final boolean _isFloorlet;
	
	/**
	 * Constructs the class and sets all product specific details.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param notional The notional of the product.
	 * @param isFloorlet A flag thatspecifies if the product is a floorlet instead of a caplet.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength, final double notional, final boolean isFloorlet) {
		super();
		_strikeRate = strikeRate;
		_fixingTime = fixingTime;
		_periodLength = periodLength;
		_notional = notional;
		_isFloorlet = isFloorlet;
	}
	
	/**
	 * Constructs the class and sets all product specific details.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity time where the interest rate is fixed and the interest period begins.
	 * @param paymentTime The time where the interest period ends and the product is paid.
	 * @param notional The notional of the product.
	 * @param isFloorlet A flag that specifies if the product is a floorlet instead of a caplet.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double paymentTime, final boolean isFloorlet) {
		this(strikeRate, fixingTime, paymentTime - fixingTime, 1.0, isFloorlet);
	}
	
	/**
	 * Constructs the class and sets all product specific details.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength) {
		this(strikeRate, fixingTime, periodLength, 1.0, false);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
		// Get random variables
		final double paymentTime = _fixingTime + _periodLength;
		final RandomVariable	forwardRate				= model.getForwardRate(_fixingTime, _fixingTime, paymentTime);
		final RandomVariable	numeraire				= model.getNumeraire(paymentTime);
		final RandomVariable	monteCarloWeights	= model.getMonteCarloWeights(paymentTime);

		/*
		 * Calculate the payoff, which is
		 *    max(L-K,0) * periodLength * notional        for caplet or
		 *   -min(L-K,0) * periodLength * notional        for floorlet.
		 */
		RandomVariable values = forwardRate;
		if(!_isFloorlet) {
			values = values.sub(_strikeRate).floor(0.0).mult(_periodLength * _notional);
		} else {
			values = values.sub(_strikeRate).cap(0.0).mult(-1.0 * _periodLength * _notional);
		}

		values = values.div(numeraire).mult(monteCarloWeights);

		// Adjust for possibility of default between evaluationTime and paymentDate:
		if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
			values = values.mult(defaultableModel.getSurvivalProbability(model.getProcess(), evaluationTime, paymentTime));
		}
		
		final RandomVariable	numeraireAtValuationTime				= model.getNumeraire(evaluationTime);
		final RandomVariable	monteCarloWeightsAtValuationTime	= model.getMonteCarloWeights(evaluationTime);
		values = values.mult(numeraireAtValuationTime).div(monteCarloWeightsAtValuationTime);

		return values;
	}

}
