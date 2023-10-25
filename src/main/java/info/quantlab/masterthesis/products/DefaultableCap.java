package info.quantlab.masterthesis.products;

import java.util.Arrays;

import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This class implements the valuation of a defaultable Cap.
 * The payoff of the product is ( <i> L(T<sub>i</sub>) - K<sub>i</sub></i> )<sup>+</sup>, for <i>i = 0, 1, ... n-1 </i>
 * and <i> n = </i><code>fixingTimes.length</code>.<p/>
 * If the valuation TermStructureModel is non defaultable, this class behaves like a normal cap.
 *
 * @author Markus Parhofer
 * @version 1.0
 */
public class DefaultableCap extends AbstractSellerDefaultableTermStructureProduct {

	final double[] _fixingTimes;
	final double[] _paymentTimes;
	final double[] _strikes;
	final double _notional;
	
	/**
	 * Constructs the class and sets all relevant product specific values.
	 * 
	 * @param fixingTimes The array of fixing times
	 * @param paymentTimes The array of payment times
	 * @param strikes The array of strikes
	 * @param notional The notional
	 */
	public DefaultableCap(double[] fixingTimes, double[] paymentTimes, double[] strikes, double notional) {
		_fixingTimes = fixingTimes;
		_paymentTimes = paymentTimes;
		_strikes = strikes;
		_notional = notional;
	}
	
	/**
	 * Constructs the class and sets all relevant product specific values.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strikes The array of strikes.
	 * @param notional The notional.
	 */
	public DefaultableCap(TimeDiscretization tenor, double[] strikes, double notional) {
		this(Arrays.copyOfRange(tenor.getAsDoubleArray(), 0, tenor.getNumberOfTimeSteps()), 
				Arrays.copyOfRange(tenor.getAsDoubleArray(), 1, tenor.getNumberOfTimes()),
				strikes, notional);
	}
	
	/**
	 * Constructs the class and sets all relevant product specific values.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strike The strike. Sets all strikes to the same value.
	 * @param notional The notional.
	 */
	public DefaultableCap(TimeDiscretization tenor, double strike, double notional) {
		this(tenor, new double[tenor.getNumberOfTimeSteps()], notional);
		Arrays.fill(_strikes, strike);
	}
	
	/**
	 * Constructs the class and sets all relevant product specific values.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strikes The array of strikes.
	 */
	public DefaultableCap(TimeDiscretization tenor, double[] strikes) {
		this(tenor, strikes, 1.0);
	}
	
	/**
	 * Constructs the class and sets all relevant product specific values.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strike The strike. Sets all strikes to the same value.
	 */
	public DefaultableCap(TimeDiscretization tenor, double strike) {
		this(tenor, strike, 1.0);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model)	throws CalculationException {

		// Allocate accumulator for values
		RandomVariable values = new RandomVariableFromDoubleArray(0.0);

		for(int period=0; period<_fixingTimes.length; period++)
		{
			final double fixingDate	= _fixingTimes[period];
			final double paymentDate	= _paymentTimes[period];

			// evaluationTime > fixingDate is allowed. Negative fixing date is allowed too (but likely not supported by the model)
			if(evaluationTime > paymentDate) {
				continue;
			}

			final double strike	 	= _strikes[period];
			final double periodLength	= paymentDate - fixingDate;

			// Get random variables
			final RandomVariable	libor					= model.getForwardRate(fixingDate, fixingDate, paymentDate);
			final RandomVariable	numeraire				= model.getNumeraire(paymentDate);
			final RandomVariable	monteCarloWeights	= model.getMonteCarloWeights(model.getTimeIndex(paymentDate));

			// Calculate payout
			RandomVariable payoff = libor.sub(strike).mult(periodLength).floor(0.0);

			// Adjust for possibility of default between evaluationTime and paymentDate:
			if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
				payoff = payoff.mult(defaultableModel.getSurvivalProbability(model.getProcess(), evaluationTime, paymentDate));
			}
			
			payoff = payoff.div(numeraire).mult(monteCarloWeights);

			// Accumulate numeraire relative values
			values = values.add(payoff);
		}

		final RandomVariable	numeraireAtEvaluationTime				= model.getNumeraire(evaluationTime);
		final RandomVariable	monteCarloWeightsAtEvaluationTime	= model.getMonteCarloWeights(evaluationTime);
		values = values.mult(numeraireAtEvaluationTime).div(monteCarloWeightsAtEvaluationTime);

		return values;
	}

}
