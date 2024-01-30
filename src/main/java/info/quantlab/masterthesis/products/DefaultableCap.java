package info.quantlab.masterthesis.products;

import java.util.Arrays;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This class implements the valuation of a Cap where a possibility for default is accounted for.
 * The payoff of the product is ( <i> L(T<sub>i</sub>) - K<sub>i</sub></i> )<sup>+</sup>, for <i>i = 0, 1, ... n-1 </i>
 * and <i> n = </i><code>fixingTimes.length</code>.<p/>
 * The model allows for a defaultable underlying and/or a defaultable issuer. If the underlying- resp. issuerModelIndex is non-negative the 
 * corresponding model used is defaultable. This results in:
 * <ul>
 * <li>If the underlying is defaultable, the underlying LIBOR rate is assumed to be defaultable, but the numeraire stays the undefaultable one.</li>
 * <li>If the issuer is defaultable the conditional expectation of the survival probability is multiplied with the original price.</li>
 * </ul>
 * The model behaves differently depending on the ProcessModel of the valuation TermStructureMonteCarloSimulationModel:
 * <ul>
 * <li>If it is a non defaultable ProcessModel, this class behaves like a normal cap.</li>
 * <li>If it is a MultiLIBORVectorModel, the underlying- resp. issuerModelIndex correspond to either the zero-based index of 
 * {@link MultiLIBORVectorModel#getArrayOfDefaultableModels()} (if non-negative) or the {@link MultiLIBORVectorModel#getNonDefaultableModel()}
 * (if negative).</li>
 * <li>If it is a DefaultableLIBORMarketModel, and the underlying- resp. issuerModelIndex is non-negative the purely defaultable model is used, 
 * else the {@link DefaultableLIBORMarketModel#getNonDefaultableLIBORModel()}.</li>
 * </ul>
 * 
 * @author Markus Parhofer
 * @version 1.0
 */
public class DefaultableCap extends AbstractDefaultableTermStructureProduct {

	final double[] _fixingTimes;
	final double[] _paymentTimes;
	final double[] _strikes;
	final int _underlyingModelIndex; /* if <0 underlying is undefaultable Model else it is a defaultable model with the specified index. */
	final int _issuerModelIndex; /* if <0 issuer has undefaultable Model else it is a defaultable model with the specified index. */
	final double _notional;
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer and 
	 * if so, which one. 
	 * See Java Doc on {@link DefaultableCap} for more information.
	 * 
	 * @param fixingTimes The array of fixing times
	 * @param paymentTimes The array of payment times
	 * @param strikes The array of strikes
	 * @param modelIndexOfUnderlying The index for the defaultable model to use as the underlying. Set -1 for non-defaultable model.
	 * @param modelIndexOfIssuer The index for the defaultable model to use as the issuer. Set -1 for non-defaultable model.
	 * @param notional The notional
	 */
	public DefaultableCap(double[] fixingTimes, double[] paymentTimes, double[] strikes, int modelIndexOfUnderlying, int modelIndexOfIssuer, double notional) {
		_fixingTimes = fixingTimes;
		_paymentTimes = paymentTimes;
		_strikes = strikes;
		_underlyingModelIndex = modelIndexOfUnderlying;
		_issuerModelIndex = modelIndexOfIssuer;
		_notional = notional;
	}
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer.
	 * If a defaultable model should be used it is always the default one with model Index 0.
	 * 
	 * @param fixingTimes The array of fixing times
	 * @param paymentTimes The array of payment times
	 * @param strikes The array of strikes
	 * @param useDefaultableUnderlying Specifies, if the product has a defaultable underlying. If so, the default model index (0) is used.
	 * @param useDefaultableIssuer Specifies, if the issuer is defaultable. If so, the default model index (0) is used.
	 * @param notional The notional
	 */
	public DefaultableCap(double[] fixingTimes, double[] paymentTimes, double[] strikes, boolean useDefaultableUnderlying, boolean useDefaultableIssuer, double notional) {
		this(fixingTimes, paymentTimes, strikes, useDefaultableUnderlying ? 0 : -1, useDefaultableIssuer ? 0 : -1, notional);
	}
	
	/**
	 * Constructs a Cap on a defaultable LIBOR rate, but without issuer risk (i.e. modelIndexOfUnderlying is 0 and modelIndexOfIssuer is -1).
	 * 
	 * @param fixingTimes The array of fixing times
	 * @param paymentTimes The array of payment times
	 * @param strikes The array of strikes
	 * @param notional The notional
	 */
	public DefaultableCap(double[] fixingTimes, double[] paymentTimes, double[] strikes, double notional) {
		this(fixingTimes, paymentTimes, strikes, 0, -1, notional);
	}
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer and 
	 * if so, which one. 
	 * See Java Doc on {@link DefaultableCap} for more information.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strikes The array of strikes.
	 * @param modelIndexOfUnderlying The index for the defaultable model to use as the underlying. Set -1 for non-defaultable model.
	 * @param modelIndexOfIssuer The index for the defaultable model to use as the issuer. Set -1 for non-defaultable model.
	 */
	public DefaultableCap(TimeDiscretization tenor, double[] strikes, int modelIndexOfUnderlying, int modelIndexOfIssuer) {
		this(Arrays.copyOfRange(tenor.getAsDoubleArray(), 0, tenor.getNumberOfTimeSteps()), 
				Arrays.copyOfRange(tenor.getAsDoubleArray(), 1, tenor.getNumberOfTimes()),
				strikes, modelIndexOfIssuer, modelIndexOfUnderlying, 1.0);
	}
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer.
	 * If a defaultable model should be used it is always the default one with model Index 0.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strikes The array of strikes.
	 * @param useDefaultableUnderlying Specifies, if the product has a defaultable underlying. If so, the default model index (0) is used.
	 * @param useDefaultableIssuer Specifies, if the issuer is defaultable. If so, the default model index (0) is used.
	 */
	public DefaultableCap(TimeDiscretization tenor, double[] strikes, boolean useDefaultableUnderlying, boolean useDefaultableIssuer) {
		this(tenor, strikes, useDefaultableUnderlying? 0 : -1, useDefaultableIssuer ? 0 : -1);
	}
	
	/**
	 * Constructs a Cap on a defaultable LIBOR rate, but without issuer risk (i.e. modelIndexOfUnderlying is 0 and modelIndexOfIssuer is -1).
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strikes The array of strikes.
	 */
	public DefaultableCap(TimeDiscretization tenor, double[] strikes) {
		this(tenor, strikes, 0, -1);
	}
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer and 
	 * if so, which one. 
	 * See Java Doc on {@link DefaultableCap} for more information.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strike The strike. Sets all strikes to the same value.
	 * @param modelIndexOfUnderlying The index for the defaultable model to use as the underlying. Set -1 for non-defaultable model.
	 * @param modelIndexOfIssuer The index for the defaultable model to use as the issuer. Set -1 for non-defaultable model.
	 */
	public DefaultableCap(TimeDiscretization tenor, double strike, int modelIndexOfUnderlying, int modelIndexOfIssuer) {
		this(tenor, new double[tenor.getNumberOfTimeSteps()], modelIndexOfUnderlying, modelIndexOfIssuer);
		Arrays.fill(_strikes, strike);
	}
	
	/**
	 * Constructs a Cap on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer.
	 * If a defaultable model should be used it is always the default one with model Index 0.
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strike The strike. Sets all strikes to the same value.
	 * @param useDefaultableUnderlying Specifies, if the product has a defaultable underlying. If so, the default model index (0) is used.
	 * @param useDefaultableIssuer Specifies, if the issuer is defaultable. If so, the default model index (0) is used.
	 */
	public DefaultableCap(TimeDiscretization tenor, double strike, boolean useDefaultableUnderlying, boolean useDefaultableIssuer) {
		this(tenor, strike, useDefaultableUnderlying ? 0 : -1, useDefaultableIssuer ? 0 : -1);
	}
	
	/**
	 * Constructs a Cap on a defaultable LIBOR rate, but without issuer risk (i.e. modelIndexOfUnderlying is 0 and modelIndexOfIssuer is -1).
	 * 
	 * @param tenor The fixing and payment times tenor, i.e. <code>fixingTimes</code>=(tenor<sub>i</sub>)<sub>i=0,...,n-1</sub>
	 * and <code>paymentTimes</code>=(tenor<sub>i</sub>)<sub>i=1,...,n</sub>.
	 * @param strike The strike. Sets all strikes to the same value.
	 */
	public DefaultableCap(TimeDiscretization tenor, double strike) {
		this(tenor, strike, 0, -1);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model)	throws CalculationException {

		// Allocate accumulator for values
		RandomVariable values = new RandomVariableFromDoubleArray(0.0);

		for(int period = 0; period < _fixingTimes.length; period++)
		{
			final double fixingDate	= _fixingTimes[period];
			final double paymentTime	= _paymentTimes[period];

			// evaluationTime > fixingDate is allowed. Negative fixing date is allowed too (but likely not supported by the model)
			if(evaluationTime > paymentTime) {
				continue;
			}

			final double strike	 	= _strikes[period];
			final double periodLength	= paymentTime - fixingDate;

			// Get random variables
			RandomVariable	forwardRate					= model.getForwardRate(fixingDate, fixingDate, paymentTime);
			RandomVariable survivalProbability			= model.getRandomVariableForConstant(1.0);
			final RandomVariable	numeraire			= model.getNumeraire(paymentTime);
			final RandomVariable	monteCarloWeights	= model.getMonteCarloWeights(paymentTime);

			
			if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
				if(_issuerModelIndex >= 0)
					survivalProbability = defaultableModel.getSurvivalProbability(model.getProcess(), paymentTime);
				if(_underlyingModelIndex < 0)
					forwardRate = defaultableModel.getNonDefaultableForwardRate(model.getProcess(), fixingDate, fixingDate, paymentTime);
			}
			else if(model.getModel() instanceof MultiLIBORVectorModel multiModel) {
				if(_issuerModelIndex >= 0)
					survivalProbability = multiModel.getSurvivalProbability(model.getProcess(), paymentTime, _issuerModelIndex);
				if(_underlyingModelIndex >= 0)
					forwardRate = multiModel.getDefaultableForwardRate(model.getProcess(), fixingDate, fixingDate, paymentTime, _underlyingModelIndex);
			}
			
			// Calculate payout
			RandomVariable payoff = forwardRate.sub(strike).mult(periodLength).floor(0.0).mult(_notional);

			// Adjust for possibility of default between evaluationTime and paymentDate:
			payoff = payoff.mult(survivalProbability);

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
