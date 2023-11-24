package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
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
 * {@link MultiLIBORVectorModel#getArrayOfDefaultableModels()} (if non-negative) or the {@link MultiLIBORVectorModel#getUndefaultableModel()} 
 * (if negative).</li>
 * <li>If it is a DefaultableLIBORMarketModel, and the underlying- resp. issuerModelIndex is non-negative the purely defaultable model is used, 
 * else the {@link DefaultableLIBORMarketModel#getUndefaultableLIBORModel()}.</li>
 * </ul>
 * 
 * @author Markus Parhofer
 * @version 1.0
 */
public class DefaultableCaplet extends AbstractSellerDefaultableTermStructureProduct {

	final double _strikeRate;
	final double _fixingTime;
	final double _periodLength;
	final double _notional;
	final boolean _isFloorlet;
	final int _underlyingModelIndex; /* if <0 underlying is undefaultable Model else it is a defaultable model with the specified index. */
	final int _issuerModelIndex; /* if <0 issuer has undefaultable Model else it is a defaultable model with the specified index. */
	
	/**
	 * Constructs a Caplet/Floorlet on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer and 
	 * if so, which one. 
	 * See Java Doc on {@link DefaultableCaplet} for more information.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param modelIndexOfUnderlying The index for the defaultable model to use as the underlying. Set -1 for non-defaultable model.
	 * @param modelIndexOfIssuer The index for the defaultable model to use as the issuer. Set -1 for non-defaultable model.
	 * @param notional The notional of the product.
	 * @param isFloorlet A flag that specifies if the product is a floorlet instead of a caplet.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength, 
			final int modelIndexOfUnderlying, final int modelIndexOfIssuer, final double notional, final boolean isFloorlet) {
		super();
		_strikeRate = strikeRate;
		_fixingTime = fixingTime;
		_periodLength = periodLength;
		_underlyingModelIndex = modelIndexOfUnderlying;
		_issuerModelIndex = modelIndexOfIssuer;
		_notional = notional;
		_isFloorlet = isFloorlet;
	}
	
	/**
	 * Constructs a Caplet/Floorlet on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer.
	 * If a defaultable model should be used it is always the default one with model Index 0.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param useDefaultableUnderlying Specifies, if the product has a defaultable underlying. If so, the default model index (0) is used.
	 * @param useDefaultableIssuer Specifies, if the issuer is defaultable. If so, the default model index (0) is used.
	 * @param isFloorlet A flag that specifies if the product is a floorlet instead of a caplet.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength, 
			final boolean useDefaultableUnderlying, final boolean useDefaultableIssuer, final boolean isFloorlet) {
		this(strikeRate, fixingTime, periodLength, useDefaultableUnderlying ? 0 : -1, useDefaultableIssuer ? 0 : -1, 1.0, isFloorlet);
	}
	
	/**
	 * Constructs a Caplet on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer and 
	 * if so, which one. 
	 * See Java Doc on {@link DefaultableCaplet} for more information.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param modelIndexOfUnderlying The index for the defaultable model to use as the underlying. Set -1 for non-defaultable model.
	 * @param modelIndexOfIssuer The index for the defaultable model to use as the issuer. Set -1 for non-defaultable model.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength, 
			final int modelIndexOfUnderlying, final int modelIndexOfIssuer) {
		this(strikeRate, fixingTime, periodLength, modelIndexOfUnderlying, modelIndexOfIssuer, 1.0, false);
	}
	
	/**
	 * Constructs a Caplet on a LIBOR rate, where one can specify if a defaultable Model should be used for the underlying or for the issuer.
	 * If a defaultable model should be used it is always the default one with model Index 0.
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param useDefaultableUnderlying Specifies, if the product has a defaultable underlying. If so, the default model index (0) is used.
	 * @param useDefaultableIssuer Specifies, if the issuer is defaultable. If so, the default model index (0) is used.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength, 
			final boolean useDefaultableUnderlying, final boolean useDefaultableIssuer) {
		this(strikeRate, fixingTime, periodLength, useDefaultableUnderlying ? 0 : -1, useDefaultableIssuer ? 0 : -1, 1.0, false);
	}
	
	/**
	 * Constructs a Caplet/Floorlet on a defaultable LIBOR rate, but without issuer risk (i.e. modelIndexOfUnderlying is 0 and modelIndexOfIssuer is -1).
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 * @param isFloorlet A flag that specifies if the product is a floorlet instead of a caplet.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double paymentTime, final boolean isFloorlet) {
		this(strikeRate, fixingTime, paymentTime - fixingTime, 0, -1, 1.0, isFloorlet);
	}
	
	/**
	 * Constructs a Caplet on a defaultable LIBOR rate, but without issuer risk (i.e. modelIndexOfUnderlying is 0 and modelIndexOfIssuer is -1).
	 * 
	 * @param strikeRate The fixed strike rate.
	 * @param fixingTime The maturity date where the interest rate is fixed and the interest period begins.
	 * @param periodLength The length of the interest period.
	 */
	public DefaultableCaplet(final double strikeRate, final double fixingTime, final double periodLength) {
		this(strikeRate, fixingTime, periodLength, 0, -1, 1.0, false);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
		// Get random variables
		final double paymentTime = _fixingTime + _periodLength;
		RandomVariable forwardRate						= model.getForwardRate(_fixingTime, _fixingTime, paymentTime);
		RandomVariable survivalProbability 				= model.getRandomVariableForConstant(1.0);
		RandomVariable survivalProbabilityUnderlying 	= model.getRandomVariableForConstant(1.0);
		final RandomVariable	numeraire				= model.getNumeraire(paymentTime);
		final RandomVariable	monteCarloWeights		= model.getMonteCarloWeights(paymentTime);

		if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
			if(_issuerModelIndex >= 0)
				survivalProbability = defaultableModel.getSurvivalProbability(model.getProcess(), evaluationTime, paymentTime);
			if(_underlyingModelIndex < 0)
				forwardRate = defaultableModel.getUndefaultableForwardRate(model.getProcess(), _fixingTime, _fixingTime, paymentTime);
			else
				survivalProbabilityUnderlying = defaultableModel.getSurvivalProbability(model.getProcess(), evaluationTime, _fixingTime);
		}
		else if(model.getModel() instanceof MultiLIBORVectorModel multiModel) {
			if(_issuerModelIndex >= 0)
				survivalProbability = multiModel.getSurvivalProbability(model.getProcess(), evaluationTime, paymentTime, _issuerModelIndex);
			if(_underlyingModelIndex >= 0) {
				forwardRate = multiModel.getDefaultableForwardRate(model.getProcess(), _fixingTime, _fixingTime, paymentTime, _underlyingModelIndex);
				survivalProbabilityUnderlying = multiModel.getSurvivalProbability(model.getProcess(), evaluationTime, _fixingTime, _underlyingModelIndex);
			}
		}
		/*
		 * Calculate the payoff, which is
		 *    max(L-K,0) * periodLength * notional        for caplet or
		 *   -min(L-K,0) * periodLength * notional        for floorlet.
		 */
		RandomVariable values = forwardRate;
		if(!_isFloorlet) {
			// We only get the payoff if the underlying survives until Fixingtime!
			values = values.sub(_strikeRate).floor(0.0).mult(survivalProbabilityUnderlying).mult(_periodLength * _notional);
		} else {
			// If the underlying does not survive until Fixingtime we get the maximum payoff (K).
			values = values.sub(_strikeRate).cap(0.0).mult(survivalProbabilityUnderlying).mult(-1.0 * _periodLength * _notional)
					.add(survivalProbabilityUnderlying.bus(1.0).mult(_strikeRate));
		}

		values = values.div(numeraire).mult(monteCarloWeights);
		
		// Adjust for possibility of default between evaluationTime and paymentDate:
		values = values.mult(survivalProbability);
		
		final RandomVariable	numeraireAtValuationTime				= model.getNumeraire(evaluationTime);
		final RandomVariable	monteCarloWeightsAtValuationTime	= model.getMonteCarloWeights(evaluationTime);
		values = values.mult(numeraireAtValuationTime).div(monteCarloWeightsAtValuationTime);

		return values;
	}

}
