/**
 * 
 */
package info.quantlab.masterthesis.defaultablecovariancemodels;

import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.stochastic.RandomVariable;

/**
 * The abstract class for the implementation of a Covariance Model for defaultable LIBOR Market Models.
 * 
 * @author Markus Parhofer
 */
public abstract class AbstractDefaultableLIBORCovarianceModel implements DefaultableLIBORCovarianceModel {
	
	private final LIBORMarketModel _underlyingUndefaultableModel;
	
	/**
	 * Initializes necessary members of the abstract model
	 * @param underlyingUndefaultableModel The undefaultable reference Model
	 */
	public AbstractDefaultableLIBORCovarianceModel(LIBORMarketModel underlyingUndefaultableModel) {
		_underlyingUndefaultableModel = underlyingUndefaultableModel;
	}

	@Override
	public	RandomVariable[] getFactorLoading(final double time, final double component, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization) {
		int componentIndex = getLiborPeriodDiscretization().getTimeIndex(component);
		if(componentIndex < 0) {
			componentIndex = -componentIndex - 2;
		}
		return getFactorLoading(time, componentIndex, defaultableRealization, undefaultableRealization);
	}

	@Override
	public	RandomVariable[] getFactorLoading(final double time, final int component, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = -timeIndex - 2;
		}
		return getFactorLoading(timeIndex, component, defaultableRealization, undefaultableRealization);
	}

	@Override
	public RandomVariable getCovariance(double time, int component1, int component2, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = -timeIndex - 2;
		}
		return getCovariance(timeIndex, component1, component2, defaultableRealization, undefaultableRealization);
	}

	@Override
	public RandomVariable getCovariance(int timeIndex, int component1, int component2, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization) {
		final RandomVariable[] factorLoadingC1 = getFactorLoading(timeIndex, component1, defaultableRealization, undefaultableRealization);
		final RandomVariable[] factorLoadingC2 = getFactorLoading(timeIndex, component2, defaultableRealization, undefaultableRealization);
		
		RandomVariable covariance = factorLoadingC1[0].mult(factorLoadingC2[0]);
		for(int factor = 1; factor < getNumberOfFactors(); factor++) {
			covariance = covariance.add(factorLoadingC1[factor].mult(factorLoadingC2[factor]));
		}
		
		return covariance;
	}
	
	@Override
	public LIBORMarketModel getUnderlyingUndefaultableModel() {
		return _underlyingUndefaultableModel;
	}

}
