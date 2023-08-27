/**
 * 
 */
package info.quantlab.masterthesis.defaultablecovariancemodels;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * The abstract class for the implementation of a Covariance Model for defaultable LIBOR Market Models.
 * 
 * @author Markus Parhofer
 */
public abstract class AbstractDefaultableLIBORCovarianceModel implements LIBORCovarianceModel {
	
	private final LIBORMarketModel _underlyingUndefaultableModel;
	
	/**
	 * Initializes necessary members of the abstract model
	 * @param underlyingUndefaultableModel The undefaultable reference Model
	 */
	public AbstractDefaultableLIBORCovarianceModel(LIBORMarketModel underlyingUndefaultableModel) {
		_underlyingUndefaultableModel = underlyingUndefaultableModel;
	}

	@Override
	public	RandomVariable[] getFactorLoading(final double time, final double component, final RandomVariable[] realizationAtTimeIndex) {
		int componentIndex = getLiborPeriodDiscretization().getTimeIndex(component);
		if(componentIndex < 0) {
			componentIndex = -componentIndex - 2;
		}
		return getFactorLoading(time, componentIndex, realizationAtTimeIndex);
	}

	@Override
	public	RandomVariable[] getFactorLoading(final double time, final int component, final RandomVariable[] realizationAtTimeIndex) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = -timeIndex - 2;
		}
		return getFactorLoading(timeIndex, component, realizationAtTimeIndex);
	}

	@Override
	public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
	}

	@Override
	public RandomVariable getCovariance(double time, int component1, int component2, RandomVariable[] realizationAtTimeIndex) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = -timeIndex - 2;
		}
		return getCovariance(timeIndex, component1, component2, realizationAtTimeIndex);
	}

	@Override
	public RandomVariable getCovariance(int timeIndex, int component1, int component2, RandomVariable[] realizationAtTimeIndex) {
		final RandomVariable[] factorLoadingC1 = getFactorLoading(timeIndex, component1, realizationAtTimeIndex);
		final RandomVariable[] factorLoadingC2 = getFactorLoading(timeIndex, component2, realizationAtTimeIndex);
		
		RandomVariable covariance = factorLoadingC1[0].mult(factorLoadingC2[0]);
		for(int factor = 1; factor < getNumberOfFactors(); factor++) {
			covariance = covariance.add(factorLoadingC1[factor].mult(factorLoadingC2[factor]));
		}
		
		return covariance;
	}

	/**
	 * Gets the number of factors (i.e. the number of Factor Loadings per component).
	 * @return Number of Factors
	 */
	public abstract int getNumberOfFactors();
	
	/**
	 * Gets the underlying undefaultable Model. Hence a reference LIBOR Market Model.
	 * @return Undefaultable LIBOR Market Model
	 */
	public LIBORMarketModel getUnderlyingUndefaultableModel() {
		return _underlyingUndefaultableModel;
	}

}
