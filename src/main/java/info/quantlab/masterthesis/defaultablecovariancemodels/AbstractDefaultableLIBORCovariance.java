package info.quantlab.masterthesis.defaultablecovariancemodels;

import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

public abstract class AbstractDefaultableLIBORCovariance extends AbstractLIBORCovarianceModelParametric implements DefaultableLIBORCovarianceModel {

	/**
	 * Default Serial UID
	 */
	private static final long serialVersionUID = 1L;
	
	private final LIBORCovarianceModel _nonDefaultableModel;
	
	public AbstractDefaultableLIBORCovariance(LIBORCovarianceModel nonDefaultableModel, TimeDiscretization timeDiscretization, TimeDiscretization liborPeriodDiscretization, int numberOfFactors) {
		super(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
		_nonDefaultableModel = nonDefaultableModel;
	}

	
	@Override
	public abstract AbstractDefaultableLIBORCovariance getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException;

	@Override
	public LIBORCovarianceModel getNonDefaultableCovarianceModel() {
		return _nonDefaultableModel;
	}

	@Override
	public abstract AbstractDefaultableLIBORCovariance getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newUndefaultableCovarianceModel);

	@Override
	public abstract double[] getParameterAsDouble();

	@Override
	public abstract int getNumberOfParameters();

	@Override
	public abstract AbstractDefaultableLIBORCovariance getCloneWithModifiedParameters(double[] parameters);

	@Override
	public AbstractDefaultableLIBORCovariance getCloneWithModifiedParameters(final RandomVariable[] parameters) {
		final double[] parameterAsDouble = new double[parameters.length];
		for(int i=0; i<parameterAsDouble.length; i++) {
			parameterAsDouble[i] = parameters[i].doubleValue();
		}
		return getCloneWithModifiedParameters(parameterAsDouble);
	}

	protected RandomVariable[] getZeroFactorLoading() {
		final RandomVariable zero = Scalar.of(0.0);
		RandomVariable[] fl = new RandomVariable[getNumberOfFactors()];
		Arrays.fill(fl, zero);
		return fl;
	}
	
	@Override
	public	RandomVariable[] getFactorLoading(final double time, final double component, final RandomVariable[] realizationAtTimeIndex) {
		
		if(time >= component) {
			return getZeroFactorLoading();
		}
		
		int componentIndex = getLiborPeriodDiscretization().getTimeIndex(component);
		if(componentIndex < 0) {
			// Get the closest libor Index that is smaller than component
			componentIndex = -componentIndex - 2;
		}
		
		// Always return the defaultable Version!
		return getFactorLoading(time, componentIndex + getLiborPeriodDiscretization().getNumberOfTimeSteps(), realizationAtTimeIndex);
	}

	@Override
	public RandomVariable[] getFactorLoadingOfSpread(double time, int component, RandomVariable[] realizationAtTimeIndex) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0)
			timeIndex = - timeIndex - 2;
		
		return getFactorLoadingOfSpread(timeIndex, component, realizationAtTimeIndex);
	}

	@Override
	public int getNumberOfLIBORPeriods() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}

	@Override
	public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException("This Method is not implemented for this model!");
	}
}
