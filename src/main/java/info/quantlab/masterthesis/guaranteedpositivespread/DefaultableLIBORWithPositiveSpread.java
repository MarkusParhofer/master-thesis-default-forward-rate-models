/**
 * 
 */
package info.quantlab.masterthesis.guaranteedpositivespread;

import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel.Measure;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel.StateSpace;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelCalibrateable;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

/**
 * This class implements a defaultable LIBOR Market Model, that guarantees a positive spread. We also implement the ProcessModel, 
 * such that the generation of a MonteCarloProcess specific for this class can be provided by using EulerSchemeFromProcessModel
 * @author Markus Parhofer
 */
public class DefaultableLIBORWithPositiveSpread extends AbstractProcessModel implements DefaultableLIBORMarketModel {

	public enum Measure { SPOT, TERMINAL }
	
	
	private final LIBORMarketModel 	_underlyingLIBORModel;
	private final ForwardCurve 		_initialCurve;
	private final AbstractDefaultableLIBORCovarianceModel _covarianceModel;
	private final Measure _measure;
	
	private final RandomVariableFactory	randomVariableFactory = new RandomVariableFromArrayFactory();
	
	/**
	 * 
	 */
	public DefaultableLIBORWithPositiveSpread(LIBORMarketModel underlyingLIBORModel, ForwardCurve initialCurve, AbstractDefaultableLIBORCovarianceModel covarianceModel, Measure measure) {
		_underlyingLIBORModel = underlyingLIBORModel;
		_initialCurve = initialCurve;
		_covarianceModel = covarianceModel;
		_measure = measure;
	}

	@Override
	public RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex)
			throws CalculationException {
		return process.getProcessValue(timeIndex, liborIndex);
	}

	@Override
	public RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd)
			throws CalculationException {
		
		final int periodStartIndex    = getLiborPeriodIndex(periodStart);
		final int periodEndIndex      = getLiborPeriodIndex(periodEnd);
		
		time = Math.min(time, periodStart);
		int timeIndex = process.getTimeIndex(time);
		
		if(periodStartIndex == periodEndIndex - 1)
			return getLIBOR(process, timeIndex, periodStartIndex);
		
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		
		final int numberOfLIBORPeriods = getNumberOfLibors();
		final RandomVariable[] initialStateRandomVariable = new RandomVariable[numberOfLIBORPeriods];
		
		for(int componentIndex = 0; componentIndex < numberOfLIBORPeriods; componentIndex++) {
			final double rate = _initialCurve.getForward(null, _underlyingLIBORModel.getLiborPeriod(componentIndex), _underlyingLIBORModel.getLiborPeriodDiscretization().getTimeStep(componentIndex));
			initialStateRandomVariable[componentIndex] = getRandomVariableForConstant(rate);
		}

		return initialStateRandomVariable;
	}
	
	@Override
	public double getLiborPeriod(int timeIndex) {
		return _underlyingLIBORModel.getLiborPeriod(timeIndex);
	}
	
	@Override
	public int getLiborPeriodIndex(double time) {
		return _underlyingLIBORModel.getLiborPeriodIndex(time);
	}

	@Override
	public RandomVariable getNumeraire(MonteCarloProcess process, double time) throws CalculationException {
		// Same Numeraire as underlying non defaultable model
		return _underlyingLIBORModel.getNumeraire(process, time);
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		final double	time					= process.getTime(timeIndex);			// t - current simulation time
		int				firstForwardRateIndex	= this.getLiborPeriodIndex(time)+1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex<0) {
			firstForwardRateIndex = - firstForwardRateIndex - 1 + 1;
		}

		final RandomVariable zero = Scalar.of(0.0);

		// Allocate drift vector and initialize to zero (will be used to sum up drift components)
		final RandomVariable[]	drift = new RandomVariable[getNumberOfComponents()];
		for(int componentIndex = firstForwardRateIndex; componentIndex < getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = zero;
		}

		// Allocate array (for each k) for the sums of delta_{i}/(1+L_{i} \delta_i) f_{i,k} (+ for spot measure, - for terminal measure)
		final RandomVariable[]	factorLoadingsSums	= new RandomVariable[process.getNumberOfFactors()];		
		Arrays.fill(factorLoadingsSums, zero);

		if(_measure == Measure.SPOT) {
			// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
			for(int componentIndex=firstForwardRateIndex; componentIndex<getNumberOfComponents(); componentIndex++) {
				final double			periodLength	= getLiborPeriodDiscretization().getTimeStep(componentIndex);
				final RandomVariable	forwardRate		= realizationAtTimeIndex[componentIndex];
				RandomVariable			oneStepMeasureTransform = Scalar.of(periodLength).discount(forwardRate, periodLength);

				final RandomVariable[]	factorLoading   	= getFactorLoading(process, timeIndex, componentIndex, realizationAtTimeIndex);
				for(int factorIndex=0; factorIndex<factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
				drift[componentIndex] = drift[componentIndex].addSumProduct(factorLoadingsSums, factorLoading);
			}
		}
		else if(_measure == Measure.TERMINAL) {
			// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
			for(int componentIndex=getNumberOfComponents()-1; componentIndex>=firstForwardRateIndex; componentIndex--) {
				final double			periodLength	= getLiborPeriodDiscretization().getTimeStep(componentIndex);
				final RandomVariable	forwardRate		= realizationAtTimeIndex[componentIndex];
				RandomVariable			oneStepMeasureTransform = Scalar.of(-periodLength).discount(forwardRate, periodLength);

				final RandomVariable[]	factorLoading   	= getFactorLoading(process, timeIndex, componentIndex, realizationAtTimeIndex);
				
				for(int factorIndex=0; factorIndex<factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
				drift[componentIndex] = drift[componentIndex].addSumProduct(factorLoadingsSums, factorLoading);
			}
		}
		else {
			throw new IllegalArgumentException("Drift not implemented for specified measure.");
		}

		return drift;
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex,
			RandomVariable[] realizationAtTimeIndex) {
		return _covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
	}

	@Override
	public RandomVariable getSpreadAtGivenTime(double evalTime, int liborIndex, double time) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable getSpreadAtGivenTimeIndex(double evalTime, int liborIndex, int timeIndex) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public AnalyticModel getAnalyticModel() {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public DiscountCurve getDiscountCurve() {
		// Same Discount curve as Non-Defaultable model
		return _underlyingLIBORModel.getDiscountCurve();
	}
	
	@Override
	public ForwardCurve getForwardRateCurve() {
		return _initialCurve;
	}

	@Override
	public LIBORMarketModel getUnderlyingNonDefaultableModel() {
		return _underlyingLIBORModel;
	}
	
	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return _underlyingLIBORModel.getLiborPeriodDiscretization();
	}
	
	@Override
	public int getNumberOfLibors() {
		return _underlyingLIBORModel.getNumberOfLibors();
	}
	
	@Override
	public int getNumberOfComponents() {
		return getNumberOfLibors();
	}
	
	@Override
	public int getNumberOfFactors() {
		return _covarianceModel.getNumberOfFactors();
	}
	
	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex,
			RandomVariable randomVariable) {
		return randomVariable;
	}
	
	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return randomVariableFactory.createRandomVariable(value);
	}
	
	@Override
	public LocalDateTime getReferenceDate() {
		throw new UnsupportedOperationException();
	}

	@Override
	public LIBORCovarianceModel getCovarianceModel() {
		return _covarianceModel;
	}

	@Override
	public LIBORMarketModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		AbstractDefaultableLIBORCovarianceModel  castedCovarianceModel = null;
		try {
			castedCovarianceModel = (AbstractDefaultableLIBORCovarianceModel)calibrationCovarianceModel;
		}
		catch(final Exception e) {
			throw new ClassCastException("The new Covariance Model must be of type AbstractDefaultableLIBORCovarianceModel.");
		}
		return new DefaultableLIBORWithPositiveSpread(_underlyingLIBORModel, _initialCurve, castedCovarianceModel, _measure);
	}

	@Override
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		// TODO Auto-generated method stub
		return null;
	}
}
