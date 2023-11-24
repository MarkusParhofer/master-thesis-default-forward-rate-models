/**
 * 
 */
package info.quantlab.masterthesis.legacy;

import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
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
	private final DefaultableLIBORCovarianceModel _covarianceModel;
	private final Measure _measure;
	
	private final RandomVariableFactory	randomVariableFactory = new RandomVariableFromArrayFactory();
	
	/**
	 * 
	 */
	public DefaultableLIBORWithPositiveSpread(LIBORMarketModel underlyingLIBORModel, ForwardCurve initialCurve, DefaultableLIBORCovarianceModel covarianceModel, Measure measure) {
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

		// If time is beyond fixing, use the fixing time.
		time = Math.min(time, periodStart);
		int timeIndex           = process.getTimeIndex(time);
		// If time is not part of the discretization, use the nearest available point.
		if(timeIndex < 0) {
			timeIndex = -timeIndex-2;
			// For now this will always get the timeIndex lower than the time
		}

		// The forward rates are provided on fractional tenor discretization points using linear interpolation. See ISBN 0470047224.
		// TODO: Implement Interpolation
		if(periodStartIndex < 0 || periodEndIndex < 0) {
			throw new AssertionError("LIBOR requested outside libor discretization points and interpolation was not performed.");
		}

		// If this is a model primitive then return it
		if(periodStartIndex+1 == periodEndIndex) {
			return getLIBOR(process, timeIndex, periodStartIndex);
		}

		// The requested LIBOR is not a model primitive. We need to calculate it (slow!)
		RandomVariable accrualAccount = null;

		// Calculate the value of the forward bond
		for(int periodIndex = periodStartIndex; periodIndex<periodEndIndex; periodIndex++)
		{
			final double subPeriodLength = getLiborPeriod(periodIndex+1) - getLiborPeriod(periodIndex);
			final RandomVariable liborOverSubPeriod = getLIBOR(process, timeIndex, periodIndex);

			accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
		}

		final RandomVariable libor = accrualAccount.sub(1.0).div(periodEnd - periodStart);

		return libor;
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
		MonteCarloProcessWithDependency castedProcess;
		try {
			castedProcess = (MonteCarloProcessWithDependency)process;
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException("process must be of type MonteCarloProcessWithDependency to get the undefaultable process values.");
		}
		return _underlyingLIBORModel.getNumeraire(castedProcess.getDependencyProcess(), time);
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
				drift[componentIndex] = drift[componentIndex].addSumProduct(factorLoadingsSums, factorLoading);
				
				for(int factorIndex=0; factorIndex<factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
			}
		}
		else {
			throw new IllegalArgumentException("Drift not implemented for specified measure.");
		}

		return drift;
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		RandomVariable[] realizationsOfUndefaultableModel = new RandomVariable[getNumberOfComponents()];
		MonteCarloProcessWithDependency castedProcess;
		try {
			castedProcess = (MonteCarloProcessWithDependency)process;
		} catch (Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException("process must be MonteCarloProcessWithDependency to get the undefaultable ProcessValues");
		}
		for(int component = 0; component < getNumberOfComponents(); component++) {
			try {
				realizationsOfUndefaultableModel[component] = castedProcess.getDependencyProcessValue(timeIndex, componentIndex);
			} catch (CalculationException e) {
				e.printStackTrace();
			}
		}
		return _covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex, realizationsOfUndefaultableModel);
	}

	@Override
	public RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcessWithDependency process, int timeIndex, int liborIndex) throws CalculationException {
		return getLIBOR(process, timeIndex, liborIndex).sub(process.getDependencyProcessValue(timeIndex, liborIndex));
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
		throw new UnsupportedOperationException("Method not yet im of statespace transform not set");
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
	public DefaultableLIBORCovarianceModel getCovarianceModel() {
		return _covarianceModel;
	}

	@Override
	public DefaultableLIBORMarketModel getCloneWithModifiedCovarianceModel(DefaultableLIBORCovarianceModel newCovarianceModel) {
		if(newCovarianceModel.getUnderlyingUndefaultableModel() != _underlyingLIBORModel) {
			throw new IllegalArgumentException("New Covariance Model must still have the same underlying undefaultable Model");
		}
		return new DefaultableLIBORWithPositiveSpread(_underlyingLIBORModel, _initialCurve, newCovarianceModel, _measure);
	}


	@Override
	public DefaultableLIBORMarketModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel) {
		DefaultableLIBORCovarianceModel newCovarianceModel = _covarianceModel.getCloneWithModifiedUndefaultableModel(newUndefaultableModel);
		return new DefaultableLIBORWithPositiveSpread(newUndefaultableModel, _initialCurve, newCovarianceModel, _measure);
	}

	

}
