package info.quantlab.masterthesis.defaultablelibormodels;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

public interface DefaultableLIBORMarketModel extends LIBORMarketModel, ProcessModel {

	/**
	 * Returns the underlying model for the non-defaultable rate. E.g. Rates through AAA Government Bonds.
	 * @return the non-defaultable model.
	 */
	LIBORMarketModel getNonDefaultableLIBORModel();

	/**
	 * Gets the drift part of the non-defaultable model. This will return an array of size <code>{@link LIBORModel#getNumberOfLibors()}</code>.
	 * @param process The discretization process generating the whole model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with MonteCarloProcess.getTimeDiscretization()).
	 * @param realizationAtTimeIndex A vector of the realizations of the model.
	 * @param realizationPredictor The given realization at <code>timeIndex + 1</code> or null if no predictor is available.
	 * @return the drift vector of the defaultable model assosciated with the realizations given.
	 */
	RandomVariable[] getDriftOfNonDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);

	/**
	 * Returns the drift vector of the model. Might give performance benefits, if the non defaultable
	 * drift vector is known. All input must be exactly as in {@link #getDrift(MonteCarloProcess, int, RandomVariable[], RandomVariable[])}.
	 * @param process The process simulating the SDE
	 * @param timeIndex The time index (associated with MonteCarloProcess.getTimeDiscretization()).
	 * @param realizationAtTimeIndex A vector of the realizations of the model.
	 * @param realizationPredictor The given realization at <code>timeIndex + 1</code> or null if no predictor is available.
	 * @param nonDefaultableDrift The drift vector of the underlying non defaultable model.
	 *                               Specifying this might give performance benefits.
	 * @return The drift vector.
	 */
	default RandomVariable[] getDriftFast(final MonteCarloProcess process, int timeIndex, final RandomVariable[] realizationAtTimeIndex, final RandomVariable[] realizationPredictor, final RandomVariable[] nonDefaultableDrift) {
		return getDrift(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
	}

	/**
	 * Gets the drift part of the defaultable model. This will return an array of size <code>{@link LIBORModel#getNumberOfLibors()}</code>. However the 
	 * realizations and realizationPredictor must still be an array of double the size, where the first half are the realizations of the non-defaultable model, 
	 * and the second half are those of the defaultable model.
	 * @param process The discretization process generating the whole model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with MonteCarloProcess.getTimeDiscretization().
	 * @param realizationAtTimeIndex A vector of the realizations of the whole model. This means the size must be 
	 * <code>2 * {@link LIBORModel#getNumberOfLibors()}</code>, where the first half are the realizations of the non-defaultable, 
	 * and the second half are those of the defaultable model.
	 * @param realizationPredictor The given realization at <code>timeIndex + 1</code> or null if no predictor is available.
	 * @return the drift vector of the defaultable model assosciated with the realizations given.
	 */
	RandomVariable[] getDriftOfDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);
	
	@Override
	DefaultableLIBORCovarianceModel getCovarianceModel();
	
	/**
	 * This method will only get either the non- or the defaultable LIBOR Rate depending on the liborIndex. 
	 */
	@Override
	default RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		return getDefaultableLIBOR(process, timeIndex, liborIndex);
	}

	/**
	 * Return the forward rate of the defaultable model at a given timeIndex and for a given liborIndex.
	 * 
	 * @param process The discretization process generating this model. The process provides call backs for 
	 * TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with {@link MonteCarloProcess#getTimeDiscretization()}).
	 * @param liborIndex The forward rate index (associated with {@link LIBORModel#getLiborPeriodDiscretization()})
	 * @return The forward rate of the defaultable model.
	 * @throws CalculationException - Thrown if calculation failed.
	 */
	RandomVariable getDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;
	
	/**
	 * Return the forward rate of the non-defaultable model at a given timeIndex and for a given liborIndex.
	 * 
	 * @param process The discretization process generating this model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with {@link MonteCarloProcess#getTimeDiscretization()}).
	 * @param liborIndex The forward rate index (associated with {@link LIBORModel#getLiborPeriodDiscretization()})
	 * @return The forward rate of the non-defaultable model.
	 * @throws CalculationException Thrown if calculation failed.
	 */
	RandomVariable getNonDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;

	/**
	 * Same as {@link #getDefaultableForwardRate(MonteCarloProcess, double, double, double)}.
	 */
	@Override
	default RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getDefaultableForwardRate(process, time, periodStart, periodEnd);
	}
	
	/**
	 * Gets the forward rate of the defaultable model. This Forward rate is <math>F<sub>t</sub></math>-measurable. Furthermore it is given 
	 * pre-default: hence default has not yet happened at time t.
	 * Hence it returns: E[L<sup>d</sup>(S,T) | F<sub>t</sub> &cap; {&tau; > t } ]
	 * 
	 * @param process The discretization process generating this model.
	 * @param time The evaluation time.
	 * @param periodStart The period start of the forward rate.
	 * @param periodEnd The period end of the forward rate.
	 * @return The defaultable forward rate.
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException;
	
	/**
	 * Gets the forward rate asscociated with the underlying non-defaultable model.
	 * 
	 * @param process The discretization process generating this model.
	 * @param time The evaluation time.
	 * @param periodStart The period start of the forward rate.
	 * @param periodEnd The period end of the forward rate.
	 * @return The non-defaultable forward rate.
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getNonDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException;
	
	/**
	 * The number of components. Here this method will return double the number of components of the non-defaultable model.
	 */
	@Override
	int getNumberOfComponents();
	
	/**
	 * The number of LIBOR Periods. Same as {@link #getNumberOfLIBORPeriods()}.
	 */
	@Override
	int getNumberOfLibors();

	/**
	 * The number of LIBOR Periods.
	 */
	int getNumberOfLIBORPeriods();
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and non-defaultable rate) at a given time index.
	 * @param process The simulation of the Model
	 * @param timeIndex The index of the time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcess process, final int timeIndex, final int liborIndex) throws CalculationException;
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and non-defaultable rate) at a given time index.
	 * @param process The simulation of the Model
	 * @param time The evaluation time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	default RandomVariable getLIBORSpreadAtGivenTime(MonteCarloProcess process, final double time, final int liborIndex) throws CalculationException {
		int timeIndex = process.getTimeIndex(time);
		if(timeIndex < 0) {
			timeIndex = - timeIndex - 2;
		}
		return getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborIndex);
	}

	/**
	 * Gets the defaultable zero coupon Bond (i.e. <i>P<sup>d</sup>(t;T)</i>) associated with this model. We assume non-default until valuation time.
	 * @param process The simulation process of the model.
	 * @param time The evaluation time of the bond with given maturity.
	 * @param maturity The maturity of the bond
	 * @return The price of a defaultable Zero Coupon Bond conditional on Pre-default (No default until <code>time</code>)
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getDefaultableBond(MonteCarloProcess process, final double time, final double maturity) throws CalculationException;
	
	/**
	 * Gets a defaultable Value associated with this product, that can be seen as a "defaultable Numeraire". Taking the Expectation of this
	 * @param process
	 * @param time
	 * @return
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getDefaultableNumeraire(final MonteCarloProcess process, final double time) throws CalculationException;
	
	/**
	 * Gets the probability of survival until maturity time. Note that this is a conditional probability conditioned
	 * on not knowing the default state, while knowing the paths of L^d until maturity time.
	 * @param process The simulation process of the model.
	 * @param maturity The time until which to get the probability of survival for
	 * @return The probability of survival
	 */
	RandomVariable getSurvivalProbability(MonteCarloProcess process, final double maturity) throws CalculationException;
	
	/**
	 * Gets the spread from the defaultable and the non-defaultable Forward rates.
	 * @param process The simulation of the Model
	 * @param time The evaluation time for which we get the spread.
	 * @param periodStart The period start of the forward rate for which to get the spread.
	 * @param periodEnd The period end of the forward rate for which to get the spread.
	 * @return The spread between the defaultable and the non-defaultable forward rate.
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getSpread(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException ;
	
	@Override
	DefaultableLIBORMarketModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel newCovarianceModel);
	
	/**
	 * Returns a clone with a modified non-defaultable state. This method must implement a way to also use a modified covariance model, 
	 * where the non-defaultable covariance structure has changed.
	 * @param newNonDefaultableModel The new non-defaultable model used as basic for this defaultable model.
	 * @return A clone with a different non-defaultable model.
	 */
	DefaultableLIBORMarketModel getCloneWithModifiedNonDefaultableModel(LIBORMarketModel newNonDefaultableModel);
	
	
}
