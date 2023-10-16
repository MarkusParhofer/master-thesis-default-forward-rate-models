package info.quantlab.masterthesis.multilibormodels;

import info.quantlab.masterthesis.defaultableliborsimulation.MonteCarloProcessWithDependency;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

public interface DefaultableLIBORMarketModel extends LIBORMarketModel {

	/**
	 * Returns the underlying model for the undefaultable rate. E.g. Rates through AAA Government Bonds.
	 * @return the undefaultable model.
	 */
	LIBORMarketModel getUndefaultableLIBORModel();
	
	/**
	 * Gets the drift part of the undefaultable model. This will return an array of size <code>{@link LIBORModel#getNumberOfLibors()}</code>.
	 * @param process The discretization process generating the whole model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with MonteCarloProcess.getTimeDiscretization().
	 * @param realizationAtTimeIndex A vector of the realizations of the model.
	 * @param realizationPredictor The given realization at <code>timeIndex + 1</code> or null if no predictor is available.
	 * @return the drift vector of the defaultable model assosciated with the realizations given.
	 */
	RandomVariable[] getDriftOfUndefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);
	
	/**
	 * Gets the drift part of the defaultable model. This will return an array of size <code>{@link LIBORModel#getNumberOfLibors()}</code>. However the 
	 * realizations and realizationPredictor must still be an array of double the size, where the first half are the realizations of the undefaultable model, 
	 * and the second half are those of the defaultable model.
	 * @param process The discretization process generating the whole model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with MonteCarloProcess.getTimeDiscretization().
	 * @param realizationAtTimeIndex A vector of the realizations of the whole model. This means the size must be 
	 * <code>2 * {@link LIBORModel#getNumberOfLibors()}</code>, where the first half are the realizations of the undefaultable, 
	 * and the second half are those of the defaultable model.
	 * @param realizationPredictor The given realization at <code>timeIndex + 1</code> or null if no predictor is available.
	 * @return the drift vector of the defaultable model assosciated with the realizations given.
	 */
	RandomVariable[] getDriftOfDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);
	
	@Override
	DefaultableLIBORCovarianceModel getCovarianceModel();
	
	/**
	 * This method will only get the defaultable LIBOR Rate. 
	 */
	@Override
	RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;

	/**
	 * Return the forward rate of the undefaultable model at a given timeIndex and for a given liborIndex.
	 * 
	 * @param process The discretization process generating this model. The process provides call backs for 
	 * 	TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (associated with {@link MonteCarloProcess#getTimeDiscretization()}).
	 * @param liborIndex The forward rate index (associated with {@link LIBORModel#getLiborPeriodDiscretization()})
	 * @return The forward rate of the undefaultable model.
	 * @throws CalculationException Thrown if calculation failed.
	 */
	RandomVariable getUndefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;

	/**
	 * This method will get the defaultable Forward Rate. 
	 */
	@Override
	RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException;
	
	/**
	 * Gets the forward rate asscociated with the underlying undefaultable model.
	 * 
	 * @param process The discretization process generating this model.
	 * @param time The evaluation time.
	 * @param periodStart The period start of the forward rate.
	 * @param periodEnd The period end of the forward rate.
	 * @return The undefaultable forward rate.
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getUndefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException;
	
	/**
	 * The number of components. Here this method will return double the number of components of the undefaultable model.
	 */
	@Override
	int getNumberOfComponents();
	
	/**
	 * The number of LIBOR Periods.
	 */
	@Override
	int getNumberOfLibors();
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and undefaultable rate) at a given time index.
	 * @param process The simulation of the Model
	 * @param timeIndex The index of the time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcessWithDependency process, final int timeIndex, final int liborIndex) throws CalculationException;
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and undefaultable rate) at a given time index.
	 * @param process The simulation of the Model
	 * @param time The evaluation time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	default RandomVariable getLIBORSpreadAtGivenTime(MonteCarloProcessWithDependency process, final double time, final int liborIndex) throws CalculationException {
		int timeIndex = process.getTimeIndex(time);
		return getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborIndex);
	}

	/**
	 * Gets the spread from the defaultable and the undefaultable Forward rates.
	 * @param process The simulation of the Model
	 * @param time The evaluation time for which we get the spread.
	 * @param periodStart The period start of the forward rate for which to get the spread.
	 * @param periodEnd The period end of the forward rate for which to get the spread.
	 * @return The spread between the defaultable and the undefaultable forward rate.
	 * @throws CalculationException - Thrown if model fails to calculate the random variable.
	 */
	RandomVariable getSpread(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException ;
	
	@Override
	DefaultableLIBORMarketModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel newCovarianceModel);
	
	/**
	 * Returns a clone with a modified undefaultable state. This method must implement a way to also use a modified covariance model, 
	 * where the undefaultable covariance structure has changed.
	 * @param newUndefaultableModel The new undefaultable model used as basic for this defaultable model.
	 * @return A clone with a different undefaultable model.
	 */
	DefaultableLIBORMarketModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel);
	
	
}
