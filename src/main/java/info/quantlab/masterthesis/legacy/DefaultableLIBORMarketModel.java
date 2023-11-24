package info.quantlab.masterthesis.legacy;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.stochastic.RandomVariable;
/**
 * Interface for a Defaultable LIBOR Market Model.
 * Each Model must have an underlying LIBORMarketModel which takes the role of the undefaultable Market Situation.
 * 
 * @author Markus Parhofer
 */
public interface DefaultableLIBORMarketModel extends LIBORModel {
	
	
	/**
	 * Return the defaultable forward rate (LIBOR) covariance model.
	 *
	 * @return The covariance model.
	 */
	DefaultableLIBORCovarianceModel getCovarianceModel();

	/**
	 * Create a new object implementing DefaultableLIBORMarketModel, using the new covariance model.
	 *
	 * @param newCovarianceModel The new covariance model.
	 * @return A new object implementing LIBORMarketModel, using the new covariance model.
	 */
	DefaultableLIBORMarketModel getCloneWithModifiedCovarianceModel(DefaultableLIBORCovarianceModel newCovarianceModel);
	
	/**
	 * Create a new object implementing DefaultableLIBORMarketModel, using the new undefaultable model as underlying undefaultable Market Situation.
	 *
	 * @param newCovarianceModel The new covariance model.
	 * @return A new object implementing LIBORMarketModel, using the new covariance model.
	 */
	DefaultableLIBORMarketModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel);
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and undefaultable rate) at a given Time.
	 * @param process The simulation of the Model
	 * @param time The time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	default RandomVariable getLIBORSpreadAtGivenTime(MonteCarloProcessWithDependency process, final double time, final int liborIndex) throws CalculationException {
		int timeIndex = process.getTimeIndex(time);
		return getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborIndex);
	}
	
	/**
	 * Gets the Spread of the LIBOR rate (i.e. the difference between the defaultable and undefaultable rate) at a given time index.
	 * @param process The simulation of the Model
	 * @param timeIndex The time for which we get the spread
	 * @param liborIndex The index of the LIBOR rate to get
	 * @return The spread
	 * @throws CalculationException
	 */
	RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcessWithDependency process, final int timeIndex, final int liborIndex) throws CalculationException;
	
	/**
	 * Returns the underlying model for the undefaultable rate. E.g. Rates through AAA Government Bonds.
	 * @return the undefaultable model.
	 */
	LIBORMarketModel getUnderlyingNonDefaultableModel();
	
}
