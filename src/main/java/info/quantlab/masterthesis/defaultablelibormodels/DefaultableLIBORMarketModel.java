package info.quantlab.masterthesis.defaultablelibormodels;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultableliborsimulation.MonteCarloProcessWithDependency;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
/**
 * Interface for a Defaultable LIBOR Market Model.
 * Each Model must have an underlying LIBORMarketModel which takes the role of the undefaultable Market Situation.
 * 
 * @author Markus Parhofer
 */
public interface DefaultableLIBORMarketModel extends LIBORModel {
	
	/* TODO: LIBORModel is maybe not suited as Interface for this Model: All methods concerning Factor Loadings need realizations from both 
	* defaultable and undefaultable Methods. Furthermore the MonteCarloProcess instances of the defaultable and the undefaultable model
	* should be forced to use the same stochastic driver.
	*/
	
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
	
	default RandomVariable getLIBORSpreadAtGivenTime(MonteCarloProcessWithDependency process, final double time, final int liborIndex) throws CalculationException {
		int timeIndex = process.getTimeIndex(time);
		return getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborIndex);
	}
	
	RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcessWithDependency process, final int timeIndex, final int liborIndex) throws CalculationException;
	
	LIBORMarketModel getUnderlyingNonDefaultableModel();
	
}
