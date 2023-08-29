package info.quantlab.masterthesis.defaultablelibormodels;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
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

	
	RandomVariable getSpreadAtGivenTime(final double evalTime, final int liborIndex, final double time);
	
	RandomVariable getSpreadAtGivenTimeIndex(final double evalTime, final int liborIndex, final int timeIndex);
	
	LIBORMarketModel getUnderlyingNonDefaultableModel();
	
}
