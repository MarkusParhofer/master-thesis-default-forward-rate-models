/**
 * 
 */
package info.quantlab.masterthesis.legacy;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * Interface for covariance models providing a vector of (possibly stochastic) factor loadings for defaultable LIBOR Models.
 * 
 * The factor loading is the vector <i>f<sub>i</sub></i> such that the scalar product <br>
 * <i>f<sub>j</sub>f<sub>k</sub> = f<sub>j,1</sub>f<sub>k,1</sub> + ... + f<sub>j,m</sub>f<sub>k,m</sub></i> <br>
 * is the instantaneous covariance of the component <i>j</i> and <i>k</i>.
 *
 * Classes implementing this interface can be used as "plug ins" for {@link DefaultableLIBORWithPositiveSpread}.
 *
 * @author Markus Parhofer
 * @version 0.1
 */
public interface DefaultableLIBORCovarianceModel {
	/**
	 * Return the factor loading for a given time and a given component.
	 *
	 * The factor loading is the vector <i>f<sub>i</sub></i> such that the scalar product <br>
	 * <i>f<sub>j</sub>f<sub>k</sub> = f<sub>j,1</sub>f<sub>k,1</sub> + ... + f<sub>j,m</sub>f<sub>k,m</sub></i> <br>
	 * is the instantaneous covariance of the component <i>j</i> and <i>k</i>.
	 *
	 * @param time The time <i>t</i> at which factor loading is requested.
	 * @param component The component time (as a double associated with the fixing of the forward rate)  <i>T<sub>i</sub></i>.
	 * @param defaultableRealization The realization of the stochastic (defaultable) process.
	 * @param undefaultableRealization The realizations of the stochastic undefaultable process.
	 * @return The factor loading <i>f<sub>i</sub>(t)</i>.
	 */
	RandomVariable[] getFactorLoading(double time, double component, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization);
	
	/**
	 * Return the factor loading for a given time and a given component.
	 *
	 * The factor loading is the vector <i>f<sub>i</sub></i> such that the scalar product <br>
	 * <i>f<sub>j</sub>f<sub>k</sub> = f<sub>j,1</sub>f<sub>k,1</sub> + ... + f<sub>j,m</sub>f<sub>k,m</sub></i> <br>
	 * is the instantaneous covariance of the component <i>j</i> and <i>k</i>.
	 *
	 * @param time The time <i>t</i> at which factor loading is requested.
	 * @param componentIndex The component index.
	 * @param defaultableRealization The realization of the stochastic (defaultable) process.
	 * @param undefaultableRealization The realizations of the stochastic undefaultable process.
	 * @return The factor loading <i>f<sub>i</sub>(t)</i>.
	 */
	RandomVariable[] getFactorLoading(double time, int componentIndex, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization);
	
	/**
	 * Return the factor loading for a given time and a given component.
	 *
	 * The factor loading is the vector <i>f<sub>i</sub></i> such that the scalar product <br>
	 * <i>f<sub>j</sub>f<sub>k</sub> = f<sub>j,1</sub>f<sub>k,1</sub> + ... + f<sub>j,m</sub>f<sub>k,m</sub></i> <br>
	 * is the instantaneous covariance of the component <i>j</i> and <i>k</i>.
	 *
	 * @param timeIndex The time index <i>t</i><sub>i</sub> at which factor loading is requested.
	 * @param componentIndex The component index.
	 * @param defaultableRealization The realization of the stochastic (defaultable) process.
	 * @param undefaultableRealization The realizations of the stochastic undefaultable process.
	 * @return The factor loading <i>f<sub>i</sub>(t)</i>.
	 */
	RandomVariable[] getFactorLoading(int timeIndex, int componentIndex, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization);
	
	/**
	 * Returns the instantaneous covariance calculated from factor loadings.
	 *
	 * @param time The time <i>t</i> at which covariance is requested.
	 * @param component1 Index of component <i>i</i>.
	 * @param component2  Index of component <i>j</i>.
	 * @param defaultableRealization The realization of the stochastic (defaultable) process.
	 * @param undefaultableRealization The realizations of the stochastic undefaultable process.
	 * @return The instantaneous covariance between component <i>i</i> and  <i>j</i>.
	 */
	RandomVariable getCovariance(double time, int component1, int component2, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization);
	
	/**
	 * Returns the instantaneous covariance calculated from factor loadings.
	 *
	 * @param timeIndex The time index <i>t</i><sub>i</sub> at which covariance is requested.
	 * @param component1 Index of component <i>i</i>.
	 * @param component2  Index of component <i>j</i>.
	 * @param defaultableRealization The realization of the stochastic (defaultable) process.
	 * @param undefaultableRealization The realizations of the stochastic undefaultable process.
	 * @return The instantaneous covariance between component <i>i</i> and  <i>j</i>.
	 */
	RandomVariable getCovariance(int timeIndex, int component1, int component2, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization);

	/**
	 * Gets the underlying undefaultable Model. Hence a reference LIBOR Market Model.
	 * @return Undefaultable LIBOR Market Model
	 */
	LIBORMarketModel getUnderlyingUndefaultableModel();
	
	/**
	 * The simulation time discretization associated with this model.
	 *
	 * @return the timeDiscretizationFromArray
	 */
	TimeDiscretization getTimeDiscretization();

	/**
	 * The forward rate time discretization associated with this model (defines the components).
	 *
	 * @return the forward rate time discretization associated with this model.
	 */
	TimeDiscretization getLiborPeriodDiscretization();

	/**
	 * @return the numberOfFactors
	 */
	int getNumberOfFactors();

	/**
	 * Returns a clone of this model where the underlying undefaultable Market Situation has changed.
	 * @param newUndefaultableModel The new underlying Market Situation.
	 * @return clone with modified data.
	 */
	DefaultableLIBORCovarianceModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel);
	
	/**
	 * Returns a clone of this model where the specified properties have been modified.
	 *
	 * Note that there is no guarantee that a model reacts on a specification of a properties in the
	 * parameter map <code>dataModified</code>. If data is provided which is ignored by the model
	 * no exception may be thrown.
	 *
	 * Furthermore the structure of the covariance model has to match changed data.
	 * A change of the time discretizations may requires a change in the parameters
	 * but this function will just insert the new time discretization without
	 * changing the parameters. An exception may not be thrown.
	 *
	 * @param dataModified Key-value-map of parameters to modify.
	 * @return A clone of this model (or a new instance of this model if no parameter was modified).
	 * @throws CalculationException Thrown when the model could not be created.
	 */
	DefaultableLIBORCovarianceModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException;
}
