package info.quantlab.masterthesis.defaultablecovariancemodels;

import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
	/**
	 * This interface represents a covariance structure for a defaultable LIBOR model.
	 * It assumes the using LIBOR model (and hence possibly given realizations) to be a vector of the form:
	 * <br>
 	 * <i>L<sup>v</sup><sub>i</sub> = </i>
	 * <br>
 	 * <ul>
 	 * 	<li>
 	 * 		<i>L<sub>i</sub></i> : The non-defaultable model, if <i>i</i> &lt; <code>liborPeriodDiscretization.getNumberOfTimeSteps()</code>
	 * 	</li>
	 * 	<li>
 	 * 		<i>L<sup>d</sup><sub>i</sub></i> : The defaultable model, else
	 * 	</li>
	 * </ul>
	 */
public interface DefaultableLIBORCovarianceModel extends LIBORCovarianceModel {
	
	LIBORCovarianceModel getNonDefaultableCovarianceModel();
	
	DefaultableLIBORCovarianceModel getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNonDefaultableCovarianceModel);
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model as double array. The parameters are usually free parameters
	 * which may be used in calibration.
	 *
	 * @return Parameter vector.
	 */
	double[] getParameterAsDouble();
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model. The parameters are usually free parameters
	 * which may be used in calibration.
	 *
	 * @return Parameter vector.
	 */
	default RandomVariable[]	getParameter() {
		final double[] parameterAsDouble = this.getParameterAsDouble();
		final RandomVariable[] parameter = new RandomVariable[parameterAsDouble.length];
		for(int i=0; i<parameter.length; i++) {
			parameter[i] = new Scalar(parameterAsDouble[i]);
		}
		return parameter;
	}
	
	/**
	 * Gets the number of Parameters.
	 */
	int getNumberOfParameters();
	
	/**
	 * Returns an instance of this model using a new set of parameters.
	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
	 */
	DefaultableLIBORCovarianceModel getCloneWithModifiedParameters(double[] parameters);

	/**
	 * Return an instance of this model using a new set of parameters.
	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
	 */
	DefaultableLIBORCovarianceModel getCloneWithModifiedParameters(final RandomVariable[] parameters);

	/**
	 * Returns the Factor Loading for the specified time and component index. While the realizations are separated here, the component index is still
	 * the index for the whole model i.e. if <code>component</code> &lt; <code>getNumberOfLIBORPeriods()</code> this function will return 
	 * the factor loading of the non-defaultable model, otherwise the ones from the defaultable model.
	 * 
	 * @param timeIndex The time index at which the factor loading is requested.
	 * @param component The component index for which the factor loading is requested. If &lt; <code>getNumberOfLIBORPeriods()</code> this will give 
	 * the factor loading of the non-defaultable model.
	 * @param realizationAtTimeIndex The realization of the defaultable model (and only the defaultable model).
	 * @param nonDefaultableRealization The realization of the non defaultable model.
	 * @return The factor loading.
	 */
	RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex, RandomVariable[] nonDefaultableRealization);
	
	int getNumberOfLIBORPeriods();

	/**
	 * Gets a factor-loading vector for simulating the spreads associated with the DefaultableLIBORMarketModel and the
	 * non defaultable model.
	 * @param timeIndex index of the time for which the fl are requested, associated with {@link #getTimeDiscretization()}.
	 * @param component component for which the fl are requested. If component &le; N this method returns the factor loadings of the non defaultable model.
	 * @param realizationAtTimeIndex A vector of realizations of the form (L_0, ..., L_N, S_0, ..., S_N), where L is the non defaultable forward rate and S is the spread.
	 * @return A vector of factor loadings to simulate the spread.
	 * @throws UnsupportedOperationException if this method is not implemented by the derived class.
	 */
	default RandomVariable[] getFactorLoadingOfSpread(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException("Method getFactorLoadingOfSpread not implemented by " + this.getClass());
	}


	default boolean isSpreadModelLogNormal() {
		return false;
	}
	
	RandomVariable[] getFactorLoadingOfSpread(double time, int component, RandomVariable[] realizationAtTimeIndex);
}
