package info.quantlab.masterthesis.multilibormodels;

import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public interface DefaultableLIBORCovarianceModel extends LIBORCovarianceModel {

	public LIBORCovarianceModel getUndefaultableCovarianceModel();
	
	public DefaultableLIBORCovarianceModel getCloneWithModifiedUndefaultableCovariance(LIBORCovarianceModel newUndefaultableCovarianceModel);
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model as double array. The parameters are usually free parameters
	 * which may be used in calibration.
	 *
	 * @return Parameter vector.
	 */
	public double[] getParameterAsDouble();
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model. The parameters are usually free parameters
	 * which may be used in calibration.
	 *
	 * @return Parameter vector.
	 */
	public default RandomVariable[]	getParameter() {
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
	public int getNumberOfParameter();
	
	/**
	 * Return an instance of this model using a new set of parameters.
	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
	 */
	public DefaultableLIBORCovarianceModel getCloneWithModifiedParameters(double[] parameters);

	/**
	 * Return an instance of this model using a new set of parameters.
	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
	 */
	public default DefaultableLIBORCovarianceModel getCloneWithModifiedParameters(final RandomVariable[] parameters) {
		final double[] parameterAsDouble = new double[parameters.length];
		for(int i=0; i<parameterAsDouble.length; i++) {
			parameterAsDouble[i] = parameters[i].doubleValue();
		}
		return getCloneWithModifiedParameters(parameterAsDouble);
	}

	/**
	 * Returns the Factor Loading for the specified time and component index. While the realizations are separated here, the component index is still
	 * the index for the whole model i.e. if <code>component</code> &lt; <code>getNumberOfLIBORPeriods()</code> this function will return 
	 * the factor loading of the non-defaultable model, otherwise the ones from the defaultable model.
	 * 
	 * @param timeIndex The time index at which the factor loading is requested.
	 * @param component The component index for which the factor loading is requested. If &lt; <code>getNumberOfLIBORPeriods()</code> this will give 
	 * the factor loading of the non-defaultable model.
	 * @param realizationAtTimeIndex The realization of the defaultable model (and only the defaultable model).
	 * @param undefaultableRealization The realization of the non defaultable model.
	 * @return The factor loading.
	 */
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex, RandomVariable[] undefaultableRealization);
	
	public int getNumberOfLIBORPeriods();
	
	public default RandomVariable[] getFactorLoadingOfSpread(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException("Spread dynamic is not implemented for this covariance Model!");
	}
	
	public default RandomVariable[] getFactorLoadingOfSpread(double time, int component, RandomVariable[] realizationAtTimeIndex) {
		int timeIndex = getTimeDiscretization().getTimeIndex(time);
		if(timeIndex < 0)
			timeIndex = - timeIndex - 2;
		
		return getFactorLoadingOfSpread(timeIndex, component, realizationAtTimeIndex);
	}
}
