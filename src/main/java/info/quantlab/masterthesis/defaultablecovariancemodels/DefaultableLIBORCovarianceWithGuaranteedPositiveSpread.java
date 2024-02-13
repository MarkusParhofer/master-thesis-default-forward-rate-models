package info.quantlab.masterthesis.defaultablecovariancemodels;

import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * This class implements the factor Loadings of a defaultable LIBOR model.
 * It assumes the using LIBOR model (and hence possibly given realizations) to be a vector of the form:
 * <br>
 * <i>L<sup>v</sup><sub>i</sub> = </i> 
 * <br>
 * <ul>
 * 	<li>
 * 		<i>L<sub>i</sub></i> : The undefaultable model, if <i>i</i> &lt; <code>liborPeriodDiscretization.getNumberOfTimeSteps()</code>
 * 	</li>
 * 	<li>
 * 		<i>L<sup>d</sup><sub>i</sub></i> : The defaultable model, else
 * 	</li>
 * </ul>
 * 
 * @author Markus Parhofer
 * @version 1.1
 *
 */
public class DefaultableLIBORCovarianceWithGuaranteedPositiveSpread extends AbstractDefaultableLIBORCovFromFreeParam {

	/**
	 * Default One
	 */
	private static final long serialVersionUID = 1L;

	private final double[][] _freeParameterMatrix;
	
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(LIBORCovarianceModel undefaultableCovarianceModel, final double[][] freeParameterMatrix) {
		super(undefaultableCovarianceModel, freeParameterMatrix[0].length);
		_freeParameterMatrix = freeParameterMatrix;
		if(freeParameterMatrix.length != getLiborPeriodDiscretization().getNumberOfTimeSteps())
			throw new IllegalArgumentException("Free Parameter Matrix must have as many rows as there are LIBOR periods!");
	}
	
	public double[][] getFreeParameterMatrix() {
		return _freeParameterMatrix;
	}

	@Override
	public RandomVariable getFreeParameter(int timeIndex, int componentIndex, int factor, RandomVariable defRealization, RandomVariable nonDefRealization) {
		return Scalar.of(_freeParameterMatrix[componentIndex][factor]);
	}

	@Override
	public double[] getParameterAsDouble() {
		double[] allParams = new double[getNumberOfParameters()];
		
		// Parameters of free Parameter Matrix
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusIndex = row * _freeParameterMatrix[0].length;
            System.arraycopy(_freeParameterMatrix[row], 0, allParams, plusIndex, _freeParameterMatrix[0].length);
		}
		
		return allParams;
	}

	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread clone() {
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(getNonDefaultableCovarianceModel(), getFreeParameterMatrix());
	}
	
	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedParameters(double[] parameters) {
		
		double[][] newFreeParameterMatrix = new double[_freeParameterMatrix.length][_freeParameterMatrix[0].length];
		
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusLength = row * _freeParameterMatrix[0].length;
            System.arraycopy(parameters, plusLength, newFreeParameterMatrix[row], 0, _freeParameterMatrix[0].length);
		}
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(getNonDefaultableCovarianceModel(), newFreeParameterMatrix);
	}

	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNonDefaultableCovarianceModel) {
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(newNonDefaultableCovarianceModel, getFreeParameterMatrix());
	}

	@Override
	public int getNumberOfParameters() {
		return getFreeParameterMatrix().length * getFreeParameterMatrix()[0].length;
	}

}
