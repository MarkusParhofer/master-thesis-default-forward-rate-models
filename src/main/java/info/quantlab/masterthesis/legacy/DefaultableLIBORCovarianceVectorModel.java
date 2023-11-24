package info.quantlab.masterthesis.legacy;

import java.util.Map;

import net.finmath.exception.CalculationException;import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * This class implements the factor Loadings of an undefaultable and a defaultable LIBOR model. 
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
 *
 */
public class DefaultableLIBORCovarianceVectorModel extends AbstractLIBORCovarianceModelParametric {

	/**
	 * Default One
	 */
	private static final long serialVersionUID = 1L;

	private final LIBORCovarianceModel _undefaultableCovarianceModel;
	
	private final double[][] _freeParameterMatrix;
	
	public DefaultableLIBORCovarianceVectorModel(LIBORCovarianceModel undefaultableCovarianceModel, final double[][] freeParameterMatrix) {
		super(undefaultableCovarianceModel.getTimeDiscretization(), 
				undefaultableCovarianceModel.getLiborPeriodDiscretization(), 
				undefaultableCovarianceModel.getNumberOfFactors() + freeParameterMatrix[0].length);
		_undefaultableCovarianceModel = undefaultableCovarianceModel;
		_freeParameterMatrix = freeParameterMatrix;
	}
	
	public LIBORCovarianceModel getUndefaultableCovarianceModel() {
		return _undefaultableCovarianceModel;
	}

	public double[][] getFreeParameterMatrix() {
		return _freeParameterMatrix;
	}
	
	@Override
	public double[] getParameterAsDouble() {
		double[] extraParams = null;
		int extraLength = 0;
		
		if(_undefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric) {
			extraParams = ((AbstractLIBORCovarianceModelParametric)getUndefaultableCovarianceModel()).getParameterAsDouble();
			extraLength = extraParams.length;
		}
		
		double[] allParams = new double[extraLength + _freeParameterMatrix.length * _freeParameterMatrix[0].length];
		
		// Parameters of undefaultable model
		for(int index = 0; index < extraLength; index++) {
			allParams[index] = extraParams[index];
		}
		
		// Parameters of free Parameter Matrix
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusIndex = extraLength + row * _freeParameterMatrix[0].length;
			for(int col = 0; col < _freeParameterMatrix[0].length; col++) {
				allParams[plusIndex + col] = _freeParameterMatrix[row][col];
			}
		}
		
		return allParams;
	}
	
	public boolean undefaultableModelIsCalibrateable() {
		return _undefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric;
	}
		
	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		final int componentNumberReal = getLiborPeriodDiscretization().getNumberOfTimeSteps();
		
		RandomVariable[] undefaultableRealizations = new RandomVariable[componentNumberReal];
		
		for(int index = 0; index < componentNumberReal; index++) {
			undefaultableRealizations[index] = realizationAtTimeIndex[index];
		}
		
		final RandomVariable[] undefaultableFactorLoading = getUndefaultableCovarianceModel().getFactorLoading(timeIndex, component, undefaultableRealizations);
		
		if(component < componentNumberReal) {
			
			// Return undefaultable Factor Loadings with higher number of Factors. Set extra factors to zero.
			
			final RandomVariable zero = new Scalar(0.0);
			
			RandomVariable[] result = new RandomVariable[getNumberOfFactors()];
			
			for(int index = 0; index < getNumberOfFactors(); index++) {
				result[index] = index < getUndefaultableCovarianceModel().getNumberOfFactors() ? undefaultableFactorLoading[index] : zero;
			}
			
			return result;
		}
		
		final RandomVariable componentRealizationUndefaultable = realizationAtTimeIndex[component - componentNumberReal];
		final RandomVariable componentRealizationDefaultable = realizationAtTimeIndex[component];
		
		final int undefaultableFactors = getUndefaultableCovarianceModel().getNumberOfFactors();
		
		final double periodLength = getLiborPeriodDiscretization().getTimeStep(component - componentNumberReal);
		
		RandomVariable[] factorLoading = new RandomVariable[getNumberOfFactors()];
		
		// Calculate from underlying undefaultable Model
		RandomVariable relationFactor = componentRealizationDefaultable.mult(periodLength).add(1.0).discount(componentRealizationUndefaultable, periodLength);
		for(int k = 0; k < undefaultableFactors; k++) {
			factorLoading[k] = relationFactor.mult(undefaultableFactorLoading[k]);
		}
		
		// Calculate from free parameters
		relationFactor = componentRealizationDefaultable.sub(componentRealizationUndefaultable);
		for(int k = undefaultableFactors; k < getNumberOfFactors(); k++) {
			factorLoading[k] = relationFactor.mult(getFreeParameterMatrix()[component - componentNumberReal][k - undefaultableFactors]);
		}
		return factorLoading;
	}

	@Override
	public	RandomVariable[] getFactorLoading(final double time, final double component, final RandomVariable[] realizationAtTimeIndex) {
		int componentIndex = getLiborPeriodDiscretization().getTimeIndex(component);
		if(componentIndex < 0) {
			componentIndex = -componentIndex - 2;
		}
		
		// Always return the defaultable Version!
		return getFactorLoading(time, componentIndex + getLiborPeriodDiscretization().getNumberOfTimeSteps(), realizationAtTimeIndex);
	}
	
	@Override
	public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public DefaultableLIBORCovarianceVectorModel clone() {
		return new DefaultableLIBORCovarianceVectorModel(getUndefaultableCovarianceModel(), getFreeParameterMatrix());
	}
	
	@Override
	public DefaultableLIBORCovarianceVectorModel getCloneWithModifiedParameters(double[] parameters) {
		LIBORCovarianceModel newUndefaultableCovariance = getUndefaultableCovarianceModel();
		int extraLength = 0;
		if(undefaultableModelIsCalibrateable()) {
			double[] newUndefParams = new double[((AbstractLIBORCovarianceModelParametric)getUndefaultableCovarianceModel()).getParameterAsDouble().length];
			extraLength = newUndefParams.length;
			for(int index = 0; index < newUndefParams.length; index++) {
				newUndefParams[index] = parameters[index];
			}
			newUndefaultableCovariance = ((AbstractLIBORCovarianceModelParametric)newUndefaultableCovariance).getCloneWithModifiedParameters(newUndefParams);
		}
		
		double[][] newFreeParameterMatrix = new double[_freeParameterMatrix.length][_freeParameterMatrix[0].length];
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusLength = extraLength + row * _freeParameterMatrix[0].length;
			for(int col=0; col < _freeParameterMatrix[0].length; col++) {
				newFreeParameterMatrix[row][col] = parameters[plusLength + col];
			}
		}
		return new DefaultableLIBORCovarianceVectorModel(newUndefaultableCovariance, newFreeParameterMatrix);
	}

	@Override
	public DefaultableLIBORCovarianceVectorModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		throw new UnsupportedOperationException();
	}

}
