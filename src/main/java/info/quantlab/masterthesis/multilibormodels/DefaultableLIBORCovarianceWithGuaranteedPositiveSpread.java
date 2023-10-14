package info.quantlab.masterthesis.multilibormodels;

import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
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
public class DefaultableLIBORCovarianceWithGuaranteedPositiveSpread extends AbstractLIBORCovarianceModelParametric implements DefaultableLIBORCovarianceModel {

	/**
	 * Default One
	 */
	private static final long serialVersionUID = 1L;

	private final LIBORCovarianceModel _undefaultableCovarianceModel;
	
	private final double[][] _freeParameterMatrix;
	
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(LIBORCovarianceModel undefaultableCovarianceModel, final double[][] freeParameterMatrix) {
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
		double[] allParams = new double[getNumberOfParameter()];
		
		// Parameters of free Parameter Matrix
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusIndex = row * _freeParameterMatrix[0].length;
			for(int col = 0; col < _freeParameterMatrix[0].length; col++) {
				allParams[plusIndex + col] = _freeParameterMatrix[row][col];
			}
		}
		
		return allParams;
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
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread clone() {
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(getUndefaultableCovarianceModel(), getFreeParameterMatrix());
	}
	
	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedParameters(double[] parameters) {
		
		double[][] newFreeParameterMatrix = new double[_freeParameterMatrix.length][_freeParameterMatrix[0].length];
		
		for(int row=0; row < _freeParameterMatrix.length; row++) {
			final int plusLength = row * _freeParameterMatrix[0].length;
			for(int col=0; col < _freeParameterMatrix[0].length; col++) {
				newFreeParameterMatrix[row][col] = parameters[plusLength + col];
			}
		}
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(getUndefaultableCovarianceModel(), newFreeParameterMatrix);
	}

	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedParameters(RandomVariable[] parameters) {
		return (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)super.getCloneWithModifiedParameters(parameters);
	}
	
	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread getCloneWithModifiedUndefaultableCovariance(LIBORCovarianceModel newUndefaultableCovarianceModel) {
		return new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(newUndefaultableCovarianceModel, getFreeParameterMatrix());
	}

	@Override
	public int getNumberOfParameter() {
		return getFreeParameterMatrix().length * getFreeParameterMatrix()[0].length;
	}
	
	public int getNumberOfLIBORPeriods() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}

	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex, RandomVariable[] undefaultableRealization) {

		RandomVariable[] undefaultableFactorLoading = getUndefaultableFactorLoading(timeIndex, component, undefaultableRealization);
		
		if(component < getNumberOfLIBORPeriods()) {
			return undefaultableFactorLoading;
		}
		return getDefaultableFactorLoading(timeIndex, component, undefaultableFactorLoading, realizationAtTimeIndex[getLIBORIndexFromComponent(component)], undefaultableRealization[getLIBORIndexFromComponent(component)]);
	}

	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		RandomVariable[] undefaultableFactorLoading = getUndefaultableFactorLoading(timeIndex, component, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));
		
		if(component < getNumberOfLIBORPeriods()) {
			return undefaultableFactorLoading;
		}
		return getDefaultableFactorLoading(timeIndex, component, undefaultableFactorLoading, realizationAtTimeIndex[component], realizationAtTimeIndex[getLIBORIndexFromComponent(component)]);
	
	}
	
	private RandomVariable[] getDefaultableFactorLoading(int timeIndex, int liborPeriodIndex, RandomVariable[] undefaultableFactorLoading, RandomVariable componentRealizationDefaultable, RandomVariable componentRealizationUndefaultable) {
		final int undefaultableFactors = getUndefaultableCovarianceModel().getNumberOfFactors();
		final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborPeriodIndex);
		
		RandomVariable[] factorLoading = new RandomVariable[getNumberOfFactors()];
		
		// Calculate from underlying undefaultable Model
		RandomVariable relationFactor = componentRealizationDefaultable.mult(periodLength).add(1.0).discount(componentRealizationUndefaultable, periodLength);
		for(int k = 0; k < undefaultableFactors; k++) {
			factorLoading[k] = relationFactor.mult(undefaultableFactorLoading[k]);
		}
		
		// Calculate from free parameters
		relationFactor = componentRealizationDefaultable.sub(componentRealizationUndefaultable);
		for(int k = undefaultableFactors; k < getNumberOfFactors(); k++) {
			factorLoading[k] = relationFactor.mult(getFreeParameterMatrix()[liborPeriodIndex][k - undefaultableFactors]);
		}
		return factorLoading;
	}
	
	private RandomVariable[] getUndefaultableFactorLoading(int timeIndex, int liborPeriodIndex, RandomVariable[] undefaultableRealization) {
		// Return undefaultable Factor Loadings with higher number of Factors. Set extra factors to zero.
		
		final RandomVariable[] undefaultableFactorLoading = getUndefaultableCovarianceModel().getFactorLoading(timeIndex, liborPeriodIndex, undefaultableRealization);
		
		final RandomVariable zero = new Scalar(0.0);
		
		RandomVariable[] result = Arrays.copyOf(undefaultableFactorLoading, getNumberOfFactors());
		
		Arrays.fill(result, getUndefaultableCovarianceModel().getNumberOfFactors(), getNumberOfFactors(), zero);

		return result;
	}
	
	public int getLIBORIndexFromComponent(int component) {
		return component % getNumberOfLIBORPeriods();
	}
}
