package info.quantlab.masterthesis.multilibormodels;

import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
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
public class DefaultableLIBORCovarianceWithGuaranteedPositiveSpread extends AbstractDefaultableLIBORCovariance implements DefaultableLIBORCovarianceModel {

	/**
	 * Default One
	 */
	private static final long serialVersionUID = 1L;

	private final double[][] _freeParameterMatrix;
	
	public DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(LIBORCovarianceModel undefaultableCovarianceModel, final double[][] freeParameterMatrix) {
		super(undefaultableCovarianceModel, undefaultableCovarianceModel.getTimeDiscretization(), 
				undefaultableCovarianceModel.getLiborPeriodDiscretization(), 
				undefaultableCovarianceModel.getNumberOfFactors() + freeParameterMatrix[0].length);
		_freeParameterMatrix = freeParameterMatrix;
		if(freeParameterMatrix.length != getLiborPeriodDiscretization().getNumberOfTimeSteps())
			throw new IllegalArgumentException("Free Parameter Matrix must have as many rows as there are LIBOR periods!");
	}
	
	public double[][] getFreeParameterMatrix() {
		return _freeParameterMatrix;
	}
	
	@Override
	public double[] getParameterAsDouble() {
		double[] allParams = new double[getNumberOfParameters()];
		
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
	public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
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
			for(int col=0; col < _freeParameterMatrix[0].length; col++) {
				newFreeParameterMatrix[row][col] = parameters[plusLength + col];
			}
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
		final int liborPeriodIndex = getLIBORIndexFromComponent(component);
		RandomVariable[] undefaultableFactorLoading = getUndefaultableFactorLoading(timeIndex, liborPeriodIndex, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));
		
		if(component < getNumberOfLIBORPeriods()) {
			return undefaultableFactorLoading;
		}
		return getDefaultableFactorLoading(timeIndex, liborPeriodIndex, undefaultableFactorLoading, realizationAtTimeIndex[component], realizationAtTimeIndex[getLIBORIndexFromComponent(component)]);
	
	}
	
	private RandomVariable[] getDefaultableFactorLoading(int timeIndex, int liborPeriodIndex, RandomVariable[] undefaultableFactorLoading, RandomVariable componentRealizationDefaultable, RandomVariable componentRealizationUndefaultable) {
		if(getTimeDiscretization().getTime(timeIndex) >= getLiborPeriodDiscretization().getTime(liborPeriodIndex))
			return getZeroFactorLoading();

		final int undefaultableFactors = getNonDefaultableCovarianceModel().getNumberOfFactors();
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
		
		final RandomVariable[] undefaultableFactorLoading = getNonDefaultableCovarianceModel().getFactorLoading(timeIndex, liborPeriodIndex, undefaultableRealization);
		
		final RandomVariable zero = new Scalar(0.0);
		
		RandomVariable[] result = Arrays.copyOf(undefaultableFactorLoading, getNumberOfFactors());
		
		Arrays.fill(result, getNonDefaultableCovarianceModel().getNumberOfFactors(), getNumberOfFactors(), zero);

		return result;
	}
	
	@Override
	public RandomVariable[] getFactorLoadingOfSpread(int timeIndex, int liborPeriodIndex, RandomVariable[] realizationAtTimeIndex) {
		if(liborPeriodIndex >= getNumberOfLIBORPeriods())
			throw new ArrayIndexOutOfBoundsException("Spread model is a model of " + getNumberOfLIBORPeriods() + " Components. Index " + liborPeriodIndex + " out of Bounds for Spread Model");
		
		if(getTimeDiscretization().getTime(timeIndex) >= getLiborPeriodDiscretization().getTime(liborPeriodIndex))
			return getZeroFactorLoading();
		final RandomVariable[] allFactorLoadings = getUndefaultableFactorLoading(timeIndex, liborPeriodIndex, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));
		
		final int nonDefNumberOfFactors = getNonDefaultableCovarianceModel().getNumberOfFactors();
		for(int k = 0; k < nonDefNumberOfFactors; k++) {
			final double deltaT = getLiborPeriodDiscretization().getTimeStep(liborPeriodIndex);
			allFactorLoadings[k] = allFactorLoadings[k].mult(deltaT).div(realizationAtTimeIndex[liborPeriodIndex].mult(deltaT).add(1.0));
		}
		
		
		for(int k = nonDefNumberOfFactors; k < getNumberOfFactors(); k++) {
			allFactorLoadings[k] = Scalar.of(getFreeParameterMatrix()[liborPeriodIndex][k - nonDefNumberOfFactors]);
		}
		
		return allFactorLoadings;
	}
	
	@Override
	public boolean isSpreadModelLogNormal() {
		return true;
	}
	
	public int getLIBORIndexFromComponent(int component) {
		return component % getNumberOfLIBORPeriods();
	}

}
