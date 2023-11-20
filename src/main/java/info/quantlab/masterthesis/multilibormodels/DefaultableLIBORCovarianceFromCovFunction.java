package info.quantlab.masterthesis.multilibormodels;

import java.util.Arrays;
import java.util.Map;
import java.util.function.BinaryOperator;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public class DefaultableLIBORCovarianceFromCovFunction extends AbstractDefaultableLIBORCovariance implements DefaultableLIBORCovarianceModel {

	/**
	 * Default Serial UID
	 */
	private static final long serialVersionUID = 1L;
	
	private final double[][] _freeParameterMatrix;
	
	// TODO: Adjust cov Function such that if (S_10 -> 0 => _covFunc(L_1, L^d_1) -> 0, ... _covFunc(L_9, L^d_9) -> 0)
	// BiFunction<RandomVariable[], RandomVariable[], RandomVariable> _covFunction;
	
	private final BinaryOperator<RandomVariable>[] _covFunction;

	public DefaultableLIBORCovarianceFromCovFunction(BinaryOperator<RandomVariable>[] covFunction, double[][] freeParameterMatrix, LIBORCovarianceModel nonDefaultableModel) {
		super(nonDefaultableModel, nonDefaultableModel.getTimeDiscretization(), nonDefaultableModel.getLiborPeriodDiscretization(), freeParameterMatrix[0].length);
		_covFunction = covFunction;
		_freeParameterMatrix = freeParameterMatrix;
		if(_freeParameterMatrix.length != getNumberOfLIBORPeriods())
			throw new IllegalArgumentException("Free Parameters' number of rows must be equal to the number of LIBOR periods!");
		
		if(_freeParameterMatrix[0].length < getNonDefaultableCovarianceModel().getNumberOfFactors())
			throw new IllegalArgumentException("Free Parameters' number of columns must be equal to the number of factors in the new model and "
					+ "the number of factors must be greater or equal to the number of factors in the non defaultable model!");
		
		
		if(_covFunction.length != 1 && _covFunction.length != _freeParameterMatrix.length)
			throw new IllegalArgumentException("Length of CovarianceFunction must be either 1 or equal to the number of LIBOR periods!");
		// Check for |sigma(L_i, L^d_i)| = O(L^d_i - L_i)
		RandomVariable zero = new Scalar(0.0);
		for(int i=0; i < _covFunction.length; i++) {
			RandomVariable result = _covFunction[i].apply(zero, zero);
			if(result.doubleValue() != 0.0)
				throw new IllegalArgumentException("If Spread is zero, also the covariance Functions must be zero.");
		}
		
	}
	
	
	public double[][] getFreeParameterMatrix() {
		return _freeParameterMatrix;
	}

	public BinaryOperator<RandomVariable>[] getCovarianceFunctions() {
		return _covFunction;
	}

	public BinaryOperator<RandomVariable> getCovarianceFunction(int liborPeriodIndex) {
		if(_covFunction.length == 1)
			return _covFunction[0];
		
		return _covFunction[liborPeriodIndex];
	}
	
	
	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		final int liborPeriodIndex = component % getNumberOfLIBORPeriods();
		RandomVariable[] nonDefaultableFL = getFactorLoadingNonDefaultable(timeIndex, liborPeriodIndex, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));
		if(component == liborPeriodIndex)
			return nonDefaultableFL;
		else
			return getFactorLoadingDefaultable(timeIndex, liborPeriodIndex, nonDefaultableFL, realizationAtTimeIndex[component], realizationAtTimeIndex[liborPeriodIndex]);
	}

	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex, RandomVariable[] nonDefaultableRealization) {
		final int liborPeriodIndex = component % getNumberOfLIBORPeriods();
		RandomVariable[] nonDefaultableFL = getFactorLoadingNonDefaultable(timeIndex, liborPeriodIndex, nonDefaultableRealization);
		if(component == liborPeriodIndex)
			return nonDefaultableFL;
		else
			return getFactorLoadingDefaultable(timeIndex, liborPeriodIndex, nonDefaultableFL, realizationAtTimeIndex[liborPeriodIndex], nonDefaultableRealization[liborPeriodIndex]);
	}
	
	public RandomVariable[] getFactorLoadingNonDefaultable(int timeIndex, int liborPeriodIndex, RandomVariable[] nonDefaultableRealization) {
		RandomVariable[] flAllFactors = Arrays.copyOf(getNonDefaultableCovarianceModel().getFactorLoading(timeIndex, liborPeriodIndex, nonDefaultableRealization), getNumberOfFactors());
		RandomVariable zero = new Scalar(0.0);
		Arrays.fill(flAllFactors, getNonDefaultableCovarianceModel().getNumberOfFactors(), getNumberOfFactors(), zero);
		return flAllFactors;
	}
	
	public RandomVariable[] getFactorLoadingDefaultable(int timeIndex, int liborPeriodIndex, RandomVariable[] nonDefaultableFactorLoadings, RandomVariable defaultableRealization, RandomVariable nonDefaultableRealization) {
		if(getTimeDiscretization().getTime(timeIndex) >= getLiborPeriodDiscretization().getTime(liborPeriodIndex))
			return getZeroFactorLoading();
		
		final double liborPeriodLength = getLiborPeriodDiscretization().getTimeStep(liborPeriodIndex);
		RandomVariable[] factorLoading = new RandomVariable[getNumberOfFactors()];
		RandomVariable relationFactor = (new Scalar(1.0)).accrue(defaultableRealization, liborPeriodLength).discount(nonDefaultableRealization, liborPeriodLength);
		RandomVariable spread = defaultableRealization.sub(nonDefaultableRealization);
		for(int k = 0; k < getNumberOfFactors(); k++) {
			factorLoading[k] = new Scalar(getFreeParameterMatrix()[liborPeriodIndex][k]);
			if(k < getNonDefaultableCovarianceModel().getNumberOfFactors()) {
				factorLoading[k] = factorLoading[k].mult(getCovarianceFunction(liborPeriodIndex).apply(nonDefaultableRealization, defaultableRealization));
				factorLoading[k] = factorLoading[k].add(relationFactor.mult(nonDefaultableFactorLoadings[k]));
			} else {
				factorLoading[k] = factorLoading[k].mult(spread);
			}
		}
		return factorLoading;
	}

	@Override
	public double[] getParameterAsDouble() {
		double[][] freeParams = getFreeParameterMatrix();
		double[] returnVec = new double[getNumberOfParameters()];
		final int numberOfCols = freeParams[0].length;
		for(int i = 0; i < freeParams.length; i++) {
			for(int j = 0; j < numberOfCols; j++) {
				returnVec[i * numberOfCols + j] = freeParams[i][j];
			}
		}
		return returnVec;
	}

	@Override
	public int getNumberOfParameters() {
		// We can not give the cov functions as Parameter :(
		return getFreeParameterMatrix().length * getFreeParameterMatrix()[0].length;
	}
	
	@Override
	public boolean isSpreadModelLogNormal() {
		return false;
	}
	
	@Override
	public DefaultableLIBORCovarianceFromCovFunction getCloneWithModifiedData(Map<String, Object> dataModified)	throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public DefaultableLIBORCovarianceFromCovFunction getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNondefaultableCovarianceModel) {
		return new DefaultableLIBORCovarianceFromCovFunction(getCovarianceFunctions(), getFreeParameterMatrix(), newNondefaultableCovarianceModel);
	}
	
	@Override
	public DefaultableLIBORCovarianceFromCovFunction getCloneWithModifiedParameters(double[] parameters) {
		double[][] freeParams = new double[getFreeParameterMatrix().length][getFreeParameterMatrix()[0].length];
		final int numberOfCols = freeParams[0].length;
		for(int i = 0; i < freeParams.length; i++) {
			for(int j = 0; j < numberOfCols; j++) {
				freeParams[i][j] = parameters[i * numberOfCols + j];
			}
		}
		return new DefaultableLIBORCovarianceFromCovFunction(getCovarianceFunctions(), freeParams, getNonDefaultableCovarianceModel());
	}
	
	
	@Override
	public DefaultableLIBORCovarianceFromCovFunction clone() {
		return new DefaultableLIBORCovarianceFromCovFunction(getCovarianceFunctions(), getFreeParameterMatrix(), getNonDefaultableCovarianceModel());
	}
	
}
