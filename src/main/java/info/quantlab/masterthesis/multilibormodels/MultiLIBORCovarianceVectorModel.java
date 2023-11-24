package info.quantlab.masterthesis.multilibormodels;

import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.DoubleStream;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public class MultiLIBORCovarianceVectorModel extends AbstractLIBORCovarianceModelParametric {

	/**
	 * Default ID
	 */
	private static final long serialVersionUID = 1L;

	private final DefaultableLIBORCovarianceModel[] _defaultableCovarianceModels;
	
	private final LIBORCovarianceModel _undefaultableCovarianceModel;
	
	private final boolean _calibrateUndefaultableModel;
	
	/**
	 * Constructs a LIBORCovarianceModel with a single defaultable Covariance model.
	 * @param defaultableCovarianceModel The (single) defaultable Covariance model.
	 * @param undefaultableCovarianceModel The covariance model of the underlying undefaultable LIBOR Model.
	 */
	public MultiLIBORCovarianceVectorModel(DefaultableLIBORCovarianceModel defaultableCovarianceModel, LIBORCovarianceModel undefaultableCovarianceModel) {
		this(new DefaultableLIBORCovarianceModel[] {defaultableCovarianceModel}, undefaultableCovarianceModel);
	}
	
	/**
	 * Constructs a LIBORCovarianceModel with a single defaultable Covariance model.
	 * @param defaultableCovarianceModel The (single) defaultable Covariance model.
	 * @param undefaultableCovarianceModel The covariance model of the underlying undefaultable LIBOR Model.
	 * @param calibrateUndefaultableModel flag determining if the undefaultable model is also to be calibrated if possible.
	 */
	public MultiLIBORCovarianceVectorModel(DefaultableLIBORCovarianceModel defaultableCovarianceModel, LIBORCovarianceModel undefaultableCovarianceModel, boolean calibrateUndefaultableModel) {
		this(new DefaultableLIBORCovarianceModel[] {defaultableCovarianceModel}, undefaultableCovarianceModel, calibrateUndefaultableModel);
	}
	
	/**
	 * Constructs a LIBORCovarianceModel with multiple defaultable Covariance models
	 * @param defaultableCovarianceModel The defaultable Covariance models as array.
	 * @param undefaultableCovarianceModel The covariance model of the underlying undefaultable LIBOR Model.
	 */
	public MultiLIBORCovarianceVectorModel(DefaultableLIBORCovarianceModel[] defaultableCovarianceModel, LIBORCovarianceModel undefaultableCovarianceModel) {
		this(defaultableCovarianceModel, undefaultableCovarianceModel, true);
	}
	
	/**
	 * Constructs a LIBORCovarianceModel with multiple defaultable Covariance models
	 * @param defaultableCovarianceModel The defaultable Covariance models as array.
	 * @param undefaultableCovarianceModel The covariance model of the underlying undefaultable LIBOR Model.
	 * @param calibrateUndefaultableModel flag determining if the undefaultable model is also to be calibrated if possible.
	 */
	public MultiLIBORCovarianceVectorModel(DefaultableLIBORCovarianceModel[] defaultableCovarianceModel, LIBORCovarianceModel undefaultableCovarianceModel, boolean calibrateUndefaultableModel) {
		super(undefaultableCovarianceModel.getTimeDiscretization(), 
				undefaultableCovarianceModel.getLiborPeriodDiscretization(), 
				undefaultableCovarianceModel.getNumberOfFactors() + 
				Arrays.stream(defaultableCovarianceModel).mapToInt(model -> model.getNumberOfFactors() - undefaultableCovarianceModel.getNumberOfFactors()).sum());
		
		_undefaultableCovarianceModel = undefaultableCovarianceModel;
		_defaultableCovarianceModels = defaultableCovarianceModel;
		_calibrateUndefaultableModel = calibrateUndefaultableModel;
		
		for(int i = 0; i < _defaultableCovarianceModels.length; i++) {
			if(!_undefaultableCovarianceModel.equals(_defaultableCovarianceModels[i].getNonDefaultableCovarianceModel())) {
				throw new InvalidParameterException("Undefaultable model is not equal to that of at least one defaultable model. Problem discovered at index " + i);
			}
		}
	}

	/**
	 * Checks if the undefaultable Model is calibrateable and if it should be calibrated.
	 * @return boolean value determining if the undefaultable model can and should be calibrated.	
	 */
	public boolean isUndefaultableModelCalibrateable() {
		return _calibrateUndefaultableModel && (_undefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric);
	}
	
	public DefaultableLIBORCovarianceModel[] getArrayOfDefaultableCovarianceModels() {
		return _defaultableCovarianceModels;
	}
	
	public DefaultableLIBORCovarianceModel getDefaultableCovarianceModel(int index) {
		return _defaultableCovarianceModels[index];
	}
	
	public LIBORCovarianceModel getUndefaultableLiborCovarianceModel() {
		return _undefaultableCovarianceModel;
	}
	
 	public int getNumberOfDefaultableModels() {
		return _defaultableCovarianceModels.length;
	}
	
	public int getNumberOfLIBORPeriods() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}
	
	@Override
	public double[] getParameterAsDouble() {
		
		double[] allParams = null;
		
		if(isUndefaultableModelCalibrateable()) {
			allParams = ((AbstractLIBORCovarianceModelParametric)_undefaultableCovarianceModel).getParameterAsDouble();
		}
		
		DoubleStream parameterStream = Arrays.stream(allParams);
		for(int i = 0; i < getNumberOfDefaultableModels(); i++) {
			parameterStream = DoubleStream.concat(parameterStream, Arrays.stream(_defaultableCovarianceModels[i].getParameterAsDouble()));
		}
		allParams = parameterStream.toArray();
		
		return allParams;
	}

	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
		final RandomVariable[] realizationOfUndefaultableModel = Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods());
		
		if(component < getNumberOfLIBORPeriods()) {
			final RandomVariable[] factorLoadingLowFactors = getUndefaultableLiborCovarianceModel().getFactorLoading(timeIndex, component, realizationOfUndefaultableModel);
			final RandomVariable[] factorLoading = Arrays.copyOf(factorLoadingLowFactors, getNumberOfFactors());
			Arrays.fill(factorLoading, getUndefaultableLiborCovarianceModel().getNumberOfFactors(), getNumberOfFactors(), Scalar.of(0.0));
			return factorLoading;
		}
		
		final int defaultableModelIndex = Math.floorDiv(component, getNumberOfLIBORPeriods()) - 1;
		final int componentStartIndex = (defaultableModelIndex + 1) * getNumberOfLIBORPeriods();
		final int componentEndIndex = componentStartIndex + getNumberOfLIBORPeriods();
		final int modelComponent = component % getNumberOfLIBORPeriods() + getNumberOfLIBORPeriods();
		final int timeIndexModified = getDefaultableCovarianceModel(defaultableModelIndex).getTimeDiscretization().
				getTimeIndex(getTimeDiscretization().getTime(timeIndex));
		
		final RandomVariable[] realizationOfDefaultableModel = Arrays.copyOfRange(realizationAtTimeIndex, componentStartIndex, componentEndIndex);
		
		final RandomVariable[] factorLoadingLowFactors = getDefaultableCovarianceModel(defaultableModelIndex).
				getFactorLoading(timeIndexModified, modelComponent, realizationOfDefaultableModel, realizationOfUndefaultableModel);
		final RandomVariable zero = Scalar.of(0.0);
		final RandomVariable[] factorLoading = Arrays.copyOf(factorLoadingLowFactors, getNumberOfFactors());
		int firstIndexZero = getDefaultableCovarianceModel(defaultableModelIndex).getNumberOfFactors();
		
		if(defaultableModelIndex != 0) {
			final int undefFactors = getUndefaultableLiborCovarianceModel().getNumberOfFactors();
			int beginningOfExtraFactors = undefFactors;
			for(int i = 0; i < defaultableModelIndex; i++) {
				beginningOfExtraFactors += getDefaultableCovarianceModel(i).getNumberOfFactors() - undefFactors;
			}
			
			Arrays.fill(factorLoading, undefFactors, beginningOfExtraFactors, zero);
			
			for(int index = undefFactors; index < firstIndexZero; index++) {
				factorLoading[beginningOfExtraFactors + index] = factorLoadingLowFactors[index];
			}
			firstIndexZero += beginningOfExtraFactors - undefFactors;
		}
		
		Arrays.fill(factorLoading, firstIndexZero, factorLoading.length, zero);
		return factorLoading;
	}

	@Override
	public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
	}

	@Override
	public MultiLIBORCovarianceVectorModel getCloneWithModifiedParameters(double[] parameters) {
		int lastIndexExclusive = parameters.length;
		final DefaultableLIBORCovarianceModel[] newCovarianceModels = new DefaultableLIBORCovarianceModel[getNumberOfDefaultableModels()];
		for(int modelIndex = getNumberOfDefaultableModels() - 1; modelIndex >= 0; modelIndex--) {
			final int firstIndex = lastIndexExclusive - _defaultableCovarianceModels[modelIndex].getNumberOfParameters();
			newCovarianceModels[modelIndex] = _defaultableCovarianceModels[modelIndex].getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, firstIndex, lastIndexExclusive));
			lastIndexExclusive = firstIndex;
		}
		
		LIBORCovarianceModel newUndefaultableCovarianceModel = _undefaultableCovarianceModel;
		if(isUndefaultableModelCalibrateable()) {
			newUndefaultableCovarianceModel = ((AbstractLIBORCovarianceModelParametric)_undefaultableCovarianceModel).getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, 0, lastIndexExclusive));
			for(int modelIndex = 0; modelIndex < getNumberOfDefaultableModels(); modelIndex++) {
				newCovarianceModels[modelIndex] = newCovarianceModels[modelIndex].getCloneWithModifiedNonDefaultableCovariance(newUndefaultableCovarianceModel);
			}
		}

		return new MultiLIBORCovarianceVectorModel(newCovarianceModels, newUndefaultableCovarianceModel, isUndefaultableModelCalibrateable());
	}

	@Override
	public MultiLIBORCovarianceVectorModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MultiLIBORCovarianceVectorModel clone() {
		return new MultiLIBORCovarianceVectorModel(_defaultableCovarianceModels, _undefaultableCovarianceModel, isUndefaultableModelCalibrateable());
	}


}
