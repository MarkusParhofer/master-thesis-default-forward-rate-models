package info.quantlab.masterthesis.multilibormodels;

import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.DoubleStream;

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
			if(!_undefaultableCovarianceModel.equals(_defaultableCovarianceModels[i].getCovarianceStructureOfUndefaultableModel())) {
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
		
		double[] allParams = new double[0];
		
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
		final int defaultableModelIndex = Math.floorDiv(component, getNumberOfLIBORPeriods()) - 1;
		final int componentStartIndex = (defaultableModelIndex + 1) * getNumberOfLIBORPeriods();
		final int componentEndIndex = componentStartIndex + getNumberOfLIBORPeriods();
		final int modelComponent = component % getNumberOfLIBORPeriods();
		
		RandomVariable[] realizationOfUndefaultableModel = Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods());
		RandomVariable[] realizationOfDefaultableModel = Arrays.copyOfRange(realizationAtTimeIndex, componentStartIndex, componentEndIndex);
		
		
		RandomVariable[] factorLoadingLowFactors = _defaultableCovarianceModels[defaultableModelIndex].getFactorLoading(timeIndex, modelComponent, realizationOfDefaultableModel, realizationOfUndefaultableModel);
		RandomVariable zero = Scalar.of(0.0);
		RandomVariable[] factorLoading = Arrays.copyOf(factorLoadingLowFactors, getNumberOfFactors());
		int firstIndexZero = _defaultableCovarianceModels[defaultableModelIndex].getNumberOfFactors();
		
		if(defaultableModelIndex != 0) {
			int beginningOfExtraFactors = _undefaultableCovarianceModel.getNumberOfFactors();
			for(int i = 0; i < defaultableModelIndex; i++) {
				beginningOfExtraFactors += _defaultableCovarianceModels[i].getNumberOfFactors() - _undefaultableCovarianceModel.getNumberOfFactors();
			}
			
			Arrays.fill(factorLoading, _undefaultableCovarianceModel.getNumberOfFactors(), beginningOfExtraFactors, zero);
			
			for(int index = _undefaultableCovarianceModel.getNumberOfFactors(); index < _defaultableCovarianceModels[defaultableModelIndex].getNumberOfFactors(); index++) {
				factorLoading[beginningOfExtraFactors + index] = factorLoadingLowFactors[index];
			}
			firstIndexZero = beginningOfExtraFactors + _defaultableCovarianceModels[defaultableModelIndex].getNumberOfFactors() - _undefaultableCovarianceModel.getNumberOfFactors();
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
		int lastIndex = parameters.length;
		DefaultableLIBORCovarianceModel[] newCovarianceModels = new DefaultableLIBORCovarianceModel[getNumberOfDefaultableModels()];
		for(int modelIndex = getNumberOfDefaultableModels() - 1; modelIndex >= 0; modelIndex--) {
			newCovarianceModels[modelIndex] = _defaultableCovarianceModels[modelIndex].getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, lastIndex - _defaultableCovarianceModels[modelIndex].getNumberOfParameter(), lastIndex));
			lastIndex -= _defaultableCovarianceModels[modelIndex].getNumberOfParameter();
		}
		
		LIBORCovarianceModel newUndefaultableCovarianceModel = _undefaultableCovarianceModel;
		if(isUndefaultableModelCalibrateable()) {
			newUndefaultableCovarianceModel = ((AbstractLIBORCovarianceModelParametric)_undefaultableCovarianceModel).getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, 0, lastIndex));
			for(int modelIndex = 0; modelIndex < getNumberOfDefaultableModels(); modelIndex++) {
				newCovarianceModels[modelIndex] = newCovarianceModels[modelIndex].getCloneWithModifiedUndefaultableCovariance(newUndefaultableCovarianceModel);
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
