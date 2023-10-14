package info.quantlab.masterthesis.multilibormodels;

import java.security.InvalidParameterException;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

public class MultiLIBORVectorModel implements LIBORMarketModel {

	private final DefaultableLIBORMarketModel[] _defaultableLIBORModels;
	
	private final LIBORMarketModel _undefaultableLIBORModel;
	
	private final MultiLIBORCovarianceVectorModel _covarianceModel;
	
	public MultiLIBORVectorModel(DefaultableLIBORMarketModel[] defaultableLIBORModels, LIBORMarketModel undefaultableLIBORModel) {
		_defaultableLIBORModels = defaultableLIBORModels;
		_undefaultableLIBORModel = undefaultableLIBORModel;
		DefaultableLIBORCovarianceModel[] defCovarianceModels = new DefaultableLIBORCovarianceModel[_defaultableLIBORModels.length];
		for(int index = 0; index < _defaultableLIBORModels.length; index++) {
			if(!_undefaultableLIBORModel.equals(_defaultableLIBORModels[index].getUndefaultableLIBORModel())) {
				throw new InvalidParameterException("Undefaultable model is not equal to that of at least one defaultable model. Problem discovered at index " + index);
			}
			defCovarianceModels[index] = _defaultableLIBORModels[index].getCovarianceModel();
		}
		_covarianceModel = new MultiLIBORCovarianceVectorModel(defCovarianceModels, _undefaultableLIBORModel.getCovarianceModel(), true);
	}

	@Override
	public RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		return process.getProcessValue(timeIndex, liborIndex);
	}

	public RandomVariable getUnDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex) throws CalculationException {
		return process.getProcessValue(timeIndex, liborPeriodIndex);
	}
	
	public RandomVariable getDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex, int liborModelIndex) throws CalculationException {
		return process.getProcessValue(timeIndex, (liborModelIndex + 1) * getNumberOfLiborPeriods() + liborPeriodIndex);
	}
	
	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return _covarianceModel.getLiborPeriodDiscretization();
	}
	
	public int getNumberOfLiborPeriods() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}
	
	@Override
	public int getNumberOfLibors() {
		return getNumberOfComponents();
	}

	@Override
	public double getLiborPeriod(int timeIndex) {
		return getLiborPeriodDiscretization().getTime(timeIndex);
	}

	@Override
	public int getLiborPeriodIndex(double time) {
		return getLiborPeriodDiscretization().getTimeIndex(time);
	}

	public int getLiborPeriodIndexFromComponent(int componentIndex) {
		return componentIndex % getNumberOfLiborPeriods();
	}

	
	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public AnalyticModel getAnalyticModel() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DiscountCurve getDiscountCurve() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ForwardCurve getForwardRateCurve() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public LocalDateTime getReferenceDate() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getNumberOfComponents() {
		return (1 + getNumberOfLiborPeriods()) * _defaultableLIBORModels.length;
	}

	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable getNumeraire(MonteCarloProcess process, double time) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,
			RandomVariable[] realizationPredictor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getNumberOfFactors() {
		return _covarianceModel.getNumberOfFactors();
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		double time = process.getTime(timeIndex);
		return _covarianceModel.getFactorLoading(time, componentIndex, realizationAtTimeIndex);
	}

	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return Scalar.of(value);
	}

	@Override
	public LIBORCovarianceModel getCovarianceModel() {
		return _covarianceModel;
	}

	@Override
	public MultiLIBORVectorModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		if(!(calibrationCovarianceModel instanceof MultiLIBORCovarianceVectorModel)) {
			throw new InvalidParameterException("Covariance model is not of type MultiLIBORCovarianceVectorModel.");
		}
		MultiLIBORCovarianceVectorModel castedCovarianceModel = (MultiLIBORCovarianceVectorModel)calibrationCovarianceModel;
		LIBORMarketModel newUndefaultableModel = _undefaultableLIBORModel.getCloneWithModifiedCovarianceModel(castedCovarianceModel.getUndefaultableLiborCovarianceModel());
		DefaultableLIBORMarketModel[] newDefaultableModels = new DefaultableLIBORMarketModel[_defaultableLIBORModels.length];
		for(int i=0; i< newDefaultableModels.length; i++) {
			newDefaultableModels[i] = _defaultableLIBORModels[i].getCloneWithModifiedCovarianceModel(castedCovarianceModel.getDefaultableCovarianceModel(i));
			newDefaultableModels[i] = newDefaultableModels[i].getCloneWithModifiedUndefaultableModel(newUndefaultableModel);
		}
		return new MultiLIBORVectorModel(newDefaultableModels, newUndefaultableModel);
	}

	@Override
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		// TODO Auto-generated method stub
		return null;
	}

}
