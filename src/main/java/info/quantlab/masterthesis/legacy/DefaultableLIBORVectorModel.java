package info.quantlab.masterthesis.legacy;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

public class DefaultableLIBORVectorModel extends AbstractProcessModel implements LIBORMarketModel {

	private final DefaultableLIBORCovarianceVectorModel _covarianceModel;
	
	private final ForwardCurve _initialRates;
	
	private final TimeDiscretization _liborPeriodDiscretization;
	
	
	public DefaultableLIBORVectorModel(DefaultableLIBORCovarianceVectorModel covarianceModel, ForwardCurve initialRates, TimeDiscretization liborPeriodDiscretization) {
		_covarianceModel = covarianceModel;
		_initialRates = initialRates;
		_liborPeriodDiscretization = liborPeriodDiscretization;
	}
	
	@Override
	public RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		return process.getProcessValue(timeIndex, liborIndex);
	}

	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return _liborPeriodDiscretization;
	}

	@Override
	public int getNumberOfLibors() {
		return getNumberOfComponents();
	}

	@Override
	public double getLiborPeriod(int timeIndex) {
		if(timeIndex >= _liborPeriodDiscretization.getNumberOfTimes() || timeIndex < 0) {
			throw new ArrayIndexOutOfBoundsException("Index for LIBOR period discretization out of bounds: " + timeIndex + ".");
		}
		return _liborPeriodDiscretization.getTime(timeIndex);
	}

	@Override
	public int getLiborPeriodIndex(double time) {
		return _liborPeriodDiscretization.getTimeIndex(time);
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
		return _initialRates;
	}

	@Override
	public int getNumberOfComponents() {
		return _liborPeriodDiscretization.getNumberOfTimeSteps() * 2;
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

	public RandomVariable[] getDriftOfUndefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		// TODO Auto-generated method stub
		return null;
	}
	
	
	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getNumberOfFactors() {
		return _covarianceModel.getNumberOfFactors();
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		return _covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
	}

	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return Scalar.of(value);
	}

	@Override
	public DefaultableLIBORCovarianceVectorModel getCovarianceModel() {
		return _covarianceModel;
	}

	@Override
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public DefaultableLIBORVectorModel clone() {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public LIBORMarketModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		// TODO Auto-generated method stub
		return null;
	}
}
