package info.quantlab.masterthesis.multilibormodels;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

import info.quantlab.masterthesis.functional.Functional;
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
				throw new IllegalArgumentException("Undefaultable model is not equal to that of at least one defaultable model. Problem discovered at index " + index);
			}
			defCovarianceModels[index] = _defaultableLIBORModels[index].getCovarianceModel();
		}
		_covarianceModel = new MultiLIBORCovarianceVectorModel(defCovarianceModels, _undefaultableLIBORModel.getCovarianceModel(), true);
	}

	/**
	 * Returns the undefaultable LIBORMarketModel that is the underlying for all defaultable Models.
	 * @return The underlying non-defaultable LIBOR model.
	 */
	public LIBORMarketModel getUndefaultableModel() {
		return _undefaultableLIBORModel;
	}
	
	/**
	 * Returns the defaultable model with the specified zero based index in this multi LIBOR Model.
	 * @param liborModelIndex The zero based index of the Defaultable Model to get.
	 * @return The defaultable model with the specified zero based index in this multi LIBOR Model
	 */
	public DefaultableLIBORMarketModel getDefaultableModel(int liborModelIndex) {
		return _defaultableLIBORModels[liborModelIndex];
	}
	
	/**
	 * Returns the full array of defaultable models of this multi LIBOR Model.
	 * @return The full array of defaultable models of this multi LIBOR Model.
	 */
	public DefaultableLIBORMarketModel[] getArrayOfDefaultableModels() {
		return _defaultableLIBORModels;
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
		return getCovarianceModel().getLiborPeriodDiscretization();
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

	public int getModelIndex(int componentIndex) {
		return Math.floorDiv(componentIndex, getNumberOfLiborPeriods());
	}
	/**
	 * This will return the non-defaultable forward rate. For the defaultable forward rate resolve to {@link#getDefaultableForwardRate(process, time, periodStart, periodEnd, liborModelIndex)}
	 */
	@Override
	public RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getUndefaultableModel().getForwardRate(getUndefaultableProcess(process), time, periodStart, periodEnd);
	}
	
	/**
	 * Return the defaultable forwad rate for a specified model index.
	 * 
	 * @param process  The discretization process generating this model. The process provides call backs for TimeDiscretization 
	 * and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param time The evaluation time.
	 * @param periodStart The period start of the forward rate.
	 * @param periodEnd The period end of the forward rate.
	 * @param liborModelIndex The zero based index of the defaultable Model to use for the calculation of the forward rate
	 * @return The defaultable Forward rate using the model specifications of {@link #getDefaultableModel(liborModelIndex)}
	 * @throws CalculationException
	 */
	public RandomVariable getDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getDefaultableForwardRate(getDefaultableProcess(process, liborModelIndex), time, periodStart, periodEnd);
	}

	@Override
	public AnalyticModel getAnalyticModel() {
		return null;
	}

	@Override
	public DiscountCurve getDiscountCurve() {
		return getUndefaultableModel().getDiscountCurve();
	}

	/**
	 * Returns the Forward Rate Curve of the undefaultable model. For the forward curves of the defaultable models see {@link#getDefaultableForwardRateCurve(liborModelIndex)}
	 */
	@Override
	public ForwardCurve getForwardRateCurve() {
		return getUndefaultableModel().getForwardRateCurve();
	}
	
	public ForwardCurve getDefaultableForwardRateCurve(int liborModelIndex) {
		return getDefaultableModel(liborModelIndex).getForwardRateCurve();
	}

	@Override
	public LocalDateTime getReferenceDate() {
		return getUndefaultableModel().getReferenceDate();
	}

	public int getNumberOfDefaultableModels() {
		return getArrayOfDefaultableModels().length;
	}
	@Override
	public int getNumberOfComponents() {
		return (1 + getNumberOfLiborPeriods()) * getNumberOfDefaultableModels();
	}

	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLiborPeriods())
			return getUndefaultableModel().applyStateSpaceTransform(getUndefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		else {
			final int modelIndex = getModelIndex(componentIndex);
			return getDefaultableModel(modelIndex).applyStateSpaceTransform(getDefaultableProcess(process, modelIndex), timeIndex, getLiborPeriodIndexFromComponent(componentIndex), randomVariable);
		}
	}

	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		RandomVariable[] initialStates = Arrays.copyOf(getDefaultableModel(0).getInitialState(getDefaultableProcess(process, 0)), getNumberOfComponents());
		int indexInitialStates = getDefaultableModel(0).getNumberOfComponents();
		for(int i = 1; i < getNumberOfDefaultableModels(); i++) {
			RandomVariable[] modelInitialStates = getDefaultableModel(i).getInitialState(getDefaultableProcess(process, i));
			for(int k=0; k < getNumberOfLiborPeriods(); k++) {
				initialStates[indexInitialStates++] = modelInitialStates[k + getNumberOfLiborPeriods()];
			}
		}
		return initialStates;
	}

	@Override
	public RandomVariable getNumeraire(MonteCarloProcess process, double time) throws CalculationException {
		return getUndefaultableModel().getNumeraire(process, time);
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		final RandomVariable[] undefaultableRealizations = Arrays.copyOf(realizationAtTimeIndex, getNumberOfLiborPeriods());
		final RandomVariable[] undefaultableRealizationPred = Arrays.copyOf(realizationPredictor, getNumberOfLiborPeriods());
		

		final int numberOfThreadsForProductValuation = Runtime.getRuntime().availableProcessors();
		final ExecutorService executor = Executors.newFixedThreadPool(numberOfThreadsForProductValuation);

		
		ArrayList<Future<RandomVariable[]>> driftFutures = new ArrayList<>(getNumberOfDefaultableModels());
		for(int i = 0; i < getNumberOfDefaultableModels(); i++) {
			final int modelIndex = i;
			final int modelComponents = getDefaultableModel(modelIndex).getNumberOfComponents();
			final RandomVariable[] realizationForModel = Arrays.copyOf(undefaultableRealizations, modelComponents);
			final RandomVariable[] realizationPredForModel = Arrays.copyOf(undefaultableRealizationPred, modelComponents);
			for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
				realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
				realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
			}
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);
			
			final Callable<RandomVariable[]> worker = new  Callable<RandomVariable[]>() {
				@Override
				public RandomVariable[] call() {
					return model.getDriftOfDefaultableModel(getDefaultableProcess(process, modelIndex), timeIndex, realizationForModel, realizationPredForModel);
				}
			};
			if(executor != null) {
				final Future<RandomVariable[]> driftFuture = executor.submit(worker);
				driftFutures.add(modelIndex, driftFuture);
			}
			else {
				final FutureTask<RandomVariable[]> driftFutureTask = new FutureTask<>(worker);
				driftFutureTask.run();
				driftFutures.add(modelIndex, driftFutureTask);
			}
		}
		
		RandomVariable[] undefaultableDrift = getUndefaultableModel().getDrift(getUndefaultableProcess(process), timeIndex, undefaultableRealizations, undefaultableRealizationPred);
		RandomVariable[] driftVector = Arrays.copyOf(undefaultableDrift, getNumberOfComponents());
		
		for(int modelIndex=0; modelIndex < getNumberOfDefaultableModels(); modelIndex++) {
			RandomVariable[] defaultableDrift = null;
			try {
				defaultableDrift = driftFutures.get(modelIndex).get();
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				throw new UnsupportedOperationException("Multi Threading did not work! Try again later.");
			}
			for(int component = 0; component < getNumberOfLiborPeriods(); component++) {
				driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = defaultableDrift[component];
			}
		}
		return driftVector;
	}

	
	public RandomVariable[] getDriftSlow(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		
		final RandomVariable[] realizationForModel = Arrays.copyOf(realizationAtTimeIndex, getDefaultableModel(0).getNumberOfComponents());
		final RandomVariable[] realizationPredForModel = Arrays.copyOf(realizationPredictor,  getDefaultableModel(0).getNumberOfComponents());
		
		RandomVariable[] firstDrift = getDefaultableModel(0).getDrift(getDefaultableProcess(process, 0), timeIndex, realizationForModel, realizationPredForModel);
		RandomVariable[] driftVector = Arrays.copyOf(firstDrift, getNumberOfComponents());
		
		for(int i = 1; i < getNumberOfDefaultableModels(); i++) {
			final int modelIndex = i;
			for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
				realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
				realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
			}
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);
			
			RandomVariable[] modelDrift = model.getDriftOfDefaultableModel(getDefaultableProcess(process, modelIndex), timeIndex, realizationForModel, realizationPredForModel);
			for(int component = 0; component < getNumberOfLiborPeriods(); component++) {
				driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = modelDrift[component];
			}
		}
		
		return driftVector;
	}
	
	@Override
	public int getNumberOfFactors() {
		return getCovarianceModel().getNumberOfFactors();
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		return getCovarianceModel().getFactorLoading(process.getTime(timeIndex), componentIndex, realizationAtTimeIndex);
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
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		// TODO Auto-generated method stub
		return null;
	}

	public int getFirstComponentOfDefaultableModel(int liborModelIndex) {
		return getNumberOfLiborPeriods() * (liborModelIndex + 1);
	}
	
	public int getLastComponentOfDefaultableModel(int liborModelIndex) {
		return getNumberOfLiborPeriods() * (liborModelIndex + 2);
	}
	
	private MonteCarloProcess getUndefaultableProcess(MonteCarloProcess process) {
		return Functional.getComponentReducedMCProcess(process, 0, getNumberOfLiborPeriods());
	}
	
	private MonteCarloProcess getDefaultableProcess(MonteCarloProcess process, int liborModelIndex) {
		MonteCarloProcess defaultableProcess = Functional.getComponentReducedMCProcess(
				process, 
				getFirstComponentOfDefaultableModel(liborModelIndex), 
				getLastComponentOfDefaultableModel(liborModelIndex));
		return Functional.getCombinedMCProcess(getUndefaultableProcess(process), defaultableProcess);
	}

	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MultiLIBORVectorModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		if(!(calibrationCovarianceModel instanceof MultiLIBORCovarianceVectorModel)) {
			throw new IllegalArgumentException("Covariance model is not of type MultiLIBORCovarianceVectorModel.");
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


}
