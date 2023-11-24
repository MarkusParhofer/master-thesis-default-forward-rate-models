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

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.functional.FunctionsOnMCProcess;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This is a class for the modelling of several Defaultable LIBOR Market models as well as one Non-Defaultable LIBOR 
 * market model.<p/>
 * <b>Caution</b>, this class assumes some basic characteristics of the underlying defaultable models:
 * <ul>
 * <li> All defaultable LIBOR Market Models are dependent on the same Non-Defaultable model.</li>
 * <li> The dependency on the non-defaultable Model is expressed such that for all SDE-methods/methods called by the 
 * MonteCarloProcess the first N components are used for modelling the non-defaultable Process (N=<code>getNumberOfLIBORPeriods()</code>)</li>
 * <li> There are no other components modelled than that of the non-defaultable process and that, which is used for the computation of the defaultable LIBORs</li>
 * <li> The statespace of the non-defaultable model must be the same in all models </li>
 * <li> The measure must be the same for all models including the non-defaultable one.</li>
 * </ul>
 * 
 * @author Markus Parhofer
 * @version 0.9
 */
public class MultiLIBORVectorModel implements LIBORMarketModel {

	private final DefaultableLIBORMarketModel[] _defaultableLIBORModels;
	
	private final LIBORMarketModel _undefaultableLIBORModel;
	
	public MultiLIBORVectorModel(DefaultableLIBORMarketModel[] defaultableLIBORModels, LIBORMarketModel undefaultableLIBORModel) {
		_defaultableLIBORModels = defaultableLIBORModels;
		_undefaultableLIBORModel = undefaultableLIBORModel;
		for(int index = 0; index < _defaultableLIBORModels.length; index++) {
			if(!_undefaultableLIBORModel.equals(_defaultableLIBORModels[index].getUndefaultableLIBORModel())) {
				throw new IllegalArgumentException("Undefaultable model is not equal to that of at least one defaultable model. Problem discovered at index " + index);
			}
		}
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
		return getUnDefaultableLIBOR(process, timeIndex, liborIndex);
	}

	public RandomVariable getUnDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex) throws CalculationException {
		return getUndefaultableModel().getLIBOR(getUndefaultableProcess(process), timeIndex, liborPeriodIndex);
	}
	
	public RandomVariable getDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getDefaultableLIBOR(getDefaultableProcess(process, liborModelIndex), timeIndex, liborPeriodIndex);
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

	public int getDefaultableModelIndex(int componentIndex) {
		return Math.floorDiv(componentIndex, getNumberOfLiborPeriods()) - 1;
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

	/**
	 * Gets the probability of survival until maturity time, given default has not yet happened at evaluation time. Note that this 
	 * is a conditional probability conditioned on F<sub>t</sub> <b>and</b> on {&tau; &gt; t}
	 * @param process The simulation process of the model.
	 * @param evaluationTime The evaluation time of the probability.
	 * @param maturity The time until which to get the probability of survival for
	 * @param liborModelIndex The zero based index of the defaultable Model to use for the calculation of the survival probability
	 * @return The probability of survival
	 */
	public RandomVariable getSurvivalProbability(MonteCarloProcess process, final double evaluationTime, final double maturity, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getSurvivalProbability(getDefaultableProcess(process, liborModelIndex), evaluationTime, maturity);
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
		return getNonDefaultableForwardRateCurve();
	}
	
	public ForwardCurve getNonDefaultableForwardRateCurve() {
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
		return getNumberOfLiborPeriods() * (getNumberOfDefaultableModels() + 1);
	}

	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLiborPeriods())
			return getUndefaultableModel().applyStateSpaceTransform(getUndefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		else {
			final int modelIndex = getDefaultableModelIndex(componentIndex);
			return getDefaultableModel(modelIndex).applyStateSpaceTransform(getDefaultableProcess(process, modelIndex), timeIndex, getLiborPeriodIndexFromComponent(componentIndex), randomVariable);
		}
	}

	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		// Initialize the Array to be non defaultable Initial states as well as the first Defaultable model Initial states:
		RandomVariable[] initialStates = Arrays.copyOf(getDefaultableModel(0).getInitialState(getDefaultableProcess(process, 0)), getNumberOfComponents());
		
		int indexInitialStates = getDefaultableModel(0).getNumberOfComponents();
		for(int i = 1; i < getNumberOfDefaultableModels(); i++) {
			RandomVariable[] modelInitialStates = getDefaultableModel(i).getInitialState(getDefaultableProcess(process, i));
			for(int k = 0; k < getNumberOfLiborPeriods(); k++) {
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
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);			
			
			final Callable<RandomVariable[]> worker = new  Callable<RandomVariable[]>() {
				@Override
				public RandomVariable[] call() {
					final int modelComponents = getDefaultableModel(modelIndex).getNumberOfComponents();
					RandomVariable[] realizationForModel = null;
					if(undefaultableRealizations != null) {
						realizationForModel = Arrays.copyOf(undefaultableRealizations, modelComponents);
						for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
							realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
						}
					}
					RandomVariable[] realizationPredForModel = null;
					if(undefaultableRealizationPred != null) {
						realizationPredForModel = Arrays.copyOf(undefaultableRealizationPred, modelComponents);
						for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
							realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
						}
					}
					
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
				// First Indices of defaultable Drift are Non defaultable Drift
				driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = defaultableDrift[component + getNumberOfLiborPeriods()];
			}
		}
		return driftVector;
	}

	
	public RandomVariable[] getDriftSlow(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		RandomVariable[] realizationForModel = null;
		if(realizationAtTimeIndex != null)
			realizationForModel = Arrays.copyOf(realizationAtTimeIndex, getDefaultableModel(0).getNumberOfComponents());
		
		RandomVariable[] realizationPredForModel = null;
		if(realizationPredictor != null)
			realizationPredForModel = Arrays.copyOf(realizationPredictor,  getDefaultableModel(0).getNumberOfComponents());
		
		RandomVariable[] firstDrift = getDefaultableModel(0).getDrift(getDefaultableProcess(process, 0), timeIndex, realizationForModel, realizationPredForModel);
		RandomVariable[] driftVector = Arrays.copyOf(firstDrift, getNumberOfComponents());
		
		for(int i = 1; i < getNumberOfDefaultableModels(); i++) {
			final int modelIndex = i;
			for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
				if(realizationForModel != null)
					realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
				if(realizationPredForModel != null)
					realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
			}
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);
			
			RandomVariable[] modelDrift = model.getDriftOfDefaultableModel(getDefaultableProcess(process, modelIndex), timeIndex, realizationForModel, realizationPredForModel);
			for(int component = 0; component < getNumberOfLiborPeriods(); component++) {
				driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = modelDrift[component + getNumberOfLiborPeriods()];
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
		final int modelIndex = getDefaultableModelIndex(componentIndex);
		if(modelIndex < 0) {
			RandomVariable[] realizations = Arrays.copyOf(realizationAtTimeIndex, getUndefaultableModel().getNumberOfComponents());
			return bringFactorLoadingsInRightPosition(getUndefaultableModel().getFactorLoading(getUndefaultableProcess(process), timeIndex, componentIndex, realizations), -1);
		}
		RandomVariable[] realizations = new RandomVariable[getDefaultableModel(modelIndex).getNumberOfComponents()];
		for(int i =0; i < getNumberOfLiborPeriods(); i++) {
			realizations[i] = realizationAtTimeIndex[i];
			realizations[i + getNumberOfLiborPeriods()] = realizationAtTimeIndex[i + getFirstComponentOfDefaultableModel(modelIndex)];
		}
		return bringFactorLoadingsInRightPosition(getDefaultableModel(modelIndex).getFactorLoading(getDefaultableProcess(process, modelIndex), timeIndex, componentIndex, realizations), modelIndex);
	}

	private RandomVariable[] bringFactorLoadingsInRightPosition(RandomVariable[] flOfModel, int modelIndex) {
		final int factorsUndefaultable = getUndefaultableModel().getNumberOfFactors();
		final int factorsModel = getDefaultableModel(modelIndex).getNumberOfFactors();
		
		RandomVariable[] allFactors = new RandomVariable[getNumberOfFactors()];
		int endInsertionIndex = modelIndex == 0? factorsModel : factorsUndefaultable;
		for(int i =0; i < endInsertionIndex; i++) {
			allFactors[i] = flOfModel[i];
		}
		
		final RandomVariable zero = getRandomVariableForConstant(0.0);
		
		if(modelIndex < 1) {
			// If we are at the non-defaultable model or the index 0 defaultable model we have an early out
			Arrays.fill(allFactors, endInsertionIndex, getNumberOfFactors(), zero);
			return allFactors;
		}
		
		int startInsertionIndex = endInsertionIndex;
		int modelIterateIndex = 0;
		while(modelIterateIndex < modelIndex)
			startInsertionIndex += getDefaultableModel(modelIterateIndex++).getNumberOfFactors() - factorsUndefaultable;
		
		Arrays.fill(allFactors, endInsertionIndex, startInsertionIndex, zero);
		endInsertionIndex = startInsertionIndex;
		for(int i = factorsUndefaultable; i < factorsModel; i++) {
			allFactors[endInsertionIndex++] = flOfModel[i];
		}
		Arrays.fill(allFactors, endInsertionIndex, getNumberOfFactors(), zero);
		
		return allFactors;
	}
	
	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return getUndefaultableModel().getRandomVariableForConstant(value);
	}

	@Override
	public MultiLIBORCovarianceVectorModel getCovarianceModel() {
		// This Model exists only for Calibration. We do not call its' methods from this class!
		DefaultableLIBORCovarianceModel[] defCovarianceModels = new DefaultableLIBORCovarianceModel[_defaultableLIBORModels.length];
		for(int index = 0; index < _defaultableLIBORModels.length; index++) {
			defCovarianceModels[index] = _defaultableLIBORModels[index].getCovarianceModel();
		}
		return new MultiLIBORCovarianceVectorModel(defCovarianceModels, _undefaultableLIBORModel.getCovarianceModel(), true);
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
		return getNumberOfLiborPeriods() * (liborModelIndex + 2) - 1;
	}
	
	private MonteCarloProcess getUndefaultableProcess(MonteCarloProcess process) {
		return FunctionsOnMCProcess.getComponentReducedMCProcess(process, 0, getNumberOfLiborPeriods());
	}
	
	private MonteCarloProcess getDefaultableProcess(MonteCarloProcess process, int liborModelIndex) {
		MonteCarloProcess defaultableProcess = FunctionsOnMCProcess.getComponentReducedMCProcess(
				process, 
				getFirstComponentOfDefaultableModel(liborModelIndex), 
				getLastComponentOfDefaultableModel(liborModelIndex));
		return FunctionsOnMCProcess.getCombinedMCProcess(getUndefaultableProcess(process), defaultableProcess);
	}

	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MultiLIBORVectorModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		if(calibrationCovarianceModel instanceof MultiLIBORCovarianceVectorModel castedModel) {
			LIBORMarketModel newUndefaultableModel = getUndefaultableModel();
			boolean undefaultableModelChanged = getUndefaultableModel().getCovarianceModel() != castedModel.getUndefaultableLiborCovarianceModel();
			if(undefaultableModelChanged)
				newUndefaultableModel = getUndefaultableModel().getCloneWithModifiedCovarianceModel(castedModel.getUndefaultableLiborCovarianceModel());
			
			DefaultableLIBORMarketModel[] newDefaultableModels = new DefaultableLIBORMarketModel[_defaultableLIBORModels.length];
			for(int i=0; i< newDefaultableModels.length; i++) {
				newDefaultableModels[i] = getDefaultableModel(i);
				if(getDefaultableModel(i).getCovarianceModel() != castedModel.getDefaultableCovarianceModel(i))
					newDefaultableModels[i] = getDefaultableModel(i).getCloneWithModifiedCovarianceModel(castedModel.getDefaultableCovarianceModel(i));
				if(undefaultableModelChanged)
					newDefaultableModels[i] = newDefaultableModels[i].getCloneWithModifiedUndefaultableModel(newUndefaultableModel);
			}
			
			return new MultiLIBORVectorModel(newDefaultableModels, newUndefaultableModel);
		}
		throw new IllegalArgumentException("Covariance model is not of type MultiLIBORCovarianceVectorModel.");
	}


}
