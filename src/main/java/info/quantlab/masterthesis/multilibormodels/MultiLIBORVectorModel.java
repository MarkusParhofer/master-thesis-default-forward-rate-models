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
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.functional.FunctionsOnMCProcess;
import info.quantlab.masterthesis.process.MonteCarloProcessView;
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

	private transient MonteCarloProcess[] _processesOfAllModel = null;
	
	private final LIBORMarketModel _nonDefaultableLIBORModel;

	private int _numberOfFactors;
	
	public MultiLIBORVectorModel(DefaultableLIBORMarketModel[] defaultableLIBORModels, LIBORMarketModel nonDefaultableLIBORModel) {
		_defaultableLIBORModels = defaultableLIBORModels;
		_nonDefaultableLIBORModel = nonDefaultableLIBORModel;
		_numberOfFactors = _nonDefaultableLIBORModel.getNumberOfFactors();
		for(int index = 0; index < _defaultableLIBORModels.length; index++) {
			if(!_nonDefaultableLIBORModel.equals(_defaultableLIBORModels[index].getNonDefaultableLIBORModel())) {
				throw new IllegalArgumentException("NonDefaultable model is not equal to that of at least one defaultable model. Problem discovered at index " + index);
			}
			_numberOfFactors += _defaultableLIBORModels[index].getNumberOfFactors() - _nonDefaultableLIBORModel.getNumberOfFactors();
		}
	}

	/**
	 * Returns the nonDefaultable LIBORMarketModel that is the underlying for all defaultable Models.
	 * @return The underlying non-defaultable LIBOR model.
	 */
	public LIBORMarketModel getNonDefaultableModel() {
		return _nonDefaultableLIBORModel;
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
		return getNonDefaultableLIBOR(process, timeIndex, liborIndex);
	}

	public RandomVariable getNonDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex) throws CalculationException {
		return getNonDefaultableModel().getLIBOR(getNonDefaultableProcess(process), timeIndex, liborPeriodIndex);
	}
	
	public RandomVariable getDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborPeriodIndex, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getDefaultableLIBOR(getDefaultableProcess(process, liborModelIndex), timeIndex, liborPeriodIndex);
	}

	/**
	 * Gets the spread for a defaultable model for a LIBOR period (model primitive or almost): L^{d m}_i(t) - L_i(t).
	 * @param process 			The Process simulating this model.
	 * @param timeIndex			The evaluation time index associated with the process.
	 * @param liborPeriodIndex	The index of the LIBOR period, which the spread should correspond to.
	 * @param liborModelIndex	The index of the defaultable model to get the spread for.
	 * @return					The spread of the defaultable model for the specified LIBOR period.
	 * @throws CalculationException If the calculation failed.
	 */
	public RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcess process, int timeIndex, int liborPeriodIndex, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getLIBORSpreadAtGivenTimeIndex(getDefaultableProcess(process, liborModelIndex), timeIndex, liborPeriodIndex);
	}

	/**
	 * Gets the spread for a defaultable model to another defaultable model: L^{d m1}_i(t) - L^{d m2}_i(t).
	 * @param process 			The Process simulating this model.
	 * @param timeIndex			The evaluation time index associated with the process.
	 * @param liborPeriodIndex	The index of the LIBOR period, which the spread should correspond to.
	 * @param liborModelIndex1	The index of the defaultable model to get the spread for.
	 * @return					The spread of the defaultable model for the specified LIBOR period.
	 * @throws CalculationException If the calculation failed.
	 */
	public RandomVariable getLIBORSpreadAtGivenTimeIndexForDefModels(MonteCarloProcess process, int timeIndex, int liborPeriodIndex, int liborModelIndex1, int liborModelIndex2) throws CalculationException {
		final RandomVariable spread1 = getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborPeriodIndex, liborModelIndex1);
		final RandomVariable spread2 = getLIBORSpreadAtGivenTimeIndex(process, timeIndex, liborPeriodIndex, liborModelIndex2);
		return spread1.sub(spread2);
	}


	/**
	 * Gets the spread for a defaultable model.
	 * @param process 			The Process simulating this model.
	 * @param time 				The evaluation time.
	 * @param periodStart 		The start of the period to get the spread for.
	 * @param periodEnd			The end of the period to get the spread for.
	 * @param liborModelIndex	The index of the defaultable model to get the spread for.
	 * @return					The spread of the defaultable model for the specified period.
	 * @throws CalculationException If the calculation failed.
	 */
	public RandomVariable getSpread(MonteCarloProcess process, double time, double periodStart, double periodEnd, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getSpread(getDefaultableProcess(process, liborModelIndex), time, periodStart, periodEnd);
	}

	/**
	 * Gets the spread between two defaultable models. Note that this might be below zero! The specification is (L^d^m1(t;T_s, T_e) - L^d^m2(t;T_s, T_e)),
	 * where m1=liborModelIndex1, m2=liborModelIndex2, T_s=periodStart, T_e=periodEnd
	 * @param process 			The Process simulating this model
	 * @param time 				The evaluation time
	 * @param periodStart 		The start of the period to get the spread for
	 * @param periodEnd			The end of the period to get the spread for
	 * @param liborModelIndex1	The index of the first defaultable model to get the spread for
	 * @param liborModelIndex2	The index of the defaultable model to get the spread for
	 * @return					The spread of the defaultable model for the specified period
	 * @throws CalculationException If the calculation failed
	 */
	public RandomVariable getSpreadBetweenDefaultableModels(MonteCarloProcess process, double time, double periodStart, double periodEnd, int liborModelIndex1, int liborModelIndex2) throws CalculationException {
		final RandomVariable fr1 = getDefaultableModel(liborModelIndex1).getForwardRate(getDefaultableProcess(process, liborModelIndex1), time, periodStart, periodEnd);
		final RandomVariable fr2 = getDefaultableModel(liborModelIndex2).getForwardRate(getDefaultableProcess(process, liborModelIndex2), time, periodStart, periodEnd);
		return fr1.sub(fr2);
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
	 * This will return the non-defaultable forward rate. For the defaultable forward rate resolve to {@link #getDefaultableForwardRate(MonteCarloProcess, double, double, double, int)}
	 */
	@Override
	public RandomVariable getForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getNonDefaultableModel().getForwardRate(getNonDefaultableProcess(process), time, periodStart, periodEnd);
	}
	
	/**
	 * Return the defaultable forward rate for a specified model index.
	 * 
	 * @param process  The discretization process generating this model. The process provides call backs for TimeDiscretization 
	 * and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param time The evaluation time.
	 * @param periodStart The period start of the forward rate.
	 * @param periodEnd The period end of the forward rate.
	 * @param liborModelIndex The zero based index of the defaultable Model to use for the calculation of the forward rate
	 * @return The defaultable Forward rate using the model specifications of {@link #getDefaultableModel(int)}
	 * @throws CalculationException If calculation fails
	 */
	public RandomVariable getDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getDefaultableForwardRate(getDefaultableProcess(process, liborModelIndex), time, periodStart, periodEnd);
	}

	/**
	 * Gets the probability of survival until maturity time. Note that this is a conditional probability conditioned
	 * on not knowing the default state, while knowing the paths of L^d until maturity time.
	 * @param process The simulation process of the model.
	 * @param maturity The time until which to get the probability of survival for
	 * @param liborModelIndex The zero based index of the defaultable Model to use for the calculation of the survival probability
	 * @return The probability of survival
	 */
	public RandomVariable getSurvivalProbability(MonteCarloProcess process, final double maturity, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getSurvivalProbability(getDefaultableProcess(process, liborModelIndex), maturity);
	}
	
	@Override
	public AnalyticModel getAnalyticModel() {
		return null;
	}

	@Override
	public DiscountCurve getDiscountCurve() {
		return getNonDefaultableModel().getDiscountCurve();
	}

	/**
	 * Returns the Forward Rate Curve of the nonDefaultable model. For the forward curves of the
	 * defaultable models see {@link #getDefaultableForwardRateCurve(int)}
	 */
	@Override
	public ForwardCurve getForwardRateCurve() {
		return getNonDefaultableForwardRateCurve();
	}
	
	public ForwardCurve getNonDefaultableForwardRateCurve() {
		return getNonDefaultableModel().getForwardRateCurve();
	}

	@SuppressWarnings("unused")
	public ForwardCurve getDefaultableForwardRateCurve(int liborModelIndex) {
		return getDefaultableModel(liborModelIndex).getForwardRateCurve();
	}

	@Override
	public LocalDateTime getReferenceDate() {
		return getNonDefaultableModel().getReferenceDate();
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
			return getNonDefaultableModel().applyStateSpaceTransform(getNonDefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		else {
			final int modelIndex = getDefaultableModelIndex(componentIndex);
			RandomVariable result = getDefaultableModel(modelIndex).applyStateSpaceTransform(getDefaultableProcess(process, modelIndex), timeIndex, getNumberOfLiborPeriods() + getLiborPeriodIndexFromComponent(componentIndex), randomVariable);
			return result;
		}
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLiborPeriods())
			return getNonDefaultableModel().applyStateSpaceTransformInverse(getNonDefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		else {
			final int modelIndex = getDefaultableModelIndex(componentIndex);
			return getDefaultableModel(modelIndex).applyStateSpaceTransformInverse(getDefaultableProcess(process, modelIndex), timeIndex, getNumberOfLiborPeriods() + getLiborPeriodIndexFromComponent(componentIndex), randomVariable);
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
		return getNonDefaultableModel().getNumeraire(getNonDefaultableProcess(process), time);
	}

	public RandomVariable getDefaultableNumeraire(MonteCarloProcess process, double time, int liborModelIndex) throws CalculationException {
		return getDefaultableModel(liborModelIndex).getDefaultableNumeraire(getDefaultableProcess(process, liborModelIndex), time);
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		final RandomVariable[] nonDefaultableRealizations = Arrays.copyOf(realizationAtTimeIndex, getNumberOfLiborPeriods());
		final RandomVariable[] nonDefaultableRealizationPred = realizationPredictor != null? Arrays.copyOf(realizationPredictor, getNumberOfLiborPeriods()) : null;

		final int numberOfThreadsForProductValuation = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = null;
		if(numberOfThreadsForProductValuation == 0) {
			executor = Executors.newFixedThreadPool(numberOfThreadsForProductValuation);
		}
		final RandomVariable[] nonDefDrift = getNonDefaultableModel().getDrift(getNonDefaultableProcess(process), timeIndex, nonDefaultableRealizations, nonDefaultableRealizationPred);
		ArrayList<Future<RandomVariable[]>> driftFutures = new ArrayList<>(getNumberOfDefaultableModels());
		for(int i = 0; i < getNumberOfDefaultableModels(); i++) {
			final int modelIndex = i;
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);			
			
			final Callable<RandomVariable[]> worker = () -> {
                final int modelComponents = getDefaultableModel(modelIndex).getNumberOfComponents();
                RandomVariable[] realizationForModel = Arrays.copyOf(nonDefaultableRealizations, modelComponents);
                for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
                    realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
                }
                RandomVariable[] realizationPredForModel = null;
                if(nonDefaultableRealizationPred != null) {
                    realizationPredForModel = Arrays.copyOf(nonDefaultableRealizationPred, modelComponents);
                    for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
                        realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
                    }
                }

                return model.getDriftFast(getDefaultableProcess(process, modelIndex), timeIndex, realizationForModel, realizationPredForModel, nonDefDrift);
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
		
		RandomVariable[] driftVector = null;
		
		for(int modelIndex=0; modelIndex < getNumberOfDefaultableModels(); modelIndex++) {
			RandomVariable[] defaultableDrift;
			try {
				defaultableDrift = driftFutures.get(modelIndex).get();
			} catch (InterruptedException | ExecutionException e) {
				System.out.println("In Drift calculation: Multithreading did not work. Calling getDriftSlow(...)...");
				return getDriftSlow(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
			}
			if(modelIndex == 0) {
				driftVector = Arrays.copyOf(defaultableDrift, getNumberOfComponents());
			}
			else {
				for (int component = 0; component < getNumberOfLiborPeriods(); component++) {
					// First Indices of defaultable Drift are Non defaultable Drift
					driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = defaultableDrift[component + getNumberOfLiborPeriods()];
				}
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
		
		RandomVariable[] nonDefDrift = getNonDefaultableModel().getDrift(getNonDefaultableProcess(process), timeIndex, realizationForModel, realizationPredForModel);
		RandomVariable[] driftVector = Arrays.copyOf(nonDefDrift, getNumberOfComponents());
		
		for(int modelIndex = 0; modelIndex < getNumberOfDefaultableModels(); modelIndex++) {
            for(int comp = 0; comp < getNumberOfLiborPeriods(); comp++) {
				if(realizationForModel != null)
					realizationForModel[comp + getNumberOfLiborPeriods()] = realizationAtTimeIndex[getFirstComponentOfDefaultableModel(modelIndex) + comp];
				if(realizationPredForModel != null)
					realizationPredForModel[comp + getNumberOfLiborPeriods()] = realizationPredictor[getFirstComponentOfDefaultableModel(modelIndex) + comp];
			}
			final DefaultableLIBORMarketModel model = getDefaultableModel(modelIndex);
			
			RandomVariable[] modelDrift = model.getDriftFast(getDefaultableProcess(process, modelIndex), timeIndex, realizationForModel, realizationPredForModel, nonDefDrift);
			for(int component = 0; component < getNumberOfLiborPeriods(); component++) {
				driftVector[getFirstComponentOfDefaultableModel(modelIndex) + component] = modelDrift[component + getNumberOfLiborPeriods()];
			}
		}
		
		return driftVector;
	}
	
	@Override
	public int getNumberOfFactors() {
		return _numberOfFactors;
	}

	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		final int modelIndex = getDefaultableModelIndex(componentIndex);
		if(modelIndex < 0) {
			RandomVariable[] realizations = Arrays.copyOf(realizationAtTimeIndex, getNonDefaultableModel().getNumberOfComponents());
			return bringFactorLoadingsInRightPosition(getNonDefaultableModel().getFactorLoading(getNonDefaultableProcess(process), timeIndex, componentIndex, realizations), -1);
		}
		RandomVariable[] realizations = new RandomVariable[getDefaultableModel(modelIndex).getNumberOfComponents()];
		for(int i =0; i < getNumberOfLiborPeriods(); i++) {
			realizations[i] = realizationAtTimeIndex[i];
			final int FirstComponent = getFirstComponentOfDefaultableModel(modelIndex);
			realizations[i + getNumberOfLiborPeriods()] = realizationAtTimeIndex[i + FirstComponent];
		}
		return bringFactorLoadingsInRightPosition(getDefaultableModel(modelIndex).getFactorLoading(getDefaultableProcess(process, modelIndex), timeIndex, componentIndex, realizations), modelIndex);
	}

	private RandomVariable[] bringFactorLoadingsInRightPosition(RandomVariable[] flOfModel, int modelIndex) {
		final int factorsNonDefaultable = getNonDefaultableModel().getNumberOfFactors();

		RandomVariable[] allFactors = new RandomVariable[getNumberOfFactors()];
		int endInsertionIndex = modelIndex == 0? getDefaultableModel(modelIndex).getNumberOfFactors() : factorsNonDefaultable;
        if (endInsertionIndex >= 0) System.arraycopy(flOfModel, 0, allFactors, 0, endInsertionIndex);
		
		final RandomVariable zero = getRandomVariableForConstant(0.0);
		
		if(modelIndex < 1) {
			// If we are at the non-defaultable model or the index 0 defaultable model we have an early out
			Arrays.fill(allFactors, endInsertionIndex, getNumberOfFactors(), zero);
			return allFactors;
		}
		
		int startInsertionIndex = endInsertionIndex;
		int modelIterateIndex = 0;
		while(modelIterateIndex < modelIndex)
			startInsertionIndex += getDefaultableModel(modelIterateIndex++).getNumberOfFactors() - factorsNonDefaultable;
		
		Arrays.fill(allFactors, endInsertionIndex, startInsertionIndex, zero);
		endInsertionIndex = startInsertionIndex;
		for(int i = factorsNonDefaultable; i < getDefaultableModel(modelIndex).getNumberOfFactors(); i++) {
			allFactors[endInsertionIndex++] = flOfModel[i];
		}
		Arrays.fill(allFactors, endInsertionIndex, getNumberOfFactors(), zero);
		
		return allFactors;
	}
	
	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return getNonDefaultableModel().getRandomVariableForConstant(value);
	}

	@Override
	public MultiLIBORCovarianceVectorModel getCovarianceModel() {
		// This Model exists only for Calibration. We do not call its methods from this class!
		DefaultableLIBORCovarianceModel[] defCovarianceModels = new DefaultableLIBORCovarianceModel[_defaultableLIBORModels.length];
		for(int index = 0; index < _defaultableLIBORModels.length; index++) {
			defCovarianceModels[index] = _defaultableLIBORModels[index].getCovarianceModel();
		}
		return new MultiLIBORCovarianceVectorModel(defCovarianceModels, _nonDefaultableLIBORModel.getCovarianceModel(), true);
	}

	@Override
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		// TODO Auto-generated method stub
		return null;
	}

	private int getFirstComponentOfDefaultableModel(int liborModelIndex) {
		return getNumberOfLiborPeriods() * (liborModelIndex + 1);
	}

	private int getLastComponentOfDefaultableModel(int liborModelIndex) {
		return getNumberOfLiborPeriods() * (liborModelIndex + 2) - 1;
	}

	public MonteCarloProcess getNonDefaultableProcess(MonteCarloProcess process) {
		if(_processesOfAllModel == null) {
			_processesOfAllModel = new MonteCarloProcess[getNumberOfDefaultableModels() + 1];
		}
		if(_processesOfAllModel[0] == null) {
			_processesOfAllModel[0] = FunctionsOnMCProcess.getComponentReducedMCProcess(process, 0, getNumberOfLiborPeriods());
		}
		return _processesOfAllModel[0];
	}
	
	public MonteCarloProcess getDefaultableProcess(MonteCarloProcess process, int liborModelIndex) {
		if(_processesOfAllModel == null) {
			_processesOfAllModel = new MonteCarloProcess[getNumberOfDefaultableModels() + 1];
		}
		if(_processesOfAllModel[liborModelIndex + 1] == null) {
			int[] components = new int[getDefaultableModel(liborModelIndex).getNumberOfComponents()];
			int i=0;
			for(; i < getNonDefaultableModel().getNumberOfComponents(); i++) components[i] = i;
			for(int j=0; i < components.length; i++, j++)
				components[i] = getFirstComponentOfDefaultableModel(liborModelIndex) + j;
			_processesOfAllModel[liborModelIndex + 1] = new MonteCarloProcessView(process, components);
		}
		return _processesOfAllModel[liborModelIndex + 1];
	}

	@Override
	public LIBORModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MultiLIBORVectorModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel calibrationCovarianceModel) {
		if(calibrationCovarianceModel instanceof MultiLIBORCovarianceVectorModel castedModel) {
			LIBORMarketModel newNonDefaultableModel = getNonDefaultableModel();
			boolean nonDefaultableModelChanged = getNonDefaultableModel().getCovarianceModel() != castedModel.getNonDefaultableLiborCovarianceModel();
			if(nonDefaultableModelChanged)
				newNonDefaultableModel = getNonDefaultableModel().getCloneWithModifiedCovarianceModel(castedModel.getNonDefaultableLiborCovarianceModel());
			
			DefaultableLIBORMarketModel[] newDefaultableModels = new DefaultableLIBORMarketModel[_defaultableLIBORModels.length];
			for(int i=0; i< newDefaultableModels.length; i++) {
				newDefaultableModels[i] = getDefaultableModel(i);
				if(getDefaultableModel(i).getCovarianceModel() != castedModel.getDefaultableCovarianceModel(i))
					newDefaultableModels[i] = getDefaultableModel(i).getCloneWithModifiedCovarianceModel(castedModel.getDefaultableCovarianceModel(i));
				if(nonDefaultableModelChanged)
					newDefaultableModels[i] = newDefaultableModels[i].getCloneWithModifiedNonDefaultableModel(newNonDefaultableModel);
			}
			
			return new MultiLIBORVectorModel(newDefaultableModels, newNonDefaultableModel);
		}
		throw new IllegalArgumentException("Covariance model is not of type MultiLIBORCovarianceVectorModel.");
	}


}
