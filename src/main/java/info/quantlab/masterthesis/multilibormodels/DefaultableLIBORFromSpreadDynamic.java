package info.quantlab.masterthesis.multilibormodels;

import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import info.quantlab.masterthesis.functional.Functional;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;



/**
 * This class models a Defaultable LIBOR Market Model on the basis of a non-defaultable Model. While the Model in itself is for getting 
 * the defaultable LIBOR rates, we model only the non-defaultable rates and the spread and get the defaultable rate, by adding the drift 
 * to the non-defaultable rates.
 * Consistent with this policy all methods from ProcessModel (getDrift(), getFactorLoadings(), getNumberOfComponents(), ... ) give values for 
 * the non-defaultable model and the Spread.
 * The Process used with this model (namely EulerSchemeFromProcessModel) will hence also only get the values corresponding to these Models.
 * 
 * @author Markus Parhofer
 *
 */
public class DefaultableLIBORFromSpreadDynamic  extends AbstractProcessModel implements DefaultableLIBORMarketModel, ProcessModel {

	public enum Measure					{ SPOT, TERMINAL }
	public enum StateSpace				{ NORMAL, @Deprecated LOGNORMAL } // No support for the LOGNORMAL version yet (drift must be adjusted)
	public enum InterpolationMethod		{ LINEAR, LOG_LINEAR_UNCORRECTED, LOG_LINEAR_CORRECTED }
	public enum HandleSimulationTime 	{ ROUND_DOWN, ROUND_NEAREST }
	public enum SimulationModel			{ SPREADS, LIBORS }
	
	// Flags
	private Measure	measure = Measure.SPOT;
	private StateSpace stateSpace = StateSpace.NORMAL;
	private HandleSimulationTime handleSimulationTime = HandleSimulationTime.ROUND_NEAREST;
	private InterpolationMethod	interpolationMethod	= InterpolationMethod.LOG_LINEAR_UNCORRECTED;
	private SimulationModel simulationModel = SimulationModel.SPREADS;

	
	private final LIBORMarketModel _nonDefaultableModel;
	private final DefaultableLIBORCovarianceModel _covarianceModel;
	private final ForwardCurve _forwardRateCurve;
	private final AnalyticModel _analyticModel;
	
	private double valueCap = 1E3d;
	private double valueFloor = Double.MIN_NORMAL;
	
	
	// Mutual/Transient Fields
	private transient Vector<RandomVariable> interpolationDriftAdjustmentsTerminal = new Vector<>();
	private transient MonteCarloProcess numerairesProcess;
	private transient double[][][] integratedLIBORCovariance;
	private transient TimeDiscretization integrationTimeDiscretizationCache;
	private transient Object integratedLIBORCovarianceLazyInitLock = new Object();


	// ------------------------------------------------------------------------- Constructors --------------------------------------------------------------------------
	
	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel undefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, AnalyticModel analyticModel, Map<String, String> propertyMap) {
		this(undefaultableModel, covarianceModel, forwardRateCurve, analyticModel);
		if(propertyMap != null) {
			Set<String> keys = propertyMap.keySet();
			for(String key: keys) {
				switch(key.toLowerCase()) {
				case "measure":
					try {
						measure = Measure.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"measure\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					break;
				case "statespace":
					try {
						stateSpace = StateSpace.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"statespace\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					break;
				case "interpolationmethod":
					try {
						interpolationMethod = InterpolationMethod.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"interpolationmethod\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					break;
				case "handlesimulationtime":
					try {
						handleSimulationTime = HandleSimulationTime.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"handlesimulationtime\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					break;
				case "simulationmodel":
					try {
						simulationModel = SimulationModel.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"simulationmodel\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					break;
				default:
					continue;

				}
			}
		}
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel undefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, Map<String, String> propertyMap) {
		this(undefaultableModel, covarianceModel, forwardRateCurve, null, propertyMap);
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel undefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, AnalyticModel analyticModel) {
		_nonDefaultableModel = undefaultableModel;
		_covarianceModel = covarianceModel;
		_forwardRateCurve = forwardRateCurve;
		_analyticModel = analyticModel;
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel undefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve) {
		this(undefaultableModel, covarianceModel, forwardRateCurve, null, null);
	}
	
	
	
	// ------------------------------------------------------------------------ Getters --------------------------------------------------------------------------------
	
	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return getUndefaultableLIBORModel().getLiborPeriodDiscretization();
	}

	public Measure getMeasure() {
		return measure;
	}
	

	public StateSpace getStateSpace() {
		return stateSpace;
	}
	

	public InterpolationMethod getInterpolationMethod() {
		return interpolationMethod;
	}
	

	public SimulationModel getSimulationModel() {
		return simulationModel;
	}
	
	public HandleSimulationTime getHandleSimulationTime() {
		return handleSimulationTime;
	}

	

	@Override
	public AnalyticModel getAnalyticModel() {
		return _analyticModel;
	}

	@Override
	public DiscountCurve getDiscountCurve() {
		return getUndefaultableLIBORModel().getDiscountCurve();
	}

	/**
	 * Returns the initial forward rate curve of the defaultable model.
	 */
	@Override
	public ForwardCurve getForwardRateCurve() {
		return _forwardRateCurve;
	}

	@Override
	public LIBORMarketModel getUndefaultableLIBORModel() {
		return _nonDefaultableModel;
	}

	@Override
	public DefaultableLIBORCovarianceModel getCovarianceModel() {
		return _covarianceModel;
	}

	

	// ------------------------------------------------------------------------ Dimension Methods ----------------------------------------------------------------------

	public int getNumberOfLibors() {
		return getNumberOfLIBORPeriods();
	}
	
	@Override
	public int getNumberOfComponents() {
		return 2 * getNumberOfLIBORPeriods();
	}

	public int getNumberOfLIBORPeriods() {
		return getLiborPeriodDiscretization().getNumberOfTimeSteps();
	}

	@Override
	public int getNumberOfFactors() {
		return getCovarianceModel().getNumberOfFactors();
	}

	@Override
	public double getLiborPeriod(int componentIndex) {
		if(componentIndex > getNumberOfLIBORPeriods())
			return getUndefaultableLIBORModel().getLiborPeriod(componentIndex - getNumberOfLIBORPeriods());
		else
			return getUndefaultableLIBORModel().getLiborPeriod(componentIndex);
	}

	@Override
	public int getLiborPeriodIndex(double time) {
		return getUndefaultableLIBORModel().getLiborPeriodIndex(time);
	}

	
	
	// --------------------------------------------------------------------- SDE Methods ------------------------------------------------------------------------
	
	@Override
	public RandomVariable[] getInitialState(MonteCarloProcess process) {
		final double[] initialStatesD = new double[getNumberOfLIBORPeriods()];
		
		// Fetch double values of Defaultable model
		for(int liborIndex = 0; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
			final double rate = getForwardRateCurve().getForward(getAnalyticModel(), getLiborPeriod(liborIndex));
			initialStatesD[liborIndex] = rate;
		}

		final RandomVariable[] initialStateRV = Arrays.copyOf(getUndefaultableLIBORModel().getInitialState(null), getNumberOfComponents());
		if(simulationModel == SimulationModel.SPREADS) {
			for(int liborIndex = 0; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				initialStateRV[getSpreadComponentIndex(liborIndex)] = initialStateRV[liborIndex].bus(initialStatesD[liborIndex]);
				initialStateRV[getSpreadComponentIndex(liborIndex)] = initialStateRV[getSpreadComponentIndex(liborIndex)].log();
			}
		}
		else {
			for(int liborIndex = 0; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				initialStateRV[getDefaultableComponentIndex(liborIndex)] = getRandomVariableForConstant(initialStatesD[liborIndex]);
				if(stateSpace == StateSpace.LOGNORMAL)
					initialStateRV[getDefaultableComponentIndex(liborIndex)] = initialStateRV[getDefaultableComponentIndex(liborIndex)].log();
			}
		}
		return initialStateRV;
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor) {
		final RandomVariable[] nonDefaultableDriftVector = getDriftOfUndefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
		RandomVariable[] fullDriftVector = Arrays.copyOf(nonDefaultableDriftVector, getNumberOfComponents());
		if(simulationModel == SimulationModel.SPREADS) {
			RandomVariable[] spreadDrift = getDriftOfSpread(process, timeIndex, realizationAtTimeIndex, realizationPredictor, nonDefaultableDriftVector);
			for(int i=getNumberOfLIBORPeriods(); i < getNumberOfComponents(); i++) {
				fullDriftVector[i] = spreadDrift[i - getNumberOfLIBORPeriods()];
			}
		}
		else {			
			RandomVariable[] defaultableDrift = getDriftOfDefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
			for(int i=getNumberOfLIBORPeriods(); i < getNumberOfComponents(); i++) {
				fullDriftVector[i] = defaultableDrift[i - getNumberOfLIBORPeriods()];
			}
		}
		return fullDriftVector;
	}

	@Override
	public RandomVariable[] getDriftOfUndefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		RandomVariable[] sRealization = null;
		if(realizationAtTimeIndex != null) {
			sRealization = Arrays.copyOf(realizationAtTimeIndex, getUndefaultableLIBORModel().getNumberOfComponents());
		}
		RandomVariable[] sRealizationPredictor = null; 
		if(realizationPredictor != null) {
			sRealizationPredictor = Arrays.copyOf(realizationPredictor, getUndefaultableLIBORModel().getNumberOfComponents());
		}
		return getUndefaultableLIBORModel().getDrift(getUndefaultableProcess(process), timeIndex, sRealization, sRealizationPredictor);
	}

	@Override
	public RandomVariable[] getDriftOfDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		if(simulationModel != SimulationModel.LIBORS) {
			throw new UnsupportedOperationException("This method only makes sense for Models modeling the defaultable LIBOR rates!");
		}
		return getDriftDefaultableInternally(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
	}

	public RandomVariable[] getDriftOfSpread(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor, RandomVariable[] nonDefaultableDrifts) {
		
		if(simulationModel != SimulationModel.SPREADS) {
			throw new UnsupportedOperationException("This method only makes sense for Models simulating the Spreads!");
		}
		
		RandomVariable[] drift = null;
		
		final double time = process.getTime(timeIndex);
		if(time >= getLiborPeriod(getNumberOfLIBORPeriods() - 1)) {
			// If the non-defaultable Drift is not null, it is filled with nulls, as we are above the last fixing time
			drift = nonDefaultableDrifts;
			if(nonDefaultableDrifts == null) {
				drift = new RandomVariable[getNumberOfLIBORPeriods()];
			}
			return drift;
		}
		
		// If non-defaultable drift is null we fetch it. 
		if(nonDefaultableDrifts == null)
			nonDefaultableDrifts = getDriftOfUndefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
		
		RandomVariable[] realizationAllLIBORS = null;
		if(realizationAtTimeIndex != null) {
			realizationAllLIBORS = new RandomVariable[getNumberOfComponents()];
			for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
				realizationAllLIBORS[libor] = realizationAtTimeIndex[libor];
				realizationAllLIBORS[libor + getNumberOfLIBORPeriods()] = realizationAtTimeIndex[libor].add(realizationAtTimeIndex[libor + getNumberOfLIBORPeriods()]);
				
			}
		}
		
		RandomVariable[] realizationPredAllLIBORS = null;
		if(realizationPredictor != null) {
			realizationPredAllLIBORS = new RandomVariable[getNumberOfComponents()];
			for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
				realizationPredAllLIBORS[libor] = realizationPredictor[libor];
				realizationPredAllLIBORS[libor + getNumberOfLIBORPeriods()] = realizationPredictor[libor].add(realizationPredictor[libor + getNumberOfLIBORPeriods()]);
				
			}
		}
		
		drift = getDriftDefaultableInternally(process, timeIndex, realizationAllLIBORS, realizationPredAllLIBORS);
		
		// TODO Add functionality of lognormal. I.e. add 0.5*sigma^2 to both drifts and multiply with realization
		
		for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
			if(drift[libor] == null || nonDefaultableDrifts[libor] == null) {
				drift[libor] = null;
				continue;
			}
			drift[libor] = drift[libor].sub(nonDefaultableDrifts[libor]);
			drift[libor] = drift[libor].div(realizationAtTimeIndex[getSpreadComponentIndex(libor)]);
			RandomVariable[] factorLoadingSpread = getFactorLoadingOfSpread(process, timeIndex, libor, realizationAtTimeIndex);
			for(int factor = 0; factor < getNumberOfFactors(); factor++) {
				drift[libor] = drift[libor].sub(factorLoadingSpread[factor].squared().mult(0.5));
			}
		}

		// This implementation only holds for one type of drift:
		/*
		RandomVariable[] realizationOfAllLIBORS = new RandomVariable[getNumberOfComponents()];
		int component = 0;
		for(; component < getNumberOfLIBORPeriods(); component++) {
			realizationOfAllLIBORS[component] = realizationAtTimeIndex[component];
		}
		for(; component < getNumberOfComponents(); component++) {
			realizationOfAllLIBORS[component] = realizationAtTimeIndex[component - getNumberOfLIBORPeriods()]
					.add(realizationAtTimeIndex[component]); // Add Spread (we have simulationModel == SPREADS)
		}
		
		// Fetch all factor
		RandomVariable[][] factorLoadingsDef = new RandomVariable[getNumberOfLIBORPeriods()][];
		
		
		for(int factor = getUndefaultableLIBORModel().getNumberOfFactors(); factor < getNumberOfFactors(); factor++) {
			RandomVariable innerSum = zero;
			for(int liborIndex = firstNonZeroComponent; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				
				final double deltaT = getLiborPeriod(liborIndex + 1) - getLiborPeriod(liborIndex);
				
				if(factor == 0) {
					drift[liborIndex] = drift[liborIndex].mult(deltaT).div(realizationAtTimeIndex[liborIndex].mult(deltaT).add(1.0));
					factorLoadingsDef[liborIndex] = getCovarianceModel().getFactorLoading(time, getDefaultableComponentIndex(liborIndex), realizationOfAllLIBORS);
				}
				
				innerSum = innerSum.add(factorLoadingsDef[liborIndex][factor].mult(deltaT).div(realizationOfAllLIBORS[getSpreadComponentIndex(liborIndex)].mult(deltaT).add(1.0)));
				
				// Ideally this is just a Constant
				RandomVariable freeParam = getFactorLoadingOfSpread(process, timeIndex, timeIndex, realizationAtTimeIndex)[factor];
				
				drift[liborIndex] = drift[liborIndex].addProduct(innerSum, freeParam);
			}
		}
		
		
		// Spread is lognormal hence we apply drift adjustment: µ - 0.5*(sigma^2) = µ - 0.5*sum(lambda*lambda)
		for(int componentIndex = firstNonZeroComponent; componentIndex < getNumberOfLIBORPeriods(); componentIndex++) {
			RandomVariable[] factorLoadingSpread = getFactorLoadingOfSpread(process, timeIndex, timeIndex, realizationAtTimeIndex);
			final RandomVariable variance = getCovarianceModel().getCovariance(time, componentIndex + getNumberOfLIBORPeriods(), componentIndex + getNumberOfLIBORPeriods(), realizationAtTimeIndex);
			drift[componentIndex] = drift[componentIndex].addProduct(variance, -0.5);
		}
		*/
		return drift;
	}
	
	@Override
	public RandomVariable[] getFactorLoading(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		if(componentIndex < getNumberOfLIBORPeriods()) {
			return getFactorLoadingOfNonDefaultableModel(process, timeIndex, componentIndex, realizationAtTimeIndex);
		}
		final int liborPeriodIndex = componentIndex % getNumberOfLIBORPeriods();
		if (simulationModel == SimulationModel.SPREADS) {
			return getFactorLoadingOfSpread(process, timeIndex, liborPeriodIndex, realizationAtTimeIndex);
		}
		else {			
			return getFactorLoadingOfDefaultableModel(process, timeIndex, liborPeriodIndex, realizationAtTimeIndex);
		}
	}

	public RandomVariable[] getFactorLoadingOfNonDefaultableModel(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		if(componentIndex >= getNumberOfLIBORPeriods())
			throw new ArrayIndexOutOfBoundsException("Non defaultable model is a model of " + getNumberOfLIBORPeriods() + " Components. Index " + componentIndex + " out of Bounds for Non-Defaultable Model");
		
		return getCovarianceModel().getFactorLoading(process.getTime(timeIndex), componentIndex, realizationAtTimeIndex);
	}
	
	public RandomVariable[] getFactorLoadingOfDefaultableModel(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		if(componentIndex >= getNumberOfLIBORPeriods())
			throw new ArrayIndexOutOfBoundsException("Defaultable model is a model of " + getNumberOfLIBORPeriods() + " Components. Index " + componentIndex + " out of Bounds for Defaultable Model");
		
		return getCovarianceModel().getFactorLoading(process.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), realizationAtTimeIndex);
	}
	
	public RandomVariable[] getFactorLoadingOfSpread(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex) {
		if(simulationModel != SimulationModel.SPREADS)
			throw new UnsupportedOperationException("This method only makes sense for Models simulating the Spreads!");
		
		if(componentIndex >= getNumberOfLIBORPeriods())
			throw new ArrayIndexOutOfBoundsException("Spread model is a model of " + getNumberOfLIBORPeriods() + " Components. Index " + componentIndex + " out of Bounds for Spread");
		
		// For now this is not necessary. The Factor Loading of the Spread depends solely on the FL of the non-defaultable model and the free parameter Matrix!
		// Maybe cause of Heap Overflow.
		RandomVariable[] realizationForCovariance = Arrays.copyOf(realizationAtTimeIndex, getNumberOfComponents());
		// Add Spread to non-defaultable LIBORs to get defaultable LIBORS
		for(int component = getNumberOfLIBORPeriods(); component < getNumberOfComponents(); component++) {
			realizationForCovariance[component] = realizationForCovariance[component - getNumberOfLIBORPeriods()].add(realizationAtTimeIndex[component]);
		}
		RandomVariable[] factorLoading = null;
		try {
			factorLoading = getCovarianceModel().getFactorLoadingOfSpread(process.getTime(timeIndex), componentIndex, realizationForCovariance);
		} catch(UnsupportedOperationException ex) {
			// Sub defaultable FactorLoading from factorLoading and divide by Spread
			factorLoading = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), realizationForCovariance);
			
			RandomVariable[] factorLoadingNonDef = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), componentIndex, realizationForCovariance);
			for(int i =0; i < getNumberOfFactors(); i++) {
				factorLoading[i] = factorLoading[i].sub(factorLoadingNonDef[i]).div(realizationAtTimeIndex[getSpreadComponentIndex(componentIndex)]);
			}
		}
		return factorLoading;
	}
	
	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLIBORPeriods()) {
			return getUndefaultableLIBORModel().applyStateSpaceTransform(getUndefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		}
		
		RandomVariable value = randomVariable;

		if(simulationModel == SimulationModel.SPREADS) {
			value = value.exp().floor(1E-5);
		}
		else if(stateSpace == StateSpace.LOGNORMAL){
			value = value.exp().floor(1E-5);
		}
		
		if(!Double.isInfinite(valueCap)) {
			value = value.cap(valueCap);
		}
		
		if(!Double.isInfinite(-valueFloor)) {
			value = value.floor(valueFloor);
		}
		return value;
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(final MonteCarloProcess process, final int timeIndex, final int componentIndex, final RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLIBORPeriods()) {
			return getUndefaultableLIBORModel().applyStateSpaceTransformInverse(getUndefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		}
		
		RandomVariable value = randomVariable;
		
		if(simulationModel == SimulationModel.SPREADS) {
			value = value.log();			
		}
		else if(stateSpace == StateSpace.LOGNORMAL){
			value = value.log();
		}
		
		return value;
	}

	
	
	// --------------------------------------------------------------------- MC Process Value Methods ------------------------------------------------------------------------
	
	@Override
	public RandomVariable getDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		if(timeIndex == 0 || getLiborPeriod(liborIndex) == 0d)
			return getRandomVariableForConstant(getForwardRateCurve().getForward(getAnalyticModel(), getLiborPeriod(liborIndex)));
		
		if (simulationModel == SimulationModel.SPREADS)
			return process.getProcessValue(timeIndex, liborIndex).add(process.getProcessValue(timeIndex, getSpreadComponentIndex(liborIndex)));
		else
			return process.getProcessValue(timeIndex, getDefaultableComponentIndex(liborIndex));
	}
	

	@Override
	public RandomVariable getDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		final int periodStartIndex    = getLiborPeriodIndex(periodStart);
		final int periodEndIndex      = getLiborPeriodIndex(periodEnd);

		// If time is beyond fixing, use the fixing time.
		time = Math.min(time, periodStart);
		int timeIndex           = process.getTimeIndex(time);
		// If time is not part of the discretization, use the nearest available point.
		if(timeIndex < 0) {
			timeIndex = -timeIndex-2;
			if(handleSimulationTime == HandleSimulationTime.ROUND_NEAREST && time-process.getTime(timeIndex) > process.getTime(timeIndex+1)-time) {
				timeIndex++;
			}
		}

		// The forward rates are provided on fractional tenor discretization points using linear interpolation. See ISBN 0470047224.

		// Interpolation on tenor using interpolationMethod
		if(periodEndIndex < 0) {
			final int		previousEndIndex	= -periodEndIndex - 2; // Closest LIBOR Index for which the LIBOR Fixing is smaller than periodEnd
			final double	nextEndTime			= getLiborPeriod(previousEndIndex+1);
			// Interpolate libor from periodStart to periodEnd on periodEnd
			final RandomVariable onePlusLongLIBORdt         = getForwardRate(process, time, periodStart, nextEndTime).mult(nextEndTime - periodStart).add(1.0);
			final RandomVariable onePlusInterpolatedLIBORDt = getOnePlusInterpolatedLIBORDt(process, timeIndex, periodEnd, previousEndIndex);
			return onePlusLongLIBORdt.div(onePlusInterpolatedLIBORDt).sub(1.0).div(periodEnd - periodStart);
		}

		// Interpolation on tenor using interpolationMethod
		if(periodStartIndex < 0) {
			final int	previousStartIndex   = -periodStartIndex - 2; // Closest LIBOR Index for which the LIBOR Fixing is smaller than periodStart
			final double nextStartTime	 = getLiborPeriod(previousStartIndex+1);
			if(nextStartTime > periodEnd) {
				throw new AssertionError("Interpolation not possible.");
			}
			if(nextStartTime == periodEnd) {
				return getOnePlusInterpolatedLIBORDt(process, timeIndex, periodStart, previousStartIndex).sub(1.0).div(periodEnd - periodStart);
			}
			//			RandomVariable onePlusLongLIBORdt         = getLIBOR(Math.min(prevStartTime, time), nextStartTime, periodEnd).mult(periodEnd - nextStartTime).add(1.0);
			final RandomVariable onePlusLongLIBORdt         = getForwardRate(process, time, nextStartTime, periodEnd).mult(periodEnd - nextStartTime).add(1.0);
			final RandomVariable onePlusInterpolatedLIBORDt = getOnePlusInterpolatedLIBORDt(process, timeIndex, periodStart, previousStartIndex);
			return onePlusLongLIBORdt.mult(onePlusInterpolatedLIBORDt).sub(1.0).div(periodEnd - periodStart);
		}

		if(periodStartIndex < 0 || periodEndIndex < 0) {
			throw new AssertionError("LIBOR requested outside libor discretization points and interpolation was not performed.");
		}

		// If this is a model primitive then return it
		if(periodStartIndex+1==periodEndIndex) {
			return getLIBOR(process, timeIndex, periodStartIndex);
		}

		// The requested LIBOR is not a model primitive. We need to calculate it (slow!)
		RandomVariable accrualAccount = null;

		// Calculate the value of the forward bond
		for(int periodIndex = periodStartIndex; periodIndex < periodEndIndex; periodIndex++)
		{
			final double subPeriodLength = getLiborPeriod(periodIndex+1) - getLiborPeriod(periodIndex);
			final RandomVariable liborOverSubPeriod = getLIBOR(process, timeIndex, periodIndex);

			accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
		}

		final RandomVariable libor = accrualAccount.sub(1.0).div(periodEnd - periodStart);

		return libor;
	}
	
	
	public RandomVariable getDefaultableBond(MonteCarloProcess process, final double time, final double maturity) throws CalculationException {
		RandomVariable forwardRate = getForwardRate(process, time, time, maturity);
		return (new Scalar(1.0)).discount(forwardRate, maturity - time);
	}
	
	
	@Override
	public RandomVariable getUndefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		if(timeIndex == 0 || getLiborPeriod(liborIndex) == 0d)
			return getRandomVariableForConstant(getUndefaultableLIBORModel().getForwardRateCurve().getForward(_analyticModel, getLiborPeriod(liborIndex)));
		return process.getProcessValue(timeIndex, liborIndex);
	}
	

	@Override
	public RandomVariable getUndefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getUndefaultableLIBORModel().getForwardRate(getUndefaultableProcess(process), time, periodStart, periodEnd);
	}

	
	@Override
	public RandomVariable getForwardDiscountBond(MonteCarloProcess process, final double time, final double maturity) throws CalculationException {
		return getUndefaultableLIBORModel().getForwardDiscountBond(getUndefaultableProcess(process), time, maturity);
	}
	
	
	
	// --------------------------------------------------------------------- Spread Methods ------------------------------------------------------------------------

	@Override
	public RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		if(timeIndex == 0 || getLiborPeriod(liborIndex) == 0d) {
			return getDefaultableLIBOR(process, timeIndex, liborIndex).sub(getUndefaultableLIBOR(process, timeIndex, liborIndex));
		}
		
		if(simulationModel == SimulationModel.SPREADS) 
			return process.getProcessValue(timeIndex, getSpreadComponentIndex(liborIndex));
		else
			return getLIBOR(process, timeIndex, liborIndex).sub(getUndefaultableLIBOR(process, timeIndex, liborIndex));
	}
	
	
	@Override
	public RandomVariable getSpread(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getForwardRate(process, time, periodStart, periodEnd).sub(getUndefaultableForwardRate(process, time, periodStart, periodEnd));
	}

	
	@Override
	public RandomVariable getSurvivalProbability(MonteCarloProcess process, double evaluationTime, double maturity) throws CalculationException {
		return getDefaultableBond(process, evaluationTime, maturity).div(getForwardDiscountBond(process, evaluationTime, maturity));
	}

	
	
	// -------------------------------------------------------------------- Clones --------------------------------------------------------------------
	
	@Override
	public DefaultableLIBORMarketModelFromCovarianceModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel) {
		Map<String, String> propertyMap = new HashMap<String, String>();
		propertyMap.put("measure", measure.toString());
		propertyMap.put("statespace", stateSpace.toString());
		propertyMap.put("handlesimulationtime", handleSimulationTime.toString());
		propertyMap.put("interpolationmethod", interpolationMethod.toString());
		
		DefaultableLIBORCovarianceModel newCovarianceModel = getCovarianceModel();
		if(newUndefaultableModel.getCovarianceModel() != newCovarianceModel.getUndefaultableCovarianceModel()) {
			newCovarianceModel = newCovarianceModel.getCloneWithModifiedUndefaultableCovariance(newUndefaultableModel.getCovarianceModel());
		}
		return new DefaultableLIBORMarketModelFromCovarianceModel(newUndefaultableModel, newCovarianceModel, getForwardRateCurve(), getAnalyticModel(), propertyMap);
	}
	
	@Override
	public DefaultableLIBORMarketModelFromCovarianceModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel newCovarianceModel) {
		Map<String, String> propertyMap = new HashMap<String, String>();
		propertyMap.put("measure", measure.toString());
		propertyMap.put("statespace", stateSpace.toString());
		propertyMap.put("handlesimulationtime", handleSimulationTime.toString());
		propertyMap.put("interpolationmethod", interpolationMethod.toString());
		
		if(!(newCovarianceModel instanceof DefaultableLIBORCovarianceModel)) {
			throw new IllegalArgumentException("newCovarianceModel must be of type DefaultableLIBORCovarianceModel!");
		}
		return new DefaultableLIBORMarketModelFromCovarianceModel(getUndefaultableLIBORModel(), (DefaultableLIBORCovarianceModel)newCovarianceModel, getForwardRateCurve(), getAnalyticModel(), propertyMap);
	}

	@Override
	public DefaultableLIBORMarketModelFromCovarianceModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		LIBORMarketModel undefaultableModel = (LIBORMarketModel)dataModified.getOrDefault("undefaultableModel", getUndefaultableLIBORModel());
		DefaultableLIBORCovarianceModel covarianceModel = (DefaultableLIBORCovarianceModel)dataModified.getOrDefault("covarianceModel", getCovarianceModel());
		ForwardCurve forwardRateCurve = (ForwardCurve)dataModified.getOrDefault("forwardRateCurve", getForwardRateCurve());
		AnalyticModel analyticModel = (AnalyticModel)dataModified.getOrDefault("analyticModel", getAnalyticModel());
		
		Map<String, String> propertyMap= new HashMap<String, String>();
		propertyMap.put("measure", measure.toString());
		propertyMap.put("statespace", stateSpace.toString());
		propertyMap.put("handlesimulationtime", handleSimulationTime.toString());
		propertyMap.put("interpolationmethod", interpolationMethod.toString());
		
		for(String key: dataModified.keySet()) {
			switch(key.toLowerCase()) {
			case "measure":
				propertyMap.put(key, dataModified.get(key).toString());
				break;
			case "statespace":
				propertyMap.put(key, dataModified.get(key).toString());
				break;
			case "handlesimulationtime":
				propertyMap.put(key, dataModified.get(key).toString());
				break;
			case "interpolationmethod":
				propertyMap.put(key, dataModified.get(key).toString());
				break;
			}
		}
		return new DefaultableLIBORMarketModelFromCovarianceModel(undefaultableModel, covarianceModel, forwardRateCurve, analyticModel, propertyMap);
	}
	
	
	
	// --------------------------------------------------------------------- Valuation Methods ----------------------------------------------------------------
	
	@Override
	public double[][][] getIntegratedLIBORCovariance(TimeDiscretization timeDiscretization) {
		synchronized (integratedLIBORCovarianceLazyInitLock) {
			if(integrationTimeDiscretizationCache != timeDiscretization) {
				integratedLIBORCovariance = null;
			}
			if(integratedLIBORCovariance == null) {
				
				integrationTimeDiscretizationCache = timeDiscretization;
				
				final TimeDiscretization liborPeriodDiscretization = getLiborPeriodDiscretization();

				integratedLIBORCovariance = new double[timeDiscretization.getNumberOfTimeSteps()][getNumberOfLIBORPeriods()][getNumberOfLIBORPeriods()];
				
				for(int timeIndex = 0; timeIndex < timeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					final double dt = timeDiscretization.getTimeStep(timeIndex);
					
					final RandomVariable[][] factorLoadings = new RandomVariable[getNumberOfLIBORPeriods()][];
					
					// Prefetch factor loadings. If they are stochastic, fetch approximated factor loadings by initial values
					{
						RandomVariable[] initialLIBORValues = null;
						for(int componentIndex = 0; componentIndex < getNumberOfLIBORPeriods(); componentIndex++) {
							try {
								
								factorLoadings[componentIndex] = getCovarianceModel().getFactorLoading(timeDiscretization.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), initialLIBORValues);
								
								if(factorLoadings[componentIndex] == null) {
									// Need Approximation
									throw new Exception();
								}
								if(factorLoadings[componentIndex][0] == null) {
									// Need Approximation
									throw new Exception();
								}
								
							} 
							catch (Exception ex) {
								if(initialLIBORValues == null) {
									// Get Approximation
									initialLIBORValues = new RandomVariable[getNumberOfComponents()];
									// Factor Loading definately needs Initial values of the form (L_0, ..., L_n, L^d_0, ..., L^d_n)
									for(int i=0; i < getNumberOfComponents(); i++) {
										try {
											initialLIBORValues[i] = i < getNumberOfLIBORPeriods() ? getUndefaultableLIBOR(null, 0, i) : getDefaultableLIBOR(null, 0, i);
										} catch (CalculationException e) {
											e.printStackTrace();
											throw new UnsupportedOperationException("Cannot get Initial LIBOR rates as Estimators!");
										}
									}
								}
								factorLoadings[componentIndex] = getCovarianceModel().getFactorLoading(timeDiscretization.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), initialLIBORValues);
							}
						}
					}
					
					
					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						final RandomVariable[] factorLoadingOfComponent1 = factorLoadings[componentIndex1];

						// Use symmetry:
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							double integratedLIBORCovarianceValue = 0.0; // Why not use getCovarianceModel().getCovariance(...)
							
							// sigma(c1, c2, t) = 0 if t >= fixing time of component
							if(getLiborPeriod(componentIndex1) > timeDiscretization.getTime(timeIndex)) {
								final RandomVariable[] factorLoadingOfComponent2 = factorLoadings[componentIndex2];
								for(int factorIndex = 0; factorIndex < factorLoadingOfComponent2.length; factorIndex++) {
									integratedLIBORCovarianceValue += factorLoadingOfComponent1[factorIndex].doubleValue() * factorLoadingOfComponent2[factorIndex].doubleValue() * dt;
								}
							}
							integratedLIBORCovariance[timeIndex][componentIndex1][componentIndex2] = integratedLIBORCovarianceValue;
						}
					}
				}

				// Integrate over time (i.e. sum up).
				for(int timeIndex = 1; timeIndex < timeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					final double[][] prevIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex-1];
					final double[][] thisIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex];
					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							thisIntegratedLIBORCovariance[componentIndex1][componentIndex2] = prevIntegratedLIBORCovariance[componentIndex1][componentIndex2] + thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
							thisIntegratedLIBORCovariance[componentIndex2][componentIndex1] = thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
						}
					}
				}
			}
		}

		return integratedLIBORCovariance;
	}

	
	@Override
	public RandomVariable getNumeraire(MonteCarloProcess process, double time) throws CalculationException {
		return getUndefaultableLIBORModel().getNumeraire(getUndefaultableProcess(process), time);
	}
	
	
	
	// --------------------------------------------------------------------- Other Public Methods ----------------------------------------------------------------
	
	@Override
	public LocalDateTime getReferenceDate() {
		return getUndefaultableLIBORModel().getReferenceDate();
	}

	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return getUndefaultableLIBORModel().getRandomVariableForConstant(value);
	}

	
	
	// -------------------------------------------------------------------- Private Methods: --------------------------------------------------------------------------

	private int getSpreadComponentIndex(int liborIndex) {
		return liborIndex + getNumberOfLIBORPeriods();
	}
	
	private int getDefaultableComponentIndex(int liborIndex) {
		return liborIndex + getNumberOfLIBORPeriods();
	}
		
	private MonteCarloProcess getUndefaultableProcess(MonteCarloProcess fullProcess) {
		if(fullProcess == null)
			return null;
		return Functional.getComponentReducedMCProcess(fullProcess, 0, getNumberOfLIBORPeriods() - 1);
		
	}
	
	private void ensureCacheConsistency(final MonteCarloProcess process) {
		/*
		 * Check if caches are valid (i.e. process did not change)
		 */
		if (process != numerairesProcess) {
			numerairesProcess = process;
			interpolationDriftAdjustmentsTerminal.clear();
		}
	}
	
	
	private RandomVariable[] getDriftDefaultableInternally(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		
		final double time = process.getTime(timeIndex);			// t - current simulation time
		int	firstForwardRateIndex = this.getLiborPeriodIndex(time) + 1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = - firstForwardRateIndex - 1 + 1;
		}

		final RandomVariable zero = Scalar.of(0.0);

		// Allocate drift vector and initialize to zero (will be used to sum up drift components)
		final RandomVariable[]	drift = new RandomVariable[getNumberOfLIBORPeriods()];
		for(int liborIndex = firstForwardRateIndex; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
			drift[liborIndex] = zero;
		}

		// Allocate array (for each k) for the sums of delta_{i}/(1+L_{i} \delta_i) f_{i,k} (+ for spot measure, - for terminal measure)
		final RandomVariable[]	factorLoadingsSums	= new RandomVariable[getNumberOfFactors()];		
		Arrays.fill(factorLoadingsSums, zero);

		if(measure == Measure.SPOT) {
			// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
			for(int liborIndex=firstForwardRateIndex; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				final double			periodLength	= getLiborPeriodDiscretization().getTimeStep(liborIndex);
				final RandomVariable	forwardRate		= realizationAtTimeIndex[getDefaultableComponentIndex(liborIndex)];
				RandomVariable			oneStepMeasureTransform = Scalar.of(periodLength).discount(forwardRate, periodLength);
				
				if(stateSpace == StateSpace.LOGNORMAL) {	// The drift has an additional forward rate factor
					oneStepMeasureTransform = oneStepMeasureTransform.mult(forwardRate);
				}
				
				// Here we are not definately in the simulationModel == LIBORS!
				final RandomVariable[]	factorLoading = getFactorLoadingOfDefaultableModel(process, timeIndex, liborIndex, realizationAtTimeIndex);
				for(int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
				
				drift[liborIndex] = drift[liborIndex].addSumProduct(factorLoadingsSums, factorLoading);
			}
		}
		else if(measure == Measure.TERMINAL) {
			// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
			for(int liborIndex=getNumberOfLIBORPeriods()-1; liborIndex>=firstForwardRateIndex; liborIndex--) {
				final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborIndex);
				final RandomVariable forwardRate = realizationAtTimeIndex[getDefaultableComponentIndex(liborIndex)];
				
				RandomVariable oneStepMeasureTransform = Scalar.of(-periodLength).discount(forwardRate, periodLength);
				if(stateSpace == StateSpace.LOGNORMAL) {
					oneStepMeasureTransform = oneStepMeasureTransform.mult(forwardRate);
				}

				final RandomVariable[] factorLoading = getFactorLoadingOfDefaultableModel(process, timeIndex, liborIndex, realizationAtTimeIndex);
				drift[liborIndex] = drift[liborIndex].addSumProduct(factorLoadingsSums, factorLoading);
				for(int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
			}
		}
		else {
			throw new IllegalArgumentException("Drift not implemented for specified measure.");
		}
		if(stateSpace == StateSpace.LOGNORMAL) {
			// Drift adjustment for log-coordinate in each component
			for(int componentIndex = firstForwardRateIndex; componentIndex < getNumberOfLIBORPeriods(); componentIndex++) {
				final RandomVariable variance = getCovarianceModel().getCovariance(time, componentIndex + getNumberOfLIBORPeriods(), componentIndex + getNumberOfLIBORPeriods(), realizationAtTimeIndex);
				drift[componentIndex] = drift[componentIndex].addProduct(variance, -0.5);
			}
		}
		return drift;
	}
	
	/**
	 * Copied from LIBORMarketModelFromCovarianceModel
	 */
	private RandomVariable getOnePlusInterpolatedLIBORDt(final MonteCarloProcess process, int timeIndex, final double periodStartTime, final int liborPeriodIndex) throws CalculationException
	{
		final double tenorPeriodStartTime       = getLiborPeriod(liborPeriodIndex);
		final double tenorPeriodEndTime         = getLiborPeriod(liborPeriodIndex + 1);
		final double tenorDt                    = tenorPeriodEndTime - tenorPeriodStartTime;
		if(tenorPeriodStartTime < process.getTime(timeIndex)) {
			// Fixed at Long LIBOR period Start.
			timeIndex  = Math.min(timeIndex, process.getTimeIndex(tenorPeriodStartTime));
			if(timeIndex < 0) {
				//				timeIndex = -timeIndex-2;			// mapping to last known fixing.
				throw new IllegalArgumentException("Tenor discretization not part of time discretization.");
			}
		}
		final RandomVariable onePlusLongLIBORDt = getLIBOR(process, timeIndex , liborPeriodIndex).mult(tenorDt).add(1.0);

		final double smallDt                    = tenorPeriodEndTime - periodStartTime;
		final double alpha                      = smallDt / tenorDt;

		RandomVariable onePlusInterpolatedLIBORDt;
		switch(interpolationMethod)
		{
		case LINEAR:
			onePlusInterpolatedLIBORDt = onePlusLongLIBORDt.mult(alpha).add(1 - alpha);
			break;
		case LOG_LINEAR_UNCORRECTED:
			onePlusInterpolatedLIBORDt = onePlusLongLIBORDt.log().mult(alpha).exp();
			break;
		case LOG_LINEAR_CORRECTED:
			final double adjustmentCoefficient     = 0.5 * smallDt * (tenorPeriodStartTime - periodStartTime);
			RandomVariable adjustment        = getInterpolationDriftAdjustment(process, timeIndex, liborPeriodIndex);
			adjustment = adjustment.mult(adjustmentCoefficient);
			onePlusInterpolatedLIBORDt = onePlusLongLIBORDt.log().mult(alpha).sub(adjustment).exp();
			break;
		default: throw new IllegalArgumentException("Method for enum " + interpolationMethod.name() + " not implemented!");
		}

		// Analytic adjustment for the interpolation
		// @TODO reference to AnalyticModelFromCuvesAndVols must not be null
		// @TODO This adjustment only applies if the corresponding adjustment in getNumeraire is enabled
		final double analyticOnePlusLongLIBORDt   = 1 + getForwardRateCurve().getForward(getAnalyticModel(), tenorPeriodStartTime, tenorDt) * tenorDt;
		final double analyticOnePlusShortLIBORDt	= 1 + getForwardRateCurve().getForward(getAnalyticModel(), periodStartTime, smallDt) * smallDt;

		double analyticOnePlusInterpolatedLIBORDt;
		switch(interpolationMethod)
		{
		case LINEAR:
			analyticOnePlusInterpolatedLIBORDt = analyticOnePlusLongLIBORDt * alpha + (1-alpha);
			break;
		case LOG_LINEAR_UNCORRECTED:
		case LOG_LINEAR_CORRECTED:
			analyticOnePlusInterpolatedLIBORDt = Math.exp(Math.log(analyticOnePlusLongLIBORDt) * alpha);
			break;
		default: throw new IllegalArgumentException("Method for enum " + interpolationMethod.name() + " not implemented!");
		}
		onePlusInterpolatedLIBORDt = onePlusInterpolatedLIBORDt.mult(analyticOnePlusShortLIBORDt / analyticOnePlusInterpolatedLIBORDt);

		return onePlusInterpolatedLIBORDt;
	}
	
	/**
	 * Copied from LIBORMarketModelFromCovarianceModel
	 */
	private RandomVariable getInterpolationDriftAdjustment(final MonteCarloProcess process, final int evaluationTimeIndex, final int liborIndex) throws CalculationException
	{
		switch(interpolationMethod)
		{
		case LINEAR:
		case LOG_LINEAR_UNCORRECTED:
			return null;

		case LOG_LINEAR_CORRECTED:
			final double tenorPeriodStartTime  = getLiborPeriod(liborIndex);
			int    tenorPeriodStartIndex = process.getTimeIndex(tenorPeriodStartTime);
			if(tenorPeriodStartIndex < 0)
			{
				tenorPeriodStartIndex = - tenorPeriodStartIndex - 2;
			}
			// Lazy init of interpolationDriftAdjustmentsTerminal
			if(evaluationTimeIndex == tenorPeriodStartIndex) {
				synchronized(interpolationDriftAdjustmentsTerminal) {
					// Invalidate cache if process has changed
					ensureCacheConsistency(process);

					// Check if value is cached
					if(interpolationDriftAdjustmentsTerminal.size() <= liborIndex) {
						interpolationDriftAdjustmentsTerminal.setSize(getNumberOfLIBORPeriods());
					}

					RandomVariable interpolationDriftAdjustment = interpolationDriftAdjustmentsTerminal.get(liborIndex);
					if(interpolationDriftAdjustment == null) {
						interpolationDriftAdjustment = getInterpolationDriftAdjustmentEvaluated(process, evaluationTimeIndex, liborIndex);
						interpolationDriftAdjustmentsTerminal.set(liborIndex, interpolationDriftAdjustment);
					}

					return interpolationDriftAdjustment;
				}
			}
			else {
				return getInterpolationDriftAdjustmentEvaluated(process, evaluationTimeIndex, liborIndex);
			}
		default: throw new IllegalArgumentException("Method for enum " + interpolationMethod.name() + " not implemented!");
		}
	}
	
	/**
	 * Copied from LIBORMarketModelFromCovarianceModel
	 */
	private RandomVariable getInterpolationDriftAdjustmentEvaluated(final MonteCarloProcess process, final int evaluationTimeIndex, final int liborIndex) throws CalculationException
	{

		final double tenorPeriodStartTime  = getLiborPeriod(liborIndex);
		final double tenorPeriodEndTime    = getLiborPeriod(liborIndex + 1);
		final double tenorDt               = tenorPeriodEndTime - tenorPeriodStartTime;

		RandomVariable driftAdjustment = getRandomVariableForConstant(0.0);

		/*
		 * Integral approximation with trapezoid method.
		 */
		RandomVariable previousIntegrand = getRandomVariableForConstant(0.0);

		/*
		 * Value in 0
		 */
		final RandomVariable[] realizationsAtZero = new RandomVariable[getNumberOfLibors()];
		final RandomVariable[] realizationsForFL = new RandomVariable[getNumberOfLibors()];
		
		for(int liborIndexForRealization = 0; liborIndexForRealization < getNumberOfLIBORPeriods(); liborIndexForRealization++)
		{
			realizationsForFL[liborIndexForRealization] = getUndefaultableLIBOR(process, 0, liborIndexForRealization);
			realizationsForFL[getDefaultableComponentIndex(liborIndexForRealization)] = getLIBOR(process, 0, liborIndexForRealization);
			
			realizationsAtZero[liborIndexForRealization] = getLIBOR(process, 0, liborIndexForRealization);
		}
		final RandomVariable[] factorLoading = getCovarianceModel().getFactorLoading(0, getDefaultableComponentIndex(liborIndex), realizationsForFL);
		//o_{Li}(t)
		for(final RandomVariable oneFactor : factorLoading)
		{
			previousIntegrand = previousIntegrand.add(oneFactor.squared());
		}
		previousIntegrand = previousIntegrand.div( (realizationsAtZero[liborIndex].mult(tenorDt).add(1.0)).squared() );
		
		if(stateSpace == StateSpace.LOGNORMAL)
		{
			previousIntegrand = previousIntegrand.mult( realizationsAtZero[liborIndex].squared() );
		}

		/*
		 * Integration
		 */
		for(int sumTimeIndex = 1; sumTimeIndex <= evaluationTimeIndex; sumTimeIndex++)
		{
			final RandomVariable[] realizationsAtTimeIndex = new RandomVariable[getNumberOfLIBORPeriods()];
			
			for(int liborIndexForRealization = 0; liborIndexForRealization < getNumberOfLIBORPeriods(); liborIndexForRealization++)
			{
				int evaluationTimeIndexForRealizations = Math.min(sumTimeIndex, process.getTimeIndex(getLiborPeriod(liborIndexForRealization)));
				if(evaluationTimeIndexForRealizations < 0)
				{
					evaluationTimeIndexForRealizations = - evaluationTimeIndexForRealizations - 2;
				}
				realizationsForFL[liborIndexForRealization] = getUndefaultableLIBOR(process, evaluationTimeIndexForRealizations, liborIndexForRealization);
				realizationsForFL[getDefaultableComponentIndex(liborIndexForRealization)] = getLIBOR(process, evaluationTimeIndexForRealizations, liborIndexForRealization);
				
				realizationsAtTimeIndex[liborIndexForRealization] = getLIBOR(process, evaluationTimeIndexForRealizations, liborIndexForRealization);
			}
			final RandomVariable[] factorLoadingAtTimeIndex = getCovarianceModel().getFactorLoading(sumTimeIndex, getDefaultableComponentIndex(liborIndex), realizationsAtTimeIndex);
			//o_{Li}(t)
			RandomVariable   integrand = getRandomVariableForConstant(0.0);
			for ( final RandomVariable oneFactor: factorLoadingAtTimeIndex)
			{
				integrand = integrand.add(oneFactor.squared());
			}
			integrand = integrand.div( (realizationsAtTimeIndex[liborIndex].mult(tenorDt).add(1.0)).squared() );
			if(stateSpace == StateSpace.LOGNORMAL)
			{
				integrand = integrand.mult( realizationsAtTimeIndex[liborIndex].squared() );
			}
			final double integralDt = 0.5 * (process.getTime(sumTimeIndex) - process.getTime(sumTimeIndex - 1));
			driftAdjustment = driftAdjustment.add( (integrand.add(previousIntegrand)).mult(integralDt) );
			previousIntegrand = integrand;
		}

		return driftAdjustment;
	}
	
}
