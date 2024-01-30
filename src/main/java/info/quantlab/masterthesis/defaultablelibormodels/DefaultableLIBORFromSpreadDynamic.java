package info.quantlab.masterthesis.defaultablelibormodels;


import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.functional.FunctionsOnMCProcess;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.model.AbstractProcessModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;



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

	// TODO: Always adjust non-defaultable drift (add 0.5*sigma^2 and multiply L^d) if stateSpace of non-defaultable LIBORMarketModel is lognormal
	// TODO: Always adjust non-defaultable factor loading (multiply L) if non defaultable LIBORMarketModel is lognormal
	
	public enum Measure					{ SPOT, TERMINAL }
	public enum InterpolationMethod		{ LINEAR, LOG_LINEAR_UNCORRECTED, LOG_LINEAR_CORRECTED }
	public enum HandleSimulationTime 	{ ROUND_DOWN, ROUND_NEAREST }
	public enum StateSpace				{ 
		/**
		 * For now the only stable StateSpace is NORMAL, because if we model the LIBORS we need to adjust the non-defaultable Factor Loadings 
		 * accordingly (but we can't actually tell if the stateSpace of the non-defaultable model is lognormal) and if we model the spreads we need 
		 * to adjust both the drift and the factor loadings.
		 */
		NORMAL, 
		/**
		 * For now the only stable StateSpace is NORMAL, because if we model the LIBORS we need to adjust the non-defaultable Factor Loadings 
		 * accordingly (but we can't actually tell if the stateSpace of the non-defaultable model is lognormal) and if we model the spreads we need 
		 * to adjust both the drift and the factor loadings.
		 */
		@Deprecated LOGNORMAL } // No support for the LOGNORMAL version yet (drift must be adjusted)
	public enum SimulationModel			{ 
		/**
		 * For now the only stable simulation model is LIBORS, or modelling spread under NORMAL statespace, because deviding by the spread 
		 * (e.g. in getDriftOfSpread(...) we use the formula drift=(drift_D - drift_ND)/Spread) poses a stability problem, because the spread 
		 * might be very small
		 */
		SPREADS,
		/**
		 * For now the only stable simulation model is LIBORS, or modelling spread under NORMAL statespace, because deviding by the spread 
		 * (e.g. in getDriftOfSpread(...) we use the formula drift=(drift_D - drift_ND)/Spread) poses a stability problem, because the spread 
		 * might be very small
		 */
		LIBORS }
	
	// Flags
	private Measure	measure = Measure.SPOT;
	private StateSpace stateSpace = StateSpace.NORMAL;
	private HandleSimulationTime handleSimulationTime = HandleSimulationTime.ROUND_NEAREST;
	private InterpolationMethod	interpolationMethod	= InterpolationMethod.LOG_LINEAR_UNCORRECTED;
	
	private SimulationModel simulationModel = SimulationModel.LIBORS;
	private StateSpace stateSpaceOfSpread = StateSpace.NORMAL;
	private boolean nondefaultableModelIsLogNormal = false;
	
	private final LIBORMarketModel _nonDefaultableModel;
	private final DefaultableLIBORCovarianceModel _covarianceModel;
	private final ForwardCurve _forwardRateCurve;
	private final DiscountCurve _defaultableDiscountCurve;
	private final AnalyticModel _analyticModel;
	
	private double valueCap = 1E3d;
	private double valueFloor = Double.MIN_NORMAL;
	
	
	// Mutable/Transient Fields
	private transient Vector<RandomVariable> interpolationDriftAdjustmentsTerminal = new Vector<>();
	private transient MonteCarloProcess numerairesProcess;
	private transient double[][][] integratedLIBORCovariance;
	private transient TimeDiscretization integrationTimeDiscretizationCache;
	private transient Object integratedLIBORCovarianceLazyInitLock = new Object();


	// ------------------------------------------------------------------------- Constructors --------------------------------------------------------------------------
	
	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel nonDefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, AnalyticModel analyticModel, Map<String, String> propertyMap) {
		this(nonDefaultableModel, covarianceModel, forwardRateCurve, analyticModel);
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
						if(stateSpace == StateSpace.LOGNORMAL)
							System.out.println("<[Warning]: No support for the lognormal statespace yet!>");
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
				case "statespaceofspread":
					try {
						stateSpaceOfSpread = StateSpace.valueOf(propertyMap.get(key).toUpperCase());
					} catch(Exception ex) {
						System.out.println("<[Warning]:\n"
								+ "Source: Constructor of DefaultableLIBORMarketModelFromCovarianceModel\n"
								+ "Warning: Key \"simulationmodel\" ignored, value unknown. Value: " + propertyMap.get(key).toUpperCase() + "\n>");
					}
					if(stateSpaceOfSpread == StateSpace.LOGNORMAL) {
						stateSpaceOfSpread = covarianceModel.isSpreadModelLogNormal() ? StateSpace.LOGNORMAL : StateSpace.NORMAL;
						try {
							for(int i = 0; i < getNumberOfLIBORPeriods(); i++) {
								final double spread = getDefaultableLIBOR(null,  0, i).sub(getNonDefaultableLIBOR(null, 0, i)).doubleValue();
								if(spread == 0.0) {
									// Cannot have lognormal spread dynamic if starting Value of Spread is zero
									stateSpaceOfSpread = StateSpace.NORMAL;
									break;
								}
							}
						} catch(CalculationException ex) {
							stateSpaceOfSpread = StateSpace.NORMAL;
						}
					}
					if(!(stateSpaceOfSpread == StateSpace.LOGNORMAL))
						valueFloor = Double.NEGATIVE_INFINITY;
					break;
				default:
					continue;

				}
			}
		}
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel nonDefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, Map<String, String> propertyMap) {
		this(nonDefaultableModel, covarianceModel, forwardRateCurve, null, propertyMap);
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel nonDefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve, AnalyticModel analyticModel) {
		_nonDefaultableModel = nonDefaultableModel;
		_covarianceModel = covarianceModel;
		_forwardRateCurve = forwardRateCurve;
		_defaultableDiscountCurve = new DiscountCurveFromForwardCurve(forwardRateCurve);
		_analyticModel = analyticModel;
		// Test for Lognormality of underlying LIBOR model:
		
	}

	public DefaultableLIBORFromSpreadDynamic(LIBORMarketModel nonDefaultableModel, DefaultableLIBORCovarianceModel covarianceModel, ForwardCurve forwardRateCurve) {
		this(nonDefaultableModel, covarianceModel, forwardRateCurve, null, null);
	}
	
	
	
	// ------------------------------------------------------------------------ Getters --------------------------------------------------------------------------------
	
	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return getNonDefaultableLIBORModel().getLiborPeriodDiscretization();
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
		return getNonDefaultableLIBORModel().getDiscountCurve();
	}

	/**
	 * Returns the initial forward rate curve of the defaultable model.
	 */
	@Override
	public ForwardCurve getForwardRateCurve() {
		return _forwardRateCurve;
	}

	@Override
	public LIBORMarketModel getNonDefaultableLIBORModel() {
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
			return getNonDefaultableLIBORModel().getLiborPeriod(componentIndex - getNumberOfLIBORPeriods());
		else
			return getNonDefaultableLIBORModel().getLiborPeriod(componentIndex);
	}

	@Override
	public int getLiborPeriodIndex(double time) {
		return getNonDefaultableLIBORModel().getLiborPeriodIndex(time);
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

		final RandomVariable[] initialStateRV = Arrays.copyOf(getNonDefaultableLIBORModel().getInitialState(null), getNumberOfComponents());
		switch(simulationModel) {
		case SPREADS:
			for(int liborIndex = 0; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				initialStateRV[getSpreadComponentIndex(liborIndex)] = initialStateRV[liborIndex].bus(initialStatesD[liborIndex]);
				if(stateSpaceOfSpread == StateSpace.LOGNORMAL)
					initialStateRV[getSpreadComponentIndex(liborIndex)] = initialStateRV[getSpreadComponentIndex(liborIndex)].log();
			}
			break;
		case LIBORS:
			for(int liborIndex = 0; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				initialStateRV[getDefaultableComponentIndex(liborIndex)] = getRandomVariableForConstant(initialStatesD[liborIndex]);
				if(stateSpace == StateSpace.LOGNORMAL)
					initialStateRV[getDefaultableComponentIndex(liborIndex)] = initialStateRV[getDefaultableComponentIndex(liborIndex)].log();
			}
			break;
		default:
			throw new IllegalArgumentException("Method not implemented for specified simulation model.");
		}
		return initialStateRV;
	}

	@Override
	public RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor) {
		final RandomVariable[] nonDefaultableDriftVector = getDriftOfNonDefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
		RandomVariable[] fullDriftVector = Arrays.copyOf(nonDefaultableDriftVector, getNumberOfComponents());
		
		switch(simulationModel) {
		case SPREADS:
			RandomVariable[] spreadDrift = getDriftOfSpread(process, timeIndex, realizationAtTimeIndex, realizationPredictor, nonDefaultableDriftVector);
			for(int i=getNumberOfLIBORPeriods(); i < getNumberOfComponents(); i++) {
				fullDriftVector[i] = spreadDrift[i - getNumberOfLIBORPeriods()];
			}
			break;
		case LIBORS:			
			RandomVariable[] defaultableDrift = getDriftOfDefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
			for(int i=getNumberOfLIBORPeriods(); i < getNumberOfComponents(); i++) {
				fullDriftVector[i] = defaultableDrift[i - getNumberOfLIBORPeriods()];
			}
			break;
		default:
			throw new IllegalArgumentException("Method not implemented for specified simulation model.");
		}
		return fullDriftVector;
	}

	@Override
	public RandomVariable[] getDriftOfNonDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		RandomVariable[] sRealization = null;
		if(realizationAtTimeIndex != null) {
			sRealization = Arrays.copyOf(realizationAtTimeIndex, getNonDefaultableLIBORModel().getNumberOfComponents());
		}
		RandomVariable[] sRealizationPredictor = null; 
		if(realizationPredictor != null) {
			sRealizationPredictor = Arrays.copyOf(realizationPredictor, getNonDefaultableLIBORModel().getNumberOfComponents());
		}
		return getNonDefaultableLIBORModel().getDrift(getNonDefaultableProcess(process), timeIndex, sRealization, sRealizationPredictor);
	}

	@Override
	public RandomVariable[] getDriftOfDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {
		if(simulationModel != SimulationModel.LIBORS) {
			throw new UnsupportedOperationException("This method only makes sense for Models modeling the defaultable LIBOR rates!");
		}
		RandomVariable[] drift = getDriftDefaultableInternally(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
		if(stateSpace == StateSpace.LOGNORMAL) {
			// Drift adjustment for log-coordinate in each component
			for(int componentIndex = 0; componentIndex < getNumberOfLIBORPeriods(); componentIndex++) {
				if(drift[componentIndex] == null)
					continue;
				final RandomVariable variance = getCovarianceModel().getCovariance(process.getTime(timeIndex), componentIndex + getNumberOfLIBORPeriods(), componentIndex + getNumberOfLIBORPeriods(), realizationAtTimeIndex);
				drift[componentIndex] = drift[componentIndex].addProduct(variance, -0.5);
			}
		}
		return drift;
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
			nonDefaultableDrifts = getDriftOfNonDefaultableModel(process, timeIndex, realizationAtTimeIndex, realizationPredictor);
	
		RandomVariable[] realizationAllLIBORS = null;
		if(realizationAtTimeIndex != null) {
			realizationAllLIBORS = new RandomVariable[getNumberOfComponents()];
			for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
				realizationAllLIBORS[libor] = realizationAtTimeIndex[libor];
				realizationAllLIBORS[getDefaultableComponentIndex(libor)] = realizationAtTimeIndex[libor].add(realizationAtTimeIndex[getSpreadComponentIndex(libor)]);
				
			}
		}
		
		RandomVariable[] realizationPredAllLIBORS = null;
		if(realizationPredictor != null) {
			realizationPredAllLIBORS = new RandomVariable[getNumberOfComponents()];
			for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
				realizationPredAllLIBORS[libor] = realizationPredictor[libor];
				realizationPredAllLIBORS[getDefaultableComponentIndex(libor)] = realizationPredictor[libor].add(realizationPredictor[getSpreadComponentIndex(libor)]);
				
			}
		}
		
		drift = getDriftDefaultableInternally(process, timeIndex, realizationAllLIBORS, realizationPredAllLIBORS);
				
		// TODO Add functionality of lognormal. I.e. add 0.5*sigma^2 to both drifts and multiply with realization
		
		
		
		
		for(int libor=0; libor < getNumberOfLIBORPeriods(); libor++) {
			if(drift[libor] == null || nonDefaultableDrifts[libor] == null) {
				drift[libor] = null;
				continue;
			}
			if(nondefaultableModelIsLogNormal) {
				// Add 0.5*sigma^2 to nonDefaultable drift
				RandomVariable[] nonDefaultableRealizations = Arrays.copyOf(realizationAtTimeIndex, getNonDefaultableLIBORModel().getNumberOfLibors());
				nonDefaultableDrifts[libor] = nonDefaultableDrifts[libor].add(getNonDefaultableLIBORModel().getCovarianceModel().getCovariance(timeIndex, libor, libor, nonDefaultableRealizations).mult(0.5));
				// Multiply L to nonDefaultable drift
				nonDefaultableDrifts[libor] = nonDefaultableDrifts[libor].mult(realizationAtTimeIndex[libor]);
			}
			if(stateSpace == StateSpace.LOGNORMAL) {
				// Multiply L^d to drift (getDefaultableDriftInternally(...) gets purely mu)
				drift[libor] = drift[libor].mult(realizationAllLIBORS[getDefaultableComponentIndex(libor)]);
			}
			drift[libor] = drift[libor].sub(nonDefaultableDrifts[libor]);
			if(stateSpaceOfSpread == StateSpace.LOGNORMAL) {
				drift[libor] = drift[libor].div(realizationAtTimeIndex[getSpreadComponentIndex(libor)]);
				RandomVariable[] factorLoadingSpread = getFactorLoadingOfSpread(process, timeIndex, libor, realizationAtTimeIndex);
				for(int factor = 0; factor < getNumberOfFactors(); factor++) {
					drift[libor] = drift[libor].sub(factorLoadingSpread[factor].squared().mult(0.5));
				}
			}
		}
		
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
		RandomVariable[] realizationAllLIBORs = Arrays.copyOf(realizationAtTimeIndex, getNumberOfComponents());
		// Add Spread to non-defaultable LIBORs to get defaultable LIBORS
		for(int component = getNumberOfLIBORPeriods(); component < getNumberOfComponents(); component++) {
			realizationAllLIBORs[component] = realizationAtTimeIndex[component - getNumberOfLIBORPeriods()].add(realizationAtTimeIndex[component]);
		}
		RandomVariable[] factorLoading = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), realizationAllLIBORs);
		RandomVariable[] factorLoadingNonDef = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), componentIndex, realizationAllLIBORs);
		
		
		for(int i =0; i < getNumberOfFactors(); i++) {
			if(stateSpace == StateSpace.LOGNORMAL) {
				// Multiply L^d to factorLoading
				factorLoading[i] = factorLoading[i].mult(realizationAllLIBORs[getDefaultableComponentIndex(componentIndex)]);
			}
			if(nondefaultableModelIsLogNormal) {
				// Multiply L to nonDefaultable factorLoading
				factorLoadingNonDef[i] = factorLoadingNonDef[i].mult(realizationAllLIBORs[componentIndex]);
			}
			factorLoading[i] = factorLoading[i].sub(factorLoadingNonDef[i]);
			if(stateSpaceOfSpread == StateSpace.LOGNORMAL) {
				factorLoading[i] = factorLoading[i].div(realizationAtTimeIndex[getSpreadComponentIndex(componentIndex)]);
			}
		}
		/*
		try {
			factorLoading = getCovarianceModel().getFactorLoadingOfSpread(process.getTime(timeIndex), componentIndex, realizationForCovariance);
		} catch(UnsupportedOperationException ex) {
			// Sub defaultable FactorLoading from factorLoading and divide by Spread
			factorLoading = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), getDefaultableComponentIndex(componentIndex), realizationForCovariance);
			
			RandomVariable[] factorLoadingNonDef = getCovarianceModel().getFactorLoading(process.getTime(timeIndex), componentIndex, realizationForCovariance);
			for(int i =0; i < getNumberOfFactors(); i++) {
				factorLoading[i] = factorLoading[i].sub(factorLoadingNonDef[i]);
				if(spreadIsLogNormal)
					factorLoading[i] = factorLoading[i].div(realizationAtTimeIndex[getSpreadComponentIndex(componentIndex)]);
			}
		}*/
		
		return factorLoading;
	}
	
	@Override
	public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		if(componentIndex < getNumberOfLIBORPeriods()) {
			return getNonDefaultableLIBORModel().applyStateSpaceTransform(getNonDefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		}
		
		RandomVariable value = randomVariable;

		if((simulationModel == SimulationModel.SPREADS && stateSpaceOfSpread == StateSpace.LOGNORMAL) || 
				(simulationModel == SimulationModel.LIBORS && stateSpace == StateSpace.LOGNORMAL)) {
			value = value.exp(); //.floor(1E-5);
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
			return getNonDefaultableLIBORModel().applyStateSpaceTransformInverse(getNonDefaultableProcess(process), timeIndex, componentIndex, randomVariable);
		}
		
		RandomVariable value = randomVariable;
		
		if((simulationModel == SimulationModel.SPREADS && stateSpaceOfSpread == StateSpace.LOGNORMAL) || 
				(simulationModel == SimulationModel.LIBORS && stateSpace == StateSpace.LOGNORMAL)) {
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
		if(periodStart == periodEnd)
			return Scalar.of(0.0);
		
		final int periodStartIndex    = getLiborPeriodIndex(periodStart);
		final int periodEndIndex      = getLiborPeriodIndex(periodEnd);

		// If time is beyond fixing, use the fixing time.
		time = Math.min(time, periodStart);
		int timeIndex = 0;
		if(time != 0.0)
			timeIndex = process.getTimeIndex(time);
			
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
		// time might be out of bounds, because we only evaluate until one Tenor Time before maturity
		// Determine if maturity is on Tenor grid:
		if(time >= maturity){
			return new Scalar(1.0);
		}
		int nextTenorIndex = getLiborPeriodIndex(maturity);
		if(nextTenorIndex >= 0) {
			double tenorStep = getLiborPeriod(nextTenorIndex) - getLiborPeriod(nextTenorIndex - 1);
			double fixingTime = maturity - time > tenorStep? time : getLiborPeriod(nextTenorIndex - 1);
			RandomVariable forwardRate = getForwardRate(process, fixingTime, fixingTime, maturity);
			return (new Scalar(1.0)).discount(forwardRate, maturity - time);
		}
		else {
			nextTenorIndex = - nextTenorIndex - 1;
			int nextTimeTenorIndex = getLiborPeriodIndex(time);
			if(nextTimeTenorIndex < 0 && (- nextTimeTenorIndex - 1) == nextTenorIndex) {
				// Handle case T_j < time < maturity < T_j+1 (short period bond)
				return getDefaultableBond(process, time, getLiborPeriod(nextTenorIndex))
						.div(getDefaultableBond(process, maturity, getLiborPeriod(nextTenorIndex)));
			} else {
				// Handle case T_i <= time =< T_i+1; T_j < maturity < T_j+1
				RandomVariable forwardRate = getForwardRate(process, time, time, maturity);
				return (new Scalar(1.0)).discount(forwardRate, maturity - time);
			}
		}
	}
	
	
	@Override
	public RandomVariable getNonDefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		if(timeIndex == 0 || getLiborPeriod(liborIndex) == 0d)
			return getRandomVariableForConstant(getNonDefaultableLIBORModel().getForwardRateCurve().getForward(_analyticModel, getLiborPeriod(liborIndex)));
		return process.getProcessValue(timeIndex, liborIndex);
	}
	

	@Override
	public RandomVariable getNonDefaultableForwardRate(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getNonDefaultableLIBORModel().getForwardRate(getNonDefaultableProcess(process), time, periodStart, periodEnd);
	}

	
	@Override
	public RandomVariable getForwardDiscountBond(MonteCarloProcess process, final double time, final double maturity) throws CalculationException {
		return getNonDefaultableLIBORModel().getForwardDiscountBond(getNonDefaultableProcess(process), time, maturity);
	}
	
	
	
	// --------------------------------------------------------------------- Spread Methods ------------------------------------------------------------------------

	@Override
	public RandomVariable getLIBORSpreadAtGivenTimeIndex(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException {
		if(timeIndex == 0 || getLiborPeriod(liborIndex) == 0d) {
			return getDefaultableLIBOR(process, timeIndex, liborIndex).sub(getNonDefaultableLIBOR(process, timeIndex, liborIndex));
		}
		
		if(simulationModel == SimulationModel.SPREADS) 
			return process.getProcessValue(timeIndex, getSpreadComponentIndex(liborIndex));
		else
			return getLIBOR(process, timeIndex, liborIndex).sub(getNonDefaultableLIBOR(process, timeIndex, liborIndex));
	}
	
	
	@Override
	public RandomVariable getSpread(MonteCarloProcess process, double time, double periodStart, double periodEnd) throws CalculationException {
		return getForwardRate(process, time, periodStart, periodEnd).sub(getNonDefaultableForwardRate(process, time, periodStart, periodEnd));
	}

	
	@Override
	public RandomVariable getSurvivalProbability(MonteCarloProcess process, double maturity) throws CalculationException {
		// We calculate the Survival Probability per path. We do not measure default other than by its probability
		RandomVariable defaultableNumeraire = getDefaultableNumeraire(process, maturity);
		RandomVariable normalNumeraire = getNumeraire(process, maturity);
		return normalNumeraire.div(defaultableNumeraire);
	}

	
	
	// -------------------------------------------------------------------- Clones --------------------------------------------------------------------
	
	@Override
	public DefaultableLIBORMarketModelFromCovarianceModel getCloneWithModifiedNonDefaultableModel(LIBORMarketModel newNonDefaultableModel) {
		Map<String, String> propertyMap = new HashMap<String, String>();
		propertyMap.put("measure", measure.toString());
		propertyMap.put("statespace", stateSpace.toString());
		propertyMap.put("handlesimulationtime", handleSimulationTime.toString());
		propertyMap.put("interpolationmethod", interpolationMethod.toString());
		
		DefaultableLIBORCovarianceModel newCovarianceModel = getCovarianceModel();
		if(newNonDefaultableModel.getCovarianceModel() != newCovarianceModel.getNonDefaultableCovarianceModel()) {
			newCovarianceModel = newCovarianceModel.getCloneWithModifiedNonDefaultableCovariance(newNonDefaultableModel.getCovarianceModel());
		}
		return new DefaultableLIBORMarketModelFromCovarianceModel(newNonDefaultableModel, newCovarianceModel, getForwardRateCurve(), getAnalyticModel(), propertyMap);
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
		return new DefaultableLIBORMarketModelFromCovarianceModel(getNonDefaultableLIBORModel(), (DefaultableLIBORCovarianceModel)newCovarianceModel, getForwardRateCurve(), getAnalyticModel(), propertyMap);
	}

	@Override
	public DefaultableLIBORMarketModelFromCovarianceModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		LIBORMarketModel nonDefaultableModel = (LIBORMarketModel)dataModified.getOrDefault("nonDefaultableModel", getNonDefaultableLIBORModel());
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
		return new DefaultableLIBORMarketModelFromCovarianceModel(nonDefaultableModel, covarianceModel, forwardRateCurve, analyticModel, propertyMap);
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
									// Factor Loading definately needs Initial values of the form (L_0, ..., L_n, L^d_0, ..., L^d_n)
									final double time = timeDiscretization.getTime(timeIndex);
									try {
										initialLIBORValues = getLiborApproximation(time, false);
									} catch (CalculationException e) {
										e.printStackTrace();
										throw new UnsupportedOperationException("Cannot get Initial LIBOR rates as Estimators!");
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
		return getNonDefaultableLIBORModel().getNumeraire(getNonDefaultableProcess(process), time);
	}
	
	@Override
	public RandomVariable getDefaultableNumeraire(final MonteCarloProcess process, final double time) throws CalculationException {
		RandomVariable numeraire = getDefaultableNumeraireUnadjusted(process, time);
		
		// Adjust:
		// final RandomVariable zeroBondAtTimeZero =  getNumeraireDefaultableZeroBondAsOfTimeZero(process, time);
		// final double expectedZeroBond = numeraire.invert().mult(getDefaultableNumeraireUnadjusted(process, 0.0)).getAverage();
		// numeraire = numeraire.mult(expectedZeroBond).div(zeroBondAtTimeZero);
		
		return numeraire;
	}
	
	
	public RandomVariable getDefaultableNumeraireUnadjusted(final MonteCarloProcess process, final double time) throws CalculationException {
		
		// Check if time is on Libor Period Tenor:
		final int liborTimeIndex = getLiborPeriodIndex(time);
		if(liborTimeIndex >= 0)
			return getDefaultableNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex);
		
		// Interpolate:
		
		RandomVariable numeraireUnadjusted;
		final int upperIndex = - liborTimeIndex - 1;
		final int lowerIndex = upperIndex - 1;
		if (lowerIndex < 0) {
			throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");
		}
		if (measure == Measure.TERMINAL) {
			/*
			 * Due to time < T_{timeIndex+1} loop is needed.
			 */
			numeraireUnadjusted = getRandomVariableForConstant(1.0);
			
			int timeIndex = process.getTimeIndex(time); // We don't need Math.min(...) because time is always smaller than getLiborPeriod(liborIndex)
			if(timeIndex < 0)
				timeIndex = - timeIndex - 2;
				
			for (int liborIndex = upperIndex; liborIndex <= getNumberOfLIBORPeriods() - 1; liborIndex++) {
				final RandomVariable forwardRate = getDefaultableLIBOR(process, timeIndex, liborIndex); 
				final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborIndex);
				numeraireUnadjusted = numeraireUnadjusted.discount(forwardRate, periodLength);
			}
		}
		else if (measure == Measure.SPOT) {
			numeraireUnadjusted = getDefaultableNumeraireUnadjustedAtLIBORIndex(process, upperIndex);
		}
		else {
			throw new IllegalArgumentException("Numeraire not implemented for specified measure.");
		}
		/*
		 * Multiply with short period bond
		 */
		numeraireUnadjusted = numeraireUnadjusted.discount(getForwardRate(process, time, time, getLiborPeriod(upperIndex)), getLiborPeriod(upperIndex) - time);

		return numeraireUnadjusted;
	}
	
	public RandomVariable getDefaultableNumeraireUnadjustedAtLIBORIndex(final MonteCarloProcess process, final int liborTimeIndex) throws CalculationException {
		
		// Prevent Resource allocation for easiest cases:
		if((liborTimeIndex == 0 && measure == Measure.SPOT) || 
				(liborTimeIndex == getNumberOfLIBORPeriods() - 1 && measure == Measure.TERMINAL)) {
			return getRandomVariableForConstant(1.0);
		}
		
		
		RandomVariable numeraireUnadjusted;
		int timeIndex = process.getTimeIndex(getLiborPeriod(liborTimeIndex));
		if(timeIndex < 0) {
			timeIndex = - timeIndex - 1; // Originally - 1 but we want to get prev. time index
		}
		
		if (measure == Measure.TERMINAL) {
			// Initialize to 1.0
			numeraireUnadjusted = getRandomVariableForConstant(1.0);

			/*
			 * Due to time < T_{timeIndex+1} loop is needed.
			 */
			for (int liborIndex = liborTimeIndex; liborIndex < getNumberOfLIBORPeriods(); liborIndex++) {
				final RandomVariable forwardRate = getDefaultableLIBOR(process, timeIndex, liborIndex);
				final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborIndex);
				numeraireUnadjusted = numeraireUnadjusted.discount(forwardRate, periodLength);
			}
		}
		else if (measure == Measure.SPOT) {
			/*
			 * If numeraire is not N(0), multiply (1 + L(Ti-1)*dt) on N(Ti-1)
			 */
			if (liborTimeIndex != 0) {
				final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborTimeIndex - 1);
				final RandomVariable forwardRate = getDefaultableLIBOR(process, timeIndex, liborTimeIndex - 1);
				// Recursive Function. Maybe better to just do for loop.
				numeraireUnadjusted = getDefaultableNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex - 1).accrue(forwardRate, periodLength);
			}
			else {
				numeraireUnadjusted = getRandomVariableForConstant(1.0);
			}
		} else {
			throw new IllegalArgumentException("Numeraire not implemented for specified measure.");
		}
		return numeraireUnadjusted;
	}
	
	
	
	// --------------------------------------------------------------------- Other Public Methods ----------------------------------------------------------------
	
	@Override
	public LocalDateTime getReferenceDate() {
		return getNonDefaultableLIBORModel().getReferenceDate();
	}

	@Override
	public RandomVariable getRandomVariableForConstant(double value) {
		return getNonDefaultableLIBORModel().getRandomVariableForConstant(value);
	}



	// -------------------------------------------------------------------- Private Methods: --------------------------------------------------------------------------

	private RandomVariable[] getLiborApproximation(double atTime, boolean withDriftApprox) throws CalculationException {
		// Get Approximation
		RandomVariable[] liborApprox = new RandomVariable[getNumberOfComponents()];
		// Factor Loading definately needs Initial values of the form (L_0, ..., L_n, L^d_0, ..., L^d_n)
		for(int i=0; i < getNumberOfComponents(); i++) {
			liborApprox[i] = i < getNumberOfLIBORPeriods() ? getNonDefaultableLIBOR(null, 0, i) : getDefaultableLIBOR(null, 0, i);
		}
		
		if(!withDriftApprox)
			return liborApprox;
		
		final TimeDiscretization timeDisc = new TimeDiscretizationFromArray(0.0, atTime);
		final int numberOfComponents = getNumberOfComponents();
		final int numberOfFactors = getNumberOfFactors();
		final MonteCarloProcess mcCallBackForTimeDiscretization = new MonteCarloProcess() {
			
			@Override
			public int getTimeIndex(double time) {
				return getTimeDiscretization().getTimeIndex(time);
			}
			
			@Override
			public TimeDiscretization getTimeDiscretization() {
				return timeDisc;
			}
			
			@Override
			public double getTime(int timeIndex) {
				return getTimeDiscretization().getTime(timeIndex);
			}
			
			@Override
			public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public int getNumberOfComponents() {
				return numberOfComponents;
			}
			
			@Override
			public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public IndependentIncrements getStochasticDriver() {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public int getNumberOfPaths() {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public int getNumberOfFactors() {
				return numberOfFactors;
			}
			
			@Override
			public MonteCarloProcess clone() {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
			
			@Override
			public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
				throw new UnsupportedOperationException("This MC Process is purely for callbacks to the timeDiscretization");
			}
		};
		
		RandomVariable[] nonDefDrift = getDriftOfNonDefaultableModel(mcCallBackForTimeDiscretization, 0, liborApprox, null);
		RandomVariable[] defDrift = getDriftDefaultableInternally(mcCallBackForTimeDiscretization, 0, liborApprox, null);
		
		// Apply Drift Approximation L_0 + mu_0 * (t - t_0):
		for(int i=0; i < getNumberOfLIBORPeriods(); i++) {
			if(nonDefDrift[i] != null)
				liborApprox[i] = liborApprox[i].add(nonDefDrift[i].mult(atTime));
			if(defDrift[i] != null)
				liborApprox[getDefaultableComponentIndex(i)] = liborApprox[getDefaultableComponentIndex(i)].add(defDrift[i].mult(atTime));
		}
		return liborApprox;
	}
	
	private int getSpreadComponentIndex(int liborIndex) {
		return liborIndex + getNumberOfLIBORPeriods();
	}
	
	private int getDefaultableComponentIndex(int liborIndex) {
		return liborIndex + getNumberOfLIBORPeriods();
	}
		
	private MonteCarloProcess getNonDefaultableProcess(MonteCarloProcess fullProcess) {
		if(fullProcess == null)
			return null;
		return FunctionsOnMCProcess.getComponentReducedMCProcess(fullProcess, 0, getNumberOfLIBORPeriods() - 1);
		
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
		int	firstForwardRateIndex = this.getLiborPeriodIndex(time);		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = - firstForwardRateIndex - 2;
		}
		firstForwardRateIndex += 1;
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
			// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
			for(int liborIndex=getNumberOfLIBORPeriods()-1; liborIndex>=firstForwardRateIndex; liborIndex--) {
				final double			periodLength	= getLiborPeriodDiscretization().getTimeStep(liborIndex);
				final RandomVariable	forwardRate		= realizationAtTimeIndex[getDefaultableComponentIndex(liborIndex)];
				RandomVariable			oneStepMeasureTransform = Scalar.of(-periodLength).discount(forwardRate, periodLength);

				if(stateSpace == StateSpace.LOGNORMAL) {
					oneStepMeasureTransform = oneStepMeasureTransform.mult(forwardRate);
				}

				final RandomVariable[]	factorLoading   	= getFactorLoadingOfDefaultableModel(process, timeIndex, liborIndex, realizationAtTimeIndex);
				drift[liborIndex] = drift[liborIndex].addSumProduct(factorLoadingsSums, factorLoading);
				for(int factorIndex=0; factorIndex<factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
			}
		}
		else {
			throw new IllegalArgumentException("Drift not implemented for specified measure.");
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
			realizationsForFL[liborIndexForRealization] = getNonDefaultableLIBOR(process, 0, liborIndexForRealization);
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
				realizationsForFL[liborIndexForRealization] = getNonDefaultableLIBOR(process, evaluationTimeIndexForRealizations, liborIndexForRealization);
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
	
	private RandomVariable getNumeraireDefaultableZeroBondAsOfTimeZero(final MonteCarloProcess process, final double maturity) {
		final boolean interpolateDFsOnLiborPeriodDiscretization = true;

		final TimeDiscretization timeDiscretizationForCurves = interpolateDFsOnLiborPeriodDiscretization ? getLiborPeriodDiscretization() : process.getTimeDiscretization();

		final int timeIndex = timeDiscretizationForCurves.getTimeIndex(maturity);
		if(timeIndex >= 0) {
			return getNumeraireDefaultableZeroBondAsOfTimeZero(process, timeIndex);
		}
		else {
			// Interpolation
			final int timeIndexPrev = Math.min(-timeIndex-2, getLiborPeriodDiscretization().getNumberOfTimes()-2);
			final int timeIndexNext = timeIndexPrev+1;
			final double timePrev = timeDiscretizationForCurves.getTime(timeIndexPrev);
			final double timeNext = timeDiscretizationForCurves.getTime(timeIndexNext);
			final RandomVariable numeraireAdjustmentPrev = getNumeraireDefaultableZeroBondAsOfTimeZero(process, timeIndexPrev);
			final RandomVariable numeraireAdjustmentNext = getNumeraireDefaultableZeroBondAsOfTimeZero(process, timeIndexNext);
			return numeraireAdjustmentPrev.mult(numeraireAdjustmentNext.div(numeraireAdjustmentPrev).pow((maturity-timePrev)/(timeNext-timePrev)));
		}
	}

	private RandomVariable getNumeraireDefaultableZeroBondAsOfTimeZero(final MonteCarloProcess process, final int maturityTimeIndex) {
		final boolean interpolateDFsOnLiborPeriodDiscretization = true;

		final TimeDiscretization timeDiscretizationForCurves = interpolateDFsOnLiborPeriodDiscretization ? getLiborPeriodDiscretization() : process.getTimeDiscretization();
		final double time = timeDiscretizationForCurves.getTime(maturityTimeIndex);
		
		RandomVariable deterministicNumeraireAdjustment = null;
		
		final double dfInitial = _defaultableDiscountCurve.getDiscountFactor(null, timeDiscretizationForCurves.getTime(0));
		deterministicNumeraireAdjustment = getRandomVariableForConstant(dfInitial);
		if(time == 0.0)
			return deterministicNumeraireAdjustment;
		for(int i=0; i<timeDiscretizationForCurves.getNumberOfTimeSteps(); i++) {
			final double dfPrev = _defaultableDiscountCurve.getDiscountFactor(null, timeDiscretizationForCurves.getTime(i));
			final double dfNext = _defaultableDiscountCurve.getDiscountFactor(null, timeDiscretizationForCurves.getTime(i+1));
			final double timeStep = timeDiscretizationForCurves.getTimeStep(i);
			final double timeNext = timeDiscretizationForCurves.getTime(i+1);
			final RandomVariable forwardRate = getRandomVariableForConstant((dfPrev / dfNext - 1.0) / timeStep);
			deterministicNumeraireAdjustment = deterministicNumeraireAdjustment.discount(forwardRate, timeStep);
			if(timeNext == time)
				return deterministicNumeraireAdjustment;
		}
		return deterministicNumeraireAdjustment;
	}
	
}
