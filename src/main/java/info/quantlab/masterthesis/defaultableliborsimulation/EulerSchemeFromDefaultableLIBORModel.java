/**
 * 
 */
package info.quantlab.masterthesis.defaultableliborsimulation;

import java.util.Map;

import org.apache.commons.lang3.Validate;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.functional.Functional;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel.Scheme;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * @author Markus Parhofer
 *
 */
public class EulerSchemeFromDefaultableLIBORModel implements MonteCarloProcessWithDependency {

	private final IndependentIncrements _stochasticDriver;
	private final DefaultableLIBORMarketModel _model;
	
	private RandomVariable[][] _defaultableProcessValues;
	private final MonteCarloProcess _undefaultableProcess;
	
	private RandomVariable[] _processWeights;
	
	public EulerSchemeFromDefaultableLIBORModel(DefaultableLIBORMarketModel defaultableLIBORModel, IndependentIncrements stochasticDriver) {
		super();
		Validate.isTrue(defaultableLIBORModel.getNumberOfFactors() == stochasticDriver.getNumberOfFactors());
		_model = defaultableLIBORModel;
		_stochasticDriver = stochasticDriver;
		
		// Monte Carlo Valuation for undefaultable Model
		LIBORMarketModel underlyingModel = defaultableLIBORModel.getUnderlyingNonDefaultableModel();
		IndependentIncrements reducedFactorBM = Functional.getFirstFactors(stochasticDriver, underlyingModel.getNumberOfFactors());
		_undefaultableProcess = new EulerSchemeFromProcessModel(underlyingModel, reducedFactorBM);
	}
	
	@Override
	public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		if(_defaultableProcessValues == null) {
			precalculateProcess();
		}
		return _defaultableProcessValues[timeIndex][componentIndex];
	}

	private void precalculateProcess() {
		
		if (_defaultableProcessValues != null && _defaultableProcessValues.length != 0) {
			return;
		}

		final int numberOfPaths			= this.getNumberOfPaths();
		final int numberOfFactors		= this.getNumberOfFactors();
		final int numberOfComponents	= this.getNumberOfComponents();
		final TimeDiscretization times = _stochasticDriver.getTimeDiscretization();
		
		// Allocate Memory
		_defaultableProcessValues = new RandomVariable[times.getNumberOfTimes()][numberOfComponents];
		_processWeights	= new RandomVariable[times.getNumberOfTimes()];

		// Given Values
		_processWeights[0] = _stochasticDriver.getRandomVariableForConstant(1.0 / numberOfPaths);

		// Set initial value
		final RandomVariable[] currentState = _model.getInitialState(this);
		for (int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
			_defaultableProcessValues[0][componentIndex] = _model.applyStateSpaceTransform(this, 0, componentIndex, currentState[componentIndex]);
		}
		
		for (int timeIndex = 1; timeIndex < times.getNumberOfTimes(); timeIndex++) {
			
			// Generate process from timeIndex-1 to timeIndex
			final double deltaT = times.getTimeStep(timeIndex - 1);

			// Fetch drift vector
			final RandomVariable[] drift;
			try {
				drift = _model.getDrift(this, timeIndex - 1, _defaultableProcessValues[timeIndex - 1], null);
			}
			catch(final Exception e) {
				throw new RuntimeException(e + " - drift calculaton failed at time index " + timeIndex + " (time=" + getTime(timeIndex - 1) + ") . See cause of this exception for details.", e);
			}

			// Fetch brownianIncrement vector
			final RandomVariable[] brownianIncrement	= _stochasticDriver.getIncrement(timeIndex - 1);

			
			for(int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
				// TODO: Implement Scheme
				// if(scheme == Scheme.EULER_FUNCTIONAL || scheme == Scheme.PREDICTOR_CORRECTOR_FUNCTIONAL) {
				//	currentState[componentIndex] = applyStateSpaceTransformInverse(timeIndex - 1, componentIndex, discreteProcess[timeIndex - 1][componentIndex]);
				//}

				if (drift[componentIndex] == null) {
					// discreteProcessAtCurrentTimeIndex.add(componentIndex, null);
					continue;
				}

				final RandomVariable[] factorLoadings = _model.getFactorLoading(this, timeIndex - 1, componentIndex, _defaultableProcessValues[timeIndex - 1]);

				// Check if the component process has stopped to evolve
				if (factorLoadings == null) {
					// discreteProcessAtCurrentTimeIndex.add(componentIndex, null);
					continue;
				}

				// Apply drift
				currentState[componentIndex] = currentState[componentIndex].addProduct(drift[componentIndex], deltaT); // mu DeltaT
				

				// Apply diffusion
				currentState[componentIndex] = currentState[componentIndex].addSumProduct(factorLoadings, brownianIncrement); // sigma DeltaW

				// Transform the state space to the value space and save it.
				_defaultableProcessValues[timeIndex][componentIndex] = _model.applyStateSpaceTransform(this, timeIndex, componentIndex, currentState[componentIndex]);
				_processWeights[timeIndex] = _processWeights[timeIndex - 1];
			}
		}
	}
	
	@Override
	public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
		if(_processWeights == null) {
			precalculateProcess();
		}
		return _processWeights[timeIndex];
	}

	@Override
	public int getNumberOfComponents() {
		return _model.getNumberOfComponents();
	}

	@Override
	public TimeDiscretization getTimeDiscretization() {
		return _stochasticDriver.getTimeDiscretization();
	}

	@Override
	public double getTime(int timeIndex) {
		return _stochasticDriver.getTimeDiscretization().getTime(timeIndex);
	}

	@Override
	public int getTimeIndex(double time) {
		return _stochasticDriver.getTimeDiscretization().getTimeIndex(time);
	}

	@Override
	public int getNumberOfPaths() {
		return _stochasticDriver.getNumberOfPaths();
	}

	@Override
	public int getNumberOfFactors() {
		return _model.getNumberOfFactors();
	}

	@Override
	public IndependentIncrements getStochasticDriver() {
		return _stochasticDriver;
	}
	
	/**
	 * Returns the model that is simulated, which is a defaultable LIBOR market model.
	 * @see info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel DefaultableLIBORMarketModel.
	 */
	public DefaultableLIBORMarketModel getModel() {
		return _model;
	}
	
	/** 
	 * Returns the undefaultable process as MonteCarloProcess. Can be used to get Values of the undefaultable Process
	 * 
	 * @return Undefaultable Process
	 */
	public MonteCarloProcess getUndefaultableProcess() {
		return _undefaultableProcess;
	}

	/**
	 * Returns the process values of the Undefaultable model.
	 * @param timeIndex The time index of the value to retrieve.
	 * @param componentIndex The component index of the value to retrieve.
	 * @return The process value of the undefaultable model.
	 * @throws CalculationException
	 */
	public RandomVariable getUndefaultableProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		return _undefaultableProcess.getProcessValue(timeIndex, componentIndex);
	}
	
	@Override
	public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
		DefaultableLIBORMarketModel castedModel;
		try {
			castedModel = (DefaultableLIBORMarketModel)model;
		} catch(final Exception e) {
			e.printStackTrace();
			throw new ClassCastException("model must be of type \"DefaultableLIBORMarketModel\"");
		}
		return new EulerSchemeFromDefaultableLIBORModel(castedModel, _stochasticDriver);
	}

	@Override
	public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MonteCarloProcess clone() {
		return new EulerSchemeFromDefaultableLIBORModel(_model, _stochasticDriver);
	}

	@Override
	public ProcessModel getDependencyModel() {
		return _undefaultableProcess.getModel();
	}

	@Override
	public MonteCarloProcess getDependencyProcess() {
		return _undefaultableProcess;
	}

	@Override
	public RandomVariable[] getDependencyProcessValue(int timeIndex) throws CalculationException {
		return _undefaultableProcess.getProcessValue(timeIndex);
	}

	@Override
	public RandomVariable getDependencyProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		return _undefaultableProcess.getProcessValue(timeIndex, componentIndex);
	}
}
