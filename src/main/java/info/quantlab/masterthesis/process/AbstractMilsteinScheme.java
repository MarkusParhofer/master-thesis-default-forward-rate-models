package info.quantlab.masterthesis.process;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.finmath.concurrency.FutureWrapper;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.montecarlo.process.MonteCarloProcessFromProcessModel;
import net.finmath.stochastic.RandomVariable;

public abstract class AbstractMilsteinScheme  extends MonteCarloProcessFromProcessModel implements MonteCarloProcess {

	private ExecutorService executor;

	private RandomVariable[][] m_Process;
	private RandomVariable[] m_ProcessWeights;
	private IndependentIncrements m_StochasticDriver;
	private boolean m_UseMultiThreadding = true;
	
	
	public AbstractMilsteinScheme(ProcessModel model, IndependentIncrements stochasticDriver) {
		super(stochasticDriver.getTimeDiscretization(), model);
		m_StochasticDriver = stochasticDriver;
		m_Process = null;
		m_ProcessWeights = null;
	}

	@Override
	public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		if(m_ProcessWeights == null)
			doPrecalculateProcess();
		
		return m_Process[timeIndex][componentIndex];
	}

	@Override
	public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
		if(m_ProcessWeights == null)
			doPrecalculateProcess();
		
		return m_ProcessWeights[timeIndex];
	}

	@Override
	public int getNumberOfPaths() {
		return getStochasticDriver().getNumberOfPaths();
	}

	@Override
	public int getNumberOfFactors() {
		return getModel().getNumberOfFactors();
	}

	@Override
	public IndependentIncrements getStochasticDriver() {
		return m_StochasticDriver;
	}

	protected RandomVariable[] getProcessValues(int timeIndex) {
		// To prevent Calculation Exception we don't precalculate the process if that is null
		if(m_Process == null)
			return null;
		
		return m_Process[timeIndex];
	}
	
	protected abstract RandomVariable[] getDifferentialOfFactorLoading(int componentIndex, int timeIndex, RandomVariable[] factorLoadings);
	
	private void doPrecalculateProcess() throws CalculationException {
		
		if (m_Process != null && m_Process.length != 0) {
			return;
		}

		final int numberOfPaths			= getNumberOfPaths();
		final int numberOfComponents	= getNumberOfComponents();

		// Allocate Memory
		m_Process	= new RandomVariable[getTimeDiscretization().getNumberOfTimeSteps() + 1][getNumberOfComponents()];
		m_ProcessWeights	= new RandomVariable[getTimeDiscretization().getNumberOfTimeSteps() + 1];

		// Set initial Monte-Carlo weights
		m_ProcessWeights[0] = m_StochasticDriver.getRandomVariableForConstant(1.0 / numberOfPaths);

		// Set initial value
		final RandomVariable[] initialState = getInitialState();
		final RandomVariable[] currentState = new RandomVariable[numberOfComponents];
		for (int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
			currentState[componentIndex] = initialState[componentIndex];
			m_Process[0][componentIndex] = applyStateSpaceTransform(0, componentIndex, currentState[componentIndex]);
		}

		/*
		 * Evolve the process using an Milstein scheme.
		 * The evolution is performed multi-threadded.
		 * Each component of the vector runs in its own thread.
		 */
		executor = Executors.newCachedThreadPool();

		// Evolve process
		for (int timeIndex2 = 1; timeIndex2 < getTimeDiscretization().getNumberOfTimeSteps()+1; timeIndex2++) {
			final int timeIndex = timeIndex2;
			// Generate process from timeIndex-1 to timeIndex
			final double deltaT = getTime(timeIndex) - getTime(timeIndex - 1);

			// Fetch drift vector
			final RandomVariable[] drift;
			try {
				drift = getDrift(timeIndex - 1, m_Process[timeIndex - 1], null);
			}
			catch(final Exception e) {
				throw new RuntimeException(e + " - drift calculaton failed at time index " + timeIndex + " (time=" + getTime(timeIndex - 1) + ") . See cause of this exception for details.", e);
			}

			// Fetch brownianIncrement vector
			final RandomVariable[] brownianIncrement	= m_StochasticDriver.getIncrement(timeIndex - 1);

			// Calculate new realization
			final ArrayList<Future<RandomVariable>> discreteProcessAtCurrentTimeIndex = new ArrayList<>(numberOfComponents);
			for (int componentIndex2 = 0; componentIndex2 < numberOfComponents; componentIndex2++) {
				final int componentIndex = componentIndex2;

				final RandomVariable	driftOfComponent	= drift[componentIndex];

				// Check if the component process has stopped to evolve
				if (driftOfComponent == null) {
					discreteProcessAtCurrentTimeIndex.add(componentIndex, null);
					continue;
				}

				final Callable<RandomVariable> worker = new  Callable<RandomVariable>() {
					@Override
					public RandomVariable call() {

						final RandomVariable[]	factorLoadings = getFactorLoading(timeIndex - 1, componentIndex, m_Process[timeIndex - 1]);

						// Check if the component process has stopped to evolve
						if (factorLoadings == null) {
							return null;
						}

						// Apply drift
						if(driftOfComponent != null) {
							currentState[componentIndex] = currentState[componentIndex].addProduct(driftOfComponent, deltaT); // mu DeltaT
						}

						// Apply diffusion
						currentState[componentIndex] = currentState[componentIndex].addSumProduct(factorLoadings, brownianIncrement); // sigma DeltaW

						// Add Milstein adjustment:
						RandomVariable[] milsteinAdjustment = getDifferentialOfFactorLoading(componentIndex, timeIndex - 1, factorLoadings);
						for(int k=0; k < getNumberOfFactors(); k++) {
							final RandomVariable milsteinIncrement = brownianIncrement[k].squared().sub(deltaT);
							milsteinAdjustment[k] = milsteinAdjustment[k].mult(factorLoadings[k]).mult(0.5d).mult(milsteinIncrement);
							currentState[componentIndex] = currentState[componentIndex].add(milsteinAdjustment[k]);
						}
						
						// Transform the state space to the value space and return it.
						return applyStateSpaceTransform(timeIndex, componentIndex, currentState[componentIndex]);
					}
				};


				/*
				 * Optional multi-threadding (asyncronous calculation of the components)
				 */
				Future<RandomVariable> result = null;
				try {
					if(m_UseMultiThreadding) {
						result = executor.submit(worker);
					} else {
						result = new FutureWrapper<>(worker.call());
					}
				} catch (final Exception e) {
					throw new RuntimeException("Euler step failed at time index " + timeIndex + " (time=" + getTime(timeIndex) + "). See cause of this exception for details.", e);
				}

				// The following line will add the result of the calculation to the vector discreteProcessAtCurrentTimeIndex
				discreteProcessAtCurrentTimeIndex.add(componentIndex, result);
			}

			// Fetch results and move to discreteProcess[timeIndex]
			for (int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
				try {
					final Future<RandomVariable> discreteProcessAtCurrentTimeIndexAndComponent = discreteProcessAtCurrentTimeIndex.get(componentIndex);
					if(discreteProcessAtCurrentTimeIndexAndComponent != null) {
						m_Process[timeIndex][componentIndex] = discreteProcessAtCurrentTimeIndexAndComponent.get().cache();
					} else {
						m_Process[timeIndex][componentIndex] = m_Process[timeIndex-1][componentIndex];
					}
				} catch (final InterruptedException | ExecutionException e) {
					throw new RuntimeException("Euler step failed at time index " + timeIndex + " (time=" + getTime(timeIndex) + "). See cause of this exception for details.", e.getCause());
				}
			}

			// Set Monte-Carlo weights
			m_ProcessWeights[timeIndex] = m_ProcessWeights[timeIndex - 1];
		} // End for(timeIndex)

		try {
			executor.shutdown();
		}
		catch(final SecurityException e) {
			// @TODO Improve exception handling here
		}
	}

	
	
}
