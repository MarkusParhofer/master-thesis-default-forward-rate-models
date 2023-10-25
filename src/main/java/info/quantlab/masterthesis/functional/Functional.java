/**
 * 
 */
package info.quantlab.masterthesis.functional;

import java.util.Map;

import org.apache.commons.lang3.Validate;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This class is used for some functional methods. They are all helpers for the master thesis.
 * 
 * @author Markus Parhofer
 * @version 1.0
 */
public class Functional {

	/**
	 * No construction of this class allowed. All static methods.
	 */
	private Functional() {}
	
	/**
	 * Creates a new instance of IndependentIncrements that is basically a copy of increments,
	 * but uses only the range factors that is specified.
	 * @param increments The original instance of independent increments
	 * @param firstFactor The first factor to consider in the new instance
	 * @param lastFactor The last factor to consider in the new instance
	 * @return A new instance with reduced factors
	 */
	public static IndependentIncrements getRangeOfFactors(IndependentIncrements increments, int firstFactor, int lastFactor) {
		Validate.isTrue(0 <= firstFactor && firstFactor < lastFactor && lastFactor < increments.getNumberOfFactors(), 
				"The following does not hold: 0 <= firstFactor (=%d) < lastFactor (=%d) < increments.numberOfFactors (=%d)",
				firstFactor, lastFactor, increments.getNumberOfFactors());

		if(firstFactor == 0 && increments.getNumberOfFactors() == lastFactor + 1) {
			return increments;			
		}
		return new IndependentIncrements() {
			private final int numberOfFactors = lastFactor - firstFactor + 1;
			
			@Override
			public RandomVariable getIncrement(int timeIndex, int factor) {
				Validate.isTrue(factor <= lastFactor - firstFactor);
				return increments.getIncrement(timeIndex, factor + firstFactor);
			}

			@Override
			public TimeDiscretization getTimeDiscretization() {
				return increments.getTimeDiscretization();
			}

			@Override
			public int getNumberOfFactors() {
				return numberOfFactors;
			}

			@Override
			public int getNumberOfPaths() {
				return increments.getNumberOfPaths();
			}

			@Override
			public RandomVariable getRandomVariableForConstant(double value) {
				return increments.getRandomVariableForConstant(value);
			}

			@Override
			public IndependentIncrements getCloneWithModifiedSeed(int seed) {
				IndependentIncrements originClone = increments.getCloneWithModifiedSeed(seed);
				return getRangeOfFactors(originClone, firstFactor, lastFactor);
			}

			@Override
			public IndependentIncrements getCloneWithModifiedTimeDiscretization(TimeDiscretization newTimeDiscretization) {
				IndependentIncrements originClone = increments.getCloneWithModifiedTimeDiscretization(newTimeDiscretization);
				return getRangeOfFactors(originClone, firstFactor, lastFactor);
			}
			
		};
	}
	
	/**
	 * Creates a new instance of IndependentIncrements that is basically a copy of increments,
	 * but uses only the first factors.
	 * @param increments The original instance of independent increments
	 * @param numberOfFactors The number factors to consider in the new instance
	 * @return A new instance with reduced factors
	 */
	public static IndependentIncrements getFirstFactors(IndependentIncrements increments, int numberOfFactors) {
		return getRangeOfFactors(increments, 0, numberOfFactors - 1);
	}

	/**
	 * Creates a MonteCarloProcess that has the same Process Values as that of <code>process</code> for a range of components.
	 * 
	 * @param process The original process.
	 * @param firstComponent The first component to include in the process (inclusive).
	 * @param lastComponent The last component to include in the new process (inclusive).
	 * @param newNumberOfFactors The number of factors, that the new process has. Will be cancelled in a future version.
	 * @return A component reduced MonteCarloProcess on the basis of the original process.
	 * @deprecated Use {@link #getComponentReducedMCProcess(process, firstComponent, lastComponent)} instead. 
	 * A new number of Factors is no longer supported.
	 */
	@Deprecated
	public static MonteCarloProcess getComponentReducedMCProcess(MonteCarloProcess process, int firstComponent, int lastComponent, int newNumberOfFactors) {
		if(firstComponent == 0 && lastComponent == process.getNumberOfComponents() - 1)
			return process;
		
		return new MonteCarloProcess() {
			private int getOriginalIndex(int componentIndex) {
				return componentIndex + firstComponent;
			}
			
			@Override
			public MonteCarloProcess clone() {
				return Functional.getComponentReducedMCProcess(process, firstComponent, lastComponent, newNumberOfFactors);
			}

			@Override
			public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
				return process.getProcessValue(timeIndex, getOriginalIndex(componentIndex));
			}

			@Override
			public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
				return process.getMonteCarloWeights(timeIndex);
			}

			@Override
			public int getNumberOfComponents() {
				return firstComponent - lastComponent + 1;
			}

			@Override
			public TimeDiscretization getTimeDiscretization() {
				return process.getTimeDiscretization();
			}

			@Override
			public double getTime(int timeIndex) {
				return process.getTime(timeIndex);
			}

			@Override
			public int getTimeIndex(double time) {
				return process.getTimeIndex(time);
			}

			@Override
			public int getNumberOfPaths() {
				return process.getNumberOfPaths();
			}

			@Override
			public int getNumberOfFactors() {
				return newNumberOfFactors == 0? process.getNumberOfFactors() : newNumberOfFactors;
			}

			@Override
			public IndependentIncrements getStochasticDriver() {
				return process.getStochasticDriver();
			}

			@Override
			public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
				return Functional.getComponentReducedMCProcess(process.getCloneWithModifiedModel(model), firstComponent, lastComponent, newNumberOfFactors);
			}

			@Override
			public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
				return Functional.getComponentReducedMCProcess(process.getCloneWithModifiedData(dataModified), firstComponent, lastComponent, newNumberOfFactors);
			}
		};
	}

	/**
	 * Creates a MonteCarloProcess that has the same Process Values as that of <code>process</code> for a range of components.
	 * 
	 * @param process The original process.
	 * @param firstComponent The first component to include in the process (inclusive).
	 * @param lastComponent The last component to include in the new process (inclusive).
	 * @return A component reduced MonteCarloProcess on the basis of the original process.
	 */
	public static MonteCarloProcess getComponentReducedMCProcess(MonteCarloProcess process, int firstComponent, int lastComponent) {
		return getComponentReducedMCProcess(process, firstComponent, lastComponent, 0);
	}
	
	/**
	 * This is a Function to combine two MC Processes that stem from the same stochastic driver and are discretized on the same tenor and 
	 * hence are compatible.
	 * 
	 * @param firstProcess The process that will be on the first component indices.
	 * @param secondProcess The process that will be on the remaining component indices.
	 * @return A combined MonteCarloProcess that gives the values of the one process for the first component indices and that of another for 
	 * higher component indices.
	 */
	public static MonteCarloProcess getCombinedMCProcess(MonteCarloProcess firstProcess, MonteCarloProcess secondProcess) {
		if(firstProcess.getStochasticDriver() != secondProcess.getStochasticDriver() && !firstProcess.getStochasticDriver().equals(secondProcess.getStochasticDriver()))
			throw new IllegalArgumentException("first and second Process must use the same stochastic driver!");
		
		if(firstProcess.getTimeDiscretization() != secondProcess.getTimeDiscretization() && !firstProcess.getTimeDiscretization().equals(secondProcess.getTimeDiscretization()))
			throw new IllegalArgumentException("first and second Process must use the same time discretization!");
		
		return new MonteCarloProcess() {
			private int getOriginalIndex(int componentIndex) {
				if(componentIndex < firstProcess.getNumberOfComponents())
					return componentIndex;
				else
					return componentIndex - firstProcess.getNumberOfComponents();
			}
			
			private MonteCarloProcess getResponsibleProcess(int componentIndex) {
				if(componentIndex < firstProcess.getNumberOfComponents())
					return firstProcess;
				else
					return secondProcess;
			}
			
			@Override
			public MonteCarloProcess clone() {
				return Functional.getCombinedMCProcess(firstProcess, secondProcess);
			}

			@Override
			public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
				return getResponsibleProcess(componentIndex).getProcessValue(timeIndex, getOriginalIndex(componentIndex));
			}

			/**
			 * Use MC Weights of first process.
			 */
			@Override
			public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
				return firstProcess.getMonteCarloWeights(timeIndex);
			}

			@Override
			public int getNumberOfComponents() {
				return firstProcess.getNumberOfComponents() + secondProcess.getNumberOfComponents();
			}

			@Override
			public TimeDiscretization getTimeDiscretization() {
				return firstProcess.getTimeDiscretization();
			}

			@Override
			public double getTime(int timeIndex) {
				return firstProcess.getTime(timeIndex);
			}

			@Override
			public int getTimeIndex(double time) {
				return firstProcess.getTimeIndex(time);
			}

			@Override
			public int getNumberOfPaths() {
				return firstProcess.getNumberOfPaths();
			}

			@Override
			public int getNumberOfFactors() {
				return Math.max(firstProcess.getNumberOfFactors(), secondProcess.getNumberOfFactors());
			}

			@Override
			public IndependentIncrements getStochasticDriver() {
				return firstProcess.getStochasticDriver();
			}

			@Override
			public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
				if(firstProcess.getModel() != secondProcess.getModel())
					throw new UnsupportedOperationException("This method makes no sense in this situation!");
				return Functional.getCombinedMCProcess(firstProcess.getCloneWithModifiedModel(model), secondProcess.getCloneWithModifiedModel(model));
			}

			@Override
			public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
				throw new UnsupportedOperationException();
			}
		};
	}
}
