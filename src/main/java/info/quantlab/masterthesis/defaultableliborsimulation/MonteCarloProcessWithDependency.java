/**
 * 
 */
package info.quantlab.masterthesis.defaultableliborsimulation;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

/**
 * 
 * @author Markus Parhofer
 */
public interface MonteCarloProcessWithDependency extends MonteCarloProcess {

	ProcessModel getDependencyModel();
	
	MonteCarloProcess getDependencyProcess();
	
	RandomVariable[] getDependencyProcessValue(int timeIndex) throws CalculationException;
	
	RandomVariable getDependencyProcessValue(int timeIndex, int componentIndex) throws CalculationException;
}
