/**
 * 
 */
package info.quantlab.masterthesis.legacy;

import java.util.Map;

import org.apache.commons.lang3.Validate;

import info.quantlab.masterthesis.functional.FunctionsOnIndependentIncrements;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

/**
 * @author Markus Parhofer
 *
 */
public class EulerSchemeWithDependencyModel extends EulerSchemeFromProcessModel implements MonteCarloProcessWithDependency {

	private final MonteCarloProcess _dependencyProcess;
	
	/**
	 * @param model
	 * @param stochasticDriver
	 */
	public EulerSchemeWithDependencyModel(ProcessModel model, ProcessModel dependencyModel, IndependentIncrements stochasticDriver) {
		super(model, stochasticDriver);
		
		Validate.isTrue(stochasticDriver.getNumberOfFactors() >= Math.max(model.getNumberOfFactors(), dependencyModel.getNumberOfFactors()), 
				"Number of Factors of the stochastic driver must be greater or equal the maximum of the number of factors of the two specified models");
		_dependencyProcess = new EulerSchemeFromProcessModel(dependencyModel, FunctionsOnIndependentIncrements.getFirstFactors(stochasticDriver, dependencyModel.getNumberOfFactors()));
	}

	@Override
	public ProcessModel getDependencyModel() {
		return getDependencyProcess().getModel();
	}
	
	public MonteCarloProcess getDependencyProcess() {
		return _dependencyProcess;
	}
	
	public RandomVariable getDependencyProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		return _dependencyProcess.getProcessValue(timeIndex, componentIndex);
	}

	public RandomVariable[] getDependencyProcessValue(int timeIndex) throws CalculationException {
		return _dependencyProcess.getProcessValue(timeIndex);
	}
	
	@Override
	public EulerSchemeWithDependencyModel clone() {
		IndependentIncrements stochasticDriver = getNumberOfFactors() >= _dependencyProcess.getNumberOfFactors()? getStochasticDriver() : _dependencyProcess.getStochasticDriver();
		return new EulerSchemeWithDependencyModel(getModel(), _dependencyProcess.getModel(), stochasticDriver);
	}
	
	@Override
	public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
		IndependentIncrements stochasticDriver = getNumberOfFactors() >= _dependencyProcess.getNumberOfFactors()? 
				getStochasticDriver() : _dependencyProcess.getStochasticDriver();
		
		return new EulerSchemeWithDependencyModel(model, _dependencyProcess.getModel(), stochasticDriver);
	}
	
	public MonteCarloProcess getCloneWithModifiedDependencyModel(ProcessModel model) {
		IndependentIncrements stochasticDriver = getNumberOfFactors() >= _dependencyProcess.getNumberOfFactors()? 
				getStochasticDriver() : _dependencyProcess.getStochasticDriver();
		
		return new EulerSchemeWithDependencyModel(getModel(), model, stochasticDriver);
	}

	public MonteCarloProcess getCloneWithModifiedStochasticDriver(IndependentIncrements stochasticDriver) {
		return new EulerSchemeWithDependencyModel(getModel(), _dependencyProcess.getModel(), stochasticDriver);
	}
	
	@Override
	public MonteCarloProcess getCloneWithModifiedSeed(final int seed) {
		IndependentIncrements stochasticDriver = getNumberOfFactors() >= _dependencyProcess.getNumberOfFactors()? 
				getStochasticDriver().getCloneWithModifiedSeed(seed) : _dependencyProcess.getStochasticDriver().getCloneWithModifiedSeed(seed);
		
		return new EulerSchemeWithDependencyModel(getModel(), _dependencyProcess.getModel(), stochasticDriver);
	}

	
	@Override
	public MonteCarloProcess getCloneWithModifiedData(final Map<String, Object> dataModified) {
		final ProcessModel newModel = (ProcessModel) dataModified.getOrDefault("model", getModel());
		final ProcessModel newDependencyModel = (ProcessModel) dataModified.getOrDefault("dependencyModel", _dependencyProcess.getModel());
		
		IndependentIncrements newStochasticDriver = getNumberOfFactors() >= _dependencyProcess.getNumberOfFactors()? 
				getStochasticDriver() : _dependencyProcess.getStochasticDriver();
		
		
		if(dataModified.containsKey("seed") && dataModified.containsKey("stochasticDriver")) {
			throw new IllegalArgumentException("Simultaneous specification of stochasticDriver and seed.");
		}

		if(dataModified.containsKey("seed")) {
			return getCloneWithModifiedSeed((int)dataModified.get("seed"));
		}
		else if(dataModified.containsKey("stochasticDriver")) {
			newStochasticDriver = (IndependentIncrements) dataModified.getOrDefault("stochasticDriver", newStochasticDriver);
		}

		return new EulerSchemeWithDependencyModel(newModel, newDependencyModel, newStochasticDriver);
	}
	
	@Override
	public String toString() {
		return "EulerSchemeFromProcessModel [stochasticDriver=" + getStochasticDriver() + ", scheme=" + getScheme() + ", " + "]";
	}

}
