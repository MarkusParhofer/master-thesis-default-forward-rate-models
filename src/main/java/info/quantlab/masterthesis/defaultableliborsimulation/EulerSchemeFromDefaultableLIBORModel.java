/**
 * 
 */
package info.quantlab.masterthesis.defaultableliborsimulation;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * @author Markus Parhofer
 *
 */
public class EulerSchemeFromDefaultableLIBORModel implements MonteCarloProcess {

	/**
	 * 
	 */
	public EulerSchemeFromDefaultableLIBORModel() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getNumberOfComponents() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public TimeDiscretization getTimeDiscretization() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getTime(int timeIndex) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getTimeIndex(double time) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getNumberOfPaths() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getNumberOfFactors() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IndependentIncrements getStochasticDriver() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
		// TODO Auto-generated method stub
		return null;
	}

}
