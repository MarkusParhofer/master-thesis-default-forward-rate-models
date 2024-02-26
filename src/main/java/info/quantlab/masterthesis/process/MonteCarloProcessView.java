package info.quantlab.masterthesis.process;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

import java.util.Map;

public class MonteCarloProcessView implements MonteCarloProcess {
    final MonteCarloProcess _mcProcess;
    final int[] _components;

    public MonteCarloProcessView(MonteCarloProcess baseProcess, int[] componentIndices) {
        _mcProcess = baseProcess;
        _components = componentIndices;
    }
    @Override
    public int getNumberOfPaths() {
        return _mcProcess.getNumberOfPaths();
    }

    @Override
    public int getNumberOfFactors() {
        return _mcProcess.getNumberOfFactors();
    }

    @Override
    public IndependentIncrements getStochasticDriver() {
        return _mcProcess.getStochasticDriver();
    }

    @Override
    public MonteCarloProcess getCloneWithModifiedModel(ProcessModel model) {
        return new MonteCarloProcessView(_mcProcess.getCloneWithModifiedModel(model), _components);
    }

    @Override
    public MonteCarloProcess getCloneWithModifiedData(Map<String, Object> dataModified) {
        return null;
    }

    @Override
    public RandomVariable getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
        return _mcProcess.getProcessValue(timeIndex, _components[componentIndex]);
    }

    @Override
    public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
        return _mcProcess.getMonteCarloWeights(timeIndex);
    }

    @Override
    public int getNumberOfComponents() {
        return _components.length;
    }

    @Override
    public TimeDiscretization getTimeDiscretization() {
        return _mcProcess.getTimeDiscretization();
    }

    @Override
    public double getTime(int timeIndex) {
        return _mcProcess.getTime(timeIndex);
    }

    @Override
    public int getTimeIndex(double time) {
        return _mcProcess.getTimeIndex(time);
    }

    @Override
    public MonteCarloProcess clone() {
        return new MonteCarloProcessView(_mcProcess, _components);
    }
}
