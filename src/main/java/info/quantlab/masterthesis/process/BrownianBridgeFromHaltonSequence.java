package info.quantlab.masterthesis.process;

import net.finmath.exception.CalculationException;
import net.finmath.functions.NormalDistribution;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.randomnumbers.HaltonSequence;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

import java.util.concurrent.Callable;

public class BrownianBridgeFromHaltonSequence {
    private final HaltonSequence _generator;
    private final RandomVariableFactory _rvFactory;
    private final int _numberOfPaths;
    private final TimeDiscretization _td;

    private final RandomVariable[][] _brownianIncrements;

    public BrownianBridgeFromHaltonSequence(HaltonSequence generator, RandomVariableFactory rvFactory, int numberOfPaths, TimeDiscretization td) {
        _generator = generator;
        _rvFactory = rvFactory;
        _numberOfPaths = numberOfPaths;
        _td = td;
        _brownianIncrements = new RandomVariable[_generator.getDimension()][];
    }

    private void calculatePaths() throws CalculationException {
        if(_brownianIncrements != null) {
            return;
        }
        final int factors = _generator.getDimension();
        for(int j = 0; j < factors; j++) {
            final int factor = j;
            Callable<Boolean> worker = () -> {
                RandomVariable[] paths = new RandomVariable[_td.getNumberOfTimes()];
                int targetIndex, rvIndex = 0, lowerIndex = 0, upperIndex = _td.getNumberOfTimes() - 1;
                paths[lowerIndex] = _rvFactory.createRandomVariable(0.0);
                if(_td.getFirstTime() != 0) {
                    paths[lowerIndex] = getStandardNormalRV(factor, rvIndex++).mult(Math.sqrt(_td.getFirstTime()));
                }
                paths[upperIndex] = getStandardNormalRV(factor, rvIndex++).mult(Math.sqrt(_td.getLastTime()));
                RandomVariable lower = paths[lowerIndex];
                RandomVariable upper = paths[upperIndex];
                for(;rvIndex < _td.getNumberOfTimes(); rvIndex++) {
                    targetIndex = (lowerIndex + upperIndex) / 2; // TODO Finish implementation
                }
                return true;
            };
        }
    }

    private RandomVariable getStandardNormalRV(int factor, int index) {
        long adjIndex = (long)_numberOfPaths * index;
        double[] states = new double[_numberOfPaths];
        for(int i= 0; i < _numberOfPaths; i++, adjIndex++) {
            final double uniform = _generator.getHaltonNumber(adjIndex, factor);
            states[i] = NormalDistribution.inverseCumulativeDistribution(uniform);
        }
        return _rvFactory.createRandomVariable(states);
    }
}
