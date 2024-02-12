package info.quantlab.easyplot;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.plots.Named;
import net.finmath.stochastic.RandomVariable;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;

@SuppressWarnings("unused")
public class PlotFactory {

    public interface Process {
        RandomVariable apply(double time) throws CalculationException;
    }

    /**
     * Specifies an operator to apply to the MonteCarloProcess values before plotting it.
     * I.e. it may not always be wanted to directly plot the paths of the MCProcess, but derivatives of these.
     * In this case one can specify a ProcessOperator.
     * The default is monteCarloProcess.getProcessValue(timeIndex, componentIndex).
     */
    public interface ProcessOperator {
        RandomVariable apply(MonteCarloProcess monteCarloProcess, int timeIndex, int componentIndex)
                throws CalculationException;
    }

    /**
     * Specifies an operator to name the Paths. The default names each path "Path {index}".
     */
    public interface PathNamer {
        String apply(int pathIndex);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime, PathNamer pathNamer) {
        final DoubleToIntFunction getTimeIndex = operand -> {
            int timeIndex = monteCarloProcess.getTimeIndex(operand);
            return timeIndex < 0 ? - timeIndex - 1 : timeIndex;
        };

        DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfPaths];
        for (int index = 0; index < numberOfPaths; index++) {
            final int path = index + firstPath;
            myFunctionArray[index] = (time) -> {
                try {
                    return processFunc.apply(monteCarloProcess, getTimeIndex.applyAsInt(time), componentIndex).get(path);
                } catch (CalculationException e) {
                    return -1.0;
                }
            };
        }
        List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
        List<Named<DoubleUnaryOperator>> myFunctionList =
                myList.stream().map(operator -> new Named<>(
                                pathNamer.apply(myList.indexOf(operator) + firstPath), operator)).toList();

        int numberOfPoints = getTimeIndex.applyAsInt(untilTime) - getTimeIndex.applyAsInt(fromTime) + 1;

        return new EasyPlot2D(fromTime, untilTime, numberOfPoints, myFunctionList);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getFirstTime(), untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime) {
        return PlotProcessPaths(monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(processFunc, monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths) {
        return PlotProcessPaths(monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime());
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, 0, numberOfPaths, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, 0, numberOfPaths, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess, componentIndex, 0, numberOfPaths, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths) {
        return PlotProcessPaths(monteCarloProcess, componentIndex, 0, numberOfPaths);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, PathNamer pathNamer) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, monteCarloProcess.getNumberOfPaths(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex) {
        return PlotProcessPaths(processFunc, monteCarloProcess, componentIndex, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex, PathNamer pathNamer) {
        return PlotProcessPaths(defaultProcessOperator, monteCarloProcess, componentIndex, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into a new plot.
     *
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(MonteCarloProcess monteCarloProcess, int componentIndex) {
        return PlotProcessPaths(monteCarloProcess, componentIndex, monteCarloProcess.getNumberOfPaths());
    }


    // Add to existing Plot

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime, PathNamer pathNamer) {
        final DoubleToIntFunction getTimeIndex = operand -> {
            int timeIndex = monteCarloProcess.getTimeIndex(operand);
            return timeIndex < 0 ? - timeIndex - 1 : timeIndex;
        };

        DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfPaths];
        for (int index = 0; index < numberOfPaths; index++) {
            final int path = index + firstPath;
            myFunctionArray[index] = (time) -> {
                try {
                    return processFunc.apply(monteCarloProcess, getTimeIndex.applyAsInt(time), componentIndex).get(path);
                } catch (CalculationException e) {
                    return -1.0;
                }
            };
        }
        List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
        List<Named<DoubleUnaryOperator>> myFunctionList =
                myList.stream().map(operator -> new Named<>(
                        pathNamer.apply(myList.indexOf(operator) + firstPath), operator)).toList();

        int numberOfPoints = getTimeIndex.applyAsInt(untilTime) - getTimeIndex.applyAsInt(fromTime) + 1;
        plot.addPlot(fromTime, untilTime, numberOfPoints, myFunctionList);
        return plot;
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess         The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param fromTime          The xMin of the plots.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double fromTime, double untilTime) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, fromTime, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getFirstTime(), untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime, PathNamer pathNamer) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param untilTime         The xMax of the plots.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex, int firstPath, int numberOfPaths, double untilTime) {
        return PlotProcessPaths(plot, monteCarloProcess, componentIndex, firstPath, numberOfPaths, untilTime, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex,  int firstPath, int numberOfPaths) {
        return PlotProcessPaths(plot, monteCarloProcess,componentIndex, firstPath, numberOfPaths, monteCarloProcess.getTimeDiscretization().getLastTime());
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, 0, numberOfPaths, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, 0, numberOfPaths, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths, PathNamer pathNamer) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess, componentIndex, 0, numberOfPaths, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex,  int numberOfPaths) {
        return PlotProcessPaths(plot, monteCarloProcess, componentIndex, 0, numberOfPaths);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex, PathNamer pathNamer) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, monteCarloProcess.getNumberOfPaths(), pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param processFunc       An operator to apply to the {@link MonteCarloProcess} before plotting them.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, ProcessOperator processFunc, MonteCarloProcess monteCarloProcess, int componentIndex) {
        return PlotProcessPaths(plot, processFunc, monteCarloProcess, componentIndex, defaultPathnamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex, PathNamer pathNamer) {
        return PlotProcessPaths(plot, defaultProcessOperator, monteCarloProcess, componentIndex, pathNamer);
    }

    /**
     * Plots the paths of a stochastic process into an existing plot.
     *
     * @param plot              The plot instance to add the paths to.
     * @param monteCarloProcess The process to plot the paths of.
     * @param componentIndex    The component of the process to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, MonteCarloProcess monteCarloProcess, int componentIndex) {
        return PlotProcessPaths(plot, monteCarloProcess, componentIndex, monteCarloProcess.getNumberOfPaths());
    }


    // Directly plot a Process given by a Double to RandomVariable operator

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The new plot.
     */
    public static EasyPlot2D PlotProcessPaths(Process process, double xMin, double xMax, int numberOfPoints, int firstPath, int numberOfPaths, PathNamer pathNamer) {
        DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfPaths];
        for (int index = 0; index < numberOfPaths; index++) {
            final int path = index + firstPath;
            myFunctionArray[index] = (time) -> {
                try {
                    return process.apply(time).get(path);
                } catch (CalculationException e) {
                    return -1.0;
                }
            };
        }
        List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
        List<Named<DoubleUnaryOperator>> myFunctionList =
                myList.stream().map(operator -> new Named<>(
                        pathNamer.apply(myList.indexOf(operator) + firstPath), operator)).toList();

        return new EasyPlot2D(xMin, xMax, numberOfPoints, myFunctionList);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The new plot.
     */
    public static EasyPlot2D PlotProcessPaths(Process process, double xMin, double xMax, int numberOfPoints, int firstPath, int numberOfPaths) {
        return PlotProcessPaths(process, xMin, xMax, numberOfPoints, firstPath, numberOfPaths, defaultPathnamer);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param numberOfPaths     The number of paths to plot.
     * @return The new plot.
     */
    public static EasyPlot2D PlotProcessPaths(Process process, double xMin, double xMax, int numberOfPoints, int numberOfPaths) {
        return PlotProcessPaths(process, xMin, xMax, numberOfPoints, 0, numberOfPaths);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @return The new plot.
     */
    public static EasyPlot2D PlotProcessPaths(Process process, double xMin, double xMax, int numberOfPoints) {
        int numberOfPaths = 0;
        try {
            numberOfPaths = Math.max(process.apply(xMax).size(), process.apply(xMin).size());
        } catch (CalculationException ex) {
            ex.printStackTrace();
        }
        return PlotProcessPaths(process, xMin, xMax, numberOfPoints, numberOfPaths);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param plot              The plot instance to add the paths to.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @param pathNamer         A String operator to apply to each path index for naming the path.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, Process process, double xMin, double xMax, int numberOfPoints, int firstPath, int numberOfPaths, PathNamer pathNamer) {
        DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfPaths];
        for (int index = 0; index < numberOfPaths; index++) {
            final int path = index + firstPath;
            myFunctionArray[index] = (time) -> {
                try {
                    return process.apply(time).get(path);
                } catch (CalculationException e) {
                    return -1.0;
                }
            };
        }
        List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
        List<Named<DoubleUnaryOperator>> myFunctionList =
                myList.stream().map(operator -> new Named<>(
                        pathNamer.apply(myList.indexOf(operator) + firstPath), operator)).toList();

        plot.addPlot(xMin, xMax, numberOfPoints, myFunctionList);
        return plot;
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param plot              The plot instance to add the paths to.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param firstPath         The index of the first path to plot.
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, Process process, double xMin, double xMax, int numberOfPoints, int firstPath, int numberOfPaths) {
        return PlotProcessPaths(plot, process, xMin, xMax, numberOfPoints, firstPath, numberOfPaths, defaultPathnamer);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param plot              The plot instance to add the paths to.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @param numberOfPaths     The number of paths to plot.
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, Process process, double xMin, double xMax, int numberOfPoints, int numberOfPaths) {
        return PlotProcessPaths(plot, process, xMin, xMax, numberOfPoints, 0, numberOfPaths);
    }

    /**
     * Directly plots a Process given by a Double to RandomVariable operator into a new Plot.
     * @param plot              The plot instance to add the paths to.
     * @param process           The Process to plot the paths of
     * @param xMin              The lower bound of x values (time) to plot
     * @param xMax              The upper bound of x values (time) to plot
     * @param numberOfPoints    The number of points to use for the plot (equally spread out on the x-axis).
     * @return The updated plot.
     */
    public static EasyPlot2D PlotProcessPaths(EasyPlot2D plot, Process process, double xMin, double xMax, int numberOfPoints) {
        int numberOfPaths = 0;
        try {
            numberOfPaths = Math.max(process.apply(xMax).size(), process.apply(xMin).size());
        } catch (CalculationException ex) {
            ex.printStackTrace();
        }
        return PlotProcessPaths(plot, process, xMin, xMax, numberOfPoints, numberOfPaths);
    }

    // Private Members:
    public static PathNamer defaultPathnamer = (index) -> "Path " + index;
    public static ProcessOperator defaultProcessOperator = net.finmath.montecarlo.process.Process::getProcessValue;


}
