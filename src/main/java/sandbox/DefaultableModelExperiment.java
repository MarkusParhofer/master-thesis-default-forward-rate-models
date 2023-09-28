package sandbox;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;

import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablecovariancemodels.SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORWithPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORWithPositiveSpread.Measure;
import info.quantlab.masterthesis.defaultableliborsimulation.EulerSchemeFromDefaultableLIBORModel;
import info.quantlab.masterthesis.defaultableliborsimulation.EulerSchemeWithDependencyModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.model.volatilities.SwaptionMarketData;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelExponentialForm5Param;
import net.finmath.plots.Named;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class DefaultableModelExperiment {

	public static void main(String[] args) throws CalculationException {
		
		System.out.println("\n");
		System.out.println("Experiment on defaultable LIBORs");
		System.out.println("_".repeat(79));

		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, 10, 0.5); // Fixing time
		
		// Set initial undefaultable LIBOR rates
		final double[] fixingTimes = liborPeriods.getAsDoubleArray();
		final double[] uForwardRates = new double[] { 0.61 / 100.0, 0.61 / 100.0, 0.67 / 100.0, 0.73 / 100.0, 0.80 / 100.0, 0.92 / 100.0, 
				1.11 / 100.0, 1.36 / 100.0, 1.60 / 100.0, 1.82 / 100.0, 2.02 / 100.0 };
		ForwardCurve uForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("forwardCurve",	fixingTimes, uForwardRates, 0.5);
		
		// Set Covariance model of undefaultable LIBORs
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, 20, 0.25); // Simulation time
		LIBORCovarianceModel uCovarianceModel = new LIBORCovarianceModelExponentialForm5Param(timeDiscretization, liborPeriods, 9);
		
		// Set some more Properties:
		SwaptionMarketData uSwaptionData = null;
		Map<String, String> properties = new HashMap<String, String>();
		properties.putIfAbsent("measure", "Spot");
		properties.putIfAbsent("stateSpace", "Normal");
		
		// Set undefaultable LIBOR model
		LIBORMarketModel uModel = new LIBORMarketModelFromCovarianceModel(liborPeriods, uForwardCurve, new DiscountCurveFromForwardCurve(uForwardCurve), 
				uCovarianceModel, uSwaptionData, properties);
		
		
		// Set initial defaultable LIBOR rates
		final double[] myForwardRates = new double[] { 0.71 / 100.0, 0.71 / 100.0, 0.74 / 100.0, 0.79 / 100.0, 0.90 / 100.0, 0.99 / 100.0, 
				1.20 / 100.0, 1.44 / 100.0, 1.71 / 100.0, 1.94 / 100.0, 2.20 / 100.0 };
		ForwardCurve myForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("defaultableForwardCurve", fixingTimes, myForwardRates, 0.5);
		
		// Set Covariance model of defaultable LIBORs
		double[][] freeParameters = new double[uModel.getNumberOfComponents()][3];
		// Will set them as random doubles between 0 and 1
		MersenneTwister randomGenerator = new MersenneTwister(1998);
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				freeParameters[componentIndex][factor] = randomGenerator.nextDouble();
			}
		}
		System.out.println("\nUndefaultable Model:\n");
		System.out.println(uModel);
		
		
		System.out.println("\nFree Parameter Matrix:\n");
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				System.out.printf("%8.6f  ", freeParameters[componentIndex][factor]);
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("_".repeat(79));
		System.out.println();
		DefaultableLIBORCovarianceModel myCovarianceModel = new SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread(uModel, freeParameters);
		
		// Set defaultable LIBOR model
		DefaultableLIBORMarketModel myModel = new DefaultableLIBORWithPositiveSpread(uModel, myForwardCurve, myCovarianceModel, Measure.SPOT);
		
		// Set stochastic driver
		int numberOfFactors = myModel.getNumberOfFactors();
		int numberOfPaths = 100000;
		int randomNumberSeed = 3141;
		BrownianMotion myBrownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, numberOfFactors, numberOfPaths, randomNumberSeed);

		// Set different MonteCarloProcesses
		EulerSchemeFromDefaultableLIBORModel myProcess = new EulerSchemeFromDefaultableLIBORModel(myModel, myBrownianMotion);
		EulerSchemeWithDependencyModel myProcess2 = new EulerSchemeWithDependencyModel(myModel, uModel, myBrownianMotion);
		
		int modelTimes = timeDiscretization.getNumberOfTimes();
		int liborTimes = liborPeriods.getNumberOfTimeSteps();
//		RandomVariable[] spreadByTimeAtLIBOR8 = new RandomVariable[modelTimes];
//		RandomVariable[] spreadByTimeAtLIBOR9 = new RandomVariable[modelTimes];
//		RandomVariable[] spreadByLIBORAtTimeIndex1 = new RandomVariable[liborTimes];
//		RandomVariable[] spreadByLIBORAtTimeIndex5 = new RandomVariable[liborTimes];
//		
//		for(int i = 0; i < Math.max(modelTimes, liborTimes); i++) {
//			if(i < modelTimes) {
//				spreadByTimeAtLIBOR8[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess, i, 8);
//				spreadByTimeAtLIBOR9[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess, i, 9);
//			}
//			if(i < liborTimes) {
//				spreadByLIBORAtTimeIndex1[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess, 1, i);
//				spreadByLIBORAtTimeIndex5[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess, 5, i);
//			}
//		}
//		
		RandomVariable[] spreadByTimeAtLIBOR0P2 = new RandomVariable[modelTimes];
		RandomVariable[] spreadByTimeAtLIBOR5P2 = new RandomVariable[modelTimes];
		RandomVariable[] spreadByTimeAtLIBOR9P2 = new RandomVariable[modelTimes];
		RandomVariable[] spreadByLIBORAtTimeIndex0P2 = new RandomVariable[liborTimes];
		RandomVariable[] spreadByLIBORAtTimeIndex1P2 = new RandomVariable[liborTimes];
		RandomVariable[] spreadByLIBORAtTimeIndex5P2 = new RandomVariable[liborTimes];
		
		for(int i = 0; i < Math.max(modelTimes, liborTimes); i++) {
			if(i < modelTimes) {
				spreadByTimeAtLIBOR0P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, i, 0);
				spreadByTimeAtLIBOR5P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, i, 5);
				spreadByTimeAtLIBOR9P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, i, 9);
			}
			if(i < liborTimes) {
				spreadByLIBORAtTimeIndex0P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, 0, i);
				spreadByLIBORAtTimeIndex1P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, 1, i);
				spreadByLIBORAtTimeIndex5P2[i] = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, 5, i);
			}
		}
		
		System.out.println();
		System.out.println("First 10 paths of independent Increments at timeIndex=2; UD = undefaultable model; DM = defaultable model");
		System.out.println();
		
		
		for(int i=-1; i<10; i++) {
			for(int j = 0; j < myModel.getNumberOfFactors(); j++) {
				if(i==-1) {
					DecimalFormat ft = new DecimalFormat("00");
					if(j < uModel.getNumberOfFactors()) {
						System.out.print("UD Factor " + ft.format(j) + " ".repeat(1) + "|");
					}
					System.out.print("DM Factor " + ft.format(j) + " ".repeat(1) + "|");
					continue;
				}
				
				final double increment = myProcess2.getStochasticDriver().getIncrement(2, j).get(i);
				if(j < uModel.getNumberOfFactors()) {
					
					double copyIncrement = myProcess2.getDependencyProcess().getStochasticDriver().getIncrement(2, j).get(i);
					if(copyIncrement == increment)
						System.out.printf("%12.6f |", copyIncrement);
					else
						System.out.printf("\033[41m%12.6f \033[m|", copyIncrement);
				    
				}
				System.out.printf("%12.6f |", increment);
			}
			System.out.println();
		}
		
		
		EasyPlot2D spreadByTime = plotPaths(spreadByTimeAtLIBOR0P2, 3, timeDiscretization, "LIBOR 0 P2", null);
		spreadByTime = plotPaths(spreadByTimeAtLIBOR5P2, 3, timeDiscretization, "LIBOR 5 P2", spreadByTime);
		// spreadByTime = plotPaths(spreadByTimeAtLIBOR8P2, 5, timeDiscretization, "LIBOR 8 P2", spreadByTime);
		// spreadByTime = plotPaths(spreadByTimeAtLIBOR9, 5, timeDiscretization, "LIBOR 9 P1", spreadByTime);
		spreadByTime = plotPaths(spreadByTimeAtLIBOR9P2, 3, timeDiscretization, "LIBOR 9 P2", spreadByTime);
		spreadByTime.setTitle("Spread by Simulation Time");
		spreadByTime.setXAxisLabel("Time Unit");
		spreadByTime.setIsLegendVisible(true);
		
		EasyPlot2D spreadByLIBOR = plotPaths(spreadByLIBORAtTimeIndex0P2, 3, liborPeriods, "TimeIndex 0 P2", null);
		spreadByLIBOR = plotPaths(spreadByLIBORAtTimeIndex1P2, 3, liborPeriods, "TimeIndex 1 P2", spreadByLIBOR);
		// spreadByLIBOR = plotPaths(spreadByLIBORAtTimeIndex1P2, 5, liborPeriods, "TimeIndex 1 P2", spreadByLIBOR);
		// spreadByLIBOR = plotPaths(spreadByLIBORAtTimeIndex5, 5, liborPeriods, "TimeIndex 5", spreadByLIBOR);
		spreadByLIBOR = plotPaths(spreadByLIBORAtTimeIndex5P2, 3, liborPeriods, "TimeIndex 5 P2", spreadByLIBOR);
		spreadByLIBOR.setTitle("Spread by LIBOR Period");
		spreadByLIBOR.setXAxisLabel("LIBOR Period Start");
		spreadByLIBOR.setIsLegendVisible(true);
		
		spreadByTime.show();
		spreadByLIBOR.show();
		
		System.out.println("\n");
		System.out.println("_".repeat(79));
		System.out.println();
		System.out.println("Defaultable/Undefaultable LIBOR Rates at Time Index 0 and their Spread");
		System.out.println();
		for(int i=0; i<liborPeriods.getNumberOfTimeSteps(); i++) {
			System.out.printf("%12.8f |", myModel.getLIBOR(myProcess2, 0, i).doubleValue());
		}
		System.out.println();
		for(int i=0; i<liborPeriods.getNumberOfTimeSteps(); i++) {
			System.out.printf("\033[4m%12.8f \033[m|", uModel.getLIBOR(myProcess.getDependencyProcess(), 0, i).doubleValue());
		}
		System.out.println();
		for(int i=0; i<liborPeriods.getNumberOfTimeSteps(); i++) {
			System.out.printf("%12.8f |", spreadByLIBORAtTimeIndex0P2[i].doubleValue());
		}
		
		
		
		
		System.out.println("\n");
		System.out.println("_".repeat(79));
		System.out.println();
		System.out.println("Defaultable/Undefaultable LIBOR Rates at LIBOR Index 9 and their Spread");
		System.out.println();
		for(int path=0; path < 50; path++) {
			System.out.println("Path " + path + ":");
			for(int i=1; i<timeDiscretization.getNumberOfTimeSteps(); i++) {
				System.out.printf("%12.8f |", myModel.getLIBOR(myProcess2, i, 9).get(path));
			}
			System.out.println();
			for(int i=1; i<timeDiscretization.getNumberOfTimeSteps(); i++) {
				System.out.printf("\033[4m%12.8f \033[m|", uModel.getLIBOR(myProcess.getDependencyProcess(), i, 9).get(path));
			}
			System.out.println();
			for(int i=1; i<timeDiscretization.getNumberOfTimeSteps(); i++) {
				final double spread = spreadByTimeAtLIBOR9P2[i].get(path);
				if(spread >= 0)
					System.out.printf("%12.8f |", spread);
				else
					System.out.printf("\033[41m%12.8f \033[m|", spread);
			}
			System.out.println("\n");
		}
		/*
		 * Caplet from undefaultableModel = cU
		 * Caplet from defaultableModel = cD
		 * E[cU] = E[cD]*P(Default)
		 */
		
		// -------------------- Control Area -----------------------------------------------
		System.out.println("\n");
		System.out.println("_".repeat(79));
		System.out.println();
		System.out.println("Control Area: Checking " + (liborTimes * modelTimes * numberOfPaths) + " numbers for positivity. Stay tuned!\n");
		// double[] allSpreads = new double[liborTimes * modelTimes * numberOfPaths];
		System.out.println("All negative Spreads:\n");
		long counter = 0;
		double min = 10.0, max = -10.0;
		long minIndex = -1;
		long maxIndex = -1;
		for(int component=0; component < liborTimes; component++) {
			long counterPerComponent = 0;
			long numberPerComp = 0;
			for(int time=0; time < modelTimes; time++) {
				if(timeDiscretization.getTime(time) > liborPeriods.getTime(component)) {
					System.out.printf("\n Jumped at TimeIndex: %3d, Time: %5.2f, LIBORIndex: %3d, LIBOR Time: %5.2f\n", time, timeDiscretization.getTime(time), component, liborPeriods.getTime(component));
					numberPerComp = (long)time * (long)numberOfPaths;
					break;
				}
				for(int path=0; path < numberOfPaths; path++) {
					
					final long index = component * modelTimes * numberOfPaths + time * numberOfPaths + path;
					/* allSpreads[index] */
					final double spread = myModel.getLIBORSpreadAtGivenTimeIndex(myProcess2, time, component).get(path);
					minIndex = Math.min(min, spread) == min? minIndex : index;
					min = Math.min(min, spread);
					maxIndex = Math.max(max, spread) == max? maxIndex : index;
					max = Math.max(max, spread); 
					if(spread < 0.0) {
						// System.out.printf("Component: %3d, Time: %4d, Path: %7d, Spread: %12.8f\n", component, time, path, spread);
						counter++;
						counterPerComponent++;
					}
				}
			}
			System.out.printf("%7d of the %15d spreads are negative at component %2d\n\n", counterPerComponent, numberPerComp, component);
		}
		
		System.out.println("Count of all numbers: " + (liborTimes * modelTimes * numberOfPaths));
		System.out.println("Count of negative numbers: " + counter);
		System.out.println("Count of positive numbers: " + (liborTimes * modelTimes * numberOfPaths - counter));
		System.out.println("Minimum:");
		long path = minIndex % numberOfPaths;
		long time = Math.floorDiv(minIndex, numberOfPaths) % modelTimes;
		long component = ((minIndex - path) / numberOfPaths - time) / modelTimes;
		System.out.printf("LIBOR Index:  %12d \n", component);
		System.out.printf("Time Index:   %12d \n", time);
		System.out.printf("Path Index:   %12d \n", path);
		System.out.printf("Min Value:    %12.8f \n\n", min);
		
		System.out.println("Maximum:");
		path = maxIndex % numberOfPaths;
		time = Math.floorDiv(maxIndex, numberOfPaths) % modelTimes;
		component = ((maxIndex - path) / numberOfPaths - time) / modelTimes;
		System.out.printf("LIBOR Index:  %12d \n", component);
		System.out.printf("Time Index:   %12d \n", time);
		System.out.printf("Path Index:   %12d \n", path);
		System.out.printf("Max Value:    %12.8f \n", max);
	}

	
	public static EasyPlot2D plotPaths(RandomVariable[] paths, int numberOfPaths, TimeDiscretization xAxis, String name, EasyPlot2D plot) {
		DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfPaths];
		for(int j = 0; j < numberOfPaths; j++) {
			final int pathIndex = j;
			myFunctionArray[j] = new DoubleUnaryOperator() {
				@Override
				public double applyAsDouble(double operand) {
					final int index = xAxis.getTimeIndex(operand);
					return paths[index].isDeterministic()? paths[index].doubleValue() : paths[index].get(pathIndex);
				}
			};
		}
		if(plot == null) {
			List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
			List<Named<DoubleUnaryOperator>> myFunctionList =
					myList.stream().map(operator -> new Named<DoubleUnaryOperator>((name==null? "Path " : (name + " ")) + myList.indexOf(operator), operator)).
					collect(Collectors.toList());
			plot = new EasyPlot2D(xAxis.getTime(0), xAxis.getTime(paths.length - 1), paths.length, myFunctionList);
		} else {
			for(int j=0; j<numberOfPaths; j++) {
				plot.addPlot(xAxis.getTime(0), xAxis.getTime(paths.length - 1), paths.length, 
						new Named<DoubleUnaryOperator>((name==null? "Path " : (name + " ")) + j, myFunctionArray[j]));
			}
		}
		return plot;
	}
}
