package sandbox;

import java.util.List;
import java.util.function.DoubleUnaryOperator;

import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablecovariancemodels.SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORWithPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORWithPositiveSpread.Measure;
import info.quantlab.masterthesis.defaultableliborsimulation.EulerSchemeFromDefaultableLIBORModel;
import info.quantlab.masterthesis.defaultableliborsimulation.EulerSchemeWithDependencyModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
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
		
		// Set undefaultable LIBOR model
		LIBORMarketModel uModel = new LIBORMarketModelFromCovarianceModel(liborPeriods, uForwardCurve, uCovarianceModel);
		
		
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
		
		
		/*
		 * Caplet from undefaultableModel = cU
		 * Caplet from defaultableModel = cD
		 * E[cU] = E[cD]*P(Default)
		 */
	}

	
	public static EasyPlot2D plotPaths(RandomVariable[] paths, int numberOfPaths, TimeDiscretization xAxis, EasyPlot2D plot) {
		DoubleUnaryOperator[] myList = new DoubleUnaryOperator[numberOfPaths];
		for(int j = 0; j < numberOfPaths; j++) {
			final int pathIndex = j;
			myList[j] = new DoubleUnaryOperator() {
				@Override
				public double applyAsDouble(double operand) {
					final int index = xAxis.getTimeIndex(operand);
					return paths[index].isDeterministic()? paths[index].doubleValue() : paths[index].get(pathIndex);
				}
			};
		}
		
		if(plot == null) {
			// TODO: Turn into List<Named<UnaryOperator>>
			List<Named<DoubleUnaryOperator>> myList;
		} else {
			for(int j=0; j<numberOfPaths; j++) {
				plot.addPlot(xAxis.getTime(0), xAxis.getTime(paths.length - 1), paths.length, 
						new Named<DoubleUnaryOperator>("Path " + j, myList[j]));
			}
		}
		return plot;
	}
}
