package sandbox;

import java.util.HashMap;
import java.util.Map;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModelFromCovarianceModel;
import info.quantlab.masterthesis.legacy.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.legacy.SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.products.DefaultableCap;
import info.quantlab.masterthesis.products.DefaultableCaplet;
import info.quantlab.masterthesis.products.DefaultableCouponBond;
import info.quantlab.masterthesis.products.DefaultableFloor;
import javafx.beans.value.ObservableDoubleValue;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.model.volatilities.SwaptionMarketData;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelExponentialForm5Param;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class NewModelsExperiment {

	public static void main(String[] args) throws CalculationException {
		System.out.println("\n");
		System.out.println("Experiment on defaultable LIBORs and their Calibration");
		System.out.println("_".repeat(79));

		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, 10, 0.5); // Fixing time
		
		// Set initial undefaultable LIBOR rates
		final double[] fixingTimes = liborPeriods.getAsDoubleArray();
		final double[] uForwardRates = new double[] { 0.61 / 100.0, 0.61 / 100.0, 0.67 / 100.0, 0.73 / 100.0, 0.80 / 100.0, 0.92 / 100.0, 
				1.11 / 100.0, 1.36 / 100.0, 1.60 / 100.0, 1.82 / 100.0, 2.02 / 100.0 };
		ForwardCurve uForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("forwardCurve",	fixingTimes, uForwardRates, 0.5);
		
		// Set Covariance model of undefaultable LIBORs
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, 250, 0.02); // Simulation time
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
		final double[] myForwardRates = new double[] { 1.2 / 100.0, 1.3 / 100.0, 1.39 / 100.0, 1.31 / 100.0, 1.54 / 100.0, 1.57 / 100.0, 
				1.83 / 100.0, 1.79 / 100.0, 1.90 / 100.0, 2.44 / 100.0, 2.90 / 100.0 };
		ForwardCurve myForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("defaultableForwardCurve", fixingTimes, myForwardRates, 0.5);
		
		
		
		// Set Covariance model of defaultable LIBORs
		double[][] originalFreeParameters = new double[uModel.getNumberOfComponents()][3];
		// Will set them as random doubles between 0 and 1
		MersenneTwister randomGenerator = new MersenneTwister(1998);
		System.out.println("\nOriginal free Parameter Matrix:\n");
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				originalFreeParameters[componentIndex][factor] = randomGenerator.nextDouble();
				System.out.printf("%8.6f  ", originalFreeParameters[componentIndex][factor]);
			}
			System.out.println();
		}
		
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread originalCovarianceModel = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(uCovarianceModel, originalFreeParameters);
		

		System.out.println();
		System.out.println("_".repeat(79));
		System.out.println();
		
		
		// Set a second covariance Model, that needs to be calibrated:
		double[][] calibrationFreeParameters = new double[uModel.getNumberOfComponents()][3];
		
		System.out.println("\nFree Parameter Matrix of different Covariance Model:\n");
		
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				calibrationFreeParameters[componentIndex][factor] = randomGenerator.nextDouble();
				System.out.printf("%8.6f  ", calibrationFreeParameters[componentIndex][factor]);
			}
			System.out.println();
		}
		
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread calibrationCovarianceModel = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(uCovarianceModel, calibrationFreeParameters);
		
		
		System.out.println();
		System.out.println("_".repeat(79));
		System.out.println();
		
		
		
		// Specify Calculation Properties
		Map<String, String> extraValues = new HashMap<String, String>();
		extraValues.put("measure", "SPOT");
		extraValues.put("statespace", "NORMAL");
		extraValues.put("handleSimulationTime", "ROUND_NEAREST");
		extraValues.put("interpolationMethod", "LOG_LINEAR_UNCORRECTED");
		
		// Original Defaultable model (for valuation of calibration Products):
		DefaultableLIBORMarketModelFromCovarianceModel originalLIBORModel = new DefaultableLIBORMarketModelFromCovarianceModel(uModel, originalCovarianceModel, myForwardCurve, null, extraValues);
		
		// Second defaultable model (for calibration):
		DefaultableLIBORMarketModelFromCovarianceModel calibrationLIBORModel = new DefaultableLIBORMarketModelFromCovarianceModel(uModel, calibrationCovarianceModel, myForwardCurve, null, extraValues);
		
		// Set stochastic driver
		TimeDiscretization simulationTenor = originalLIBORModel.getCovarianceModel().getTimeDiscretization();
		int numberOfFactors = originalLIBORModel.getNumberOfFactors();
		int numberOfPaths = 5000;
		int randomNumberSeed = 3141;
		BrownianMotion myBrownianMotion = new BrownianMotionFromMersenneRandomNumbers(simulationTenor, numberOfFactors, numberOfPaths, randomNumberSeed);

		// Set Process and SimulationModel
		MonteCarloProcess originalMCProcess = new EulerSchemeFromProcessModel(originalLIBORModel, myBrownianMotion);
		TermStructureMonteCarloSimulationModel originalSimulationModel = new LIBORMonteCarloSimulationFromLIBORModel(originalMCProcess);
		
		// Calibration Products:
		CalibrationProduct[] calibrationProducts = getCalibrationProducts(originalSimulationModel);
		System.out.println("Generated calibration Products:");
		System.out.println();
		System.out.printf("%3s | %10s | %74s | Weight %n","Nr.", "Value", "Item");
		for(int i=0; i < calibrationProducts.length; i++) {
			System.out.printf("%3d | %10.5f | %74s | %3.1f %n", i, calibrationProducts[i].getTargetValue().getAverage(), calibrationProducts[i].getName(), calibrationProducts[i].getWeight());
		}
		// Prepare Calibration:
		Map<String, Object> parameters = new HashMap<>();
		parameters.put("brownianMotion", myBrownianMotion);
		parameters.put("maxIterations", 400);
		parameters.put("accuracy", 1E-7);
		
		// Try calibrate model:
		AbstractLIBORCovarianceModelParametric calibratedCovarianceModel = calibrationCovarianceModel.getCloneCalibrated(calibrationLIBORModel, calibrationProducts, parameters);
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread calibCov = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)calibratedCovarianceModel;
		System.out.println();
		System.out.println("_".repeat(79));
		System.out.println("\nResults:\n");
		System.out.println("Original vs. Calibrated Free Parameter Matrix:\n");
		System.out.println();
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				System.out.printf("%8.6f  ", originalFreeParameters[componentIndex][factor]);
			}
			System.out.print(" ".repeat(20));
			for(int factor = 0; factor < 3; factor++) {
				System.out.printf("%8.6f  ", calibCov.getFreeParameterMatrix()[componentIndex][factor]);
			}
			System.out.println();
		}
		System.out.println();
		System.out.println();
		
		DefaultableLIBORMarketModelFromCovarianceModel calibratedLIBORModel = calibrationLIBORModel.getCloneWithModifiedCovarianceModel(calibratedCovarianceModel);
		MonteCarloProcess calibratedMCProcess = new EulerSchemeFromProcessModel(calibratedLIBORModel, myBrownianMotion);
		TermStructureMonteCarloSimulationModel calibratedSimulationModel = new LIBORMonteCarloSimulationFromLIBORModel(calibratedMCProcess);
		
		CalibrationProduct[] calibratedProducts = getCalibrationProducts(calibratedSimulationModel);
		System.out.println("Calibrated Products:");
		System.out.println();
		System.out.printf("\033[4m%3s | %14s   | %16s | %16s | %16s | Weight \033[m%n","Nr.", "Original Value", "Calibrated Value", "Delta", "Relative Delta");
		for(int i=0; i < calibrationProducts.length; i++) {
			final double origValue = calibrationProducts[i].getTargetValue().getAverage();
			final double caliValue = calibratedProducts[i].getTargetValue().getAverage();
			final double delta = Math.abs(origValue - caliValue);
			final double relDelta = delta/origValue;
			
			if(delta/origValue > 0.05) {
				System.out.printf("%3d | %16.12f | %16.12f | \033[41m%16.12f\033[m | \033[41m%16.12f\033[m | %3.1f %n", i, origValue, caliValue, delta, relDelta, calibrationProducts[i].getWeight());				
			} else {
				System.out.printf("%3d | %16.12f | %16.12f | %16.12f | %16.12f | %3.1f %n", i, origValue, caliValue, delta, relDelta, calibrationProducts[i].getWeight());
			}
		}
		
	}
	
	public static CalibrationProduct[] getCalibrationProducts(TermStructureMonteCarloSimulationModel valuationModel) throws CalculationException {
		CalibrationProduct[] results = new CalibrationProduct[12];
		int index = 0;
		final double notional = 1.0;
		final double totalWeights = 1.0;
		// Calibration item 0:
		{ 
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(0.3, 8, 0.5);
			final double strikeRate = 0.06;
			
			final DefaultableCap product = new DefaultableCap(paymentTenor, strikeRate, true, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 2.0 / totalWeights;
			
			final String name = "Cap: {Strike: 0.06, Notional: 1, Times: 0.3 + 0.5 * i, i < 8 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 1:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(0.9, 16, 0.25);
			final double strikeRate = 0.03;
			final DefaultableCap product = new DefaultableCap(paymentTenor, strikeRate, true, false);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 3.0 / totalWeights;
			
			final String name = "Cap: {Strike: 0.03, Notional: 1, Times: 0.9 + 0.25 * i, i < 16 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 2:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(4.5, 4, 0.125);
			final double strikeRate = 0.01;
			final DefaultableCap product = new DefaultableCap(paymentTenor, strikeRate, false, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.5 / totalWeights;
			
			final String name = "Cap: {Strike: 0.01, Notional: 1, Times: 4.5 + 0.125 * i, i < 4 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 3:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(0.1, 5, 0.25);
			final double strikeRate = 0.02;
			final DefaultableCap product = new DefaultableCap(paymentTenor, strikeRate, true, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.5 / totalWeights;
			
			final String name = "Cap: {Strike: 0.02, Notional: 1, Times: 0.1 + 0.25 * i, i < 5 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 4:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(0.15, 6, 0.35);
			final double strikeRate = 0.006;
			final DefaultableFloor product = new DefaultableFloor(paymentTenor, strikeRate, true, false);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.9 / totalWeights;
			
			final String name = "Floor: {Strike: 0.006, Notional: 1, Times: 0.15 + 0.35 * i, i < 6 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 5:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(1.25, 3, 1.0);
			final double strikeRate = 0.01;
			final DefaultableFloor product = new DefaultableFloor(paymentTenor, strikeRate, false, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 2.5 / totalWeights;
			
			final String name = "Floor: {Strike: 0.01, Notional: 1, Times: 1.25 + 1.0 * i, i < 3 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 6:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(4.25, 4, 0.14);
			final double strikeRate = 0.03;
			final DefaultableFloor product = new DefaultableFloor(paymentTenor, strikeRate, true, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.5 / totalWeights;
			
			final String name = "Floor: {Strike: 0.03, Notional: 1, Times: 4.25 + 0.14 * i, i < 4 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 7:
		{
			final TimeDiscretization paymentTenor = new TimeDiscretizationFromArray(0.235, 20, 0.20);
			final double strikeRate = 0.015;
			final DefaultableFloor product = new DefaultableFloor(paymentTenor, strikeRate, true, false);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 4.0 / totalWeights;
			
			final String name = "Floor: {Strike: 0.015, Notional: 1, Times: 0.235 + 0.2 * i, i < 20 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 8:
		{
			final double fixingTime = 2.0;
			final double periodLength = 2.0;
			final double strikeRate = 0.029;
			final DefaultableCaplet product = new DefaultableCaplet(strikeRate, fixingTime, periodLength, 0, 0, notional, false);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.0 / totalWeights;
			
			final String name = "Caplet: {Strike: 0.029, Fixing: 2.0, Payment: 4.0, Notional: 1 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 9:
		{
			final double fixingTime = 3.865;
			final double periodLength = 0.566;// 4.431
			final double strikeRate = 0.024;
			final DefaultableCaplet product = new DefaultableCaplet(strikeRate, fixingTime, periodLength, 0, -1, notional, false);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.0 / totalWeights;
			
			final String name = "Caplet: {Strike: 0.024, Fixing: 3.865, Payment: 4.431, Notional: 1 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 10:
		{
			final double fixingTime = 1.5;
			final double periodLength = 1.0;
			final double strikeRate = 0.014;
			final DefaultableCaplet product = new DefaultableCaplet(strikeRate, fixingTime, periodLength, -1, 0, notional, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.0 / totalWeights;
			
			final String name = "Floorlet: {Strike: 0.014, Fixing: 1.5, Payment: 2.5, Notional: 1 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		// Calibration item 11:
		{
			final double fixingTime = 4.25;
			final double periodLength = 0.15;
			final double strikeRate = 0.028;
			final DefaultableCaplet product = new DefaultableCaplet(strikeRate, fixingTime, periodLength, 0, 0, notional, true);
			
			final RandomVariable targetValue = product.getValue(0.0, valuationModel);
			
			final double weight = 1.0 / totalWeights;
			
			final String name = "Floorlet: {Strike: 0.028, Fixing: 4.25, Payment: 4.4, Notional: 1 }";
			
			results[index++] = new CalibrationProduct(name, product, targetValue, weight);
		}
		
		
		return results;
	}
}
