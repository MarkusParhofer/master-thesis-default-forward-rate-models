package info.quantlab.masterthesis.defaultablemodels.testing;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntUnaryOperator;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import info.quantlab.debug.Debug;
import info.quantlab.debug.EulerSchemeFromProcessModel;
import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.functional.Functional;
import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.products.DefaultableCaplet;
import info.quantlab.masterthesis.products.DefaultableCapletAnalyticApproximation;
import info.quantlab.masterthesis.products.DefaultableZeroCouponBond;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.TermStructureModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Bond;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.Swaption;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation.StateSpace;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.montecarlo.process.MonteCarloProcessFromProcessModel;
import net.finmath.plots.Named;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

@RunWith(Parameterized.class)
public class ModelFromSpreadTest extends info.quantlab.debug.Time{
	
	@Parameters(name="{0}")
	public static Collection<Object[]> generateData()
	{
		return Arrays.asList(new Object[][] {
			// Put here an array of arrays where each array represents input for the constructor
			{"Run 0: Baseline",						0.01, 	2, "SPREADS", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 0},
			{"Run 1: Modelling defaultable LIBORs", 0.01, 	2, "LIBORS", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 1},
			{"Run 2: Spread is 0", 					0.01, 	2, "SPREADS", 	"SPOT", 	new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 },  2},
			{"Run 3: Rougher time delta", 			0.1, 	2, "SPREADS", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 3},
			{"Run 4: Terminal Measure",			 	0.01, 	2, "SPREADS", 	"TERMINAL",	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 4},
			{"Run 5: More free parameters", 		0.01, 	4, "SPREADS", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 5},
		});
	}
	
	public ModelFromSpreadTest(final String runName, final double simulationTimeDelta, final int extraFactors, final String simulationProduct, final String measure, final double[] initialRatesNonDefaultable, final int run) throws CalculationException {
		model = getValuationModel(simulationTimeDelta, extraFactors, simulationProduct, measure, initialRatesNonDefaultable);
		runCount = run;
		runningName = runName;
	}
	
	
	
	// Simulation parameters:
	private static final int numberOfPaths = 10000;
	private static final int bmSeed = 3124;
	
	// Non Defaultable Model Parameters:
	private static final int numberOfFactors = 5; // No Factor reduction
	private static final double[] initialRatesNonDefaultable = { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 };
	
	// General Parameters:
	private static final double liborPeriodLength = 2.0;
	private static final int numberOfLiborPeriods = 5;

	private final static String stateSpace = "NORMAL"; /* No support for the lognormal model (of the defaultable LIBOR) supplied yet*/
	
	private static String runningName;
	private static int runCount = 0; /* For saving plots*/
	
	
	private final LIBORModelMonteCarloSimulationModel model; /*Simulation model*/
	

	
	private static DecimalFormat formatterMaturity	= new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue		= new DecimalFormat(" 00.000%;-00.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	// Unnecessary: private static DecimalFormat formatterMoneyness	= new DecimalFormat(" 000.0%;-000.0%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	// Unnecessary: private static DecimalFormat formatterLong	= new DecimalFormat(" ##0.000;-##0.000", new DecimalFormatSymbols(Locale.ENGLISH));

	

	@Test
	public void testGeneralSpecs() throws CalculationException, InterruptedException, IOException {
		final int numberOfLIBORs = model.getNumberOfLibors();
		DefaultableLIBORFromSpreadDynamic defModel = (DefaultableLIBORFromSpreadDynamic)model.getModel();
		
		// TODO: Plot 10 Paths of Defaultable model and of undefaultable model and the spread
		
		
		// Print Drift at t=0:
		/*RandomVariable[] initials = defModel.getInitialValue(model.getProcess());
		RandomVariable[] drift = defModel.getDrift(model.getProcess(), 0, initials, null);
		System.out.print("Drift Libor at t=0:   ");
		for (int i = 0; i < drift.length; i++) {
			if(i == numberOfLiborPeriods)
				System.out.print("\nDrift Spread at t=0:  ");
			System.out.print(formatterLong.format(drift[i] == null? 0 : drift[i].doubleValue()) + "  ");
		}
		System.out.println("\n");*/
		
		// Plot Average of Defaultable model and of undefaultable model and the spread
		DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[3 * numberOfLIBORs];
		final DoubleToIntFunction getTimeIndex = operand -> {
			int timeIndex = model.getTimeIndex(operand);
			return timeIndex < 0 ? - timeIndex - 2 : timeIndex;
		};
		
		for(int j = 0; j < numberOfLIBORs; j++) {
			final int compIndex = j;
			myFunctionArray[j] = (operand) -> {
					try {
						return defModel.getUndefaultableLIBOR(model.getProcess(), getTimeIndex.applyAsInt(operand), compIndex).getAverage();
					} catch (CalculationException e) {
						return - 1.0;
					}
			};
			
			myFunctionArray[j + numberOfLIBORs] = (operand) -> {
				try {
					return defModel.getDefaultableLIBOR(model.getProcess(), getTimeIndex.applyAsInt(operand), compIndex).getAverage();
				} catch (CalculationException e) {
					return - 1.0;
				}
			};
			
			myFunctionArray[j + 2 * numberOfLIBORs] = (operand) -> {
				try {
					return defModel.getLIBORSpreadAtGivenTime(model.getProcess(), operand, compIndex).getAverage();
				} catch (CalculationException e) {
					return - 1.0;
				}
			};
		}
		
		List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
		List<Named<DoubleUnaryOperator>> myFunctionList =
				myList.stream().map(operator -> new Named<DoubleUnaryOperator>(
						(myList.indexOf(operator) / 5d < 1.0 ? "LIBOR " : myList.indexOf(operator) / 5d < 2.0 ? "Defaultable " : "Spread ") + 
						(myList.indexOf(operator) % 5), operator)).
				collect(Collectors.toList());
		EasyPlot2D plot = new EasyPlot2D(model.getTime(0), model.getTimeDiscretization().getLastTime(), 51, myFunctionList);
		plot.setTitle("Average of LIBORs and Spreads (" + runningName + ")");
		plot.setIsLegendVisible(true);
		
		boolean savePlots = false;
		
		if(savePlots) {
			String timeStamp = new SimpleDateFormat("MM-dd-yyyy_HH-mm-ss-SSSS").format(new Date());
			String pathString = "Graphs\\" + timeStamp + "\\";
			File directory = new File(pathString);
			if(!directory.exists())
				directory.mkdirs();
			File picFile = new File(pathString + "Run_" + runCount + ".png");
			plot.saveAsPNG(picFile, 1200, 1200);
			File specsFile = new File(pathString + "Run_" + runCount + "_Specs.txt");
			if(!writeSpecsToFile(specsFile))
				System.out.println("\nSaving Specs failed.");
		}
		else
			plot.show();
		
		
		
		// Print Free Parameter Matrix
		DefaultableLIBORMarketModel defTheoModel = (DefaultableLIBORFromSpreadDynamic)model.getModel();
		/*DefaultableLIBORCovarianceWithGuaranteedPositiveSpread covModel = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)(defTheoModel.getCovarianceModel());
		double[][] freeParameters = covModel.getFreeParameterMatrix();
		System.out.println("Matrix of non defaultable FL      |      Free Parameters:");
		for(int row = 0; row < freeParameters.length; row++) {
			RandomVariable[] flMatrixNonDef = covModel.getNonDefaultableCovarianceModel().getFactorLoading(0.0, row, null);
			for(int col = 0; col < flMatrixNonDef.length; col++) {
				System.out.printf("%9.4f      ", flMatrixNonDef[col].doubleValue());
			}
			System.out.print("|      ");
			for(int col = 0; col < freeParameters[row].length; col++) {
				System.out.printf("%9.4f      ", freeParameters[row][col]);
			}
			System.out.println();
		}
		
		
		System.out.println("\nPath 0 of all LIBORs");*/
		final int maxTimeIndexToPrint = 10;
		final TimeDiscretization timesForPrinting = new TimeDiscretizationFromArray(0.0, maxTimeIndexToPrint + 1, model.getTimeDiscretization().getLastTime() / maxTimeIndexToPrint);
		final IntUnaryOperator getModelTimeIndex = (printIndex) -> {
			final int modelIndex = model.getTimeIndex(timesForPrinting.getTime(printIndex));
			if(modelIndex < 0)
				return - modelIndex - 2;
			else
				return modelIndex;
		};
		System.out.printf(" Time:%8s ", " ");
		for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
			System.out.printf("%2d: %6.4f     ", timeIndex, timesForPrinting.getTime(timeIndex));
		}
		System.out.println();
		System.out.println();/*
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "NonDefRate "+ row);
			for(int timeIndex = 0; timeIndex < maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", defTheoModel.getUndefaultableLIBOR(model.getProcess(),timeIndex, row).get(0));
			}
			System.out.println();
		}
		System.out.println();
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex < maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(timeIndex, row).get(0));
			}
			System.out.println();
		}
		System.out.println();
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "Spread Rate "+ row);
			for(int timeIndex = 0; timeIndex < maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", defTheoModel.getLIBORSpreadAtGivenTimeIndex(model.getProcess(), timeIndex, row).get(0));
			}
			System.out.println();
		}
		*/
		System.out.println("\nAverage, Minimum and Maximum (" + runningName + "):");
		
		System.out.println("Mean:");
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(getModelTimeIndex.applyAsInt(timeIndex), row).getAverage());				
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Min:");
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(getModelTimeIndex.applyAsInt(timeIndex), row).getMin());				
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Max:");
		for (int row = 0; row < numberOfLiborPeriods; row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(getModelTimeIndex.applyAsInt(timeIndex), row).getMax());				
			}
			System.out.println();
		}
		
		System.out.println();
				
		System.out.println("\n" + "_".repeat(300) + "\n");
		
		Assert.assertTrue(true);
	}
	
	@Test
 	public void testPositivity() throws CalculationException {
		
		Assert.assertTrue(true);
	}
	
	@Test
 	public void testCaplet() throws CalculationException {
		/*
		 * Value a caplet
		 */
		System.out.println("\nCaplet prices (" + runningName + "):\n");
		System.out.println("Maturity       Strike          Simulation     Analytic       Deviation      Non Def        Non Def Alt");

		double maxAbsDeviation = 0.0;
		for (double maturity = 2.0; maturity <= liborPeriodLength * (numberOfLiborPeriods - 1); maturity += liborPeriodLength) {
			System.out.print(formatterMaturity.format(maturity) + "          ");
			
			final double strikeRate = getParSwaprate(model, new double[] { maturity, maturity + liborPeriodLength});
			System.out.print(formatterValue.format(strikeRate) + "        ");
			
			DefaultableCaplet defCaplet = new DefaultableCaplet(strikeRate, maturity, liborPeriodLength, true, false);
			final double defCapletSimValue = defCaplet.getValue(model);
			System.out.print(formatterValue.format(defCapletSimValue) + "       ");
			
			DefaultableCapletAnalyticApproximation defCapletAnalytic = new DefaultableCapletAnalyticApproximation(strikeRate, maturity, liborPeriodLength, true, false);
			final double defCapletAnalyticValue = defCapletAnalytic.getValue(model);
			System.out.print(formatterValue.format(defCapletAnalyticValue) + "      ");
			
			final double absDeviation = Math.abs(defCapletAnalyticValue - defCapletSimValue);
			System.out.print(formatterDeviation.format(absDeviation) + "    ");
			
			DefaultableCaplet nonDefCaplet = new DefaultableCaplet(strikeRate, maturity, liborPeriodLength, false, false);
			final double nonDefCapletValue = nonDefCaplet.getValue(model);
			System.out.print(formatterValue.format(nonDefCapletValue) + "       ");
			
			Caplet nonDefCapletAlt = new Caplet(maturity, liborPeriodLength, strikeRate, false);
			final double nonDefCapletValueAlt = nonDefCapletAlt.getValue(getNonDefaultableValuationModel());
			System.out.print(formatterValue.format(nonDefCapletValueAlt) + "       ");
			
			DefaultableCapletAnalyticApproximation nonDefCapletAnalytic = new DefaultableCapletAnalyticApproximation(strikeRate, maturity, liborPeriodLength, false, false);
			final double nonDefCapletAnalyticValue = nonDefCapletAnalytic.getValue(model);
			System.out.println(formatterValue.format(nonDefCapletAnalyticValue));
			
			maxAbsDeviation = Math.max(maxAbsDeviation, absDeviation); 
		}

		System.out.println("Maximum abs deviation: " + formatterDeviation.format(maxAbsDeviation));
		System.out.println("_".repeat(300) + "\n");
		
		/*
		 * jUnit assertion: condition under which we consider this test successful
		 */
		Assert.assertTrue(Math.abs(maxAbsDeviation) < 8E-3);
	}
	
	
	@Test
	public void testBond() throws CalculationException {
		System.out.println("\nTesting Bond Prices (" + runningName + "):\nTime     Simulation       Analytic      Deviation      NonDefBond      NonDefReal      Deviation");
		double absDev = 0.0;
		LIBORModelMonteCarloSimulationModel nonDefModel = getNonDefaultableValuationModel();
		for(double maturity = 2.0; maturity < liborPeriodLength * (numberOfLiborPeriods - 1); maturity+= liborPeriodLength) {
			DefaultableZeroCouponBond bond = new DefaultableZeroCouponBond(maturity);
			final Bond normalBond = new Bond(maturity);
			final double value = bond.getValue(model);
			final double normalValue = normalBond.getValue(nonDefModel);
			double realValue = 1.0;
			// Bond price analytic
			double normalAnalyticValue = 1.0;

			final double lastPeriodIndex = nonDefModel.getLiborPeriodIndex(maturity) - 1;
			for(int periodIndex=0; periodIndex<=lastPeriodIndex; periodIndex++) {
				realValue /= 1.0 + model.getLIBOR(0, periodIndex).doubleValue() * (model.getLiborPeriod(periodIndex+1) - model.getLiborPeriod(periodIndex));
				normalAnalyticValue /= 1.0 + nonDefModel.getLIBOR(0, periodIndex).doubleValue() * (nonDefModel.getLiborPeriod(periodIndex+1) - nonDefModel.getLiborPeriod(periodIndex));
			}

			
			final double deviation = realValue - value;
			final double normalDev = normalAnalyticValue - normalValue;
			System.out.println(formatterMaturity.format(maturity) + "      " + formatterValue.format(value) + "      " + formatterValue.format(realValue) + "      " + formatterDeviation.format(deviation) + "      " + formatterValue.format(normalValue) + "      " + formatterValue.format(normalAnalyticValue) + "      " + formatterDeviation.format(normalDev));
			absDev = Math.max(absDev, Math.abs(deviation));
		}
		
		System.out.println("\nMaximal absolute Deviation: " + formatterDeviation.format(absDev));
		
		System.out.println("\n" + "_".repeat(300) + "\n");
		
		Assert.assertTrue(absDev < 1E-3);
	}
	
	private static LIBORModelMonteCarloSimulationModel getValuationModel(double simulationTimeDelta, int numberOfExtraFactors, String simulationProduct, String measure, double[] initialRatesDefaultable) throws CalculationException {
		
		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, numberOfLiborPeriods, liborPeriodLength); // Fixing time
		
		// Set initial forward curves
		final double[] fixingTimes = new double[] { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 };
		if(initialRatesDefaultable.length != fixingTimes.length || initialRatesNonDefaultable.length != fixingTimes.length)
			throw new IllegalArgumentException("Initial rates must have length " + fixingTimes.length);
		
		final ForwardCurve nonDefaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"nonDefaultableForwardCurve",
				fixingTimes,
				initialRatesNonDefaultable,
				liborPeriodLength);
		
		/*
		 * Create a simulation time discretization. We save space by just simulationg to the last fixing time.
		 */
		final int numberOfTimes = (int)(liborPeriods.getTime(liborPeriods.getNumberOfTimeSteps() - 1)/simulationTimeDelta);
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, numberOfTimes, simulationTimeDelta);

		/*
		 * Create volatility structure v[i][j] = sigma_j(t_i) and correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		final double a = 0.1, b = 0.0, c = 0.25, d = 0.1;
		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriods, a, b, c, d, true);
		final double correlationDecayParam = 0.2;
		final LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriods, numberOfFactors,	correlationDecayParam, true);
		
		final LIBORCovarianceModel baseCovarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriods, volatilityModel, correlationModel);
		
		// Set model properties
		final Map<String, String> properties = new HashMap<>();

		// Choose the simulation measure
		properties.put("measure", measure);

		// Choose normal model
		properties.put("stateSpace", stateSpace);

		// Empty array of calibration items - hence, model will use given covariance
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		final LIBORMarketModel baseLiborMarketModel = LIBORMarketModelFromCovarianceModel.of(
				liborPeriods, 		/* LIBORPeriodDiscretization */
				null,			 	/* analyticModel */
				nonDefaultableForwards, /* ForwardCurve */
				new DiscountCurveFromForwardCurve(nonDefaultableForwards),  /* DiscountCurve */
				new RandomVariableFromArrayFactory(true),  /* RV Factory */
				baseCovarianceModel,  /* covarianceModel */
				calibrationItems,  /* calibrationItems */
				properties);  /* properties */
		
		final ForwardCurve defaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"defaultableForwardCurve",
				fixingTimes,
				initialRatesDefaultable,
				liborPeriodLength);
		
		properties.put("simulationModel", simulationProduct);
		
		final DefaultableLIBORCovarianceModel defaultableCovariance = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(baseCovarianceModel, generateFreeParameterMatrix(numberOfExtraFactors));
		
		final DefaultableLIBORFromSpreadDynamic defaultableModel = new DefaultableLIBORFromSpreadDynamic(baseLiborMarketModel, defaultableCovariance, defaultableForwards, properties);
		
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers(timeDiscretization, defaultableModel.getNumberOfFactors(), numberOfPaths, bmSeed);
		
		EulerSchemeFromProcessModel.Scheme eulerScheme = EulerSchemeFromProcessModel.Scheme.EULER;
		
		if(simulationProduct.toUpperCase() == "SPREADS")
			eulerScheme = EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL;
		
		final MonteCarloProcessFromProcessModel process = new EulerSchemeFromProcessModel(defaultableModel, brownianMotion, eulerScheme);

		return new LIBORMonteCarloSimulationFromLIBORModel(process);
	}
	
	private static double[][] generateFreeParameterMatrix(int numberOfExtraFactors) {
		double[][] matrix = new double[numberOfLiborPeriods][numberOfExtraFactors];
		MersenneTwister randomGenerator = new MersenneTwister(1021);
		for (int row = 0; row < matrix.length; row++) {
			for (int col = 0; col < matrix[row].length; col++) {
				matrix[row][col] = randomGenerator.nextDoubleFast() - 0.5d;
			}
		}
		return matrix;
	}
	
	private static double getParSwaprate(final LIBORModelMonteCarloSimulationModel liborMarketModel, final double[] swapTenor) {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretizationFromArray(swapTenor), new TimeDiscretizationFromArray(swapTenor), liborMarketModel.getModel().getForwardRateCurve(), liborMarketModel.getModel().getDiscountCurve());
	}
	
	private LIBORModelMonteCarloSimulationModel getNonDefaultableValuationModel() {
		
		final MonteCarloProcess normalProcess = model.getProcess();
		final DefaultableLIBORMarketModel defaultableTheoModel = (DefaultableLIBORMarketModel)(model.getModel());
		
		LIBORModelMonteCarloSimulationModel nonDefSimModel = new LIBORModelMonteCarloSimulationModel() {

			@Override
			public RandomVariable getForwardRate(double time, double periodStart, double periodEnd) throws CalculationException {
				return defaultableTheoModel.getUndefaultableForwardRate(normalProcess, time, periodStart, periodEnd);
			}

			@Override
			public RandomVariable getNumeraire(double time) throws CalculationException {
				return defaultableTheoModel.getNumeraire(normalProcess, time);
			}

			@Override
			public TermStructureModel getModel() {
				return defaultableTheoModel.getUndefaultableLIBORModel();
			}

			@Override
			public MonteCarloProcess getProcess() {
				return Functional.getComponentReducedMCProcess(normalProcess, 0, getModel().getNumberOfComponents() - 1);
			}

			@Override
			public Object getCloneWithModifiedSeed(int seed) {
				return null;
			}

			@Override
			public int getNumberOfPaths() {
				return numberOfPaths;
			}

			@Override
			public TimeDiscretization getTimeDiscretization() {
				return normalProcess.getTimeDiscretization();
			}

			@Override
			public double getTime(int timeIndex) {
				return normalProcess.getTime(timeIndex);
			}

			@Override
			public int getTimeIndex(double time) {
				return normalProcess.getTimeIndex(time);
			}

			@Override
			public RandomVariable getRandomVariableForConstant(double value) {
				return defaultableTheoModel.getRandomVariableForConstant(value);
			}

			@Override
			public RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException {
				return normalProcess.getMonteCarloWeights(timeIndex);
			}

			@Override
			public RandomVariable getMonteCarloWeights(double time) throws CalculationException {
				int timeIndex = getTimeIndex(time);
				if(timeIndex < 0) {
					timeIndex = - timeIndex - 2;
				}
				return normalProcess.getMonteCarloWeights(timeIndex);
			}

			@Override
			public MonteCarloSimulationModel getCloneWithModifiedData(Map<String, Object> dataModified)
					throws CalculationException {
				return null;
			}

			@Override
			public TimeDiscretization getLiborPeriodDiscretization() {
				return null;
			}

			@Override
			public int getNumberOfLibors() {
				return defaultableTheoModel.getUndefaultableLIBORModel().getNumberOfLibors();
			}

			@Override
			public double getLiborPeriod(int timeIndex) {
				return defaultableTheoModel.getLiborPeriod(timeIndex);
			}

			@Override
			public int getLiborPeriodIndex(double time) {
				return defaultableTheoModel.getLiborPeriodIndex(time);
			}

			@Override
			public RandomVariable getLIBOR(int timeIndex, int liborIndex) throws CalculationException {
				return defaultableTheoModel.getUndefaultableLIBOR(normalProcess, timeIndex, liborIndex);
			}

			@Override
			public RandomVariable[] getLIBORs(int timeIndex) throws CalculationException {
				RandomVariable[] libors = new RandomVariable[getNumberOfLibors()];
				for(int i = 0; i < getNumberOfLibors(); i++) {
					libors[i] = getLIBOR(timeIndex, i);
				}
				return libors;
			}
			
		};
		
		return nonDefSimModel;
	}
	
	private boolean writeSpecsToFile(File file) {
		boolean result = true;
		try {
			result = file.createNewFile();
		} catch(IOException ioEx) {
			Debug.logln(ioEx.getMessage());
			return false;
		} catch(SecurityException secEx) {
			Debug.logln(secEx.getMessage());
			Debug.logln("Access to File Denied!");
			return false;
		}
		if(!result)
			return false;
		
		try {
			FileWriter myWriter = new FileWriter(file);
			myWriter.write("Model Specs:\n\n");
			DecimalFormat formatLongDouble	= new DecimalFormat("##0.00000000000000;-##0.00000000000000", new DecimalFormatSymbols(Locale.ENGLISH));
			DecimalFormat formatterDouble	= new DecimalFormat("##0.00000;-##0.00000", new DecimalFormatSymbols(Locale.ENGLISH));
			
			{
				TimeDiscretization libors = model.getLiborPeriodDiscretization();
				myWriter.write("LiborPeriods:\n" + libors.getFirstTime() + " " + libors.getNumberOfTimeSteps() + " " + libors.getTimeStep(0) + "\n\n");
				TimeDiscretization timeDisc = model.getTimeDiscretization();
				myWriter.write("TimeDiscretization:\n" + timeDisc.getFirstTime() + " " + timeDisc.getNumberOfTimeSteps() + " " + timeDisc.getTimeStep(0) + "\n\n");

			}
			{
				AbstractLIBORCovarianceModelParametric nonDefCovModel = (AbstractLIBORCovarianceModelParametric)((DefaultableLIBORFromSpreadDynamic)model.getModel()).getCovarianceModel().getNonDefaultableCovarianceModel();
				myWriter.write("DefaultableCovarianceModel: " + nonDefCovModel.getClass() + "\n");
				myWriter.write("NumberOfFactors=" + nonDefCovModel.getNumberOfFactors() + "\n");
				myWriter.write("Parameters=[");
				double[] parameters = nonDefCovModel.getParameterAsDouble();
				for(double par: parameters) {
					myWriter.write(formatterDouble.format(par) + " ");
				}
				myWriter.write("]\n\n");
			}
			{
				DefaultableLIBORCovarianceWithGuaranteedPositiveSpread defCovModel = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)(((DefaultableLIBORFromSpreadDynamic)model.getModel()).getCovarianceModel());
				myWriter.write("DefaultableCovarianceModel: " + defCovModel.getClass() + "\n");
				myWriter.write("NumberOfFactors=" + defCovModel.getNumberOfFactors() + "\n");
				myWriter.write("FreeParameters=\n");
				double[][] freeParameters = defCovModel.getFreeParameterMatrix();
				for(int row = 0; row < freeParameters.length; row++) {
					for(int col = 0; col < freeParameters[row].length; col++) {
						myWriter.write(formatLongDouble.format(freeParameters[row][col]) + " ");
					}
					myWriter.write("\n");
				}
				myWriter.write("\n");
			}
			{
				LIBORMarketModelFromCovarianceModel nonDefModel = (LIBORMarketModelFromCovarianceModel)((DefaultableLIBORFromSpreadDynamic)model.getModel()).getUndefaultableLIBORModel();
				List<ForwardCurveInterpolation.Point> forwardPoints = ((ForwardCurveInterpolation)nonDefModel.getForwardRateCurve()).getPoints();
				myWriter.write("NonDefaultableForwardCurve: " + nonDefModel.getForwardRateCurve().getClass() + "\nTimes=[");
				for(ForwardCurveInterpolation.Point point: forwardPoints) {
					myWriter.write(formatterDouble.format(point.getTime()) + " ");
				}
				myWriter.write("]\nValues=[");
				for(ForwardCurveInterpolation.Point point: forwardPoints) {
					myWriter.write(formatterDouble.format(point.getValue()) + " ");
				}
				myWriter.write("]\n\n");
				myWriter.write("NonDefaultableModel: " + nonDefModel.getClass() + "\n");
				myWriter.write("Measure=" + nonDefModel.getMeasure() + "\n");
				myWriter.write("InterpolationMethod=" + nonDefModel.getInterpolationMethod() + "\n");
				myWriter.write("StateSpace=" + stateSpace + "\n");
				myWriter.write("DriftApproximationMethod=" + nonDefModel.getDriftApproximationMethod() + "\n\n");
			}
			{
				DefaultableLIBORFromSpreadDynamic defModel = (DefaultableLIBORFromSpreadDynamic)model.getModel();
				List<ForwardCurveInterpolation.Point> forwardPoints = ((ForwardCurveInterpolation)defModel.getForwardRateCurve()).getPoints();
				myWriter.write("DefaultableForwardCurve: " + defModel.getForwardRateCurve().getClass() + "\nTimes=[");
				for(ForwardCurveInterpolation.Point point: forwardPoints) {
					myWriter.write(formatterDouble.format(point.getTime()) + " ");
				}
				myWriter.write("]\nValues=[");
				for(ForwardCurveInterpolation.Point point: forwardPoints) {
					myWriter.write(formatterDouble.format(point.getValue()) + " ");
				}
				myWriter.write("]\n\n");
				myWriter.write("ProcessModel: " + defModel.getClass() + "\n");
				myWriter.write("Measure=" + defModel.getMeasure() + "\n");
				myWriter.write("InterpolationMethod=" + defModel.getInterpolationMethod() + "\n");
				myWriter.write("StateSpace=" + defModel.getStateSpace() + "\n");
				myWriter.write("HandleSimulationTime=" + defModel.getHandleSimulationTime() + "\n");
				myWriter.write("SimulationModel=" + defModel.getSimulationModel() + "\n\n");
			}
			{
				BrownianMotion bm = model.getBrownianMotion();
				myWriter.write("BrownianMotion: " + bm.getClass() + "\n");
				myWriter.write("NumberOfPaths=" + bm.getNumberOfPaths() + "\n");
				myWriter.write("Seed=" + bmSeed + "\n\n");
			}
			{
				EulerSchemeFromProcessModel euler = (EulerSchemeFromProcessModel)model.getProcess();
				myWriter.write("MCProcess: " + euler.getClass() + "\n");
				myWriter.write("Scheme=" + euler.getScheme().name() + "\n\n");
			}
			myWriter.write("SimulationModel: " + model.getClass() + "\n\n");
			
			{
				LIBORMarketModelFromCovarianceModel nonDefModel = (LIBORMarketModelFromCovarianceModel)((DefaultableLIBORFromSpreadDynamic)model.getModel()).getUndefaultableLIBORModel();
				myWriter.write(nonDefModel.toString() + "\n\n");
				
			}
			
			myWriter.close();
			System.out.println("Successfully wrote to the file.");
		} catch (IOException e) {
			System.out.println("An error occurred.");
			return false;
		}
		
		return true;
	}
}
