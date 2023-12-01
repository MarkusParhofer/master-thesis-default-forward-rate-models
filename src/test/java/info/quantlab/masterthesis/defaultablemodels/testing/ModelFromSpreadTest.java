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
import java.util.Set;
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
import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.functional.FunctionsOnMCProcess;
import info.quantlab.masterthesis.functional.FunctionsOnRandomVariables;
import info.quantlab.masterthesis.products.DefaultableCaplet;
import info.quantlab.masterthesis.products.DefaultableCapletAnalyticApproximation;
import info.quantlab.masterthesis.products.DefaultableZeroCouponBond;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.TermStructureModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.products.Bond;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation.StateSpace;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.plots.Named;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.util.TriFunction;

@RunWith(Parameterized.class)
public class ModelFromSpreadTest extends info.quantlab.debug.Time{
	
	@Parameters(name="{0}")
	public static Collection<Object[]> generateData()
	{
		return Arrays.asList(new Object[][] {
			// Put here an array of arrays where each array represents input for the constructor
			{"Run 0: Baseline",						0.01, 		2, "SPREADS",	"EULER_FUNCTIONAL",			"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 0},
			
			// For now the most stable version:
			//{"Run 1: Modelling defaultable LIBORs", 		 0.001, 	2, "LIBORS", 	"EULER",					"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 1},
			
			{"Run 2: Modelling Spreads Milstein", 	0.01, 		2, "SPREADS", 	"MILSTEIN_FDCENTRAL",		"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 2},
			//{"Run 2: Spread with normal model",		0.001, 	2, "SPREADS",	"NORMAL",		"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 0},
			/*{"Run 2: Spread is 0", 					0.01, 	2, "SPREADS",	"LOGNORMAL", 	"SPOT", 	new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 },  2},
			{"Run 3: Rougher time delta", 			0.1, 	2, "SPREADS",	"LOGNORMAL", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 3},
			{"Run 4: Terminal Measure",			 	0.01, 	2, "SPREADS",	"LOGNORMAL", 	"TERMINAL",	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 4},
			{"Run 5: More free parameters", 		0.01, 	4, "SPREADS",	"LOGNORMAL", 	"SPOT", 	new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 5},*/
		});
	}
	
	public ModelFromSpreadTest(final String runName, final double simulationTimeDelta, final int extraFactors, final String simulationProduct, final String scheme, final String measure, final double[] initialRatesNonDefaultable, final int run) throws CalculationException {
		model = getValuationModel(simulationTimeDelta, extraFactors, simulationProduct, scheme, measure, initialRatesNonDefaultable);
		runCount = run;
		runningName = runName;
	}

	private static final boolean savePlots = true;
	private static final String stateSpace = "NORMAL";
	private static final int bmSeed = 7824;

	private static String runningName;
	private static int runCount = 0; /* For saving plots*/
	private static final int timeStepJumpForPos = 2;

	
	private static DefaultableLIBORModelFactory factory;
	private final LIBORModelMonteCarloSimulationModel model; /*Simulation model*/
	

	
	private static final DecimalFormat formatterMaturity	= new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static final DecimalFormat formatterValue		= new DecimalFormat(" 00.000%;-00.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	// Unnecessary: private static DecimalFormat formatterMoneyness	= new DecimalFormat(" 000.0%;-000.0%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static final DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	// Unnecessary: private static DecimalFormat formatterLong	= new DecimalFormat(" ##0.000;-##0.000", new DecimalFormatSymbols(Locale.ENGLISH));

	

	@Test
	public void testGeneralSpecs() throws CalculationException, IOException {
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
		
		
		
		EasyPlot2D plotAverages;
		{
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
					myList.stream().map(operator -> new Named<>(
							(myList.indexOf(operator) / 5d < 1.0 ? "LIBOR " : myList.indexOf(operator) / 5d < 2.0 ? "Defaultable " : "Spread ") + 
							(myList.indexOf(operator) % 5), operator)).
					collect(Collectors.toList());
			plotAverages = new EasyPlot2D(model.getTime(0), model.getTimeDiscretization().getLastTime(), 51, myFunctionList);
			plotAverages.setTitle("Average of LIBORs and Spreads (" + runningName + ")");
			plotAverages.setIsLegendVisible(true);
		}
		
		if(savePlots) {
			final int numberOfPathsToPlot = 100;
			// If we don't have to look at the plots now we can dump all Infos we have into different plots:
			final DoubleToIntFunction getTimeIndex = operand -> {
				int timeIndex = model.getTimeIndex(operand);
				return timeIndex < 0 ? - timeIndex - 2 : timeIndex;
			};

			EasyPlot2D[] plotPathsOfNonDefLIBOR = new EasyPlot2D[model.getNumberOfLibors()];
			EasyPlot2D[] plotPathsOfDefLIBOR = new EasyPlot2D[model.getNumberOfLibors()];
			EasyPlot2D[] plotPathsOfSpread = new EasyPlot2D[model.getNumberOfLibors()];
			for(int liborIndex = 0; liborIndex < model.getNumberOfLibors(); liborIndex++) {
				final int libor = liborIndex;
				// Plot Average of Defaultable model and of undefaultable model and the spread
				DoubleUnaryOperator[] fArrayNonDefModel = new DoubleUnaryOperator[numberOfPathsToPlot];
				DoubleUnaryOperator[] fArrayDefModel = new DoubleUnaryOperator[numberOfPathsToPlot];
				DoubleUnaryOperator[] fArraySpreadModel = new DoubleUnaryOperator[numberOfPathsToPlot];
				
				for(int path = 0; path < numberOfPathsToPlot; path++) {
					final int pathIndex = path;
					fArrayNonDefModel[path] = (operand) -> {
							try {
								return defModel.getUndefaultableLIBOR(model.getProcess(), getTimeIndex.applyAsInt(operand), libor).get(pathIndex);
							} catch (CalculationException e) {
								return - 1.0;
							}
					};
					
					fArrayDefModel[path] = (operand) -> {
						try {
							return defModel.getDefaultableLIBOR(model.getProcess(), getTimeIndex.applyAsInt(operand), libor).get(pathIndex);
						} catch (CalculationException e) {
							return - 1.0;
						}
					};
					
					fArraySpreadModel[path] = (operand) -> {
						try {
							return defModel.getLIBORSpreadAtGivenTime(model.getProcess(), operand, libor).get(pathIndex);
						} catch (CalculationException e) {
							return - 1.0;
						}
					};
				}
				List<DoubleUnaryOperator> fListNonDefModel = Arrays.asList(fArrayNonDefModel);
				List<Named<DoubleUnaryOperator>> fListNamedNonDefModel =
						fListNonDefModel.stream().map(operator -> new Named<>(
                                "Path " + fListNonDefModel.indexOf(operator), operator)).collect(Collectors.toList());
				plotPathsOfNonDefLIBOR[libor] = new EasyPlot2D(model.getTime(0), model.getTimeDiscretization().getLastTime(), 51, fListNamedNonDefModel);
				plotPathsOfNonDefLIBOR[libor].setTitle("Sample Paths non Defaultable LIBOR " + libor + " (" + runningName + ")");
				
				List<DoubleUnaryOperator> fListDefModel = Arrays.asList(fArrayDefModel);
				List<Named<DoubleUnaryOperator>> fListNamedDefModel =
						fListDefModel.stream().map(operator -> new Named<>("Path " + fListDefModel.indexOf(operator), operator)).collect(Collectors.toList());
				plotPathsOfDefLIBOR[libor] = new EasyPlot2D(model.getTime(0), model.getTimeDiscretization().getLastTime(), 51, fListNamedDefModel);
				plotPathsOfDefLIBOR[libor].setTitle("Sample Paths Defaultable LIBOR " + libor + " (" + runningName + ")");
				
				List<DoubleUnaryOperator> fListSpreadModel = Arrays.asList(fArraySpreadModel);
				List<Named<DoubleUnaryOperator>> fListNamedSpreadModel =
						fListSpreadModel.stream().map(operator -> new Named<>(
                                "Path " + fListSpreadModel.indexOf(operator), operator)).collect(Collectors.toList());
				plotPathsOfSpread[libor] = new EasyPlot2D(model.getTime(0), model.getTimeDiscretization().getLastTime(), 51, fListNamedSpreadModel);
				plotPathsOfSpread[libor].setTitle("Sample Paths Spread " + libor + " (" + runningName + ")");
			}
			
			String timeStamp = new SimpleDateFormat("MM-dd-yyyy_HH-mm-ss-SSSS").format(new Date());
			String pathString = "Graphs\\" + timeStamp + "_" + runningName.substring(7) + "\\";
			File directory = new File(pathString);
			boolean directoryExists = directory.exists();
			if(!directoryExists)
				directoryExists = directory.mkdirs();

			if(directoryExists) {
				File picFile = new File(pathString + "Run_" + runCount + ".png");
				plotAverages.saveAsPNG(picFile, 1200, 1200);
				for (int liborIndex = 0; liborIndex < model.getNumberOfLibors(); liborIndex++) {
					File pathPicFileNonDefLIBOR = new File(pathString + "Sample Paths non Defaultable LIBOR " + liborIndex + ".png");
					plotPathsOfNonDefLIBOR[liborIndex].saveAsPNG(pathPicFileNonDefLIBOR, 1200, 1200);

					File pathPicFileDefLIBOR = new File(pathString + "Sample Paths Defaultable LIBOR " + liborIndex + ".png");
					plotPathsOfDefLIBOR[liborIndex].saveAsPNG(pathPicFileDefLIBOR, 1200, 1200);

					File pathPicFileSpread = new File(pathString + "Sample Paths Spread " + liborIndex + ".png");
					plotPathsOfSpread[liborIndex].saveAsPNG(pathPicFileSpread, 1200, 1200);
				}


				File specsFile = new File(pathString + "Run_" + runCount + "_Specs.txt");
				if (!writeSpecsToFileNew(specsFile))
					System.out.println("\nSaving Specs failed.");
			}
			else
				System.out.println("Cannot create directory for Graphs");
		}
		else
			plotAverages.show();
		
		
		
		// Print Free Parameter Matrix
		/*DefaultableLIBORMarketModel defTheoModel = (DefaultableLIBORFromSpreadDynamic)model.getModel();
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread covModel = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)(defTheoModel.getCovarianceModel());
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
		for (int row = 0; row < model.getNumberOfLibors(); row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(getModelTimeIndex.applyAsInt(timeIndex), row).getAverage());				
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Min:");
		for (int row = 0; row < model.getNumberOfLibors(); row++) {
			System.out.printf(" %12s: ", "Def Rate "+ row);
			for(int timeIndex = 0; timeIndex <= maxTimeIndexToPrint; timeIndex+=1) {
				System.out.printf("%9.4f      ", model.getLIBOR(getModelTimeIndex.applyAsInt(timeIndex), row).getMin());				
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Max:");
		for (int row = 0; row < model.getNumberOfLibors(); row++) {
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
		System.out.println("\n\nTesting Spread for Positivity (" + runningName + "):\n");
		double minimum = Double.POSITIVE_INFINITY;
		int timeIndexOfLowestSpread = -1;
		int liborIndexOfLowestSpread = -1;
		int pathIndexOfLowestSpread = -1;
		System.out.println("LIBOR -index/-time | Minimum Spread      time      path");
		DefaultableLIBORMarketModel defModel = (DefaultableLIBORMarketModel)model.getModel();
		for(int libor = 0; libor < model.getNumberOfLibors(); libor++) {
			double minAtLIBOR = Double.POSITIVE_INFINITY;
			int pathAtLIBOR = 0;
			int timeAtLIBOR = 0;
			final double liborTime = model.getLiborPeriod(libor);
			for(int timeIndex = 0; timeIndex < model.getTimeDiscretization().getNumberOfTimes() / timeStepJumpForPos; timeIndex++) {
				final RandomVariable spreadAtTimeIndex = defModel.getLIBORSpreadAtGivenTimeIndex(model.getProcess(), timeIndex, libor);
				final double minAtTimeIndex = spreadAtTimeIndex.getMin();
				if(minAtTimeIndex < minAtLIBOR) {
					minAtLIBOR = minAtTimeIndex;
					pathAtLIBOR = FunctionsOnRandomVariables.findFirstPathWithValue(spreadAtTimeIndex, minAtTimeIndex);
					timeAtLIBOR = timeIndex;
				}
				if(liborTime <= model.getTime(timeIndex))
					// No need for further Evaluation, Spread stays the same from hereon
					break;
			}
			// Print result for LIBOR:
			System.out.printf("LIBOR %2d    %4.2f   | ", libor, liborTime);
			System.out.printf("%11s        ", formatterDeviation.format(minAtLIBOR));
			System.out.printf("%4d        ", timeAtLIBOR);
			System.out.printf("%6d %n", pathAtLIBOR);
			// Adjust General smallest value
			if(minAtLIBOR < minimum) {
				minimum = minAtLIBOR;
				pathIndexOfLowestSpread = pathAtLIBOR;
				timeIndexOfLowestSpread = timeAtLIBOR;
				liborIndexOfLowestSpread = libor;
			}
			
		}
		
		System.out.println("\nOverall Minimum is:");
		System.out.printf("LIBOR %2d    %4.2f   | ", liborIndexOfLowestSpread, model.getLiborPeriod(liborIndexOfLowestSpread));
		System.out.printf("%11s        ", formatterDeviation.format(minimum));
		System.out.printf("%4d        ", timeIndexOfLowestSpread);
		System.out.printf("%6d", pathIndexOfLowestSpread);
		System.out.println();
		System.out.println("\n" + "_".repeat(300) + "\n");
		Assert.assertTrue(minimum > -1E-5);
	}
	
	@Test
 	public void testCaplet() throws CalculationException {
		System.out.println("\nCaplet prices (" + runningName + "):\n");
		System.out.println("                    Defaultable Caplets                                     |                        Non-Defaultable Caplets");
		System.out.println("Maturity       Strike          Simulation     Analytic       Deviation      |  Simulation         FinMath-Caplet      Analytic      FinMath-Analytic from swaption");

		double maxAbsDeviation = 0.0;
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		for (double maturity = 2.0; maturity <= liborPeriodLength  * (model.getNumberOfLibors() - 1); maturity += liborPeriodLength) {
			System.out.print(formatterMaturity.format(maturity) + "          ");
			
			final double strikeRate = getParSwaprate(model, new double[] { maturity, maturity + liborPeriodLength});
			System.out.print(formatterValue.format(strikeRate) + "        ");
			
			DefaultableCaplet defCaplet = new DefaultableCaplet(strikeRate, maturity, liborPeriodLength, true, false);
			final double simulationValue = defCaplet.getValue(model);
			System.out.print(formatterValue.format(simulationValue) + "       ");
			
			DefaultableCapletAnalyticApproximation defCapletAnalytic = new DefaultableCapletAnalyticApproximation(strikeRate, maturity, liborPeriodLength, true, false);
			final double analyticValue = defCapletAnalytic.getValue(model);
			System.out.print(formatterValue.format(analyticValue) + "      ");
			
			final double absDeviation = Math.abs(analyticValue - simulationValue);
			System.out.print(formatterDeviation.format(absDeviation) + "    |  ");
			
			DefaultableCaplet nonDefCaplet = new DefaultableCaplet(strikeRate, maturity, liborPeriodLength, false, false);
			final double nonDefSimulationValue = nonDefCaplet.getValue(model);
			System.out.print(formatterValue.format(nonDefSimulationValue) + "            ");
			
			Caplet nonDefCapletAlt = new Caplet(maturity, liborPeriodLength, strikeRate, false);
			final double nonDefFinMathSimulationValue = nonDefCapletAlt.getValue(getNonDefaultableValuationModel());
			System.out.print(formatterValue.format(nonDefFinMathSimulationValue) + "          ");
			
			DefaultableCapletAnalyticApproximation nonDefCapletAnalytic = new DefaultableCapletAnalyticApproximation(strikeRate, maturity, liborPeriodLength, false, false);
			final double nonDefAnalyticValue = nonDefCapletAnalytic.getValue(model);
			System.out.print(formatterValue.format(nonDefAnalyticValue) + "       ");
			
			TimeDiscretization swapTenor = new TimeDiscretizationFromArray(maturity, maturity + liborPeriodLength);
			SwaptionGeneralizedAnalyticApproximation nonDefCapletAnalyticFromSwaption = new SwaptionGeneralizedAnalyticApproximation(strikeRate, swapTenor, StateSpace.NORMAL);
			final double nonDefFinMathAnalyticValue = nonDefCapletAnalyticFromSwaption.getValue(getNonDefaultableValuationModel());
			System.out.println(formatterValue.format(nonDefFinMathAnalyticValue));
			
			
			maxAbsDeviation = Math.max(maxAbsDeviation, absDeviation); 
		}

		System.out.println("Maximum abs deviation: " + formatterDeviation.format(maxAbsDeviation));
		System.out.println("_".repeat(300) + "\n");
		
		Assert.assertTrue(Math.abs(maxAbsDeviation) < 8E-3);
	}
		
	@Test
	public void testBond() throws CalculationException {
		System.out.println("\nTesting Bond Prices (" + runningName + "):\nTime     Simulation       Analytic      Deviation      NonDefBond      NonDefReal      Deviation");
		double absDev = 0.0;
		LIBORModelMonteCarloSimulationModel nonDefModel = getNonDefaultableValuationModel();
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		for(double maturity = 2.0; maturity < liborPeriodLength * (model.getNumberOfLibors() - 1); maturity+= liborPeriodLength) {
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
	
	private static LIBORModelMonteCarloSimulationModel getValuationModel(double simulationTimeDelta, int numberOfExtraFactors, String simulationProduct, String scheme, String measure, double[] initialRatesDefaultable) throws CalculationException {
		
		factory = new DefaultableLIBORModelFactory();
		Map<String, Object> properties = new HashMap<>();
		

		properties.put("fixingTimes", new double[] { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 });
		properties.put("liborPeriodLength", 2.0);
		properties.put("numberOfLiborPeriods", 5);
		properties.put("stateSpace", stateSpace);
		properties.put("measure", measure);
		
		properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
		properties.put("numberOfFactors", 5);
		properties.put("volatilityParams", new double[] {0.1, 0.0, 0.25, 0.1});
		properties.put("correlationDecayParam", 0.2);
		properties.put("displacement", 0.6);
		properties.put("initialRatesNonDefaultable", new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 });
		
		properties.put("numberOfExtraFactors", numberOfExtraFactors);
		properties.put("initialRatesDefaultable", initialRatesDefaultable);
		properties.put("simulationModel", simulationProduct);
		properties.put("stateSpaceOfSpread", "NORMAL");
		properties.put("freeParamsSeed", 2000);
		properties.put("freeParamsRange", 0.5);
		// properties.put("freeParamsGenerator", freeParamsGenerator); // None needed because the default one is great! lol...
		properties.put("simulationTimeDelta", simulationTimeDelta);
		properties.put("numberOfPaths", 10000);
		properties.put("brownianMotionSeed", bmSeed);
		properties.put("numericalScheme", scheme);
		properties.put("finiteDifferenceDelta", 1E-6);
		
		// properties.put("analyticDifferentialFactorLoadings", analyticDifferentialFactorLoadings); // None needed, for now we only use FINITE Difference Milstein Scheme

		factory.setProperties(properties);
		
		DefaultableLIBORMarketModel defaultableModel = factory.createDefaultableModel();
				
		TriFunction<Integer, Double, RandomVariable[], RandomVariable[]> analyticDifferentialFactorLoadings;
		if(scheme.equals("MILSTEIN_ANALYTIC")){
			final int numberOfLiborPeriods = (int) properties.get("numberOfLiborPeriods");
			final int numberOfFactors = (int)properties.get("numberOfFactors") + (int)properties.get("numberOfExtraFactors");
			final int numberOfNonDefFactors = (int)properties.get("numberOfFactors");
			final double liborPeriodLength = defaultableModel.getLiborPeriod(1) - defaultableModel.getLiborPeriod(0);
			final DefaultableLIBORCovarianceWithGuaranteedPositiveSpread covModel = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)defaultableModel.getCovarianceModel();
			analyticDifferentialFactorLoadings = (componentIndex, time, realizations) -> {
				
				// Does not work for Blended or displaced model
				if(componentIndex < numberOfLiborPeriods) {
					RandomVariable[] zeros = new RandomVariable[numberOfFactors];
					RandomVariable zero = Scalar.of(0.0d);
                    Arrays.fill(zeros, zero);
					return zeros;
				} else {
					RandomVariable[] differential = covModel.getFactorLoading(time, componentIndex - numberOfLiborPeriods, realizations);
					RandomVariable factor = Scalar.of(liborPeriodLength).div(realizations[componentIndex - numberOfLiborPeriods].mult(liborPeriodLength).add(1.0d));
					for(int k = 0; k < numberOfNonDefFactors; k++) {
						differential[k] = differential[k].mult(factor);
					}
					for(int k = numberOfNonDefFactors; k < numberOfFactors; k++) {
						differential[k] = Scalar.of(covModel.getFreeParameterMatrix()[componentIndex - numberOfLiborPeriods][k - numberOfNonDefFactors]);
					}
					return differential;
				}
			};
			properties.put("analyticDifferentialFactorLoadings", analyticDifferentialFactorLoadings);

		}
		
		factory.setProperties(properties);
		
		final MonteCarloProcess process = factory.createNumericalScheme(defaultableModel);

		return new LIBORMonteCarloSimulationFromLIBORModel(process);
	}
	
	private static double getParSwaprate(final LIBORModelMonteCarloSimulationModel liborMarketModel, final double[] swapTenor) {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretizationFromArray(swapTenor), new TimeDiscretizationFromArray(swapTenor), liborMarketModel.getModel().getForwardRateCurve(), liborMarketModel.getModel().getDiscountCurve());
	}
	
	private LIBORModelMonteCarloSimulationModel getNonDefaultableValuationModel() {
		
		final MonteCarloProcess normalProcess = model.getProcess();
		final DefaultableLIBORMarketModel defaultableTheoModel = (DefaultableLIBORMarketModel)(model.getModel());

        return new LIBORModelMonteCarloSimulationModel() {

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
				return FunctionsOnMCProcess.getComponentReducedMCProcess(normalProcess, 0, getModel().getNumberOfComponents() - 1);
			}

			@Override
			public Object getCloneWithModifiedSeed(int seed) {
				return null;
			}

			@Override
			public int getNumberOfPaths() {
				return model.getNumberOfPaths();
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
			public MonteCarloSimulationModel getCloneWithModifiedData(Map<String, Object> dataModified) {
				return null;
			}

			@Override
			public TimeDiscretization getLiborPeriodDiscretization() {
				return defaultableTheoModel.getLiborPeriodDiscretization();
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
	}
	
	
	private boolean writeSpecsToFileNew(File file) {
		boolean result;
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
			Map<String, Object> properties = factory.getProperties();
			Set<String> keys = properties.keySet();
			for(String key: keys) {
				// TODO: Handle Lambda exressions
				if(properties.get(key) instanceof double[] vec) {
					myWriter.write(key + ": {");
					for(int i=0; i < vec.length; i++) {
						String ender = i == vec.length - 1? "}\n" : ", ";
						myWriter.write(vec[i] + ender);
					}
				}
				else {
					myWriter.write(key + ": {" + properties.get(key) + "}\n");
				}
			}

			myWriter.close();
			System.out.println("Successfully wrote to the file.");
			
		} catch (IOException e) {
			System.out.println("An error occurred.");
			return false;
		}
		
		return true;
	}


	private boolean writeSpecsToFile(File file) {
		boolean result;
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
                for (double[] freeParameter : freeParameters) {
                    for (double v : freeParameter) {
                        myWriter.write(formatLongDouble.format(v) + " ");
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
