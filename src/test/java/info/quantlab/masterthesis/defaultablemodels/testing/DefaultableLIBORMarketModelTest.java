package info.quantlab.masterthesis.defaultablemodels.testing;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;
import java.util.function.ToDoubleBiFunction;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;


import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModelFromCovarianceModel;
import info.quantlab.masterthesis.functional.FunctionsOnMCProcess;
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
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelCalibrateable;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Swaption;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionGeneralizedAnalyticApproximation.StateSpace;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.montecarlo.process.MonteCarloProcessFromProcessModel;
import net.finmath.plots.Named;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

@RunWith(Parameterized.class)
public class DefaultableLIBORMarketModelTest {
	

	@Parameters(name="{0}")
	public static Collection<Object[]> generateData()
	{
		return Arrays.asList(new Object[][] {
			// Put here an array of arrays where each array represents input for the constructor
			{new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 2, 0.3, 0.09, "NORMAL"},
			{new double[] { 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }, 2, 0.01, 0.0001, "LOGNORMAL"}
		});
	}

	// Non defaultable model Specifications
	private static final int numberOfFactors = 6;
	private static final boolean calibrateableModel = false;
	
	// General model specifications
	private static final int numberOfPaths = 1100;
	private static final int bmSeed = 9; /* seed */
	
	// RandomGenerator for Subsets
	private static final int subsetSeed  = 754693;
	
	private final LIBORModelMonteCarloSimulationModel defaultableLiborMarketModel;
	private final String stateSpace;
	
	private static DecimalFormat formatterMaturity	= new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	// Unnecessary: private static DecimalFormat formatterMoneyness	= new DecimalFormat(" 000.0%;-000.0%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterLong	= new DecimalFormat(" ##0,000;-##0,000", new DecimalFormatSymbols(Locale.ENGLISH));

	public DefaultableLIBORMarketModelTest(final double[] initialRatesDefaultable, int extraFactors, double correlationDecayParam, double exponentialFactor, String stateSpace) throws CalculationException {
		final double[] initialRatesUndefaultableModel = { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 };
		final String measure = "SPOT";
		final double simulationTimeDelta = 0.001;
		defaultableLiborMarketModel = getValuationModel(initialRatesDefaultable, initialRatesUndefaultableModel, simulationTimeDelta, stateSpace, measure, extraFactors, correlationDecayParam, exponentialFactor);
		this.stateSpace = stateSpace;
	}
	
	@Test
	public void testGeneralSpecs() throws CalculationException {
		// Print All Specifications:
		TimeDiscretization timeDiscretization = defaultableLiborMarketModel.getTimeDiscretization();
		DefaultableLIBORMarketModelFromCovarianceModel defaultableModel = (DefaultableLIBORMarketModelFromCovarianceModel)defaultableLiborMarketModel.getModel();
		MonteCarloProcess process = defaultableLiborMarketModel.getProcess();
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread covModel = (DefaultableLIBORCovarianceWithGuaranteedPositiveSpread)(defaultableModel.getCovarianceModel());
		double[][] freeParameters = covModel.getFreeParameterMatrix();
		
		// Number of pointers to RVs
		long heapRequirement = (long)timeDiscretization.getNumberOfTimes() * (long)defaultableModel.getNumberOfComponents() * 8;
		// Number of different RandomVariables:
		long numberOfRandomVariables = (long)defaultableModel.getNumberOfComponents();
		for(int component = 0; component < defaultableModel.getNumberOfLibors(); component++) {
			long timeIndexOfFixing = (long)timeDiscretization.getTimeIndex(defaultableModel.getLiborPeriod(component));
			if(timeIndexOfFixing < 0)
				timeIndexOfFixing = - timeIndexOfFixing - 1;
			numberOfRandomVariables += timeIndexOfFixing * 2; // Same for non-defaultable and defaultable model
		}
		// Add number of Pointers for realizations in RV (might be 0, still take space)
		heapRequirement += numberOfRandomVariables * 8;
		// Add number of doubles for nonStochasticValue in RV (might be NaN, still take space)
		heapRequirement += numberOfRandomVariables * 8;
		// Number of doubles in realizations of the RVs:
		heapRequirement += (numberOfRandomVariables - defaultableModel.getNumberOfComponents()) * numberOfPaths * 8;
		
		System.out.println("\nEulerScheme Heap requirement for model: " + formatterLong.format(heapRequirement/1E6)+ " Megabytes");
		System.out.println("\nGeneral Model Specifications:\n");
		System.out.println("Measure:    " + defaultableModel.getMeasure());
		System.out.println("Statespace: " + defaultableModel.getStateSpace());
		
		System.out.println("\nStarting Rates:");
		System.out.printf("LIBOR Period %21s Non-Defaultable %5s Defaultable %9s Spread %n", " ", " ", " ");
		tic();
		for(int i=0; i < defaultableModel.getNumberOfLibors(); i++) {
			final double nonDefRate = defaultableModel.getUndefaultableLIBOR(process, 0, i).doubleValue();
			final double defRate = defaultableModel.getLIBOR(process, 0, i).doubleValue();
			final double spread = defaultableModel.getLIBORSpreadAtGivenTimeIndex(process, 0, i).doubleValue();
			System.out.printf("%2d (fixing time: %4.1f): %10s %10.8f %10s %10.8f %10s %10.8f %n", 
					i, defaultableModel.getLiborPeriod(i), " ", nonDefRate, " ", defRate, " ", spread);
		}
		long result = toc(false);
		System.out.println("\nElapsed Time for calculation of Processes is " + (result / 1E9f) + "Seconds.");
		System.out.println();
		System.out.println("Matrix of FL factors at time 0:");
		
		for(int row = 0; row < freeParameters.length; row++) {
			final RandomVariable[] factorLoadings = covModel.getNonDefaultableCovarianceModel().getFactorLoading(0.0, row, null);
			for(int col = 0; col < factorLoadings.length; col++) {
				System.out.printf("%12.8f      ", factorLoadings[col].doubleValue());
			}
			System.out.printf("| %5s ", " ");
			for(int col = 0; col < freeParameters[row].length; col++) {
				System.out.printf("%12.8f      ", freeParameters[row][col]);
			}
			System.out.println();
		}
		
		System.out.println("\nCovariance Matrix of non-Defaultable Model at time 0:");
		RandomVariable[] liborsAtZero = defaultableModel.getInitialValue(process);
		for(int row = 0; row < defaultableModel.getNumberOfLIBORPeriods(); row++) {
			for(int col = 0; col < defaultableModel.getNumberOfLIBORPeriods(); col++) {
				final RandomVariable covariance = covModel.getNonDefaultableCovarianceModel().getCovariance(0.0, row, col, null);
				System.out.printf("%10.6f      ", covariance.doubleValue());
			}
			System.out.printf("| %5s ", " ");
			for(int col = 0; col < defaultableModel.getNumberOfLIBORPeriods(); col++) {
				final RandomVariable covariance = covModel.getCovariance(0.0, row + defaultableModel.getNumberOfLIBORPeriods(), col + defaultableModel.getNumberOfLIBORPeriods(), liborsAtZero);
				System.out.printf("%10.6f      ", covariance.doubleValue());
			}
			System.out.println();
		}
		System.out.println("\n" + "_".repeat(300) + "\n");
		Assert.assertTrue(1.0 > 0.0);
	}
	
	
	@Test
 	public void testPositivity() throws CalculationException {
		/*
		 * Checking subsets of Spread for positivity
		 */
		System.out.println("Testing Positivity of Spread:\n");
		System.out.println(" ".repeat(20) + "Minimum at LIBOR");
		System.out.print("Time");
		for(int component = defaultableLiborMarketModel.getNumberOfLibors() - 1; component >= 0; component--) {
			System.out.printf(" ".repeat(18) + "%2d", component);
		}
		System.out.println();
		MersenneTwister randomGenerator = new MersenneTwister(subsetSeed);
		final int numberOfPointsToCheck = 20;
		final DefaultableLIBORMarketModel defModel = (DefaultableLIBORMarketModel)(defaultableLiborMarketModel.getModel());
		double overAllMin = 1.0;
		
		
		for(int i = 0; i < numberOfPointsToCheck; i++) {
			final double time = randomGenerator.nextDoubleFast() * defaultableLiborMarketModel.getTimeDiscretization().getLastTime();
			System.out.print(formatterMaturity.format(time) + " ".repeat(19));
			for(int component = defaultableLiborMarketModel.getNumberOfLibors() - 1; component >= 0; component--) {
				if(time > defaultableLiborMarketModel.getLiborPeriod(component)) {
					break;
				}
				final double minimumAtComponent = defModel.getLIBORSpreadAtGivenTime(defaultableLiborMarketModel.getProcess(), time, component).getMin();
				System.out.print(formatterDeviation.format(minimumAtComponent) + " ".repeat(8));
				overAllMin = Math.min(overAllMin, minimumAtComponent);
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Overall Minimum was:\t" + formatterDeviation.format(overAllMin) + "\n");
		
		// Plot Average of each component (also undefaultable ones):
		final int numberOfComponents = defModel.getNumberOfComponents();
		
		final int numberOfLIBORs = defModel.getNumberOfLibors();
		DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[numberOfComponents];
		int j = 0;
		for(; j < numberOfLIBORs; j++) {
			final int compIndex = j;
			myFunctionArray[j] = new DoubleUnaryOperator() {
				@Override
				public double applyAsDouble(double operand) {
					try {
						int timeIndex = defaultableLiborMarketModel.getTimeIndex(operand);
						if(timeIndex < 0)
							timeIndex = - timeIndex - 2;
						return defModel.getUndefaultableLIBOR(defaultableLiborMarketModel.getProcess(), timeIndex, compIndex).getAverage();
					} catch (CalculationException e) {
						return - 1.0;
					}
				}
			};
		}
		for(; j < numberOfComponents; j++) {
			final int compIndex = j;
			myFunctionArray[j] = new DoubleUnaryOperator() {
				@Override
				public double applyAsDouble(double operand) {
					try {
						int timeIndex = defaultableLiborMarketModel.getTimeIndex(operand);
						if(timeIndex < 0)
							timeIndex = - timeIndex - 2;
						return defaultableLiborMarketModel.getLIBOR(timeIndex, compIndex - defModel.getNumberOfLibors()).getAverage();
					} catch (CalculationException e) {
						return - 1.0;
					}
				}
			};
		}
		List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
		List<Named<DoubleUnaryOperator>> myFunctionList =
				myList.stream().map(operator -> new Named<DoubleUnaryOperator>("LIBOR  " + myList.indexOf(operator), operator)).
				collect(Collectors.toList());
		EasyPlot2D plot = new EasyPlot2D(defaultableLiborMarketModel.getTime(0), defaultableLiborMarketModel.getTimeDiscretization().getLastTime(), 51, myFunctionList);
		plot.show();
		
		System.out.println("_".repeat(300) + "\n");
		Assert.assertTrue(overAllMin > -1E-3); // smaller than - 0.001
	}
	
	@Test
 	public void testCaplet() throws CalculationException {
		/*
		 * Value a caplet
		 */
		System.out.println("Caplet prices:\n");
		System.out.print("Maturity      Simulation      Sim Non-Def     Analytic");
		
		System.out.println("        Deviation");
		double maxAbsDeviation = 0.0;
		for (double maturity = 2.0; maturity <= 19.0; maturity += 2.0) {

			final double exerciseDate = maturity;
			System.out.print(formatterMaturity.format(exerciseDate) + "         ");

			final int numberOfPeriods = 1;

			// Create a swaption

			final double[] fixingDates = new double[numberOfPeriods];
			final double[] paymentDates = new double[numberOfPeriods];
			final double[] swapTenor = new double[numberOfPeriods + 1];
			final double swapPeriodLength = 2.0;

			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
				paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex + 1) * swapPeriodLength;
				swapTenor[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
			}
			swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;

			// Swaptions swap rate
			final double swaprate = getParSwaprate(defaultableLiborMarketModel, swapTenor);

			// Set swap rates for each period
			final double[] swaprates = new double[numberOfPeriods];
			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				swaprates[periodStartIndex] = swaprate;
			}

			// Value with Monte Carlo
			final Swaption swaptionMonteCarlo	= new Swaption(exerciseDate, fixingDates, paymentDates, swaprates);
			final double valueSimulation = swaptionMonteCarlo.getValue(defaultableLiborMarketModel);
			System.out.print(formatterValue.format(valueSimulation) + "        ");
			final double valueNonDefSim = swaptionMonteCarlo.getValue(getNonDefaultableValuationModel(defaultableLiborMarketModel));
			System.out.print(formatterValue.format(valueNonDefSim) + "        ");
			
			// Value analytic.
			AbstractTermStructureMonteCarloProduct swaptionAnalytic = null;
			final TimeDiscretization swapTenorTD = new TimeDiscretizationFromArray(swapTenor);
			if(stateSpace.toUpperCase() == "LOGNORMAL") {
				swaptionAnalytic = new SwaptionGeneralizedAnalyticApproximation(swaprate, swapTenorTD, StateSpace.LOGNORMAL);
			}
			else {
				swaptionAnalytic = new SwaptionGeneralizedAnalyticApproximation(swaprate, swapTenorTD, StateSpace.NORMAL);
			}
			final double valueAnalytic = swaptionAnalytic.getValue(defaultableLiborMarketModel);
			System.out.print(formatterValue.format(valueAnalytic) + "          ");

			// Absolute deviation
			final double deviation = (valueSimulation - valueAnalytic);
			System.out.println(formatterDeviation.format(deviation) + "          ");

			maxAbsDeviation = Math.max(maxAbsDeviation, Math.abs(deviation));
		}

		System.out.println("Maximum abs deviation: " + formatterDeviation.format(maxAbsDeviation));
		System.out.println("_".repeat(300) + "\n");
		
		/*
		 * jUnit assertion: condition under which we consider this test successful
		 */
		Assert.assertTrue(Math.abs(maxAbsDeviation) < 8E-3);
	}
	
	private static LIBORModelMonteCarloSimulationModel getValuationModel(
			double[] initialRatesDefaultable, 
			double[] initialRatesUndefaultable, 
			double simulationTimeDelta, 
			String stateSpace, 
			String measure, 
			int numberOfExtraFactors,
			double correlationDecayParam,
			double exponentialFactorFreeParameters
			) throws CalculationException {
		
		// Set LIBOR times
		final double liborPeriodLength = 2.0;
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, 10, liborPeriodLength); // Fixing time
		
		// Set initial forward curves
		final double[] fixingTimes = new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0};
		if(initialRatesDefaultable.length != fixingTimes.length || initialRatesUndefaultable.length != fixingTimes.length)
			throw new IllegalArgumentException("Initial rates must have length " + fixingTimes.length);
		
		final ForwardCurve nonDefaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"nonDefaultableForwardCurve",
				fixingTimes,
				initialRatesUndefaultable,
				liborPeriodLength);
		
		/*
		 * Create a simulation time discretization. We save space by just simulationg to the last fixing time.
		 */
		final int numberOfTimes = (int)(liborPeriods.getTime(liborPeriods.getNumberOfTimeSteps() - 1)/simulationTimeDelta);
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, numberOfTimes, simulationTimeDelta);

		/*
		 * Create volatility structure v[i][j] = sigma_j(t_i) and correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		final double a = 0.2, b = 0.0, c = 0.25, d = 0.3;
		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriods, a, b, c, d, calibrateableModel);
		final LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriods, numberOfFactors,	correlationDecayParam, calibrateableModel);
		
		final LIBORCovarianceModelCalibrateable baseCovarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriods, volatilityModel, correlationModel);
		
		// Set model properties
		final Map<String, String> properties = new HashMap<>();

		// Choose the simulation measure
		properties.put("measure", measure);

		// Choose log normal model
		properties.put("stateSpace", stateSpace);

		// Empty array of calibration items - hence, model will use given covariance
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		final LIBORMarketModel baseLiborMarketModel = LIBORMarketModelFromCovarianceModel.of(
				liborPeriods, 		/* LIBORPeriodDiscretization */
				null,			 	/* analyticModel */
				nonDefaultableForwards, /* LIBORPeriodDiscretization */
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

		final double[][] freeParameters = new double[baseLiborMarketModel.getNumberOfComponents()][numberOfExtraFactors];
		
		ToDoubleBiFunction<Integer, Integer> freeParameterFunc = 
				(row, col) -> 1.0 - Math.exp(- Math.abs(exponentialFactorFreeParameters) * Math.sqrt(liborPeriods.getTime(row)));
				
				
		for(int i=0; i < baseLiborMarketModel.getNumberOfComponents(); i++) {
			for(int j=0; j < numberOfExtraFactors; j++) {
				freeParameters[i][j] = freeParameterFunc.applyAsDouble(i, j);
			}
		}
		
		final DefaultableLIBORCovarianceModel defaultableCovariance = 
				new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(baseCovarianceModel, freeParameters);
		
		final DefaultableLIBORMarketModel defaultableModel = new DefaultableLIBORMarketModelFromCovarianceModel(baseLiborMarketModel, defaultableCovariance, defaultableForwards, properties);
		
		
		
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers(timeDiscretization, defaultableModel.getNumberOfFactors(), numberOfPaths, bmSeed);
		EulerSchemeFromProcessModel.Scheme eulerScheme = EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR;
		
		if(stateSpace.toUpperCase() == "LOGNORMAL")
			eulerScheme = EulerSchemeFromProcessModel.Scheme.PREDICTOR_CORRECTOR_FUNCTIONAL;
		
		final MonteCarloProcessFromProcessModel process = new EulerSchemeFromProcessModel(defaultableModel, brownianMotion, eulerScheme);

		return new LIBORMonteCarloSimulationFromLIBORModel(process);
	}
	
	private static LIBORModelMonteCarloSimulationModel getNonDefaultableValuationModel(LIBORModelMonteCarloSimulationModel defaultableModel) {
		
		final MonteCarloProcess normalProcess = defaultableModel.getProcess();
		final DefaultableLIBORMarketModel defaultableTheoModel = (DefaultableLIBORMarketModel)(defaultableModel.getModel());
		
		LIBORModelMonteCarloSimulationModel newModel = new LIBORModelMonteCarloSimulationModel() {

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
				return defaultableTheoModel.getDefaultableLIBOR(normalProcess, timeIndex, liborIndex);
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
		
		return newModel;
	}
	
	private static double getParSwaprate(final LIBORModelMonteCarloSimulationModel liborMarketModel, final double[] swapTenor) {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretizationFromArray(swapTenor), new TimeDiscretizationFromArray(swapTenor), liborMarketModel.getModel().getForwardRateCurve(), liborMarketModel.getModel().getDiscountCurve());
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Timing
	
	
	private static long startTime = 0;
	
	public static void tic() {
		startTime = System.nanoTime();
	}
	
	public static long toc(boolean printResult) {
		long result = System.nanoTime();
		result -= startTime;
		if(printResult)
			printTimeResult(result);
		return result;
	}
	
	public static long toc() {
		return toc(true);
	}
	
	public static void printTimeResult(long result) {
		System.out.println("\n\nElapsed Time is " + (result / 1000) + " Microseconds\n");
	}
}
