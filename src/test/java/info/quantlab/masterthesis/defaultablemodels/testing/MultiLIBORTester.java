package info.quantlab.masterthesis.defaultablemodels.testing;


import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleBinaryOperator;

import org.junit.Assert;
import org.junit.Test;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.functional.FunctionsOnIndependentIncrements;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel.Scheme;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class MultiLIBORTester {
	
	private static final int brownianMotionSeed = 3156;
	private static final int modelGeneratorSeed = 1453;
	private static final int numberOfPaths = 1000;
	private static final String stateSpace = LIBORMarketModelFromCovarianceModel.StateSpace.NORMAL.name();
	private static final String measure = LIBORMarketModelFromCovarianceModel.Measure.SPOT.name();
	private static final int numberOfLiborPeriods = 5;
	private static final double liborPeriodLength = 2.0;
	private static final double simulationTimeDelta = 0.01;
	private static final double[] fixingTimes = new double[] { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 };
	private static final double[] initialRatesNonDefaultable = new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 };
	private static final int numberOfFactorsNonDefaultable = 5;
	private static final double freeParameterRange = 0.5;
	

	public MultiLIBORTester() throws CalculationException {
		multiModel = getModel(2, measure, stateSpace);
		multiProcess = getProcess(brownianMotionSeed, numberOfPaths);
	}
	
	// We use static variables so that the process does not have to be calculated every time.
	private static MultiLIBORVectorModel multiModel;
	private static MonteCarloProcess multiProcess;
	
	
	@Test
	public void testNonDefaultableModel() {
		IndependentIncrements nonDefStochasticDriver = FunctionsOnIndependentIncrements.getRangeOfFactors(multiProcess.getStochasticDriver(), 0, numberOfFactorsNonDefaultable - 1);
		MonteCarloProcess nonDefProcess = new EulerSchemeFromProcessModel(multiModel.getNonDefaultableModel(), nonDefStochasticDriver, Scheme.EULER_FUNCTIONAL);
		
		
		Assert.assertTrue(true);
	}
	
	
	
	private MultiLIBORVectorModel getModel(int numberOfDefaultableModels, String measureOfModel, String stateSpaceOfModel) throws CalculationException {
		if(multiModel != null)
			return multiModel;
		
		LIBORMarketModelFromCovarianceModel nonDefaultableModel = getNonDefaultableModel(measureOfModel, stateSpaceOfModel);
		
		DefaultableLIBORMarketModel[] defaultableModels = new DefaultableLIBORMarketModel[numberOfDefaultableModels];
		MersenneTwister generator = new MersenneTwister(modelGeneratorSeed);
		// Set model properties
		final Map<String, String> properties = new HashMap<>();
		properties.put("measure", measureOfModel);
		properties.put("stateSpace", stateSpaceOfModel);
		properties.put("simulationModel", "LIBORS");
		DoubleBinaryOperator rand = (lower, upper) -> {return generator.nextDoubleFast() * (upper - lower) + lower; };
		for(int i=0; i < numberOfDefaultableModels; i++) {
			final int numberOfExtraFactors = (int)Math.floor(rand.applyAsDouble(1.0, (double)(nonDefaultableModel.getNumberOfFactors() + 1) - 0.3)); 
			// Only 70% probability for highest possible number of extra factors
			double[][] freeParameter = new double[numberOfLiborPeriods][numberOfExtraFactors];
			for (int row = 0; row < freeParameter.length; row++) {
				for (int col = 0; col < freeParameter[row].length; col++) {
					freeParameter[row][col] = rand.applyAsDouble(-freeParameterRange, freeParameterRange);
				}
			}
			
			DefaultableLIBORCovarianceWithGuaranteedPositiveSpread defCovarianceModel = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(nonDefaultableModel.getCovarianceModel(), freeParameter);

			final double[] initialRatesDefaultable = new double[] { 
					0.035 + rand.applyAsDouble(0.001, 0.02), 
					0.043 + rand.applyAsDouble(0.001, 0.02), 
					0.05 + rand.applyAsDouble(0.001, 0.02), 
					0.041 + rand.applyAsDouble(0.001, 0.02), 
					0.035 + rand.applyAsDouble(0.001, 0.02), 
					0.02 + rand.applyAsDouble(0.001, 0.02)
			};
			
			final ForwardCurve defaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
					"defaultableForwardCurve",
					fixingTimes,
					initialRatesDefaultable,
					liborPeriodLength);
			
			defaultableModels[i] = new DefaultableLIBORFromSpreadDynamic(nonDefaultableModel, defCovarianceModel, defaultableForwards, properties);
		}
		return new MultiLIBORVectorModel(defaultableModels, nonDefaultableModel);
	}

	private MonteCarloProcess getProcess(int bmSeed, int paths) {
		if(multiProcess != null)
			return multiProcess;
		
		TimeDiscretization times = multiModel.getNonDefaultableModel().getCovarianceModel().getTimeDiscretization();
		BrownianMotion myBM = new BrownianMotionFromMersenneRandomNumbers(times, multiModel.getNumberOfFactors(), paths, bmSeed);
		return new EulerSchemeFromProcessModel(multiModel, myBM, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);
	}
	
	private LIBORMarketModelFromCovarianceModel getNonDefaultableModel(String measureOfModel, String stateSpaceOfModel) throws CalculationException {
		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, numberOfLiborPeriods, liborPeriodLength); // Fixing time
		
		// Set initial forward curve
		
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
		final LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriods, numberOfFactorsNonDefaultable,	correlationDecayParam, true);
		
		final LIBORCovarianceModel baseCovarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriods, volatilityModel, correlationModel);
		
		// Set model properties
		final Map<String, String> properties = new HashMap<>();

		// Choose the simulation measure
		properties.put("measure", measureOfModel);

		// Choose normal model
		properties.put("stateSpace", stateSpaceOfModel);

		// Empty array of calibration items - hence, model will use given covariance
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		final LIBORMarketModelFromCovarianceModel baseLiborMarketModel = LIBORMarketModelFromCovarianceModel.of(
				liborPeriods, 		/* LIBORPeriodDiscretization */
				null,			 	/* analyticModel */
				nonDefaultableForwards, /* ForwardCurve */
				new DiscountCurveFromForwardCurve(nonDefaultableForwards),  /* DiscountCurve */
				new RandomVariableFromArrayFactory(true),  /* RV Factory */
				baseCovarianceModel,  /* covarianceModel */
				calibrationItems,  /* calibrationItems */
				properties);  /* properties */
		
		return baseLiborMarketModel;
	}
}
