package info.quantlab.masterthesis.defaultablemodels.testing;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import org.junit.Assert;
import org.junit.jupiter.api.Test;

import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORCovarianceWithInitialUndefaultableCovariance;
import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class CovarianceTest {
	
	public static final double[] fixingTimes = new double[] { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 };
	public static double[] initialRatesDefaultable = new double[] { 0.3, 0.23, 0.19, 0.24, 0.33, 0.23 }; //{ 0.1, 0.13, 0.11, 0.14, 0.18, 0.2 }; //{ 0.04, 0.049, 0.062, 0.049, 0.044, 0.031 }; 		// needs to be of length fixingTimes
	public static double[] initialRatesNondefaultable = new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 };		// needs to be of length fixingTimes
	public static final double liborPeriodLength = 2.0;
	public static final double simulationTimeDelta = 0.1;
	
	public static final ForwardCurve nonDefaultableForwards;
	public static final ForwardCurve defaultableForwards;
	public static final LIBORCovarianceModel baseCovarianceModel;
	
	public static final int numberOfExtraFactors = 5;
	public final double squeezeParam = 1.0;
	
	
	private static DecimalFormat formatterValue	= new DecimalFormat(" 0.0000;-0.0000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	
	
	public CovarianceTest() {
		// TODO Auto-generated constructor stub
	}

	@Test
	public void testCovarianceStructureFromForwardCurve() {
		DefaultableLIBORCovarianceWithGuaranteedPositiveSpread defCovModel = 
				new DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(baseCovarianceModel, defaultableForwards, nonDefaultableForwards, numberOfExtraFactors, squeezeParam);
		
		System.out.println("\nNumber of Factors: " + defCovModel.getNumberOfFactors() + "\n");
		
		System.out.println("\nFree Parameter Matrix:\n");
		double[][] freeParams = defCovModel.getFreeParameterMatrix();
		
		for (int row = 0; row < freeParams.length; row++) {
			for (int col = 0; col < freeParams[row].length; col++) {
				System.out.printf("%9.4f   ", freeParams[row][col]);
			}
			System.out.println();
		}
		System.out.println("\n\nCovariance Matrix:\n");
		TimeDiscretization liborPeriods = baseCovarianceModel.getLiborPeriodDiscretization();
		RandomVariable[] initialValues = new RandomVariable[2 * liborPeriods.getNumberOfTimeSteps()];
		for (int component = 0; component < initialValues.length; component++) {
			final double initValue = component < liborPeriods.getNumberOfTimeSteps() ? 
					nonDefaultableForwards.getForward(null, liborPeriods.getTime(component)) :
						defaultableForwards.getForward(null, liborPeriods.getTime(component - liborPeriods.getNumberOfTimeSteps()));
			initialValues[component] = Scalar.of(initValue);
		}
		
		double maxAbsDeviation = 0.0;
		for(int row = 0; row < liborPeriods.getNumberOfTimeSteps(); row++) {
			for (int col = 0; col < liborPeriods.getNumberOfTimeSteps(); col++) {
				final double initCovNonDef = baseCovarianceModel.getCovariance(0, row, col, null).doubleValue();
				final double initCovDef = defCovModel.getCovariance(0, row+ liborPeriods.getNumberOfTimeSteps(), col+ liborPeriods.getNumberOfTimeSteps(), initialValues).doubleValue();
				maxAbsDeviation = Math.max(maxAbsDeviation, Math.abs(initCovNonDef - initCovDef));
			}
			
			for (int col = 0; col < liborPeriods.getNumberOfTimeSteps(); col++) {
				final double initCovNonDef = baseCovarianceModel.getCovariance(0, row, col, null).doubleValue();
				System.out.print(formatterValue.format(initCovNonDef) + "     ");
			}
			System.out.print("|     ");
			for (int col = 0; col < liborPeriods.getNumberOfTimeSteps(); col++) {
				final double initCovDef = defCovModel.getCovariance(0, row + liborPeriods.getNumberOfTimeSteps(), col + liborPeriods.getNumberOfTimeSteps(), initialValues).doubleValue();
				System.out.print(formatterValue.format(initCovDef) + "     ");
			}
			System.out.println();
		}
		System.out.println("\nMaximum Absolute Deviation is: " + formatterDeviation.format(maxAbsDeviation));
		System.out.println();
		double maxAbsDeviationLowFac = 0.0;
		
		for(int row = 0; row < liborPeriods.getNumberOfTimeSteps(); row++) {
			
			for (int col = 0; col < liborPeriods.getNumberOfTimeSteps(); col++) {
				double initCovDef = 0.0;
				RandomVariable[] flsRow = defCovModel.getFactorLoading(0, row + defCovModel.getLiborPeriodDiscretization().getNumberOfTimeSteps(), initialValues);
				RandomVariable[] flsCol = defCovModel.getFactorLoading(0, col + defCovModel.getLiborPeriodDiscretization().getNumberOfTimeSteps(), initialValues);
				for(int factor = 0; factor < baseCovarianceModel.getNumberOfFactors(); factor++) {
					initCovDef += flsRow[factor].doubleValue() * flsCol[factor].doubleValue();
				}
				final double initCovNonDef = baseCovarianceModel.getCovariance(0, row, col, null).doubleValue();
				maxAbsDeviationLowFac = Math.max(maxAbsDeviationLowFac, Math.abs(initCovNonDef - initCovDef));
				System.out.print(formatterValue.format(initCovDef) + "     ");
			}
			System.out.println();
		}
		System.out.println("\nMaximum Absolute Deviation is: " + formatterDeviation.format(maxAbsDeviationLowFac));
		
		
		Assert.assertTrue(maxAbsDeviation < 1E-2d);
	}
	
	static {
		
		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, 5, liborPeriodLength); // Fixing time
		
		// Set initial forward curves
		if(initialRatesDefaultable.length != fixingTimes.length || initialRatesNondefaultable.length != fixingTimes.length)
			throw new IllegalArgumentException("Initial rates must have length " + fixingTimes.length);
		
		nonDefaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"nonDefaultableForwardCurve",
				fixingTimes,
				initialRatesNondefaultable,
				liborPeriodLength);
		
		defaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"defaultableForwardCurve",
				fixingTimes,
				initialRatesDefaultable,
				liborPeriodLength);
		
		/*
		 * Create a simulation time discretization. We save space by just simulationg to the last fixing time.
		 */
		final int numberOfTimes = (int)(liborPeriods.getTime(liborPeriods.getNumberOfTimeSteps() - 1)/simulationTimeDelta);
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, numberOfTimes, simulationTimeDelta);

		/*
		 * Create volatility structure v[i][j] = sigma_j(t_i) and correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		final int numberOfFactors = 5;
		final double a = 0.2, b = 0.0, c = 0.25, d = 0.3;
		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriods, a, b, c, d, false);
		final double correlationDecayParam = 0.1;
		final LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriods, numberOfFactors,	correlationDecayParam, false);
		
		baseCovarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriods, volatilityModel, correlationModel);
		
	}
}
