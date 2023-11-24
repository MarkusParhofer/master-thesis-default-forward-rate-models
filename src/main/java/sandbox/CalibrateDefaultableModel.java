package sandbox;

import java.util.HashMap;
import java.util.Map;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModelFromCovarianceModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.model.volatilities.SwaptionMarketData;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelExponentialForm5Param;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Swaption;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class CalibrateDefaultableModel {

	public static void main(String[] args) throws CalculationException {
		// Set general Settings:
		
		// LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, 10, 0.5);
		final double[] fixingTimes = liborPeriods.getAsDoubleArray();
		// Simulation times
		TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, 20, 0.02);
		
		// Model properties of as well undefaultable as defaultable model
		Map<String, String> properties = new HashMap<String, String>();
		properties.putIfAbsent("measure", "Spot");
		properties.putIfAbsent("stateSpace", "Normal");
		
		
		
		
		
		// ------------------------------------------------------------------------------------------------------------------------------ //
		
		// Create undefaultable model
		
		// Set initial LIBOR rates
		final double[] uForwardRates = new double[] { 0.61 / 100.0, 0.61 / 100.0, 0.67 / 100.0, 0.73 / 100.0, 0.80 / 100.0, 0.92 / 100.0, 
				1.11 / 100.0, 1.36 / 100.0, 1.60 / 100.0, 1.82 / 100.0, 2.02 / 100.0 };
		ForwardCurve uForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("forwardCurve",	fixingTimes, uForwardRates, 0.5);
		
		// Set Covariance model
		LIBORCovarianceModel uCovarianceModel = new LIBORCovarianceModelExponentialForm5Param(timeDiscretization, liborPeriods, 9);
		
		SwaptionMarketData uSwaptionData = null;
		
		// Set undefaultable LIBOR model
		LIBORMarketModel uModel = new LIBORMarketModelFromCovarianceModel(liborPeriods, uForwardCurve, new DiscountCurveFromForwardCurve(uForwardCurve), uCovarianceModel, uSwaptionData, properties);
		
		
		
		
		
		// ------------------------------------------------------------------------------------------------------------------------------ //
		
		// Create original defaultable model
				
		// Set initial LIBOR rates
		final double[] myForwardRates = new double[] { 0.71 / 100.0, 0.71 / 100.0, 0.74 / 100.0, 0.79 / 100.0, 0.90 / 100.0, 0.99 / 100.0, 
				1.20 / 100.0, 1.44 / 100.0, 1.71 / 100.0, 1.94 / 100.0, 2.20 / 100.0 };
		ForwardCurve myForwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards("defaultableForwardCurve", fixingTimes, myForwardRates, 0.5);
		
		// Set Covariance model
		double[][] originalFreeParameters = new double[uModel.getNumberOfComponents()][3];
		// Free Parameters as random doubles between 0 and 1
		System.out.println("\nOriginal Free Parameter Matrix:\n");
		MersenneTwister randomGenerator = new MersenneTwister(1998);
		for(int componentIndex = 0; componentIndex < uModel.getNumberOfComponents(); componentIndex++) {
			for(int factor = 0; factor < 3; factor++) {
				originalFreeParameters[componentIndex][factor] = randomGenerator.nextDouble();
				System.out.printf("%8.6f  ", originalFreeParameters[componentIndex][factor]);
			}
			System.out.println();
		}
		
		DefaultableLIBORCovarianceModel originalCovariance = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(uCovarianceModel, originalFreeParameters);
		DefaultableLIBORMarketModelFromCovarianceModel originalModel = new DefaultableLIBORMarketModelFromCovarianceModel(uModel, originalCovariance, myForwardCurve, null, properties);
		
		
		
		
		
		// ------------------------------------------------------------------------------------------------------------------------------ //
		
		// TODO: Value some options
		AbstractTermStructureMonteCarloProduct[] products = new AbstractTermStructureMonteCarloProduct[10];
		
		// AbstractTermStructureMonteCarloProduct product = new Swaption(exerciseDate, swapTenor, swaprate);
		
		
		// TODO: Create Calibration Products
		// CalibrationProduct(final AbstractTermStructureMonteCarloProduct product, final double targetValue, final double weight);
		// TODO: Create false model
		// TODO: Calibrate false model
	}

}
