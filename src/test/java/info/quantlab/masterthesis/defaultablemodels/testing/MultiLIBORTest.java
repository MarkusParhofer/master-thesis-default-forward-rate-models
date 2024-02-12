package info.quantlab.masterthesis.defaultablemodels.testing;


import java.util.HashMap;
import java.util.Map;

import info.quantlab.debug.Time;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import org.junit.Assert;
import org.junit.Test;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.MonteCarloProcess;

public class MultiLIBORTest extends Time {

	final static int numberOfLibors = 15;
	final static double liborPeriodLength = 0.5;



	public static class ReturnType {
		public ReturnType(MonteCarloProcess process, LIBORMarketModel model) {
			p = process;
			m = model;
		}

		public final MonteCarloProcess p;
		public final LIBORMarketModel m;
	}

	@Test
	public void debtorTestTest() throws CalculationException {
		ReturnType mm = createMultiModel();
		ReturnType dm = debtorModel();
		System.out.println("_".repeat(300));
		System.out.println();
		System.out.println("Testing debtor model: ");
		double overallMaxDev = 0.0;
		for(int i=0; i < dm.m.getNumberOfLibors(); i++) {
			final int timeStepsToCheck = 5;
			final int stepsToJump = dm.p.getTimeDiscretization().getNumberOfTimeSteps() / timeStepsToCheck;
			double maxRelDev = 0.0;

			System.out.printf("Libor %2d (T=%4.2f):     ", i, dm.m.getLiborPeriod(i));
			for (int j = 0; j < timeStepsToCheck; j++) {
				double mmAv = ((MultiLIBORVectorModel)mm.m).getDefaultableLIBOR(mm.p, j * stepsToJump, i, 1).getAverage();
				double dmAv = ((DefaultableLIBORMarketModel)dm.m).getDefaultableLIBOR(dm.p, j * stepsToJump, i).getAverage();
				double e = Math.abs(mmAv - dmAv);
				maxRelDev = Math.max(e / dmAv, maxRelDev);
				System.out.printf("%8.6f     %8.6f     %10.8f  |     ", mmAv, dmAv, e);
			}
			System.out.printf("%6.2f %%", maxRelDev * 100.0);
			System.out.println();
			overallMaxDev = Math.max(maxRelDev, overallMaxDev);
		}
		System.out.printf("\nOverall maximum relative deviation was %6.2f %%\n\n", overallMaxDev * 100);
		System.out.println("_".repeat(300));

		Assert.assertTrue(overallMaxDev < 0.05); // 5 % deviation at most
	}

	@Test
	public void creditorTest() throws CalculationException {
		ReturnType mm = createMultiModel();
		ReturnType dm = creditorModel();
		System.out.println("_".repeat(300));
		System.out.println();
		System.out.println("Testing creditor model: ");
		double overallMaxDev = 0.0;

		for(int i=0; i < dm.m.getNumberOfLibors(); i++) {
			final int timeStepsToCheck = 5;
			final int stepsToJump = dm.p.getTimeDiscretization().getNumberOfTimeSteps() / timeStepsToCheck;
			double maxRelDev = 0.0;
			System.out.printf("Libor %2d (T=%4.2f):     ", i, dm.m.getLiborPeriod(i));
			for (int j = 0; j < timeStepsToCheck; j++) {
				double mmAv = ((MultiLIBORVectorModel)mm.m).getDefaultableLIBOR(mm.p, j * stepsToJump, i, 0).getAverage();
				double dmAv = ((DefaultableLIBORMarketModel)dm.m).getDefaultableLIBOR(dm.p, j * stepsToJump, i).getAverage();
				double e = Math.abs(mmAv - dmAv);
				maxRelDev = Math.max(e / dmAv, maxRelDev);
				System.out.printf("%8.6f     %8.6f     %10.8f  |     ", mmAv, dmAv, e);
			}
			System.out.printf("%6.2f %%", maxRelDev * 100.0);
			System.out.println();
			overallMaxDev = Math.max(maxRelDev, overallMaxDev);
		}

		System.out.printf("\nOverall maximum relative deviation was %6.2f %%\n\n", overallMaxDev * 100);
		System.out.println("_".repeat(300));
		Assert.assertTrue(overallMaxDev < 0.05); // 5 % deviation at most
	}

	public static ReturnType creditorModel() throws CalculationException {
		tic();
		DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
		Map<String, Object> properties = new HashMap<>();

		properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
		properties.put("numberOfLiborPeriods", numberOfLibors);
		properties.put("liborPeriodLength", liborPeriodLength);
		properties.put("stateSpace", "NORMAL");
		properties.put("measure", "SPOT");

		properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
		properties.put("numberOfFactors", 5);
		properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
		properties.put("correlationDecayParam", 0.2);
		properties.put("displacement", 0.0001);
		properties.put("initialRatesNonDefaultable", new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032});
		properties.put("simulationTimeDelta", 0.01);

		factory.setProperties(properties);

		LIBORMarketModel base = factory.createBaseModel();

		// Def Model 1: Creditor (lower Covariance and lower Spread):
		properties.put("numberOfExtraFactors", 2);
		properties.put("initialRatesDefaultable", new double[] {0.028, 0.03, 0.031, 0.033, 0.034, 0.036});
		properties.put("simulationModel", "SPREADS");
		properties.put("stateSpaceOfSpread", "LOGNORMAL");
		properties.put("freeParamsSeed", 2000);
		properties.put("freeParamsRange", 0.3);
		factory.setProperties(properties);

		DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);

		// Simulation:
		properties.put("numberOfPaths", 5000);
		properties.put("brownianMotionSeed", 3666);
		properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
		factory.setProperties(properties);

		MonteCarloProcess process = factory.createNumericalScheme(modelCreditor);

		process.getProcessValue(10, 3); // Such that calculation takes place here
		printTimeResult("Creditor Model", toc(false));
		return new ReturnType(process, modelCreditor);
	}

	public static ReturnType debtorModel() throws CalculationException {
		tic();
		DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
		Map<String, Object> properties = new HashMap<>();

		properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
		properties.put("numberOfLiborPeriods", numberOfLibors);
		properties.put("liborPeriodLength", liborPeriodLength);
		properties.put("stateSpace", "NORMAL");
		properties.put("measure", "SPOT");

		properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
		properties.put("numberOfFactors", 5);
		properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
		properties.put("correlationDecayParam", 0.2);
		properties.put("displacement", 0.0001);
		properties.put("initialRatesNonDefaultable", new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032});
		properties.put("simulationTimeDelta", 0.01);

		factory.setProperties(properties);

		LIBORMarketModel base = factory.createBaseModel();

		// Def Model 1: Debtor (lower Covariance and lower Spread):
		properties.put("numberOfExtraFactors", 3);
		properties.put("initialRatesDefaultable", new double[] {0.032, 0.036, 0.039, 0.04, 0.04, 0.043});
		properties.put("simulationModel", "SPREADS");
		properties.put("stateSpaceOfSpread", "LOGNORMAL");
		properties.put("freeParamsSeed", 1030);
		properties.put("freeParamsRange", 0.6);
		factory.setProperties(properties);

		DefaultableLIBORMarketModel modelDebtor = factory.createDefaultableModel(base);

		// Simulation:
		properties.put("numberOfPaths", 5000);
		properties.put("brownianMotionSeed", 3436);
		properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
		factory.setProperties(properties);

		MonteCarloProcess process = factory.createNumericalScheme(modelDebtor);

		process.getProcessValue(10, 3); // Such that calculation takes place here
		printTimeResult("Debtor Model", toc(false));

		return new ReturnType(process, modelDebtor);
	}

	public static ReturnType createMultiModel() throws CalculationException {
		tic();
		DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
		Map<String, Object> properties = new HashMap<>();

		properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
		properties.put("liborPeriodLength", liborPeriodLength);
		properties.put("numberOfLiborPeriods", numberOfLibors);
		properties.put("stateSpace", "NORMAL");
		properties.put("measure", "SPOT");

		properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
		properties.put("numberOfFactors", 5);
		properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
		properties.put("correlationDecayParam", 0.2);
		properties.put("displacement", 0.0001);
		properties.put("initialRatesNonDefaultable", new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032});
		properties.put("simulationTimeDelta", 0.01);

		factory.setProperties(properties);

		LIBORMarketModel base = factory.createBaseModel();

		// Def Model 1: Creditor (lower Covariance and lower Spread):
		properties.put("numberOfExtraFactors", 2);
		properties.put("initialRatesDefaultable", new double[] {0.028, 0.03, 0.031, 0.033, 0.034, 0.036});
		properties.put("simulationModel", "SPREADS");
		properties.put("stateSpaceOfSpread", "LOGNORMAL");
		properties.put("freeParamsSeed", 2000);
		properties.put("freeParamsRange", 0.3);
		factory.setProperties(properties);

		DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);


		// Def Model 2: Debtor (lower Covariance and lower Spread):
		properties.put("numberOfExtraFactors", 3);
		properties.put("initialRatesDefaultable", new double[] {0.032, 0.036, 0.039, 0.04, 0.04, 0.043});
		properties.put("simulationModel", "SPREADS");
		properties.put("stateSpaceOfSpread", "LOGNORMAL");
		properties.put("freeParamsSeed", 1030);
		properties.put("freeParamsRange", 0.6);
		factory.setProperties(properties);

		DefaultableLIBORMarketModel modelDebtor = factory.createDefaultableModel(base);


		MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditor, modelDebtor}, base);

		// Simulation:
		properties.put("numberOfPaths", 5000);
		properties.put("brownianMotionSeed", 7878);
		properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
		factory.setProperties(properties);

		MonteCarloProcess process = factory.createNumericalScheme(multiModel);

		process.getProcessValue(10, 3); // Such that calculation takes place here
		printTimeResult("Multi Model", toc(false));

		return new ReturnType(process, multiModel);
	}

}
