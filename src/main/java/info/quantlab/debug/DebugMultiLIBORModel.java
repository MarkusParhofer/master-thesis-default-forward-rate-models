package info.quantlab.debug;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionView;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;

import java.beans.PropertyEditorSupport;
import java.util.HashMap;
import java.util.Map;

public class DebugMultiLIBORModel {
    public static class Models {
        private Models() {
            baseModel = null;
            debtorModel = null;
            multiModel = null;
            baseProcess = null;
            debtorProcess = null;
            multiProcess = null;
        }
        public LIBORMarketModel baseModel;
        public DefaultableLIBORMarketModel debtorModel;
        public MultiLIBORVectorModel multiModel;
        public MonteCarloProcess baseProcess;
        public MonteCarloProcess debtorProcess;
        public MonteCarloProcess multiProcess;


        static Models create() throws CalculationException {
            Models result = new Models();
            DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
            Map<String, Object> properties = new HashMap<>();

            properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
            properties.put("liborPeriodLength", 0.5);
            properties.put("numberOfLiborPeriods", 10);
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

            result.baseModel = factory.createBaseModel();

            // Def Model 1: Creditor (lower Covariance and lower Spread):
            properties.put("numberOfExtraFactors", 2);
            properties.put("initialRatesDefaultable", new double[] {0.028, 0.03, 0.031, 0.033, 0.034, 0.036});
            properties.put("simulationModel", "SPREADS");
            properties.put("stateSpaceOfSpread", "LOGNORMAL");
            properties.put("freeParamsSeed", 2000);
            properties.put("freeParamsRange", 0.4);
            factory.setProperties(properties);

            DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(result.baseModel);


            // Def Model 2: Debtor (lower Covariance and lower Spread):
            properties.put("numberOfExtraFactors", 3);
            properties.put("initialRatesDefaultable", new double[] {0.032, 0.036, 0.039, 0.04, 0.04, 0.043});
            properties.put("simulationModel", "SPREADS");
            properties.put("stateSpaceOfSpread", "LOGNORMAL");
            properties.put("freeParamsSeed", 1030);
            properties.put("freeParamsRange", 0.6);
            factory.setProperties(properties);

            result.debtorModel = factory.createDefaultableModel(result.baseModel);


            result.multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditor, result.debtorModel}, result.baseModel);

            // Simulation:
            properties.put("numberOfPaths", 5000);
            properties.put("brownianMotionSeed", 30312);
            properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
            factory.setProperties(properties);

            result.multiProcess = factory.createNumericalScheme(result.multiModel);
            BrownianMotion motion = (BrownianMotion) result.multiProcess.getStochasticDriver();

            Integer[] baseFactors = new Integer[result.baseModel.getNumberOfFactors()];
            for(int i=0; i < baseFactors.length; i++) baseFactors[i] = i;
            BrownianMotion baseBM = new BrownianMotionView(motion, baseFactors);
            result.baseProcess = new EulerSchemeFromProcessModel(result.baseModel, baseBM, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);

            Integer[] debtorFactors = new Integer[result.debtorModel.getNumberOfFactors()];
            System.arraycopy(baseFactors, 0, debtorFactors, 0, baseFactors.length);
            for(int i=baseFactors.length; i < result.debtorModel.getNumberOfFactors(); i++)
                debtorFactors[i] = i + result.multiModel.getDefaultableModel(0).getNumberOfFactors() - result.baseModel.getNumberOfFactors();
            BrownianMotion debtorBM = new BrownianMotionView(motion, debtorFactors);
            result.debtorProcess = new net.finmath.montecarlo.process.EulerSchemeFromProcessModel(result.debtorModel, debtorBM, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);
            return result;
        }
    }


    public static void main(String[] args) throws CalculationException {
        Models models = Models.create();
        int timeIndicesToCheck = 9;
        System.out.print("mean LIBOR 9");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9s", "i="+i);
        }

        System.out.println("\n" + "_".repeat(70) + "\n");
        System.out.printf("%15s", "MLM Base:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.multiModel.getNonDefaultableLIBOR(models.multiProcess, i, 9).getAverage());
        }
        System.out.printf("\n%15s", "SA  Base:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.baseModel.getLIBOR(models.baseProcess, i, 9).getAverage());
        }

        System.out.println("\n" + "_".repeat(70) + "\n");
        System.out.printf("%15s", "MLM Debtor:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.multiModel.getDefaultableLIBOR(models.multiProcess, i, 9, 1).getAverage());
        }
        System.out.printf("\n%15s", "SA  Debtor:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.debtorModel.getLIBOR(models.debtorProcess, i, 9).getAverage());
        }
        // Mean
        System.out.println("\n" + "_".repeat(70) + "\n");
        System.out.printf("%15s", "MLM µ Debtor:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.multiModel.getDrift(models.multiProcess, i, models.multiProcess.getProcessValue(i), null)[29].getAverage());
        }
        System.out.printf("\n%15s", "SA  µ Debtor:");
        for(int i = 0; i < timeIndicesToCheck; i++) {
            System.out.printf("%9.6f", models.debtorModel.getDrift(models.debtorProcess, i, models.debtorProcess.getProcessValue(i), null)[19].getAverage());
        }

        // FL
        System.out.println("\n" + "_".repeat(70) + "\n");
        for(int j=0; j < models.debtorModel.getNumberOfFactors(); j++) {
            System.out.printf("\n%19s", "MLM fl Debtor " + j + ":");
            if(j < models.baseModel.getNumberOfFactors()) {
                for (int i = 0; i < timeIndicesToCheck; i++) {
                    System.out.printf("  %18.15f", models.multiModel.getFactorLoading(models.multiProcess, i, 29, models.multiProcess.getProcessValue(i))[j].getAverage());
                }
            } else {
                final int facIndex = j + models.multiModel.getDefaultableModel(0).getNumberOfFactors() - models.baseModel.getNumberOfFactors();
                for (int i = 0; i < timeIndicesToCheck; i++) {
                    System.out.printf("  %18.15f", models.multiModel.getFactorLoading(models.multiProcess, i, 29, models.multiProcess.getProcessValue(i))[facIndex].getAverage());
                }
            }
            System.out.printf("\n%19s", "SA  fl Debtor " + j + ":");
            for (int i = 0; i < timeIndicesToCheck; i++) {
                System.out.printf("  %18.15f", models.debtorModel.getFactorLoading(models.debtorProcess, i, 29, models.debtorProcess.getProcessValue(i))[j].getAverage());
            }
        }
    }
}
