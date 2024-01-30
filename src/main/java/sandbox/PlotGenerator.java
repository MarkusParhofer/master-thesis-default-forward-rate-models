package sandbox;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import info.quantlab.masterthesis.functional.FunctionsOnRandomVariables;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class PlotGenerator {

    public static void main(String[] args) throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();


        properties.put("fixingTimes", new double[] { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 });
        properties.put("liborPeriodLength", 2.0);
        properties.put("numberOfLiborPeriods", 5);
        properties.put("stateSpace", "NORMAL");
        properties.put("measure", "SPOT");

        properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("numberOfFactors", 5);
        properties.put("volatilityParams", new double[] {0.1, 0.0, 0.25, 0.1});
        properties.put("correlationDecayParam", 0.2);
        properties.put("displacement", 0.0001);
        properties.put("initialRatesNonDefaultable", new double[] { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 });

        properties.put("numberOfExtraFactors", 2);
        properties.put("initialRatesDefaultable", new double[] { 0.045, 0.047, 0.07, 0.050, 0.039, 0.025 });
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "NORMAL");
        properties.put("freeParamsSeed", 2000);
        properties.put("freeParamsRange", 0.5);
        // properties.put("freeParamsGenerator", freeParamsGenerator); // None needed because the default one is great! lol...
        properties.put("numberOfPaths", 10000);
        properties.put("brownianMotionSeed", 1959);

        properties.put("finiteDifferenceDelta", 1E-6);


        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER);
        factory.setProperties(properties);


        double[] timeDeltas = {1.0, 0.5, 0.1, 0.05, 0.01, 0.001};
        int[] numberOfNegativePathsEuler = new int[timeDeltas.length];
        int[] numberOfNegativePathsMilstein = new int[timeDeltas.length];

        for(int index=0; index < timeDeltas.length; index++) {
            Map<String, Object> specialProps = new HashMap<>();
            specialProps.put("simulationtimedelta", timeDeltas[index]);
            factory.setProperties(specialProps);
            DefaultableLIBORMarketModel model = factory.createDefaultableModel();
            MonteCarloProcess process = factory.createNumericalScheme(model);
            boolean[] negativePath = new boolean[10000];
            Arrays.fill(negativePath, false);
            int numberOfNegativePaths = 0;
            for(int libor = 0; libor < model.getNumberOfLibors(); libor++) {
                final double liborTime = model.getLiborPeriod(libor);
                for(int timeIndex = 1; timeIndex < process.getTimeDiscretization().getNumberOfTimes(); timeIndex++) {
                    final RandomVariable spreadAtTimeIndex = model.getLIBORSpreadAtGivenTimeIndex(process, timeIndex, libor);
                    if (spreadAtTimeIndex.getMin() > 0) {
                        continue;
                    }
                    for (int path = 0; path < spreadAtTimeIndex.size(); path++) {
                        if (negativePath[path])
                            continue;

                        if (spreadAtTimeIndex.get(path) < 0) {
                            negativePath[path] = true;
                            numberOfNegativePaths++;
                        }
                    }
                }
            }
        }

        DefaultableLIBORMarketModel defaultableModel = factory.createDefaultableModel();
    }
}
