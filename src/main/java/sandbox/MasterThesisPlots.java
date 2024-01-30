package sandbox;

import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.plots.Named;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

public class MasterThesisPlots {


    public static void main(String[] args) throws CalculationException {
        createSurvivalProbabilityPlot();
    }

    public static char subscript(int number) {
        return (char)('\u2080' + number);
    }
    public static void createSurvivalProbabilityPlot() throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();


        properties.put("fixingTimes", new double[]{0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 2.0);
        properties.put("numberOfLiborPeriods", 2);
        properties.put("stateSpace", "NORMAL");
        properties.put("measure", "SPOT");

        properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("numberOfFactors", 2);
        properties.put("volatilityParams", new double[]{0.1, 0.0, 0.25, 0.1});
        properties.put("correlationDecayParam", 0.2);
        properties.put("displacement", 0.0001);
        properties.put("initialRatesNonDefaultable", new double[]{0.035, 0.043, 0.05, 0.041, 0.035, 0.02});

        properties.put("numberOfExtraFactors", 2);
        properties.put("initialRatesDefaultable", new double[]{0.045, 0.047, 0.07, 0.050, 0.039, 0.025});
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "NORMAL");
        properties.put("freeParamsSeed", 2000);
        properties.put("freeParamsRange", 0.5);
        // properties.put("freeParamsGenerator", freeParamsGenerator); // None needed because the default one is great! lol...
        properties.put("numberOfPaths", 3); // Very little paths for efficiency
        properties.put("brownianMotionSeed", 1959);
        properties.put("simulationtimedelta", 0.01);

        properties.put("finiteDifferenceDelta", 1E-6);


        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER);
        factory.setProperties(properties);
        DefaultableLIBORMarketModel defModel = factory.createDefaultableModel();
        MonteCarloProcess defSimModel = factory.createNumericalScheme(defModel);

        final double[] defaultTimes = {0.5, 2.1, 3.4};

        DoubleUnaryOperator[] myFunctionArray = new DoubleUnaryOperator[defaultTimes.length + 1];
        final DoubleToIntFunction getTimeIndex = operand -> {
            int timeIndex = defSimModel.getTimeIndex(operand);
            return timeIndex < 0 ? - timeIndex - 1 : timeIndex;
        };

        myFunctionArray[0] = (operand) -> {
            try {
                return defModel.getDefaultableBond(defSimModel, operand, 4.0).get(1);
            } catch (CalculationException e) {
                return - 1.0;
            }
        };

        for (int i = 1; i < myFunctionArray.length; i++) {
            final int index = i - 1;
            myFunctionArray[i] = (operand) -> {
                if (operand >= defaultTimes[index])
                    return 0.0;
                try {
                    return defModel.getDefaultableBond(defSimModel, operand, 4.0).get(1);
                } catch (CalculationException e) {
                    return -1.0;
                }
            };
        }

        EasyPlot2D plotDefaultableBonds;


        List<DoubleUnaryOperator> myList = Arrays.asList(myFunctionArray);
        List<Named<DoubleUnaryOperator>> myFunctionList =
                myList.stream().map(operator -> new Named<>(
                                (myList.indexOf(operator) < 1 ? "P\u1D48*(t, 4.0)(\u03c9)" : ("P\u1D48(t, 4.0)(\u03c9" + subscript(myList.indexOf(operator)) + ")")), operator)).
                        collect(Collectors.toList());
        plotDefaultableBonds = new EasyPlot2D(defSimModel.getTime(0), 4.0, 51, myFunctionList);

        for (int i = 0; i < myFunctionArray.length; i++) {
            plotDefaultableBonds.changePlotMarker(i, null);
        }
        plotDefaultableBonds.setYAxisLabel("Bond Value");
        plotDefaultableBonds.setXAxisLabel("t");
        plotDefaultableBonds.setXRange(0.0, 4.0);
        plotDefaultableBonds.setYRange(0.5, 1.0);
        plotDefaultableBonds.setTitle("Defaultable Bonds");
        plotDefaultableBonds.setIsLegendVisible(true);
        plotDefaultableBonds.show();

        String pathString = "Graphs\\" + "ForMasterThesis\\";
        File directory = new File(pathString);
        boolean directoryExists = directory.exists();
        if(!directoryExists)
            directoryExists = directory.mkdirs();
        File picFile = new File(pathString + "SurvivalProbability.png");
        try {
            plotDefaultableBonds.saveAsPNG(picFile, 1200, 600);
        } catch(IOException ex) {
            System.out.println("Could not save Plot to file. See Stack trace for more.");
            ex.printStackTrace();
        }

    }
}
