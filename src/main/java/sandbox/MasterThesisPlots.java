package sandbox;

import info.quantlab.debug.Time;
import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.easyplot.PlotFactory;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import info.quantlab.masterthesis.functional.FunctionsOnRandomVariables;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import info.quantlab.masterthesis.process.MilsteinSchemeFiniteDifference;
import info.quantlab.masterthesis.products.DefaultableCouponBondForward;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.assetderivativevaluation.models.BlackScholesModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.plots.GraphStyle;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

@SuppressWarnings("unused")
public class MasterThesisPlots extends Time {


    public static void main(String[] args) throws CalculationException {
        // creditorModel();
        // debtorModel();
        createCreditorAndPayerDefaultablePlot();
    }

    public static void creditorModel() throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 0.5);
        properties.put("numberOfLiborPeriods", 20);
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
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        MonteCarloProcess process = factory.createNumericalScheme(modelCreditor);

        final PlotFactory.PathNamer namer = (ind) -> "L(t;T" + subscript(ind) + ",T" + subscript(ind+1) + ")";
        PlotFactory.ProcessOperator pOP_Cred = (_process, timeIndex, component) -> new Scalar(modelCreditor.getDefaultableLIBOR(_process, timeIndex, component).getAverage());
        EasyPlot2D plotCreditor = PlotFactory.PlotProcessPaths(pOP_Cred, process, 0, 1, 1, 0.0, base.getLiborPeriod(19), (ind)-> namer.apply(0));
        for(int i=1; i < base.getNumberOfComponents(); i++) {
            final int comp = i;
            PlotFactory.PlotProcessPaths(plotCreditor, pOP_Cred, process, i, 1, 1, 0.0, base.getLiborPeriod(19), (ind)->namer.apply(comp));
        }
        for(int i = 0; i < plotCreditor.getNumberOfPlots(); i++) {
            //plotCreditor.changePlotMarker(i, null);
        }
        plotCreditor.setTitle("Creditor defaultable LIBOR rates");
        plotCreditor.setYAxisLabel("L^d(t)");
        plotCreditor.setXAxisLabel("t");
        plotCreditor.show();

    }

    public static void debtorModel() throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 0.5);
        properties.put("numberOfLiborPeriods", 20);
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
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        MonteCarloProcess process = factory.createNumericalScheme(modelDebtor);

        final PlotFactory.PathNamer namer = (ind) -> "L(t;T" + subscript(ind) + ",T" + subscript(ind+1) + ")";
        PlotFactory.ProcessOperator pOP_Cred = (_process, timeIndex, component) -> new Scalar(modelDebtor.getDefaultableLIBOR(_process, timeIndex, component).getAverage());
        EasyPlot2D plotDebtor = PlotFactory.PlotProcessPaths(pOP_Cred, process, 0, 1, 1, 0.0, base.getLiborPeriod(19), (ind)-> namer.apply(0));
        for(int i=1; i < base.getNumberOfComponents(); i++) {
            final int comp = i;
            PlotFactory.PlotProcessPaths(plotDebtor, pOP_Cred, process, i, 1, 1, 0.0, base.getLiborPeriod(19), (ind)->namer.apply(comp));
        }
        for(int i = 0; i < plotDebtor.getNumberOfPlots(); i++) {
            //plotCreditor.changePlotMarker(i, null);
        }
        plotDebtor.setTitle("Debtor defaultable LIBOR rates");
        plotDebtor.setYAxisLabel("L^d(t)");
        plotDebtor.setXAxisLabel("t");
        plotDebtor.show();

    }

    public static void createCreditorAndPayerDefaultablePlot() throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 0.5);
        properties.put("numberOfLiborPeriods", 20);
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
        properties.put("freeParamsRange", 0.4);
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
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        MonteCarloProcess process = factory.createNumericalScheme(multiModel);
        MonteCarloProcess secProcess = factory.createNumericalScheme(modelDebtor);

        double[] couponRates = new double[5];
        Arrays.fill(couponRates, 0.03 * 0.5);
        double nominal = 1000;
        int maturityIndex = 3;
        final int numberOfValuations = 13;
        double[] valuesCouponForwardBothD = new double[numberOfValuations];
        double[] valuesCouponForwardDebtD = new double[numberOfValuations];
        double[] valuesCouponForwardDebtD2 = new double[numberOfValuations];
        double[] valuesCouponForwardNoD = new double[numberOfValuations];
        double[] xPoints = new double[numberOfValuations];
        System.out.printf("%10s      %10s     %10s     %10s     %10s%n",
                "Maturity", "C, D def.", "D def. 1", "D def. 2", "Non def.");
        for(int i = 0; i < numberOfValuations; i++) {
            DefaultableCouponBondForward forward = new DefaultableCouponBondForward(nominal, nominal, couponRates, i + maturityIndex);
            valuesCouponForwardBothD[i] = forward.getValue(0, multiModel, process, 1, 0).getAverage();
            valuesCouponForwardDebtD[i] = forward.getValue(0, modelDebtor, multiModel.getDefaultableProcess(process, 1)).getAverage();
            valuesCouponForwardDebtD2[i] = forward.getValue(0, modelDebtor, secProcess).getAverage();
            valuesCouponForwardNoD[i] = forward.getValue(0, base, multiModel.getNonDefaultableProcess(process)).getAverage();
            xPoints[i] = multiModel.getLiborPeriod(i+maturityIndex);
            System.out.printf("%2d (T_S=%3.1f):     %10.5f     %10.5f     %10.5f     %10.5f%n",
                    i + maturityIndex, xPoints[i],
                    valuesCouponForwardBothD[i], valuesCouponForwardDebtD[i], valuesCouponForwardDebtD2[i], valuesCouponForwardNoD[i]);
        }

        EasyPlot2D plot = new EasyPlot2D("Creditor and Debtor defaultable", xPoints, valuesCouponForwardBothD);
        plot.addPlot("Debtor defaultable", xPoints, valuesCouponForwardDebtD);
        plot.changePlotColor(1, EasyPlot2D.getDefaultColor(1));
        // plot.addPlot("Non defaultable valuation", xPoints, valuesCouponForwardNoD);

        plot.setTitle("Valuation of Coupon Bond Forwards");
        plot.setXAxisLabel("Time to maturity");
        plot.setYAxisLabel("Price");
        plot.setIsLegendVisible(true);
        plot.show();
    }


    public static void createNumericalSchemesProcessPlot() throws CalculationException {
        final int numberOfPaths = 1000;
        final int numberOfTimeSteps = 30;
        final double initialValue = 0.1;
        final double riskFreeRate = -0.1;
        final double volatility = 0.7;

        BlackScholesModel geometricBMExact = new BlackScholesModel(initialValue, riskFreeRate, volatility);
        BlackScholesModel geometricBMAppr = new BlackScholesModel(initialValue, riskFreeRate, volatility) {
            @Override
            public RandomVariable[] getInitialState(MonteCarloProcess process) {
                return super.getInitialValue(process);
            }

            @Override
            public RandomVariable[] getDrift(final MonteCarloProcess process, final int timeIndex, final RandomVariable[] realizationAtTimeIndex, final RandomVariable[] realizationPredictor) {
                RandomVariable[] drift = super.getDrift(process, timeIndex, realizationAtTimeIndex, realizationPredictor).clone();
                drift[0] = drift[0].mult(realizationAtTimeIndex[0]);
                return drift;
            }

            @Override
            public RandomVariable[] getFactorLoading(final MonteCarloProcess process, final int timeIndex, final int component, final RandomVariable[] realizationAtTimeIndex) {
                RandomVariable[] fl = super.getFactorLoading(process, timeIndex, component, realizationAtTimeIndex).clone();
                fl[0] = fl[0].mult(realizationAtTimeIndex[0]);
                return fl;
            }

            @Override
            public RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, final int componentIndex, final RandomVariable randomVariable) {
                return randomVariable;
            }

            @Override
            public RandomVariable applyStateSpaceTransformInverse(final MonteCarloProcess process, final int timeIndex, final int componentIndex, final RandomVariable randomVariable) {
                return randomVariable;
            }
        };

        TimeDiscretization timeDisc = new TimeDiscretizationFromArray(0.0, numberOfTimeSteps, 0.2);

        BrownianMotion bm = new BrownianMotionFromMersenneRandomNumbers(timeDisc, 1, numberOfPaths, 21052);

        MonteCarloProcess simulationExact = new EulerSchemeFromProcessModel(geometricBMExact, bm, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);
        MonteCarloProcess simulationMilst = new MilsteinSchemeFiniteDifference(geometricBMAppr, bm, MilsteinSchemeFiniteDifference.FiniteDifferenceMethod.CENTRAL, 0.02);
        MonteCarloProcess simulationEuler = new EulerSchemeFromProcessModel(geometricBMAppr, bm, EulerSchemeFromProcessModel.Scheme.EULER);

        int firstPath = 15;
        int numberOfPathsPlot = 5;
        int firstNeg = 0;
        for(int i=1; i < numberOfTimeSteps + 1; i++) {
            firstNeg = FunctionsOnRandomVariables.findFirstNegativePath(simulationEuler.getProcessValue(i, 0));
            if(firstNeg < 0) {
                if(i == numberOfTimeSteps) {
                    System.out.println("Negative Path not found!");
                }
            }
            else {
                System.out.printf("Negative Path found at time index %d for path %d!%n", i, firstNeg);
                if(firstNeg + numberOfPathsPlot > numberOfPaths) {
                    firstPath = firstNeg - numberOfPathsPlot;
                }
                else{
                    firstPath = firstNeg;
                }
                break;
            }

        }

        final int firstIndexN = firstNeg - 1;
        EasyPlot2D directVergleich = PlotFactory.PlotProcessPaths(simulationExact, 0, firstNeg, 5, ind -> "Exact Scheme \u03c9" + subscript(ind - firstIndexN));
        PlotFactory.PlotProcessPaths(directVergleich, simulationMilst, 0, firstNeg, 5, ind -> "Milstein Scheme  \u03c9" + subscript(ind - firstIndexN));
        PlotFactory.PlotProcessPaths(directVergleich, simulationEuler, 0, firstNeg, 5, ind -> "Euler Scheme  \u03c9" + subscript(ind - firstIndexN));
        for (int i = 0; i < directVergleich.getNumberOfPlots(); i++) {
            GraphStyle gs;
            Shape shape;
            if(Math.floorDiv(i, 2) == 0) {
                shape = new Rectangle(-3, -3, 6, 6);
            }
            else if(Math.floorDiv(i, 2) == 1) {
                shape = new Polygon(new int[]{-3, 3, 0}, new int[]{3, 3, -3}, 3);
            }
            else{
                shape = new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0);
            }

            if(i % 2 == 0) {
                double color_mult = (1.0 - Math.floorDiv(i, 2) / 3.0) + 0.3;

                final Stroke stroke = new BasicStroke(2.0f);
                final Color color = new Color(0, (int) (0.4470 * 255 * color_mult), (int) (0.7410 * 255 * color_mult));
                gs = new GraphStyle(shape, stroke, color);
            }
            else {
                double color_mult = (1.0 - Math.floorDiv(i, 2) / 5.0);
                final Stroke stroke = new BasicStroke(2.0f);
                final Color color = new Color((int) (0.8500 * 255 * color_mult), (int) (0.3250 * 255 * color_mult), (int) (0.0980 * 255 * color_mult));
                gs = new GraphStyle(shape, stroke, color);
            }
            directVergleich.changePlotStyle(i, gs);
        }
        directVergleich.setYAxisLabel("X(t)");
        directVergleich.setXAxisLabel("t");
        directVergleich.setTitle("Comparison of Schemes");
        directVergleich.setIsLegendVisible(true);
        directVergleich.show();

        EasyPlot2D plotEulerFunc = PlotFactory.PlotProcessPaths(simulationExact, 0, firstPath, numberOfPathsPlot);
        for (int i = 0; i < plotEulerFunc.getNumberOfPlots(); i++) {
            plotEulerFunc.changePlotMarker(i, null);
        }
        plotEulerFunc.setYAxisLabel("X(t)");
        plotEulerFunc.setXAxisLabel("t");
        plotEulerFunc.setTitle("Functional Euler Scheme");
        plotEulerFunc.show();

        EasyPlot2D plotMilst = PlotFactory.PlotProcessPaths(simulationMilst, 0, firstPath, numberOfPathsPlot);
        for (int i = 0; i < plotMilst.getNumberOfPlots(); i++) {
            plotMilst.changePlotMarker(i, null);
        }
        plotMilst.setYAxisLabel("X(t)");
        plotMilst.setXAxisLabel("t");
        plotMilst.setTitle("Milstein Scheme");
        plotMilst.show();

        EasyPlot2D plotEuler = PlotFactory.PlotProcessPaths(simulationEuler, 0, firstPath, numberOfPathsPlot);
        for (int i = 0; i < plotEuler.getNumberOfPlots(); i++) {
            plotEuler.changePlotMarker(i, null);
        }
        plotEuler.setYAxisLabel("X(t)");
        plotEuler.setXAxisLabel("t");
        plotEuler.setTitle("Euler Scheme");
        plotEuler.show();

        String pathString = "Graphs\\" + "ForMasterThesis\\";
        File directory = new File(pathString);
        boolean directoryExists = directory.exists();
        if(!directoryExists)
            directoryExists = directory.mkdirs();

        try {
            if(!directoryExists)
                throw new IOException("Directory could not be created!");
            final int width = 800;
            final int height = 500;
            directVergleich.saveAsPNG(new File(pathString + "Schemes_DirectCompare.png"), width, height);
            plotEulerFunc.saveAsPNG(new File(pathString + "EulerFunctional_Scheme.png"), width, height);
            plotMilst.saveAsPNG(new File(pathString + "Milstein_Scheme.png"), width, height);
            plotEuler.saveAsPNG(new File(pathString + "Euler_Scheme.png"), width, height);
            System.out.println("Successfully saved Plots");
        } catch(IOException ex) {
            System.out.println("Could not save Plot to file. See Stack trace for more.");
            ex.printStackTrace();
        }
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

        PlotFactory.Process operator = (time) -> {
            final double[] defaultState = new double[] {
                    1.0,
                    time < defaultTimes[0] ? 1.0 : 0.0,
                    time < defaultTimes[1] ? 1.0 : 0.0,
                    time < defaultTimes[2] ? 1.0 : 0.0
            };
            return (new RandomVariableFromDoubleArray(time, defaultState)).mult(defModel.getDefaultableBond(defSimModel, time, 4.0).get(1));
        };

        PlotFactory.PathNamer namer = (index) -> (index < 1 ? "P\u1D48*(t, 4.0)(\u03c9)" : ("P\u1D48(t, 4.0)(\u03c9" + subscript(index) + ")"));

        EasyPlot2D plotDefaultableBonds = PlotFactory.PlotProcessPaths(operator, 0.0, 4.0, 51, 0, 4, namer);

        for (int i = 0; i < plotDefaultableBonds.getNumberOfPlots(); i++) {
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
