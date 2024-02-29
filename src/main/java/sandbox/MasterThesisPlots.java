package sandbox;

import info.quantlab.debug.Time;
import info.quantlab.easyplot.EasyPlot2D;
import info.quantlab.easyplot.PlotFactory;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.factory.DefaultableLIBORModelFactory;
import info.quantlab.masterthesis.functional.FunctionsOnRandomVariables;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import info.quantlab.masterthesis.process.MilsteinSchemeFiniteDifference;
import info.quantlab.masterthesis.products.CancellableLoan;
import info.quantlab.masterthesis.products.DefaultableCouponBondForward;
import info.quantlab.masterthesis.products.LoanProduct;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.BrownianMotionView;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.assetderivativevaluation.models.BlackScholesModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.plots.GraphStyle;
import net.finmath.plots.Plotable2D;
import net.finmath.plots.PlotablePoints2D;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import static info.quantlab.easyplot.PlotRunner.runConfigLoop;

@SuppressWarnings("unused")
public class MasterThesisPlots extends Time {


    public static void main(String[] args) throws CalculationException {
        runConfigLoop(createCancellableLoanPlot());
    }

    public static EasyPlot2D[] createBondPriceErrorPlot() throws CalculationException {

        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();
        {
            final double initialSpreadDebt = 0.025;


            properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
            properties.put("liborPeriodLength", 1.0);
            properties.put("numberOfLiborPeriods", 15);
            properties.put("stateSpace", "NORMAL");
            properties.put("measure", "SPOT");

            final double[] initRatesNonDef = new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032};
            properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
            properties.put("numberOfFactors", 5);
            properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
            properties.put("correlationDecayParam", 0.2);
            properties.put("displacement", 0.0001);
            properties.put("initialRatesNonDefaultable", initRatesNonDef);
            properties.put("simulationTimeDelta", 0.02);

            factory.setProperties(properties);

            // Def Model 1: Creditor (lower Covariance and lower Spread):
            double[] initRatesCred = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
            for (int i = 0; i < initRatesCred.length; i++) {
                initRatesCred[i] += initialSpreadDebt;
            }

            properties.put("numberOfExtraFactors", 2);
            properties.put("initialRatesDefaultable", initRatesCred);
            properties.put("simulationModel", "SPREADS");
            properties.put("stateSpaceOfSpread", "LOGNORMAL");
            properties.put("freeParamsSeed", 33);
            properties.put("freeParamsRange", 0.4);
            factory.setProperties(properties);

            // Simulation:
            properties.put("numberOfPaths", 25000);
            properties.put("brownianMotionSeed", 4692);
            properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
            factory.setProperties(properties);
        }
        LIBORMarketModel base = factory.createBaseModel();
        DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);
        MonteCarloProcess process = factory.createNumericalScheme(modelCreditor);

        EasyPlot2D[] plots = new EasyPlot2D[6];
        // Bond Price error and Forward error
        final int jump = 5;
        final int numberOfTimes0 = process.getTimeDiscretization().getNumberOfTimeSteps() / 5 + 1;
        double[][] bondPrices = new double[2][numberOfTimes0]; // 0 Analytic, 1 Approximated
        double[] bondErrorsAbs = new double[numberOfTimes0];
        double[] bondErrorsRel = new double[numberOfTimes0];
        double[] xPoints = new double[numberOfTimes0];
        for (int i = 0; i < numberOfTimes0; i++) {
            // Price Bond:
            xPoints[i] = process.getTime(i * 5);
            bondPrices[0][i] = 1.0 / (modelCreditor.getDefaultableForwardRate(process, 0, 0, xPoints[i]).doubleValue()* xPoints[i] + 1);
            bondPrices[1][i] = modelCreditor.getDefaultableNumeraire(process, xPoints[i]).vid(1.0).getAverage();
            bondErrorsAbs[i] = bondPrices[0][i] - bondPrices[1][i];
            bondErrorsRel[i] = bondErrorsAbs[i] / bondPrices[0][i];
            bondPrices[0][i] *= 100.0;
            bondPrices[1][i] *= 100.0;
        }

        plots[0] = new EasyPlot2D(new String[] {"Analytic", "Monte Carlo"}, xPoints, bondPrices);
        plots[0].setTitle("Defaultable Bond Prices (Δt=0.02, n=50000)");
        plots[0].setXAxisLabel("Maturity");
        plots[0].setYAxisLabel("Bond Price P\u1D48(0;T)");
        plots[0].setYAxisNumberFormat(new DecimalFormat("00.00%;-00.00%", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[0].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[0].setIsLegendVisible(true);
        plots[0].show();

        plots[1] = new EasyPlot2D("", xPoints, bondErrorsAbs);
        plots[1].setTitle("Approximation Error of Defaultable Bond Prices (Δt=0.02, n=50000)");
        plots[1].setXAxisLabel("Maturity");
        plots[1].setYAxisLabel("Error");
        plots[1].setYAxisNumberFormat(new DecimalFormat("0.0E00;-0.0E00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[1].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[1].show();

        plots[2] = new EasyPlot2D("", xPoints, bondErrorsRel);
        plots[2].setTitle("Relative Approximation Error of Defaultable Bond Prices (Δt=0.02, n=50000)");
        plots[2].setXAxisLabel("Maturity");
        plots[2].setYAxisLabel("rel. Error");
        plots[2].setYAxisNumberFormat(new DecimalFormat("0.0E00;-0.0E00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[2].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[2].show();


        final int numberOfTimes1 = modelCreditor.getNumberOfLIBORPeriods() - 1;
        double[][] forwardPrices = new double[2][numberOfTimes1]; // 0 Analytic, 1 Approximated
        double[] forwardErrorsAbs = new double[numberOfTimes1];
        double[] forwardErrorsRel = new double[numberOfTimes1];
        double[] fixingTimes = new double[numberOfTimes1];
        for (int i = 0; i < numberOfTimes1; i++) {
            fixingTimes[i] = modelCreditor.getLiborPeriod(i + 1);
            int timeIndex = process.getTimeIndex(fixingTimes[i]);
            if(timeIndex < 0) {
                timeIndex = -timeIndex - 1;
            }
            final double payment = modelCreditor.getLiborPeriod(i + 2);
            final double bondLiborPlusOne = 1.0 / (modelCreditor.getDefaultableForwardRate(process, 0, 0, payment).doubleValue() * payment + 1);
            forwardPrices[0][i] = modelCreditor.getDefaultableLIBOR(process, 0, i + 1).doubleValue() * (payment - fixingTimes[i]) * bondLiborPlusOne;
            final RandomVariable libor = modelCreditor.getDefaultableLIBOR(process, timeIndex, i+1).mult(payment - fixingTimes[i]);
            forwardPrices[1][i] = libor.div(modelCreditor.getDefaultableNumeraire(process, payment)).getAverage();
            forwardErrorsAbs[i] = forwardPrices[0][i] - forwardPrices[1][i];
            forwardErrorsRel[i] = forwardErrorsAbs[i] / forwardPrices[0][i];
            forwardPrices[0][i] *= 100.0;
            forwardPrices[1][i] *= 100.0;
        }


        plots[3] = new EasyPlot2D(new String[] {"Analytic", "Monte Carlo"}, fixingTimes, forwardPrices);
        plots[3].setTitle("Approximation error of (def.) LIBOR rate (Δt=0.02, n=25000)");
        plots[3].setXAxisLabel("Fixing Time (S)");
        plots[3].setYAxisLabel("L\u1D48(0;S,S+ΔT) P\u1D48(0;S+ΔT) in %");
        plots[3].setYAxisNumberFormat(new DecimalFormat("00.00%;-00.00%", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[3].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[3].setIsLegendVisible(true);
        plots[3].show();

        plots[4] = new EasyPlot2D("", fixingTimes, forwardErrorsAbs);
        plots[4].setTitle("Approximation error of (def.) LIBOR rate (Δt=0.02, n=25000)");
        plots[4].setXAxisLabel("Fixing Time (S)");
        plots[4].setYAxisLabel("Error");
        plots[4].setYAxisNumberFormat(new DecimalFormat("0.0E00;-0.0E00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[4].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[4].show();

        plots[5] = new EasyPlot2D("", fixingTimes, forwardErrorsRel);
        plots[5].setTitle("Relative Approximation error of (def.) LIBOR rate (Δt=0.02, n=25000)");
        plots[5].setXAxisLabel("Fixing Time (S)");
        plots[5].setYAxisLabel("rel. Error");
        plots[5].setYAxisNumberFormat(new DecimalFormat("0.0E00;-0.0E00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[5].setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plots[5].show();


        return plots;
    }

    public static EasyPlot2D[] createQuickMeanPlot() throws CalculationException {
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 0.5);
        properties.put("numberOfLiborPeriods", 20);
        properties.put("stateSpace", "NORMAL");
        properties.put("measure", "SPOT");

        properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("numberOfFactors", 5);
        properties.put("volatilityParams", new double[] {0.02, 0.00, 0.25, 0.15});
        properties.put("correlationDecayParam", 0.2);
        properties.put("displacement", 0.0001);
        properties.put("initialRatesNonDefaultable", new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032});
        properties.put("simulationTimeDelta", 0.02);

        factory.setProperties(properties);

        LIBORMarketModel base = factory.createBaseModel();

        // Def Model 1: Creditor (lower Covariance and lower Spread):
        final double spread = 0.015;
        properties.put("numberOfExtraFactors", 2);
        properties.put("initialRatesDefaultable", new double[] {0.025 + spread, 0.023 + spread, 0.028 + spread, 0.031 + spread, 0.031 + spread, 0.032 + spread});
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsDesc", DefaultableLIBORModelFactory.FreeParametersModel.RANDOM);
        properties.put("freeParams", new double[] {0.04, 0.0, 0.35, 0.02, 0.1});
        properties.put("freeParamsCovModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("freeParamsDisplacement", 0.0001);
        properties.put("freeParamsSeed", 30303);
        properties.put("freeParamsRange", 0.5);
        factory.setProperties(properties);

        final DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);

        // Simulation:
        properties.put("numberOfPaths", 5000);
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        MonteCarloProcess process = factory.createNumericalScheme(modelCreditor);

        EasyPlot2D[] plots = new EasyPlot2D[4 + 3]; // Mean Variance and Correlation to non defaultable model

        final int skipTimePoints = 10;
        final int numberOfPoints = process.getTimeDiscretization().getNumberOfTimeSteps() / skipTimePoints + 1;
        double[] xPoints = new double[numberOfPoints];
        double[][] yPointsMean = new double[modelCreditor.getNumberOfLIBORPeriods() - 1][numberOfPoints]; // We exclude L_0
        double[][] yPointsVar  = new double[modelCreditor.getNumberOfLIBORPeriods() - 1][numberOfPoints]; // We exclude L_0
        double[][] yPointsCorr = new double[modelCreditor.getNumberOfLIBORPeriods() - 1][numberOfPoints]; // We exclude L_0
        double[][] yPointsCorrSL = new double[modelCreditor.getNumberOfLIBORPeriods() - 1][numberOfPoints]; // We exclude L_0
        String[] namesS = new String[modelCreditor.getNumberOfLIBORPeriods() - 1];
        for(int i = 0; i < numberOfPoints; i++) {
            xPoints[i] = process.getTime(i * skipTimePoints);
            for (int j = 0; j < modelCreditor.getNumberOfLIBORPeriods() - 1; j++) {
                if(i==0) {
                    namesS[j] = "Process for Period " + j;
                }
                yPointsMean[j][i] = modelCreditor.getLIBORSpreadAtGivenTimeIndex(process, i * skipTimePoints, j + 1).getAverage();
                yPointsVar[j][i] = modelCreditor.getLIBORSpreadAtGivenTimeIndex(process, i * skipTimePoints, j + 1).getVariance();

                final RandomVariable nonDefLIBOR = modelCreditor.getNonDefaultableLIBOR(process, i * skipTimePoints, j + 1);
                final RandomVariable defLIBOR = modelCreditor.getDefaultableLIBOR(process, i * skipTimePoints, j + 1);
                final RandomVariable spread1 = modelCreditor.getLIBORSpreadAtGivenTimeIndex(process, i * skipTimePoints, j + 1);
                yPointsCorr[j][i] = defLIBOR.covariance(nonDefLIBOR).doubleValue() / Math.sqrt(nonDefLIBOR.getVariance() * defLIBOR.getVariance());
                yPointsCorrSL[j][i] = spread1.covariance(nonDefLIBOR).doubleValue() / Math.sqrt(nonDefLIBOR.getVariance() * spread1.getVariance());
            }
        }
        plots[0] = new EasyPlot2D(namesS, xPoints, yPointsMean);
        plots[0].setTitle("Mean of Spreads");
        plots[0].setXAxisLabel("t");
        plots[0].setYAxisLabel("Spread");
        plots[0].setIsLegendVisible(true);
        plots[0].show();
        //plots[0].writeRecreationJavaFile();

        plots[1] = new EasyPlot2D(namesS, xPoints, yPointsVar);
        plots[1].setTitle("Variance of Spreads");
        plots[1].setXAxisLabel("t");
        plots[1].setYAxisLabel("Variance");
        plots[1].setIsLegendVisible(true);
        plots[1].show();
        //plots[1].writeRecreationJavaFile();

        plots[2] = new EasyPlot2D(namesS, xPoints, yPointsCorr);
        plots[2].setTitle("Correlation between non-defaultable and defaultable LIBORs");
        plots[2].setXAxisLabel("t");
        plots[2].setYAxisLabel("Correlation");
        plots[2].setIsLegendVisible(true);
        plots[2].show();
        //plots[2].writeRecreationJavaFile();

        plots[3] = new EasyPlot2D(namesS, xPoints, yPointsCorrSL);
        plots[3].setTitle("Correlation between non-defaultable LIBORs and Spreads");
        plots[3].setXAxisLabel("t");
        plots[3].setYAxisLabel("Correlation");
        plots[3].setIsLegendVisible(true);
        plots[3].show();
        //plots[3].writeRecreationJavaFile();

        PlotFactory.ProcessOperator pOp = modelCreditor::getLIBORSpreadAtGivenTimeIndex;
        PlotFactory.ProcessOperator pOp2 = modelCreditor::getNonDefaultableLIBOR;
        PlotFactory.ProcessOperator pOp3 = modelCreditor::getDefaultableLIBOR;
        plots[4] = PlotFactory.PlotProcessPaths(pOp, process, modelCreditor.getNumberOfLIBORPeriods() - 1, 0, 150);
        plots[4].show();
        plots[5] = PlotFactory.PlotProcessPaths(pOp2, process, modelCreditor.getNumberOfLIBORPeriods() - 1, 0, 150);
        plots[5].show();
        plots[6] = PlotFactory.PlotProcessPaths(pOp3, process, modelCreditor.getNumberOfLIBORPeriods() - 1, 0, 150);
        plots[6].show();


        return plots;
    }

    public static EasyPlot2D[] createDependencePlots() throws CalculationException {

        final double initialSpreadCred = 0.015;
        final double initialSpreadDebt = 0.025;

        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 1.0);
        properties.put("numberOfLiborPeriods", 10);
        properties.put("stateSpace", "NORMAL");
        properties.put("measure", "SPOT");

        final double[] initRatesNonDef = new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032};
        properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("numberOfFactors", 5);
        properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
        properties.put("correlationDecayParam", 0.2);
        properties.put("displacement", 0.0001);
        properties.put("initialRatesNonDefaultable", initRatesNonDef);
        properties.put("simulationTimeDelta", 0.02);

        factory.setProperties(properties);

        LIBORMarketModel base = factory.createBaseModel();

        // Def Model 1: Creditor (lower Covariance and lower Spread):
        double[] initRatesCred = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
        for (int i = 0; i < initRatesCred.length; i++) {
            initRatesCred[i] += initialSpreadCred;
        }

        properties.put("numberOfExtraFactors", 2);
        properties.put("initialRatesDefaultable", initRatesCred);
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 33);
        properties.put("freeParamsRange", 0.4);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);

        // Def Model 2: Debtor (lower Covariance and lower Spread):
        double[] initRatesDebt = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
        for (int i = 0; i < initRatesDebt.length; i++) {
             initRatesDebt[i] += initialSpreadDebt;
        }

        properties.put("numberOfExtraFactors", 3);
        properties.put("initialRatesDefaultable", initRatesDebt);
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 1030);
        properties.put("freeParamsRange", 0.4);
        factory.setProperties(properties);

        // Simulation:
        properties.put("numberOfPaths", 10000);
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelDebtor = factory.createDefaultableModel(base);

        final int lastTimeIndex = base.getCovarianceModel().getTimeDiscretization().getNumberOfTimeSteps();
        final int numberOfXPoints = 21;
        final int numberOfPlots = 5;
        final double[][] xPoints = new double[2][numberOfXPoints];
        EasyPlot2D[] plots = new EasyPlot2D[numberOfPlots];
        {
            final double[][] yPointsPlot0 = new double[8][numberOfXPoints]; // Covariance and Variance
            final double[][] yPointsPlot1 = new double[4][numberOfXPoints]; // Correlation

            // Plot Creditor Debtor Variance: x = freeParamsRangeD, x = freeParamsRangeC, x = spreadCreditor, x = spreadDebtor, x = spreadDebtor and Creditor
            // First X: freeParamsRangeD
            final double startX0 = 0.4;
            final double deltaX0 = (0.55 - startX0) / 20;
            for (int i = 0; i < numberOfXPoints; i++) {
                xPoints[0][i] = startX0 + i * deltaX0;
                properties.put("numberOfExtraFactors", 2);
                properties.put("initialRatesDefaultable", initRatesCred);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 33);
                properties.put("freeParamsRange", xPoints[0][i]);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelCreditorF = factory.createDefaultableModel(base);

                properties.put("numberOfExtraFactors", 3);
                properties.put("initialRatesDefaultable", initRatesDebt);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 1030);
                properties.put("freeParamsRange", xPoints[0][i]);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelDebtorF = factory.createDefaultableModel(base);
                MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditorF, modelDebtorF, modelCreditor, modelDebtor}, base);

                MonteCarloProcess process = factory.createNumericalScheme(multiModel);

                final RandomVariable credF = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 0);
                final RandomVariable debtF = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 1);
                final RandomVariable cred = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 2);
                final RandomVariable debt = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 3);


                yPointsPlot0[0][i] = credF.getVariance();
                yPointsPlot0[1][i] = debtF.getVariance();
                yPointsPlot0[2][i] = cred.getVariance();
                yPointsPlot0[3][i] = debt.getVariance();
                System.out.println("At " + i + " (fp range): Creditor Variance: " + yPointsPlot0[0][i]);
                System.out.println("At " + i + " (fp range): Debtor   Variance: " + yPointsPlot0[1][i]);

                yPointsPlot0[4][i] = credF.covariance(debtF).doubleValue();
                yPointsPlot0[5][i] = debtF.covariance(cred).doubleValue();
                yPointsPlot0[6][i] = cred.covariance(debt).doubleValue();
                yPointsPlot0[7][i] = debt.covariance(credF).doubleValue();

                yPointsPlot1[0][i] = yPointsPlot0[4][i] / Math.sqrt(yPointsPlot0[0][i] * yPointsPlot0[1][i]);
                yPointsPlot1[1][i] = yPointsPlot0[5][i] / Math.sqrt(yPointsPlot0[1][i] * yPointsPlot0[2][i]);
                yPointsPlot1[2][i] = yPointsPlot0[6][i] / Math.sqrt(yPointsPlot0[2][i] * yPointsPlot0[3][i]);
                yPointsPlot1[3][i] = yPointsPlot0[7][i] / Math.sqrt(yPointsPlot0[3][i] * yPointsPlot0[0][i]);
            }

            String[] plotNames0 = new String[] {
                    "Var(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x)) (model affected by x Axis)",
                    "Var(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x)) (model affected by x Axis)",
                    "Var(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")) (constant model)",
                    "Var(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")) (constant model)",
                    "Cov(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x))",
                    "Cov(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1d9c" + subscript(9) + "(T" + subscript(9) + "))",
                    "Cov(S\u1d9c" + subscript(9) + "(T" + subscript(9) + "), S\u1D48" + subscript(9) + "(T" + subscript(9) + "))",
                    "Cov(S\u1D48" + subscript(9) + "(T" + subscript(9) + "), S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x))",
            };

            String[] plotNames1 = new String[] {
                    "Corr(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x))",
                    "Corr(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1d9c" + subscript(9) + "(T" + subscript(9) + "))",
                    "Corr(S\u1d9c" + subscript(9) + "(T" + subscript(9) + "), S\u1D48" + subscript(9) + "(T" + subscript(9) + "))",
                    "Corr(S\u1D48" + subscript(9) + "(T" + subscript(9) + "), S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x))",
            };

            // Plot result
            plots[0] = new EasyPlot2D(plotNames0, xPoints[0], yPointsPlot0);
            plots[0].setTitle("Terminal Variance and Covariance of Spreads by range of free parameters");
            plots[0].setIsLegendVisible(true);
            plots[0].setXAxisLabel("Free Parameter range");
            plots[0].setYAxisLabel("Var, Cov");
            plots[0].show();
            plots[0].writeRecreationJavaFile();

            // Plot result
            plots[1] = new EasyPlot2D(plotNames1, xPoints[0], yPointsPlot1);
            plots[1].setTitle("Terminal Correlation of Spreads by range of the free Parameter");
            plots[1].setIsLegendVisible(true);
            plots[1].setXAxisLabel("Free Parameter range");
            plots[1].setYAxisLabel("Correlation");
            plots[1].show();
            plots[1].writeRecreationJavaFile();
        }

        {
            final double[][] yPointsPlot2 = new double[8][numberOfXPoints]; // Covariance and Variance
            final double[][] yPointsPlot3 = new double[4][numberOfXPoints]; // Correlation

            // Plot Creditor Debtor Variance: x = freeParamsRangeD, x = freeParamsRangeC, x = spreadCreditor, x = spreadDebtor, x = spreadDebtor and Creditor
            // First X: freeParamsRangeD
            final double startSpread = 0.005;
            final double deltaSpread = 0.002;
            for (int i = 0; i < numberOfXPoints; i++) {
                xPoints[1][i] = startSpread + i * deltaSpread;
                double[] initRatesCredF = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
                double[] initRatesDebtF = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
                for (int j = 0; j < initRatesCredF.length; j++) {
                    initRatesCredF[j] += xPoints[1][i];
                    initRatesDebtF[j] += xPoints[1][i];
                }

                properties.put("numberOfExtraFactors", 2);
                properties.put("initialRatesDefaultable", initRatesCredF);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 33);
                properties.put("freeParamsRange", 0.4);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelCreditorF = factory.createDefaultableModel(base);

                properties.put("numberOfExtraFactors", 3);
                properties.put("initialRatesDefaultable", initRatesDebtF);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 1030);
                properties.put("freeParamsRange", 0.4);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelDebtorF = factory.createDefaultableModel(base);
                MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditorF, modelDebtorF, modelCreditor, modelDebtor}, base);

                MonteCarloProcess process = factory.createNumericalScheme(multiModel);

                final RandomVariable credF = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 0);
                final RandomVariable debtF = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 1);
                final RandomVariable cred = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 2);
                final RandomVariable debt = multiModel.getLIBORSpreadAtGivenTimeIndex(process, lastTimeIndex, 9, 3);


                yPointsPlot2[0][i] = credF.getVariance();
                yPointsPlot2[1][i] = debtF.getVariance();
                yPointsPlot2[2][i] = cred.getVariance();
                yPointsPlot2[3][i] = debt.getVariance();

                System.out.println("At " + i + " (spread): Creditor Variance: " + yPointsPlot2[0][i]);
                System.out.println("At " + i + " (spread): Debtor   Variance: " + yPointsPlot2[1][i]);

                yPointsPlot2[4][i] = credF.covariance(debtF).doubleValue();
                yPointsPlot2[5][i] = debtF.covariance(cred).doubleValue();
                yPointsPlot2[6][i] = cred.covariance(debt).doubleValue();
                yPointsPlot2[7][i] = debt.covariance(credF).doubleValue();

                yPointsPlot3[0][i] = yPointsPlot2[4][i] / Math.sqrt(yPointsPlot2[0][i] * yPointsPlot2[1][i]);
                yPointsPlot3[1][i] = yPointsPlot2[5][i] / Math.sqrt(yPointsPlot2[1][i] * yPointsPlot2[2][i]);
                yPointsPlot3[2][i] = yPointsPlot2[6][i] / Math.sqrt(yPointsPlot2[2][i] * yPointsPlot2[3][i]);
                yPointsPlot3[3][i] = yPointsPlot2[7][i] / Math.sqrt(yPointsPlot2[3][i] * yPointsPlot2[0][i]);
            }

            String[] plotNames0 = new String[] {
                    "Var(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x)) (model affected by x Axis)",
                    "Var(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x)) (model affected by x Axis)",
                    "Var(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")) (constant model)",
                    "Var(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")) (constant model)",
                    "Cov(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x))",
                    "Cov(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1d9c" + subscript(9) + "(T" + subscript(9) + "))",
                    "Cov(S\u1d9c" + subscript(9) + "(T" + subscript(9) + "), S\u1D48" + subscript(9) + "(T" + subscript(9) + "))",
                    "Cov(S\u1D48" + subscript(9) + "(T" + subscript(9) + "), S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x))",
            };

            String[] plotNames1 = new String[] {
                    "Corr(S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x))",
                    "Corr(S\u1D48" + subscript(9) + "(T" + subscript(9) + ")(x), S\u1d9c" + subscript(9) + "(T" + subscript(9) + "))",
                    "Corr(S\u1d9c" + subscript(9) + "(T" + subscript(9) + "), S\u1D48" + subscript(9) + "(T" + subscript(9) + "))",
                    "Corr(S\u1D48" + subscript(9) + "(T" + subscript(9) + "), S\u1d9c" + subscript(9) + "(T" + subscript(9) + ")(x))",
            };

            // Plot result
            plots[2] = new EasyPlot2D(plotNames0, xPoints[1], yPointsPlot2);
            plots[2].setTitle("Terminal Variance and Covariance by initial Spread difference");
            plots[2].setIsLegendVisible(true);
            plots[2].setXAxisLabel("Initial Spread diff");
            plots[2].setYAxisLabel("Var, Cov");
            plots[2].show();
            plots[2].writeRecreationJavaFile();

            // Plot result
            plots[3] = new EasyPlot2D(plotNames1, xPoints[1], yPointsPlot3);
            plots[3].setTitle("Terminal Correlation by initial Spread difference");
            plots[3].setIsLegendVisible(true);
            plots[3].setXAxisLabel("Initial Spread diff");
            plots[3].setYAxisLabel("Correlation");
            plots[3].show();
            plots[3].writeRecreationJavaFile();
        }

        // Plot Creditor Debtor Correlation: x = freeParamsRangeD, x = spreadCreditor, x = spreadDebtor, x = spreadDebtor and Creditor


        return plots;
    }


    public static EasyPlot2D[] createCancellableLoanPlot() throws CalculationException {
        final double initialSpreadCred = 0.025;
        final double initialSpreadDebt = 0.035;

        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        Map<String, Object> properties = new HashMap<>();

        properties.put("fixingTimes", new double[] {0.5, 1.0, 2.0, 4.0, 8.0, 25.0});
        properties.put("liborPeriodLength", 1.0);
        properties.put("numberOfLiborPeriods", 16);
        properties.put("stateSpace", "NORMAL");
        properties.put("measure", "SPOT");

        final double[] initRatesNonDef = new double[] {0.025, 0.023, 0.028, 0.031, 0.031, 0.032};
        properties.put("covarianceModel", DefaultableLIBORModelFactory.CovarianceModel.BLENDED);
        properties.put("numberOfFactors", 5);
        properties.put("volatilityParams", new double[] {0.02, 0.0, 0.25, 0.15});
        properties.put("correlationDecayParam", 0.2);
        properties.put("displacement", 0.0001);
        properties.put("initialRatesNonDefaultable", initRatesNonDef);
        properties.put("simulationTimeDelta", 0.02);

        factory.setProperties(properties);

        LIBORMarketModel base = factory.createBaseModel();

        // Def Model 1: Creditor (lower Covariance and lower Spread):
        double[] initRatesCred = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
        for (int i = 0; i < initRatesCred.length; i++) {
            initRatesCred[i] += initialSpreadCred;
        }

        properties.put("numberOfExtraFactors", 2);
        properties.put("initialRatesDefaultable", initRatesCred);
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 33);
        properties.put("freeParamsRange", 0.5);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);

        // Def Model 2: Debtor (lower Covariance and lower Spread):
        double[] initRatesDebt = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
        for (int i = 0; i < initRatesDebt.length; i++) {
            initRatesDebt[i] += initialSpreadDebt;
        }

        properties.put("numberOfExtraFactors", 3);
        properties.put("initialRatesDefaultable", initRatesDebt);
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 1030);
        properties.put("freeParamsRange", 0.6);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelDebtor = factory.createDefaultableModel(base);

        // Simulation:
        properties.put("numberOfPaths", 10000);
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);


        final int lastTimeIndex = base.getCovarianceModel().getTimeDiscretization().getNumberOfTimeSteps();
        final int numberOfXPoints = 20;
        final int numberOfPlots = 6;
        final double[][] xPoints = new double[numberOfPlots][numberOfXPoints];
        String[] plotNames = new String[numberOfPlots];
        plotNames[0] = "Cost of a cancellable loan in regards to Free Parameters of the debtor";
        plotNames[1] = "Cost of a cancellable loan in regards to Free Parameters of the creditor";
        plotNames[2] = "Cost of a cancellable loan in regards to Free Parameters of both parties";
        plotNames[3] = "Cost of a cancellable loan in regards to the initial spread of the debtor";      
        plotNames[4] = "Cost of a cancellable loan in regards to the initial spread of the creditor";    
        plotNames[5] = "Cost of a cancellable loan in regards to the initial spread of both parties";

        EasyPlot2D[] plotArray = new EasyPlot2D[numberOfPlots];
        {
            final double[][] yPoints0 = new double[3][numberOfXPoints];
            final double[][] yPoints1 = new double[3][numberOfXPoints];
            final double[][] yPoints2 = new double[3][numberOfXPoints];
            final double liborPeriodLength = base.getLiborPeriod(1);
            TimeDiscretization loanTenor = new TimeDiscretizationFromArray(liborPeriodLength, 14, liborPeriodLength);

            final double startX0 = 0.4;
            final double deltaX0 = (0.6 - startX0) / 20;
            for (int i = 0; i < numberOfXPoints; i++) {
                xPoints[0][i] = startX0 + i * deltaX0;
                properties.put("numberOfExtraFactors", 2);
                properties.put("initialRatesDefaultable", initRatesCred);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 33);
                properties.put("freeParamsRange", xPoints[0][i]);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelCreditorF = factory.createDefaultableModel(base);

                properties.put("numberOfExtraFactors", 3);
                properties.put("initialRatesDefaultable", initRatesDebt);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 1030);
                properties.put("freeParamsRange", xPoints[0][i]);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelDebtorF = factory.createDefaultableModel(base);
                MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditorF, modelDebtorF, modelCreditor, modelDebtor}, base);

                MonteCarloProcess process = factory.createNumericalScheme(multiModel);

                final double[] coupons = new double[15];
                final double nominal = 100;
                for(int j=0; j < 15; j++) {
                    coupons[j] = multiModel.getDefaultableLIBOR(process, 0, j, 1).doubleValue();
                    coupons[j] *= nominal * liborPeriodLength;
                }
                coupons[14] += nominal;
                CancellableLoan loanCC0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);

                CancellableLoan loanCC1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);

                CancellableLoan loanCC2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);


                yPoints0[0][i] = loanCC0.getValue(0, multiModel, process).getAverage();
                yPoints0[1][i] = loanCD0.getValue(0, multiModel, process).getAverage();
                yPoints0[2][i] = loanDD0.getValue(0, multiModel, process).getAverage();
                yPoints1[0][i] = loanCC1.getValue(0, multiModel, process).getAverage();
                yPoints1[1][i] = loanCD1.getValue(0, multiModel, process).getAverage();
                yPoints1[2][i] = loanDD1.getValue(0, multiModel, process).getAverage();
                yPoints2[0][i] = loanCC2.getValue(0, multiModel, process).getAverage();
                yPoints2[1][i] = loanCD2.getValue(0, multiModel, process).getAverage();
                yPoints2[2][i] = loanDD2.getValue(0, multiModel, process).getAverage();

            }

            String[] plotableNames = new String[] {"Valuation of Creditor","Valuation of Creditor with Stopping as by Debtor","Valuation of Debtor"};
            // Plot result
            plotArray[0] = new EasyPlot2D(plotableNames, xPoints[0], yPoints0);
            plotArray[0].setTitle(plotNames[0]);
            plotArray[0].setIsLegendVisible(true);
            plotArray[0].setXAxisLabel("Free Parameter range");
            plotArray[0].setYAxisLabel("Value");
            plotArray[0].show();
            plotArray[0].writeRecreationJavaFile();

            plotArray[1] = new EasyPlot2D(plotableNames, xPoints[0], yPoints1);
            plotArray[1].setTitle(plotNames[1]);
            plotArray[1].setIsLegendVisible(true);
            plotArray[1].setXAxisLabel("Free Parameter range");
            plotArray[1].setYAxisLabel("Value");
            plotArray[1].show();
            plotArray[1].writeRecreationJavaFile();

            plotArray[2] = new EasyPlot2D(plotableNames, xPoints[0], yPoints2);
            plotArray[2].setTitle(plotNames[2]);
            plotArray[2].setIsLegendVisible(true);
            plotArray[2].setXAxisLabel("Free Parameter range");
            plotArray[2].setYAxisLabel("Value");
            plotArray[2].show();
            plotArray[2].writeRecreationJavaFile();
        }


        // Plot by Spread
        {
            final double[][] yPoints0 = new double[3][numberOfXPoints];
            final double[][] yPoints1 = new double[3][numberOfXPoints];
            final double[][] yPoints2 = new double[3][numberOfXPoints];
            final double liborPeriodLength = base.getLiborPeriod(1);
            TimeDiscretization loanTenor = new TimeDiscretizationFromArray(liborPeriodLength, 14, liborPeriodLength);

            final double startSpread = 0.01;
            final double deltaSpread = 0.002;
            for (int i = 0; i < numberOfXPoints; i++) {
                xPoints[1][i] = startSpread + i * deltaSpread;
                double[] initRatesCredF = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
                double[] initRatesDebtF = Arrays.copyOf(initRatesNonDef, initRatesNonDef.length);
                for (int j = 0; j < initRatesCredF.length; j++) {
                    initRatesCredF[j] += xPoints[1][i];
                    initRatesDebtF[j] += xPoints[1][i];
                }

                properties.put("numberOfExtraFactors", 2);
                properties.put("initialRatesDefaultable", initRatesCredF);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 33);
                properties.put("freeParamsRange", 0.5);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelCreditorF = factory.createDefaultableModel(base);

                properties.put("numberOfExtraFactors", 3);
                properties.put("initialRatesDefaultable", initRatesDebtF);
                properties.put("simulationModel", "SPREADS");
                properties.put("stateSpaceOfSpread", "LOGNORMAL");
                properties.put("freeParamsSeed", 1030);
                properties.put("freeParamsRange", 0.7);
                factory.setProperties(properties);

                DefaultableLIBORMarketModel modelDebtorF = factory.createDefaultableModel(base);
                MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditorF, modelDebtorF, modelCreditor, modelDebtor}, base);

                MonteCarloProcess process = factory.createNumericalScheme(multiModel);

                final double[] coupons = new double[15];
                final double nominal = 100;
                for(int j=0; j < 15; j++) {
                    coupons[j] = multiModel.getDefaultableLIBOR(process, 0, j, 1).doubleValue();
                    coupons[j] *= nominal * liborPeriodLength;
                }
                coupons[14] += nominal;
                CancellableLoan loanCC0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD0 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 2, 1, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);

                CancellableLoan loanCC1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD1 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 3, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);

                CancellableLoan loanCC2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.CREDITOR);
                CancellableLoan loanCD2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.CREDITOR, LoanProduct.Perspective.DEBTOR);
                CancellableLoan loanDD2 = new CancellableLoan(loanTenor, coupons, 4.5, nominal, 0, 1, LoanProduct.Perspective.DEBTOR, LoanProduct.Perspective.DEBTOR);


                yPoints0[0][i] = loanCC0.getValue(0, multiModel, process).getAverage();
                yPoints0[1][i] = loanCD0.getValue(0, multiModel, process).getAverage();
                yPoints0[2][i] = loanDD0.getValue(0, multiModel, process).getAverage();
                yPoints1[0][i] = loanCC1.getValue(0, multiModel, process).getAverage();
                yPoints1[1][i] = loanCD1.getValue(0, multiModel, process).getAverage();
                yPoints1[2][i] = loanDD1.getValue(0, multiModel, process).getAverage();
                yPoints2[0][i] = loanCC2.getValue(0, multiModel, process).getAverage();
                yPoints2[1][i] = loanCD2.getValue(0, multiModel, process).getAverage();
                yPoints2[2][i] = loanDD2.getValue(0, multiModel, process).getAverage();

            }

            String[] plotableNames = new String[] {"Valuation of Creditor","Valuation of Creditor with Stopping as by Debtor","Valuation of Debtor"};
            // Plot result
            plotArray[3] = new EasyPlot2D(plotableNames, xPoints[1], yPoints0);
            plotArray[3].setTitle(plotNames[3]);
            plotArray[3].setIsLegendVisible(true);
            plotArray[3].setXAxisLabel("Initial Spread");
            plotArray[3].setYAxisLabel("Value");
            plotArray[3].show();
            plotArray[3].writeRecreationJavaFile();

            plotArray[4] = new EasyPlot2D(plotableNames, xPoints[1], yPoints1);
            plotArray[4].setTitle(plotNames[4]);
            plotArray[4].setIsLegendVisible(true);
            plotArray[4].setXAxisLabel("Initial Spread");
            plotArray[4].setYAxisLabel("Value");
            plotArray[4].show();
            plotArray[4].writeRecreationJavaFile();

            plotArray[5] = new EasyPlot2D(plotableNames, xPoints[1], yPoints2);
            plotArray[5].setTitle(plotNames[5]);
            plotArray[5].setIsLegendVisible(true);
            plotArray[5].setXAxisLabel("Initial Spread");
            plotArray[5].setYAxisLabel("Value");
            plotArray[5].show();
            plotArray[5].writeRecreationJavaFile();
        }

        return plotArray;
    }

    public static void createTimingPlotDefModel() throws CalculationException {
        /*
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

        properties.put("numberOfExtraFactors", 3);
        properties.put("initialRatesDefaultable", new double[] {0.032, 0.036, 0.039, 0.04, 0.04, 0.043});
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 1030);
        properties.put("freeParamsRange", 0.6);

        // Simulation:
        properties.put("numberOfPaths", 500);
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        */
        // Plot
        final double[] paths = new double[] {500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 7500.0, 10000.0, 12500.0, 15000.0}; //, 20000.0, 50000.0};
        final int[] pathsL = new int[] {500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000, 12500, 15000, 20000, 50000 };
        final int howManyNumbers = paths.length;
        final double[] timesDefModel = { 3230.0, 2891.0, 3660.0, 5719.0, 7560.0, 11324.0, 17111.0, 19318.0, 31185.0, 47956.0, 55509.0, 76337.0}; //, 0.0, 0.0 };
        final double[] timesNonDefModel = { 585.0, 607.0, 769.0, 1032.0, 1600.0, 2435.0, 3176.0, 4062.0, 5447.0, 7435.0, 9202.0, 14842.0}; //, 0.0, 0.0 };
        /*
        Runnable Printer = () -> {
            System.out.print("timesDefModel = { ");
            for(int i=0; i < howManyNumbers - 1; i++) {
                System.out.print(timesDefModel[i] + ", ");
            }
            System.out.print(timesDefModel[howManyNumbers - 1] + " };\n");

            System.out.print("timesNonDefModel = { ");
            for(int i=0; i < howManyNumbers - 1; i++) {
                System.out.print(timesNonDefModel[i] + ", ");
            }
            System.out.print(timesNonDefModel[howManyNumbers - 1] + " };\n");
        };
        DefaultableLIBORModelFactory factory = new DefaultableLIBORModelFactory();
        factory.setProperties(properties);
        LIBORMarketModel baseModel= factory.createBaseModel();
        DefaultableLIBORMarketModel model = factory.createDefaultableModel(baseModel);
        for(int j=0; j < howManyNumbers; j++) {
            // Create defaultable model
            properties.put("numberOfPaths", pathsL[j]);
            factory.setProperties(properties);
            try {
                tic();
                MonteCarloProcess process = factory.createNumericalScheme(model);
                model.getLIBOR(process, 9, 10);
                timesDefModel[j] = (double)(toc() / 1000000L);
            } catch(Exception ex) {
                ex.printStackTrace();
                System.out.println("\n\nExited at j=" + j + " checking defaulable model!");
                Printer.run();
            }
            try {
                tic();
                MonteCarloProcess process = factory.createNumericalScheme(baseModel);
                baseModel.getLIBOR(process, 9, 10);
                timesNonDefModel[j] = (double)(toc() / 1000000L);
            } catch(Exception ex) {
                ex.printStackTrace();
                System.out.println("\n\nExited at j=" + j + " checking base model!");
                Printer.run();
            }
        }
        Printer.run();*/
        GraphStyle gsDef = new GraphStyle(new Rectangle(-3, -3, 6, 6), new BasicStroke(2.0f), EasyPlot2D.getDefaultColor(0));
        GraphStyle gsNonDef = new GraphStyle(new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0), new BasicStroke(2.0f), EasyPlot2D.getDefaultColor(1));
        Plotable2D plotableDef = PlotablePoints2D.of("Defaultable Model", paths, timesDefModel, gsDef);
        Plotable2D plotableNonDef = PlotablePoints2D.of("Non Defaultable Model", paths, timesNonDefModel, gsNonDef);
        EasyPlot2D plot = new EasyPlot2D(plotableDef);
        plot.addPlot(plotableNonDef);
        plot.setYAxisLabel("Time in Millis");
        plot.setXAxisLabel("Number of Paths");
        plot.setTitle("Comparison of Calculation Time");
        plot.setIsLegendVisible(true);
        plot.show();
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

    public static EasyPlot2D[] createCreditorAndPayerDefaultablePlotByCouponRate() throws CalculationException {
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
        properties.put("initialRatesDefaultable", new double[] {0.03, 0.03, 0.031, 0.033, 0.034, 0.036});
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 33);
        properties.put("freeParamsRange", 0.4);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelCreditor = factory.createDefaultableModel(base);


        // Def Model 2: Debtor (lower Covariance and lower Spread):
        properties.put("numberOfExtraFactors", 3);
        properties.put("initialRatesDefaultable", new double[] {0.04, 0.043, 0.042, 0.044, 0.045, 0.043});
        properties.put("simulationModel", "SPREADS");
        properties.put("stateSpaceOfSpread", "LOGNORMAL");
        properties.put("freeParamsSeed", 1030);
        properties.put("freeParamsRange", 0.7);
        factory.setProperties(properties);

        DefaultableLIBORMarketModel modelDebtor = factory.createDefaultableModel(base);


        MultiLIBORVectorModel multiModel = new MultiLIBORVectorModel(new DefaultableLIBORMarketModel[]{modelCreditor, modelDebtor}, base);

        // Simulation:
        properties.put("numberOfPaths", 5000);
        properties.put("brownianMotionSeed", 30312);
        properties.put("numericalScheme", DefaultableLIBORModelFactory.Scheme.EULER_FUNCTIONAL);
        factory.setProperties(properties);

        MonteCarloProcess process = factory.createNumericalScheme(multiModel);

        double nominal = 1;
        final int numberOfValuations = 10;
        final int fMaturityIndex = 9; // i + maturityIndex;
        final int tenorCount = 6;
        final double[] tenor = new double[tenorCount];
        for (int i = 0; i < tenor.length; i++) {
            tenor[i] = base.getLiborPeriod(fMaturityIndex + i);
        }
        double[] valuesCouponForwardBothDCredPersp = new double[numberOfValuations];
        double[] valuesCouponForwardBothDDebtPersp = new double[numberOfValuations];
        double[] valuesCouponForwardDebtD = new double[numberOfValuations];
        double[] valuesCouponForwardNoD = new double[numberOfValuations];
        double[] xPoints = new double[numberOfValuations];

        final double survivalDebtor = modelDebtor.getSurvivalProbability(multiModel.getDefaultableProcess(process, 1), multiModel.getLiborPeriod(fMaturityIndex)).getAverage();
        final double survivalCreditor = modelCreditor.getSurvivalProbability(multiModel.getDefaultableProcess(process, 0), multiModel.getLiborPeriod(fMaturityIndex)).getAverage();
        System.out.printf("Survival Probability Debtor: %5.1f %% %n", survivalDebtor * 100);
        System.out.printf("Survival Probability Creditor: %5.1f %% %n", survivalCreditor * 100);
        System.out.printf("%10s      %10s     %10s     %10s     %10s     %10s%n",
                "Coupons", "C, D def.", "D def.", "C coll.", "D,C coll.", "Non def.");
        final double ciDelta = 0.0002;
        double ci = 0.043;
        for(int i = 0; i < numberOfValuations; i++) {
            double[] couponRates = new double[tenorCount - 1];
            Arrays.fill(couponRates, ci);
            DefaultableCouponBondForward forwardBDCP = new DefaultableCouponBondForward(tenor, couponRates, nominal, nominal, 1, 0, LoanProduct.Perspective.CREDITOR);
            DefaultableCouponBondForward forwardBDDP = new DefaultableCouponBondForward(tenor, couponRates, nominal, nominal, 1, 0, LoanProduct.Perspective.DEBTOR);
            DefaultableCouponBondForward forwardDDDP = new DefaultableCouponBondForward(tenor, couponRates, nominal, nominal, -1, 0, LoanProduct.Perspective.DEBTOR);
            DefaultableCouponBondForward forwardNODP = new DefaultableCouponBondForward(tenor, couponRates, nominal, nominal, -1, -1, LoanProduct.Perspective.DEBTOR);

            valuesCouponForwardBothDCredPersp[i] = forwardBDCP.getValue(0, multiModel, process).getAverage();
            valuesCouponForwardBothDDebtPersp[i] = forwardBDDP.getValue(0, multiModel, process).getAverage();
            valuesCouponForwardDebtD[i] = forwardDDDP.getValue(0, multiModel, process).getAverage();
            valuesCouponForwardNoD[i] = forwardNODP.getValue(0, multiModel, process).getAverage();
            xPoints[i] = ci * 100;
            ci += ciDelta;

            System.out.printf("c=%6.4f:     %10.5f     %10.5f     %10.5f     %10.5f%n",
                    xPoints[i],
                    valuesCouponForwardBothDCredPersp[i], valuesCouponForwardBothDDebtPersp[i], valuesCouponForwardDebtD[i], valuesCouponForwardNoD[i]);
        }
        double[][] yPoints = new double[][]{valuesCouponForwardBothDCredPersp, valuesCouponForwardBothDDebtPersp, valuesCouponForwardDebtD, valuesCouponForwardNoD};
        String[] names = new String[] {"Both defaultable, Creditors Perspective", "Both defaultable, Debtors Perspective", "Debtor defaultable", "No default probability"};

        EasyPlot2D plot = new EasyPlot2D(names, xPoints, yPoints);
        plot.setTitle("Valuation of Coupon Bond Forwards");
        plot.setXAxisLabel("Coupons in %");
        plot.setYAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plot.setXAxisNumberFormat(new DecimalFormat("00.00;-00.00", new DecimalFormatSymbols(Locale.ENGLISH)));
        plot.setYAxisLabel("Price");
        plot.setIsLegendVisible(true);
        plot.show();
        return new EasyPlot2D[] {plot};
    }

    /*
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
        BrownianMotion motion2 = (BrownianMotion) process.getStochasticDriver();
        Integer[] factors = new Integer[modelDebtor.getNumberOfFactors()];
        for(int i=0; i < base.getNumberOfFactors(); i++) factors[i] = i;
        for(int i=base.getNumberOfFactors(); i < modelDebtor.getNumberOfFactors(); i++)
            factors[i] = i + multiModel.getDefaultableModel(0).getNumberOfFactors() - base.getNumberOfFactors();
        motion2 = new BrownianMotionView((BrownianMotion) process.getStochasticDriver(), factors);
        MonteCarloProcess secProcess = new EulerSchemeFromProcessModel(modelDebtor, motion2, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);

        double nominal = 1000000;
        int maturityIndex = 3;
        final int numberOfValuations = 13;
        double[] valuesCouponForwardBothD = new double[numberOfValuations];
        double[] valuesCouponForwardDebtD = new double[numberOfValuations];
        double[] valuesCouponForwardBothDCredC = new double[numberOfValuations];
        double[] valuesCouponForwardBothDBothC = new double[numberOfValuations];
        double[] valuesCouponForwardNoD = new double[numberOfValuations];
        double[] xPoints = new double[numberOfValuations];
        System.out.printf("%10s      %10s     %10s     %10s     %10s     %10s%n",
                "Maturity", "C, D def.", "D def.", "C coll.", "D,C coll.", "Non def.");
        final double ciDelta = 0.002;
        double ci = 0.03;
        for(int i = 0; i < numberOfValuations; i++) {
            double[] couponRates = new double[5];
            Arrays.fill(couponRates, ci * 0.5);
            // ci += ciDelta;
            final int fMaturityIndex = i + maturityIndex;
            DefaultableCouponBondForward forward = new DefaultableCouponBondForward(nominal, nominal, couponRates, i + maturityIndex);
            valuesCouponForwardBothD[i] = forward.getValue(0, multiModel, process, 1, 0).getAverage();
            valuesCouponForwardDebtD[i] = forward.getValue(0, modelDebtor, multiModel.getDefaultableProcess(process, 1)).getAverage();
            valuesCouponForwardNoD[i] = forward.getValue(0, base, multiModel.getNonDefaultableProcess(process)).getAverage();
            DefaultableCouponBondForward multiForwardCC = new DefaultableCouponBondForward(nominal, nominal, couponRates, i + maturityIndex, true, false);
            DefaultableCouponBondForward multiForwardBC = new DefaultableCouponBondForward(nominal, nominal, couponRates, i + maturityIndex, true);
            valuesCouponForwardBothDCredC[i] = multiForwardCC.getDoubleValue(0, multiModel, process, 1, 0);
            valuesCouponForwardBothDBothC[i] = multiForwardBC.getDoubleValue(0, multiModel, process, 1, 0);
            xPoints[i] = multiModel.getLiborPeriod(i+maturityIndex);

            System.out.printf("%2d (T_S=%3.1f):     %10.5f     %10.5f     %10.5f     %10.5f     %10.5f%n",
                    i + maturityIndex, xPoints[i],
                    valuesCouponForwardBothD[i], valuesCouponForwardDebtD[i], valuesCouponForwardBothDCredC[i], valuesCouponForwardBothDBothC[i], valuesCouponForwardNoD[i]);
        }

        EasyPlot2D plot = new EasyPlot2D("Creditor and Debtor defaultable", xPoints, valuesCouponForwardBothD);
        plot.addPlot("Debtor defaultable", xPoints, valuesCouponForwardDebtD);
        //plot.addPlot("Pure Cash-Settled, Both defaultable", xPoints, valuesCouponForwardBothDBothC);
        plot.addPlot("Cash-Settled, debtor does not collect, Both defaultable", xPoints, valuesCouponForwardBothDCredC);
        plot.changePlotColor(1, EasyPlot2D.getDefaultColor(1));
        plot.changePlotColor(2, EasyPlot2D.getDefaultColor(2));
        // plot.changePlotColor(1, EasyPlot2D.getDefaultColor(3));
        // plot.addPlot("Non defaultable valuation", xPoints, valuesCouponForwardNoD);

        plot.setTitle("Valuation of Coupon Bond Forwards");
        plot.setXAxisLabel("Time to maturity");
        plot.setYAxisLabel("Price");
        plot.setIsLegendVisible(true);
        plot.show();
    }*/


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


    public static EasyPlot2D[] createSurvivalProbabilityPlot() throws CalculationException {
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
        return new EasyPlot2D[] {plotDefaultableBonds};
    }

    /**
     * Finds the subscript of the specified number as unicode literal.
     * See <a href="https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts">Wikipedia</a>
     * for more subscripts.
     * @param number the number (0-9) to apply as subscript.
     * @return the literal for the subscript version of the specified number.
     */
    public static char subscript(int number) {
        return (char)('\u2080' + number);
    }

    /**
     * Finds the superscript of the specified number as unicode literal.
     * See <a href="https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts">Wikipedia</a>
     * for more superscripts.
     * @param number the number (0-9) to apply as superscript.
     * @return the literal for the superscript version of the specified number.
     */
    public static char superscript(int number) {
        if(number > 3 || number == 0)
            return (char)('\u2070' + number);
        else if(number > 1)
            return (char)('\u00b0' + number);
        else
            return '\u00b9';
    }

}
