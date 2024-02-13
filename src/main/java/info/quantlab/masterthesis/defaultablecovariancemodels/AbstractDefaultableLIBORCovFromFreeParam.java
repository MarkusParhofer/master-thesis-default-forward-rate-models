package info.quantlab.masterthesis.defaultablecovariancemodels;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

import java.util.Arrays;
import java.util.Map;

public abstract class AbstractDefaultableLIBORCovFromFreeParam extends AbstractDefaultableLIBORCovariance {
    public AbstractDefaultableLIBORCovFromFreeParam(LIBORCovarianceModel nonDefaultableCovarianceModel, int extraNumberOfFactors) {
        super(nonDefaultableCovarianceModel, nonDefaultableCovarianceModel.getTimeDiscretization(),
                nonDefaultableCovarianceModel.getLiborPeriodDiscretization(),
                nonDefaultableCovarianceModel.getNumberOfFactors() + extraNumberOfFactors);
    }

    public RandomVariable[] getFreeParameter(int timeIndex, int componentIndex, final RandomVariable defRealization, final RandomVariable nonDefRealization) {
        final RandomVariable[] result = new RandomVariable[getNumberOfFactors() - getNonDefaultableCovarianceModel().getNumberOfFactors()];
        for(int i=0; i < result.length; i++) {
            result[i] = getFreeParameter(timeIndex, componentIndex, i, defRealization, nonDefRealization);
        }
        return result;
    }

    public abstract RandomVariable getFreeParameter(int timeIndex, int componentIndex, int factor, final RandomVariable defRealization, final RandomVariable nonDefRealization);

    @Override
    public abstract double[] getParameterAsDouble();

    @Override
    public RandomVariable getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariable[] realizationAtTimeIndex) {
        throw new UnsupportedOperationException();
    }

    @Override
    public abstract AbstractDefaultableLIBORCovFromFreeParam clone();

    @Override
    public abstract AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedParameters(double[] parameters);

    @Override
    public abstract AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException;

    @Override
    public abstract AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNonDefaultableCovarianceModel);

    @Override
    public abstract int getNumberOfParameters();

    @Override
    public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex, RandomVariable[] undefaultableRealization) {

        RandomVariable[] undefaultableFactorLoading = getNonDefaultableFactorLoading(timeIndex, component, undefaultableRealization);

        if(component < getNumberOfLIBORPeriods()) {
            return undefaultableFactorLoading;
        }
        return getDefaultableFactorLoading(timeIndex, component, undefaultableFactorLoading, realizationAtTimeIndex[getLIBORIndexFromComponent(component)], undefaultableRealization[getLIBORIndexFromComponent(component)]);
    }

    @Override
    public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] realizationAtTimeIndex) {
        final int liborPeriodIndex = getLIBORIndexFromComponent(component);
        RandomVariable[] undefaultableFactorLoading = getNonDefaultableFactorLoading(timeIndex, liborPeriodIndex, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));

        if(component < getNumberOfLIBORPeriods()) {
            return undefaultableFactorLoading;
        }
        return getDefaultableFactorLoading(timeIndex, liborPeriodIndex, undefaultableFactorLoading, realizationAtTimeIndex[component], realizationAtTimeIndex[getLIBORIndexFromComponent(component)]);

    }

    private RandomVariable[] getDefaultableFactorLoading(int timeIndex, int liborPeriodIndex, RandomVariable[] nondefaultableFactorLoading, RandomVariable componentRealizationDefaultable, RandomVariable componentRealizationUndefaultable) {
        if(getTimeDiscretization().getTime(timeIndex) >= getLiborPeriodDiscretization().getTime(liborPeriodIndex))
            return getZeroFactorLoading();

        final int undefaultableFactors = getNonDefaultableCovarianceModel().getNumberOfFactors();
        final double periodLength = getLiborPeriodDiscretization().getTimeStep(liborPeriodIndex);

        RandomVariable[] factorLoading = new RandomVariable[getNumberOfFactors()];

        // Calculate from underlying undefaultable Model
        RandomVariable spread = componentRealizationDefaultable.sub(componentRealizationUndefaultable);

        RandomVariable relationFactor = Scalar.of(1.0);
        if(spread.getMin() != spread.getMax() || spread.getMax() != 0.0)
            relationFactor = relationFactor.accrue(componentRealizationDefaultable,periodLength).discount(componentRealizationUndefaultable, periodLength);

        for(int k = 0; k < undefaultableFactors; k++) {
            factorLoading[k] = relationFactor.mult(nondefaultableFactorLoading[k]);
        }

        // Calculate from free parameters
        RandomVariable[] freeParams = getFreeParameter(timeIndex, liborPeriodIndex, componentRealizationDefaultable, componentRealizationUndefaultable);
        for(int k = undefaultableFactors; k < getNumberOfFactors(); k++) {
            factorLoading[k] = spread.mult(freeParams[k - undefaultableFactors]);
        }
        return factorLoading;
    }

    private RandomVariable[] getNonDefaultableFactorLoading(int timeIndex, int liborPeriodIndex, RandomVariable[] undefaultableRealization) {
        // Return undefaultable Factor Loadings with higher number of Factors. Set extra factors to zero.

        final RandomVariable[] undefaultableFactorLoading = getNonDefaultableCovarianceModel().getFactorLoading(timeIndex, liborPeriodIndex, undefaultableRealization);

        final RandomVariable zero = new Scalar(0.0);

        RandomVariable[] result = Arrays.copyOf(undefaultableFactorLoading, getNumberOfFactors());

        Arrays.fill(result, getNonDefaultableCovarianceModel().getNumberOfFactors(), getNumberOfFactors(), zero);

        return result;
    }

    @Override
    public RandomVariable[] getFactorLoadingOfSpread(int timeIndex, int liborPeriodIndex, RandomVariable[] realizationAtTimeIndex) {
        if(liborPeriodIndex >= getNumberOfLIBORPeriods())
            throw new ArrayIndexOutOfBoundsException("Spread model is a model of " + getNumberOfLIBORPeriods() + " Components. Index " + liborPeriodIndex + " out of Bounds for Spread Model");

        if(getTimeDiscretization().getTime(timeIndex) >= getLiborPeriodDiscretization().getTime(liborPeriodIndex))
            return getZeroFactorLoading();
        final RandomVariable[] allFactorLoadings = getNonDefaultableFactorLoading(timeIndex, liborPeriodIndex, Arrays.copyOf(realizationAtTimeIndex, getNumberOfLIBORPeriods()));

        final int nonDefNumberOfFactors = getNonDefaultableCovarianceModel().getNumberOfFactors();
        for(int k = 0; k < nonDefNumberOfFactors; k++) {
            final double deltaT = getLiborPeriodDiscretization().getTimeStep(liborPeriodIndex);
            allFactorLoadings[k] = allFactorLoadings[k].mult(deltaT).div(realizationAtTimeIndex[liborPeriodIndex].mult(deltaT).add(1.0));
        }


        final RandomVariable defRealization = realizationAtTimeIndex[getNumberOfLIBORPeriods() + liborPeriodIndex].add(realizationAtTimeIndex[liborPeriodIndex]);

        RandomVariable[] freeParams = getFreeParameter(timeIndex, liborPeriodIndex, defRealization, realizationAtTimeIndex[liborPeriodIndex]);

        System.arraycopy(freeParams, 0, allFactorLoadings, nonDefNumberOfFactors, getNumberOfFactors() - nonDefNumberOfFactors);

        return allFactorLoadings;
    }

    @Override
    public boolean isSpreadModelLogNormal() {
        return true;
    }

    private int getLIBORIndexFromComponent(int component) {
        return component % getNumberOfLIBORPeriods();
    }
}
