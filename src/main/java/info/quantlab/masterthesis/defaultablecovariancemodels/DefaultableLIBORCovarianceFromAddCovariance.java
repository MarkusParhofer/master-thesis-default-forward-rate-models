package info.quantlab.masterthesis.defaultablecovariancemodels;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

import java.util.Arrays;
import java.util.Map;
/**
 * This class represents a covariance structure for defaultable LIBOR models, that generates a log-normal spread dynamic.
 * As input, it takes a covariance structure, that functions as <b>addon</b> to the
 * naturally derived volatility structure, which depends fully on the non defaultable model:
 * <p/>
 * <p style="text-align:center">
 * (&sigma;<sub>i</sub><sup>d</sup>)<sup>2</sup> = (Q(&tau; > T<sub>i</sub> | &tau; > T<sub>i-1</sub>)<sup>-1</sup>
 * &sigma;<sub>i</sub>))<sup>2</sup> + (&sigma;&#771;<sub>i</sub><sup>d</sup>)<sup>2</sup>
 * </p>
 * and
 * <p style="text-align:center">
 * &sigma;<sub>i</sub><sup>d</sup> &sigma;<sub>j</sub><sup>d</sup> &rho;<sub>i,j</sub><sup>d</sup>=
 * ((Q(&tau; > T<sub>i</sub> | &tau; > T<sub>i-1</sub>)Q(&tau; > T<sub>j</sub> | &tau; > T<sub>j-1</sub>))<sup>-1</sup>
 * &sigma;<sub>i</sub> &sigma;<sub>j</sub> &rho;<sub>i,j</sub> +
 * &sigma;&#771;<sub>i</sub><sup>d</sup> &sigma;&#771;<sub>j</sub><sup>d</sup> &rho;&#771;<sub>i,j</sub><sup>d</sup>
 * </p>
 * <p/>
 * where &sigma;&#771; and &rho;&#771; represent the inputted covariance model.
 *
 * @author Markus Parhofer
 * @version 1.0
 */
@SuppressWarnings("unused")
public class DefaultableLIBORCovarianceFromAddCovariance extends AbstractDefaultableLIBORCovFromFreeParam {
    LIBORCovarianceModel _covModel;
    int _nonDefParams, _covParams;

    public DefaultableLIBORCovarianceFromAddCovariance(LIBORCovarianceModel additionalCovModel, LIBORCovarianceModel nonDefaultableCovarianceModel) {
        this(additionalCovModel, nonDefaultableCovarianceModel, nonDefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric);
    }

    public DefaultableLIBORCovarianceFromAddCovariance(LIBORCovarianceModel additionalCovModel, LIBORCovarianceModel nonDefaultableCovarianceModel, boolean nonDefaultableModelIsCalibrateable) {
        this(additionalCovModel, nonDefaultableCovarianceModel, 0, 0);
        if(_covModel instanceof AbstractLIBORCovarianceModelParametric paramModel)
            _covParams = paramModel.getParameter().length;
        if (nonDefaultableModelIsCalibrateable)
            _nonDefParams = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getParameter().length;
    }

    private DefaultableLIBORCovarianceFromAddCovariance(LIBORCovarianceModel additionalCovModel, LIBORCovarianceModel nonDefCovModel, int nonDefParams, int covParams) {
        super(nonDefCovModel, additionalCovModel.getNumberOfFactors());
        _covModel = additionalCovModel;
        _nonDefParams = nonDefParams;
        _covParams = covParams;
    }

    @Override
    public RandomVariable[] getFreeParameter(int timeIndex, int liborPeriodIndex, RandomVariable defRealization, RandomVariable nonDefRealization) {
        RandomVariable[] realVector = new RandomVariable[getNumberOfLIBORPeriods()];
        realVector[liborPeriodIndex] = defRealization;
        RandomVariable[] resultVector = _covModel.getFactorLoading(timeIndex, liborPeriodIndex, realVector);
        for (int i = 0; i < realVector.length; i++) {
            resultVector[i] = resultVector[i].div(defRealization.sub(nonDefRealization));
        }
        return resultVector;
    }

    @Override
    public RandomVariable getFreeParameter(int timeIndex, int liborPeriodIndex, int factor, RandomVariable defRealization, RandomVariable nonDefRealization) {
        RandomVariable[] realVector = new RandomVariable[getNumberOfLIBORPeriods()];
        realVector[liborPeriodIndex] = defRealization;
        return _covModel.getFactorLoading(timeIndex, liborPeriodIndex, realVector)[factor].div(defRealization.sub(nonDefRealization));
    }

    @Override
    public double[] getParameterAsDouble() {
        int index = 0;
        double[] newParams = new double[getNumberOfParameters()];
        if(_nonDefParams > 0) {
            double[] nonDefParam = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getParameterAsDouble();
            for (; index < nonDefParam.length; index++) {
                newParams[index] = nonDefParam[index];
            }
        }
        if(_covParams > 0) {
            double[] covParam = ((AbstractLIBORCovarianceModelParametric)_covModel).getParameterAsDouble();
            for (double v : covParam) {
                newParams[index++] = v;
            }
        }
        return newParams;
    }

    @Override
    public RandomVariable[] getParameter() {
        int index = 0;
        RandomVariable[] newParams = new RandomVariable[getNumberOfParameters()];
        if(_nonDefParams > 0) {
            RandomVariable[] nonDefParam = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getParameter();
            for (; index < nonDefParam.length; index++) {
                newParams[index] = nonDefParam[index];
            }
        }
        if(_covParams > 0) {
            RandomVariable[] covParam = ((AbstractLIBORCovarianceModelParametric)_covModel).getParameter();
            for (RandomVariable randomVariable : covParam) {
                newParams[index++] = randomVariable;
            }
        }
        return newParams;
    }

    @Override
    public DefaultableLIBORCovarianceFromAddCovariance clone() {
        return new DefaultableLIBORCovarianceFromAddCovariance(_covModel, getNonDefaultableCovarianceModel(), _nonDefParams, _covParams);
    }

    @Override
    public DefaultableLIBORCovarianceFromAddCovariance getCloneWithModifiedParameters(double[] parameters) {
        RandomVariable[] paramAsRV = new RandomVariable[parameters.length];
        Arrays.setAll(paramAsRV, i -> Scalar.of(parameters[i]));
        return getCloneWithModifiedParameters(paramAsRV);
    }

    @Override
    public DefaultableLIBORCovarianceFromAddCovariance getCloneWithModifiedParameters(RandomVariable[] parameters) {
        LIBORCovarianceModel newNonDefCovModel = getNonDefaultableCovarianceModel();
        if(_nonDefParams > 0)
            newNonDefCovModel = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, 0, _nonDefParams));
        LIBORCovarianceModel newCovModel = _covModel;
        if(_covParams > 0)
            newCovModel = ((AbstractLIBORCovarianceModelParametric)_covModel).getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, _nonDefParams, _nonDefParams + _covParams));

        return new DefaultableLIBORCovarianceFromAddCovariance(newCovModel, newNonDefCovModel, _nonDefParams, _covParams);
    }

    @Override
    public DefaultableLIBORCovarianceFromAddCovariance getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
        LIBORCovarianceModel additionalCov = _covModel;
        LIBORCovarianceModel nonDefModel = getNonDefaultableCovarianceModel();
        boolean nonDefModelIsCalibrateable = _nonDefParams > 0;
        for (String key: dataModified.keySet()) {
            switch(key.toUpperCase()) {
                case "NONDEFAULTABLEMODEL":
                    nonDefModel = (LIBORCovarianceModel)dataModified.get(key);
                    break;
                case "ADDITIONALCOVARIANCEMODEL":
                    additionalCov = (LIBORCovarianceModel)dataModified.get(key);
                    break;
                case "NODEFAULTABLEMODELISCALIBRATEABLE":
                    nonDefModelIsCalibrateable = (Boolean)dataModified.get(key);
                    break;
            }
        }
        return new DefaultableLIBORCovarianceFromAddCovariance(additionalCov, nonDefModel, nonDefModelIsCalibrateable);
    }

    @Override
    public DefaultableLIBORCovarianceFromAddCovariance getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNonDefaultableCovarianceModel) {
        int nonDefParams = 0;
        if(newNonDefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric newModel) {
            nonDefParams = newModel.getParameter().length;
        }
        return new DefaultableLIBORCovarianceFromAddCovariance(_covModel, newNonDefaultableCovarianceModel, nonDefParams, _covParams);
    }

    @Override
    public int getNumberOfParameters() {
        return _nonDefParams + _covParams;
    }
}