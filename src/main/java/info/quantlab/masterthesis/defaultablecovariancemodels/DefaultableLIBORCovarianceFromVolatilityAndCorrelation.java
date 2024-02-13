package info.quantlab.masterthesis.defaultablecovariancemodels;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

import java.util.Arrays;
import java.util.Map;

/**
 * This class represents a covariance structure for defaultable LIBOR models, that generates a log-normal spread dynamic.
 * As input, it takes a volatility and correlation structure, that functions as <b>addons</b> to the
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
 * where &sigma;&#771; and &rho;&#771; represent the inputted volatility and correlation models respectively.
 *
 * @author Markus Parhofer
 * @version 1.0
 */
@SuppressWarnings("unused")
public class DefaultableLIBORCovarianceFromVolatilityAndCorrelation extends AbstractDefaultableLIBORCovFromFreeParam {
    LIBORVolatilityModel _vola;
    LIBORCorrelationModel _corr;
    int _nonDefParams, _volaParams, _corrParams;

    public DefaultableLIBORCovarianceFromVolatilityAndCorrelation(LIBORVolatilityModel volatilityModel, LIBORCorrelationModel correlationModel, LIBORCovarianceModel nonDefaultableCovarianceModel) {
        this(volatilityModel, correlationModel, nonDefaultableCovarianceModel, nonDefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric);
    }

    public DefaultableLIBORCovarianceFromVolatilityAndCorrelation(LIBORVolatilityModel volatilityModel, LIBORCorrelationModel correlationModel, LIBORCovarianceModel nonDefaultableCovarianceModel, boolean nonDefaultableModelIsCalibrateable) {
        this(volatilityModel, correlationModel, nonDefaultableCovarianceModel, 0, 0, 0);
        _volaParams = volatilityModel.getParameter().length;
        _corrParams = correlationModel.getParameter().length;
        if (nonDefaultableModelIsCalibrateable)
            _nonDefParams = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getParameter().length;
    }

    private DefaultableLIBORCovarianceFromVolatilityAndCorrelation(LIBORVolatilityModel volatilityModel, LIBORCorrelationModel correlationModel, LIBORCovarianceModel nonDefCovModel, int nonDefParams, int volaParams, int corrParams) {
        super(nonDefCovModel, correlationModel.getNumberOfFactors());
        _vola = volatilityModel;
        _corr = correlationModel;
        _nonDefParams = nonDefParams;
        _volaParams = volaParams;
        _corrParams = corrParams;
    }

    @Override
    public RandomVariable getFreeParameter(int timeIndex, int componentIndex, int factor, RandomVariable defRealization, RandomVariable nonDefRealization) {
        return _vola.getVolatility(timeIndex, componentIndex).mult(_corr.getFactorLoading(timeIndex, factor, componentIndex)).div(defRealization.sub(nonDefRealization));
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
        double[] volaParam = _vola.getParameterAsDouble();
        int otherIndex = 0;
        for (; otherIndex < volaParam.length; otherIndex++, index++) {
            newParams[index] = volaParam[otherIndex];
        }
        double[] corrParam = _corr.getParameterAsDouble();
        otherIndex = 0;
        for (; otherIndex < volaParam.length; otherIndex++, index++) {
            newParams[index] = corrParam[otherIndex];
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
        RandomVariable[] volaParam = _vola.getParameter();
        int otherIndex = 0;
        for (; otherIndex < volaParam.length; otherIndex++, index++) {
            newParams[index] = volaParam[otherIndex];
        }
        RandomVariable[] corrParam = _corr.getParameter();
        otherIndex = 0;
        for (; otherIndex < volaParam.length; otherIndex++, index++) {
            newParams[index] = corrParam[otherIndex];
        }
        return newParams;
    }

    @Override
    public AbstractDefaultableLIBORCovFromFreeParam clone() {
        return new DefaultableLIBORCovarianceFromVolatilityAndCorrelation(_vola, _corr, getNonDefaultableCovarianceModel(), _nonDefParams, _volaParams, _corrParams);
    }

    @Override
    public AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedParameters(double[] parameters) {
        RandomVariable[] paramAsRV = new RandomVariable[parameters.length];
        Arrays.setAll(paramAsRV, i -> Scalar.of(parameters[i]));
        return getCloneWithModifiedParameters(paramAsRV);
    }

    @Override
    public AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedParameters(RandomVariable[] parameters) {
        LIBORCovarianceModel newNonDefCovModel = getNonDefaultableCovarianceModel();
        if(_nonDefParams > 0)
            newNonDefCovModel = ((AbstractLIBORCovarianceModelParametric)getNonDefaultableCovarianceModel()).getCloneWithModifiedParameters(Arrays.copyOfRange(parameters, 0, _nonDefParams));
        LIBORVolatilityModel newVolaModel = _vola.getCloneWithModifiedParameter(Arrays.copyOfRange(parameters, _nonDefParams, _nonDefParams + _volaParams));
        LIBORCorrelationModel newCorrModel = _corr.getCloneWithModifiedParameter(Arrays.copyOfRange(parameters, _nonDefParams + _volaParams, _nonDefParams + _volaParams + _corrParams));
        return new DefaultableLIBORCovarianceFromVolatilityAndCorrelation(newVolaModel, newCorrModel, newNonDefCovModel, _nonDefParams, _volaParams, _corrParams);
    }

    @Override
    public AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
        LIBORVolatilityModel volatilityModel = _vola;
        LIBORCorrelationModel correlationModel = _corr;
        LIBORCovarianceModel nonDefModel = getNonDefaultableCovarianceModel();
        boolean nonDefModelIsCalibrateable = _nonDefParams > 0;
        for (String key: dataModified.keySet()) {
            switch(key.toUpperCase()) {
                case "NONDEFAULTABLEMODEL":
                    nonDefModel = (LIBORCovarianceModel)dataModified.get(key);
                    break;
                case "VOLATILITYMODEL":
                    volatilityModel = (LIBORVolatilityModel)dataModified.get(key);
                    break;
                case "CORRELATIONMODEL":
                    correlationModel = (LIBORCorrelationModel)dataModified.get(key);
                    break;
                case "NODEFAULTABLEMODELISCALIBRATEABLE":
                    nonDefModelIsCalibrateable = (Boolean)dataModified.get(key);
                    break;
            }
        }
        return new DefaultableLIBORCovarianceFromVolatilityAndCorrelation(volatilityModel, correlationModel, nonDefModel, nonDefModelIsCalibrateable);
    }

    @Override
    public AbstractDefaultableLIBORCovFromFreeParam getCloneWithModifiedNonDefaultableCovariance(LIBORCovarianceModel newNonDefaultableCovarianceModel) {
        int nonDefParams = 0;
        if(newNonDefaultableCovarianceModel instanceof AbstractLIBORCovarianceModelParametric newModel) {
            nonDefParams = newModel.getParameter().length;
        }
        return new DefaultableLIBORCovarianceFromVolatilityAndCorrelation(_vola, _corr, newNonDefaultableCovarianceModel, nonDefParams, _volaParams, _corrParams);
    }

    @Override
    public int getNumberOfParameters() {
        return _nonDefParams + _volaParams + _corrParams;
    }
}
