package info.quantlab.masterthesis.multilibormodels;


import net.finmath.functions.LinearAlgebra;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * Calculates a free Parameter Matrix such that the initial covariance Structure is the same as that of the defaultable model.
 * @author Markus Parhofer
 * @version 1.0
 */
public class DefaultableLIBORCovarianceWithInitialUndefaultableCovariance extends DefaultableLIBORCovarianceWithGuaranteedPositiveSpread {

	/**
	 * Default Version ID
	 */
	private static final long serialVersionUID = 1L;


	/**
	 * Creates a DefaultableLIBORCovarianceModel of type DefaultableLIBORCovarianceWithGuaranteedPositiveSpread. The free parameter matrix is 
	 * calculated to generate the same initial covariance as the underlying Non Defaultable model.
	 * @param undefaultableCovarianceModel The underlying non defaultable covariance Model.
	 * @param initialRatesDefaultable The initial rates of the defaultable model associated with the libor period discretization.
	 * @param initialRatesNonDefaultable The initial rates of the underlying non defaultable model associated with the libor period discretization.
	 * @param numberOfExtraFactors The number of extra factors to use (i.e the number ofcolumns of the free parameter matrix).
	 */
	public DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(LIBORCovarianceModel undefaultableCovarianceModel, double[] initialRatesDefaultable, double[] initialRatesNonDefaultable, int numberOfExtraFactors) {
		super(undefaultableCovarianceModel, calculateFreeParameterMatrix(undefaultableCovarianceModel, initialRatesNonDefaultable, initialRatesDefaultable, numberOfExtraFactors));
	}
	
	/**
	 * Creates a DefaultableLIBORCovarianceModel of type DefaultableLIBORCovarianceWithGuaranteedPositiveSpread. The free parameter matrix is 
	 * calculated to generate the same initial covariance as the underlying Non Defaultable model.
	 * @param undefaultableCovarianceModel The underlying non defaultable covariance Model.
	 * @param forwardsDefaultable 
	 * @param forwardsNonDefaultable
	 * @param numberOfExtraFactors The number of extra factors to use (i.e the number ofcolumns of the free parameter matrix).
	 */
	public DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(LIBORCovarianceModel undefaultableCovarianceModel, ForwardCurve forwardsDefaultable, ForwardCurve forwardsNonDefaultable, int numberOfExtraFactors) {
		this(undefaultableCovarianceModel, fromCurveToArray(forwardsDefaultable, undefaultableCovarianceModel.getLiborPeriodDiscretization()), fromCurveToArray(forwardsNonDefaultable, undefaultableCovarianceModel.getLiborPeriodDiscretization()), numberOfExtraFactors);
	}
	

	
	private static double[][] calculateFreeParameterMatrix(LIBORCovarianceModel undefaultableModel, double[] initialRatesDefaultable, double[] initialRatesNonDefaultable, int numberOfExtraFactors) {
		final int numberOfLIBORs = undefaultableModel.getLiborPeriodDiscretization().getNumberOfTimeSteps();
		double[][] restCovariance = new double[numberOfLIBORs][numberOfLIBORs];
		for(int i=0; i < numberOfLIBORs; i++) {
			for(int j=0; j < numberOfLIBORs; j++) {
				restCovariance[i][j] = undefaultableModel.getCovariance(0, i, j, null).doubleValue();
				RandomVariable[] factorLoadingI = undefaultableModel.getFactorLoading(0.0, i, null);
				RandomVariable[] factorLoadingJ = undefaultableModel.getFactorLoading(0.0, j, null);
				for(int factor=0; factor < undefaultableModel.getNumberOfFactors(); factor++) {
					restCovariance[i][j] -= factorLoadingI[factor].doubleValue() * initialRatesDefaultable[i] / initialRatesNonDefaultable[i] * factorLoadingJ[factor].doubleValue() * initialRatesDefaultable[j] / initialRatesNonDefaultable[j];
				}
				restCovariance[j][i] = restCovariance[i][j];
			}
		}
		double[][] freeParamerters = LinearAlgebra.factorReduction(restCovariance, numberOfExtraFactors);
		for(int i=0; i < freeParamerters.length; i++) {
			for(int j=0; j < freeParamerters[0].length; j++) {
				freeParamerters[i][j] /= initialRatesDefaultable[i] - initialRatesNonDefaultable[i];
			}
		}
		return freeParamerters;
	}

	private static double[] fromCurveToArray(ForwardCurve curve, AnalyticModel model, TimeDiscretization liborPeriods) {
		double[] forwards = new double[liborPeriods.getNumberOfTimeSteps()];
		for(int i=0; i < forwards.length; i++) {
			forwards[i] = curve.getForward(model, liborPeriods.getTime(i));
		}
		return forwards;
	}

	private static double[] fromCurveToArray(ForwardCurve curve, TimeDiscretization liborPeriods) {
		return fromCurveToArray(curve, null, liborPeriods);
	}
}
