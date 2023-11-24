package info.quantlab.masterthesis.defaultablecovariancemodels;


import net.finmath.functions.LinearAlgebra;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
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
	public DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(LIBORCovarianceModel undefaultableCovarianceModel, double[] initialRatesDefaultable, double[] initialRatesNonDefaultable, int numberOfExtraFactors, double maxFreeParam) {
		super(undefaultableCovarianceModel, calculateFreeParameterMatrix(undefaultableCovarianceModel, initialRatesNonDefaultable, initialRatesDefaultable, numberOfExtraFactors, maxFreeParam));
	}
	
	/**
	 * Creates a DefaultableLIBORCovarianceModel of type DefaultableLIBORCovarianceWithGuaranteedPositiveSpread. The free parameter matrix is 
	 * calculated to generate the same initial covariance as the underlying Non Defaultable model.
	 * @param undefaultableCovarianceModel The underlying non defaultable covariance Model.
	 * @param forwardsDefaultable 
	 * @param forwardsNonDefaultable
	 * @param numberOfExtraFactors The number of extra factors to use (i.e the number ofcolumns of the free parameter matrix).
	 */
	public DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(LIBORCovarianceModel undefaultableCovarianceModel, ForwardCurve forwardsDefaultable, ForwardCurve forwardsNonDefaultable, int numberOfExtraFactors, double maxFreeParam) {
		this(undefaultableCovarianceModel, fromCurveToArray(forwardsDefaultable, undefaultableCovarianceModel.getLiborPeriodDiscretization()), fromCurveToArray(forwardsNonDefaultable, undefaultableCovarianceModel.getLiborPeriodDiscretization()), numberOfExtraFactors, maxFreeParam);
	}
	

	private static double[][] calculateFreeParameterMatrix(LIBORCovarianceModel undefaultableModel, double[] initialRatesNonDefaultable, double[] initialRatesDefaultable, int numberOfExtraFactors, double maximumRange) {
		final int numberOfLIBORs = undefaultableModel.getLiborPeriodDiscretization().getNumberOfTimeSteps();
		double[][] restCorrelation = new double[numberOfLIBORs][numberOfLIBORs];
		double[] volatility = new double[numberOfLIBORs];
		double[] normFreeParamVecs = new double[numberOfLIBORs];
		for(int i=0; i < numberOfLIBORs; i++) {
			final double liborPeriodI = undefaultableModel.getLiborPeriodDiscretization().getTimeStep(i);
			final double relationFactorI = (1.0 + initialRatesDefaultable[i]*liborPeriodI) / (1.0 + initialRatesNonDefaultable[i]*liborPeriodI);
			final double spreadI = initialRatesDefaultable[i] - initialRatesNonDefaultable[i];
			volatility[i] = Math.sqrt(undefaultableModel.getCovariance(0, i, i, null).doubleValue());
			volatility[i] = volatility[i] == 0.0 ? 1.0 : volatility[i];
			
			restCorrelation[i][i] = 0.0;
			//normFreeParamVecs[i] = Double.MIN_VALUE;// (1.0 - relationFactorI * relationFactorI) / (spreadI * spreadI);
			//normFreeParamVecs[i] = Math.sqrt(normFreeParamVecs[i]);
			
			for(int j=0; j < i; j++) {
				final double liborPeriodJ = undefaultableModel.getLiborPeriodDiscretization().getTimeStep(j);
				final double relationFactorJ = (1.0 + initialRatesDefaultable[j]*liborPeriodJ) / (1.0 + initialRatesNonDefaultable[j]*liborPeriodJ);
				final double spreadJ = initialRatesDefaultable[j] - initialRatesNonDefaultable[j];
				restCorrelation[i][j] = undefaultableModel.getCovariance(0, i, j, null).doubleValue() / (volatility[i] * volatility[j]);
				restCorrelation[i][j] *= 1.0 - relationFactorI * relationFactorJ;
				restCorrelation[i][j] /= spreadI * spreadJ; // * normFreeParamVecs[i] * normFreeParamVecs[j];
				restCorrelation[j][i] = restCorrelation[i][j];
			}
		}
		double[][] freeParameters = LinearAlgebra.factorReduction(restCorrelation, numberOfExtraFactors);
		double rangeMax = 0.0;
		
		for(int i=0; i < freeParameters.length; i++) {
			for(int k=0; k < freeParameters[0].length; k++) {
				freeParameters[i][k] *= volatility[i] * (normFreeParamVecs[i]);
				rangeMax = Math.max(rangeMax, Math.abs(freeParameters[i][k]));
			}
		}
		
		/* Squeezing between -0.5 to 0.5
		for(int i=0; i < freeParameters.length; i++) {
			for(int j=0; j < freeParameters[0].length; j++) {
				//freeParameters[i][j] *= maximumRange/rangeMax;
			}
		}
		*/
		return freeParameters;
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
