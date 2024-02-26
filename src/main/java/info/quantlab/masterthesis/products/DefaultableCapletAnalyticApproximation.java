package info.quantlab.masterthesis.products;


import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

import java.util.Map;


/* TODO: Handle Analytic SurvivalProbability!!! */
public class DefaultableCapletAnalyticApproximation extends AbstractDefaultableTermStructureProduct{
	final double _strikeRate;
	final double _fixingTime;
	final double _periodLength;
	final double _notional;
	final int _underlyingModelIndex; /* if <0 underlying is undefaultable Model else it is a defaultable model with the specified index. */
	final int _issuerModelIndex; /* if <0 issuer has nondefaultable Model else it is a defaultable model with the specified index. */
	
	public DefaultableCapletAnalyticApproximation(final double strikeRate, final double fixingTime, final double periodLength, 
			final int modelIndexOfUnderlying, final int modelIndexOfIssuer, final double notional) {
		super();
		_strikeRate = strikeRate;
		_fixingTime = fixingTime;
		_periodLength = periodLength;
		_underlyingModelIndex = modelIndexOfUnderlying;
		_issuerModelIndex = modelIndexOfIssuer;
		_notional = notional;
	}
	
	public DefaultableCapletAnalyticApproximation(final double strikeRate, final double fixingTime, final double periodLength, 
			final boolean useDefaultableUnderlying, final boolean useDefaultableIssuer) {
		this(strikeRate, fixingTime, periodLength, useDefaultableUnderlying ? 0 : -1, useDefaultableIssuer ? 0 : -1, 1.0);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model)	throws CalculationException {
		if(evaluationTime != 0.0)
			throw new UnsupportedOperationException("Analytic value can only be approximated at time 0.0");
		try {
			return getValue((LIBORMarketModel)model.getModel());
		} catch(ClassCastException ex) {
			throw new IllegalArgumentException("ProcessModel must be of type DefaultableLIBORMarketModel!");
		}
	}

	public RandomVariable getValue(LIBORMarketModel model) throws CalculationException {
		int fixingTimeIndexOnLIBORPeriods = model.getLiborPeriodIndex(_fixingTime);
		if(fixingTimeIndexOnLIBORPeriods < 0) {
			throw new UnsupportedOperationException("For now fixing time must lie on LIBOR Tenor.");
			// TODO: Handle fixing time not on Libor Periods
		}
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(fixingTimeIndexOnLIBORPeriods);
		if(liborPeriodLength != _periodLength) {
			throw new UnsupportedOperationException("For now period length must be the same as the LIBOR period length of the model.");
			// TODO: Handle payment time not on Libor Periods
		}
		
		double liborForward = model.getForwardRateCurve().getForward(null, _fixingTime);
		double payOffUnit = model.getDiscountCurve().getDiscountFactor(_fixingTime + _periodLength);
		double integratedLIBORCov = 0.0;
		// For getting integrated LIBOR Cov:
		final TimeDiscretization simulationTenor = model.getCovarianceModel().getTimeDiscretization();
		int fixingTimeIndexOnSimTenor = simulationTenor.getTimeIndex(_fixingTime);
		if(fixingTimeIndexOnSimTenor < 0)
			fixingTimeIndexOnSimTenor = - fixingTimeIndexOnSimTenor - 2;
		fixingTimeIndexOnSimTenor--;
		
		if(model instanceof DefaultableLIBORMarketModel defModel) {
			// Adjust Payoff Unit:
			if(_issuerModelIndex >=0)
				// Issuer must survive until payment date for the buyer to get money!
			// TODO: payOffUnit *= defModel.getSurvivalProbability(null, 0.0, _fixingTime + _periodLength).doubleValue();
			{
				final double defFR = 1.0 + (_fixingTime + _periodLength)*calculateFR(defModel, fixingTimeIndexOnLIBORPeriods + 1, false);
				final double fr = 1.0 + (_fixingTime + _periodLength)*calculateFR(defModel, fixingTimeIndexOnLIBORPeriods + 1, true);
				payOffUnit *= fr / defFR;
				// System.out.println("There is no analytic approximation of the Survival Probability");
			}
			else if(_underlyingModelIndex >= 0) {
				// Underlying must survive until fixing date for the buyer to get money: (Otherwise L(fixing; fixing, payment) = 0 => (L-K)^+ = 0)
				final double defFR = 1.0 + _fixingTime*calculateFR(defModel, fixingTimeIndexOnLIBORPeriods, false);
				final double fr = 1.0 + _fixingTime*calculateFR(defModel, fixingTimeIndexOnLIBORPeriods, true);
				payOffUnit *= fr / defFR;
				// System.out.println("There is no analytic approximation of the Survival Probability");
			}
			// Adjust Initial LIBOR
			if(_underlyingModelIndex < 0) {
				liborForward = defModel.getNonDefaultableLIBORModel().getForwardRateCurve().getForward(null, _fixingTime);
			}

		
			// Calculate Covariance:
			if(_underlyingModelIndex < 0)
				integratedLIBORCov = defModel.getNonDefaultableLIBORModel().getIntegratedLIBORCovariance(simulationTenor)[fixingTimeIndexOnSimTenor][fixingTimeIndexOnLIBORPeriods][fixingTimeIndexOnLIBORPeriods];
			else
				integratedLIBORCov = defModel.getIntegratedLIBORCovariance(simulationTenor)[fixingTimeIndexOnSimTenor][fixingTimeIndexOnLIBORPeriods][fixingTimeIndexOnLIBORPeriods];
			
		}
		else if(model instanceof MultiLIBORVectorModel vecModel) {
			// TODO: Handle Implementation
		}
		final double volatility = Math.sqrt(integratedLIBORCov/_fixingTime);
		return Scalar.of(AnalyticFormulas.bachelierOptionValue(liborForward, volatility, _fixingTime, _strikeRate, payOffUnit)).mult(_periodLength*_notional);
	}

	private double calculateFR(DefaultableLIBORMarketModel model, int periodEndIndex, boolean nonDefaultableModel) throws CalculationException {
		double accrualAccount = 1.0;

		// Calculate the value of the forward bond
		for(int periodIndex = 0; periodIndex < periodEndIndex; periodIndex++)
		{
			final double subPeriodLength = model.getLiborPeriod(periodIndex+1) - model.getLiborPeriod(periodIndex);
			double liborOverSubPeriod = 0.0;
			if(nonDefaultableModel) {
				liborOverSubPeriod = model.getNonDefaultableLIBOR(null, 0, periodIndex).doubleValue();
			}
			else {
				liborOverSubPeriod = model.getDefaultableLIBOR(null, 0, periodIndex).doubleValue();
			}
			accrualAccount *= liborOverSubPeriod * subPeriodLength + 1.0;
		}

		return (accrualAccount - 1.0) / model.getLiborPeriod(periodEndIndex);
	}

}
