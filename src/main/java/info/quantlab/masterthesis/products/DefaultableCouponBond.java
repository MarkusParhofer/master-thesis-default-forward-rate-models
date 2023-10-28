package info.quantlab.masterthesis.products;

import java.util.Arrays;

import info.quantlab.masterthesis.multilibormodels.DefaultableLIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

public class DefaultableCouponBond extends AbstractSellerDefaultableTermStructureProduct {

	final private TimeDiscretization _paymentTenor;
	final private double[] _coupons;
	final private double _notional;
	
	public DefaultableCouponBond(TimeDiscretization paymentTenor, double[] couponRates, double notional) {
		_paymentTenor = paymentTenor;
		_coupons = couponRates;
		_notional = notional;
	}
	
	public DefaultableCouponBond(TimeDiscretization paymentTenor, double couponRate, double notional) {
		this(paymentTenor, new double[paymentTenor.getNumberOfTimes()], notional);
		Arrays.fill(_coupons, couponRate);
	}

	@Override
	public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
		RandomVariable values = model.getRandomVariableForConstant(0.0);
		
		for(int i = 0; i < _paymentTenor.getNumberOfTimes(); i++) {
			final RandomVariable	numeraire				= model.getNumeraire(_paymentTenor.getTime(i));
			final RandomVariable	monteCarloWeights	= model.getMonteCarloWeights(model.getTimeIndex(_paymentTenor.getTime(i)));

			RandomVariable payoff = model.getRandomVariableForConstant(_coupons[i] * _notional);
			if(i == _paymentTenor.getNumberOfTimes() - 1) {
				payoff = payoff.add(_notional);
			}
			// Adjust for possibility of default between evaluationTime and paymentDate:
			if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
				payoff = payoff.mult(defaultableModel.getSurvivalProbability(model.getProcess(), evaluationTime, _paymentTenor.getTime(i)));
			}
			
			payoff = payoff.div(numeraire).mult(monteCarloWeights);
			values = values.add(payoff);
		}
		final RandomVariable	numeraireAtEvaluationTime = model.getNumeraire(evaluationTime);
		final RandomVariable	monteCarloWeightsAtEvaluationTime = model.getMonteCarloWeights(model.getTimeIndex(evaluationTime));

		return values.mult(numeraireAtEvaluationTime).div(monteCarloWeightsAtEvaluationTime);
	}

	public static double getCouponRateForNotionalPrice(TimeDiscretization paymentTenor, TermStructureMonteCarloSimulationModel model) throws CalculationException {
		double value = 1.0;
		if(model.getModel() instanceof DefaultableLIBORMarketModel defaultableModel) {
			value -= defaultableModel.getDefaultableBond(model.getProcess(), 0.0, paymentTenor.getTime(paymentTenor.getNumberOfTimes() - 1)).getAverage();
			double nominator = 0.0;
			for(double time: paymentTenor) {
				nominator += defaultableModel.getDefaultableBond(model.getProcess(), 0.0, time).getAverage();
			}
			value /= nominator;
		}
		else {
			value -= model.getModel().getForwardDiscountBond(model.getProcess(), 0.0, paymentTenor.getTime(paymentTenor.getNumberOfTimes() - 1)).getAverage();
			double nominator = 0.0;
			for(double time: paymentTenor) {
				nominator += model.getModel().getForwardDiscountBond(model.getProcess(), 0.0, time).getAverage();
			}
			value /= nominator;
		}
		return value;
	}
}
