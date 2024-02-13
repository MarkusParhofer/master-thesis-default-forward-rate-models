package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

public class DefaultableCouponBondForward extends AbstractTermStructureMonteCarloProduct {

    public DefaultableCouponBondForward(double strike, double nominal, double[] couponRates, int maturityIndex) {
        _strike = strike;
        _nominal = nominal;
        _couponRates = couponRates;
        _maturityIndex = maturityIndex;
    }

    double _strike;
    double _nominal;
    double[] _couponRates;
    int _maturityIndex;



    @Override
    public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
        if(model.getModel() instanceof DefaultableLIBORMarketModel defModel) {
            return getValue(evaluationTime, defModel, model.getProcess(), true);
        }
        else if(model.getModel() instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, model.getProcess(), 1, 0);
        }

        // Case No Defaultable model given:
        return getValue(evaluationTime, (LIBORMarketModel) model.getModel(), model.getProcess());
    }

    public RandomVariable getValue(double evaluationTime, LIBORMarketModel model, MonteCarloProcess process) throws CalculationException {
        if(model instanceof DefaultableLIBORMarketModel defModel) {
            return getValue(evaluationTime, defModel, process, true);
        }
        else if(model instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, process, 1, 0);
        }

        // Case No Defaultable model given:
        if(evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime higher than 0 not supported");
        }
        final int terminalIndex = Math.min(model.getLiborPeriodDiscretization().getNumberOfTimes(), _couponRates.length + _maturityIndex + 1);

        final double maturityTime = model.getLiborPeriod(_maturityIndex);
        RandomVariable payoff = model.getRandomVariableForConstant(_strike);
        for(int i = _maturityIndex + 1; i < terminalIndex; i++) {
            final double couponTime = model.getLiborPeriod(i);
            RandomVariable coupon = model.getRandomVariableForConstant(_couponRates[i - _maturityIndex - 1]);
            if(i == terminalIndex - 1) {
                coupon = coupon.add(1.0); // Add redemption
            }
            RandomVariable forwardRate = model.getForwardRate(process, maturityTime, maturityTime, couponTime);
            RandomVariable bond = new Scalar(1.0).discount(forwardRate, couponTime - maturityTime); // P(T_s;T_i)
            coupon = coupon.mult(_nominal).mult(bond);
            payoff = payoff.sub(coupon);
        }
        payoff = payoff.div(model.getNumeraire(process, maturityTime));
        return payoff;
    }

    private RandomVariable getValue(double evaluationTime, DefaultableLIBORMarketModel model, MonteCarloProcess process, boolean reserved) throws CalculationException {
        if(evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime higher than 0 not supported");
        }
        final int terminalIndex = Math.min(model.getLiborPeriodDiscretization().getNumberOfTimes(), _couponRates.length + _maturityIndex + 1);
        final double maturityTime = model.getLiborPeriod(_maturityIndex);
        RandomVariable payoff = model.getRandomVariableForConstant(_strike);
        for(int i = _maturityIndex + 1; i < terminalIndex; i++) {
            final double couponTime = model.getLiborPeriod(i);
            RandomVariable coupon = model.getRandomVariableForConstant(_couponRates[i - _maturityIndex - 1]);
            if(i == terminalIndex - 1) {
                coupon = coupon.add(1.0); // Add redemption
            }
            coupon = coupon.mult(_nominal).mult(model.getDefaultableBond(process, maturityTime, couponTime));
            payoff = payoff.sub(coupon);
        }
        payoff = payoff.div(model.getDefaultableNumeraire(process, maturityTime));
        return payoff;
    }

    public RandomVariable getValue(double evaluationTime, MultiLIBORVectorModel model, MonteCarloProcess process, int debtorIndex, int creditorIndex) throws CalculationException {
        final double maturityTime = model.getLiborPeriod(_maturityIndex);
        final DefaultableLIBORMarketModel debtorModel = model.getDefaultableModel(debtorIndex);
        final MonteCarloProcess debtorProcess= model.getDefaultableProcess(process, debtorIndex);
        final DefaultableLIBORMarketModel creditorModel = model.getDefaultableModel(creditorIndex);
        final MonteCarloProcess creditorProcess= model.getDefaultableProcess(process, creditorIndex);

        return getValue(evaluationTime, debtorModel, debtorProcess).mult(creditorModel.getSurvivalProbability(creditorProcess, maturityTime));
    }

}
