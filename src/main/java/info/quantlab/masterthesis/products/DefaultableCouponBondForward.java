package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.products.AbstractTermStructureMonteCarloProduct;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

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
            return getValue(evaluationTime, defModel, model.getProcess());
        }
        else if(model.getModel() instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, model.getProcess(), 1, 0);
        }
        // TODO: Implement non defaultable version!
        return null;
    }

    public RandomVariable getValue(double evaluationTime, DefaultableLIBORMarketModel model, MonteCarloProcess process) throws CalculationException {
        if(evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime higher than 0 not supported");
        }
        final double maturityTime = model.getLiborPeriod(_maturityIndex);
        RandomVariable payoff = model.getRandomVariableForConstant(_strike);
        for(int i = _maturityIndex + 1; i < model.getLiborPeriodDiscretization().getNumberOfTimes(); i++) {
            final double couponTime = model.getLiborPeriod(i);
            RandomVariable coupon = model.getRandomVariableForConstant(_couponRates[i - _maturityIndex - 1]);
            if(i == model.getLiborPeriodDiscretization().getNumberOfTimeSteps()) {
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
