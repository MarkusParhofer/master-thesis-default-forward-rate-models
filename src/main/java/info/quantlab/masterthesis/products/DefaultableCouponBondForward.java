package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.util.TriFunction;

public class DefaultableCouponBondForward extends LoanProduct {

    public DefaultableCouponBondForward(TimeDiscretization tenor, double[] coupons, double strike, int creditorIndex, int debtorIndex, Perspective valuationPerspective) {
        super(valuationPerspective);
        _tenor = tenor;
        _strike = strike;
        _coupons = coupons;
        _creditorIndex = creditorIndex;
        _debtorIndex = debtorIndex;
    }

    public DefaultableCouponBondForward(double[] tenor, double[] coupons, double strike, int creditorIndex, int debtorIndex, Perspective valuationPerspective) {
        this(new TimeDiscretizationFromArray(tenor), coupons, strike, creditorIndex, debtorIndex, valuationPerspective);
    }

    public DefaultableCouponBondForward(double[] tenor, double[] couponRates, double nominal, double strike, int creditorIndex, int debtorIndex, Perspective valuationPerspective) {
        this(tenor, couponRates.clone(), strike, creditorIndex, debtorIndex, valuationPerspective);
        for(int i =0; i < _coupons.length; i++) {
            _coupons[i] *= (_tenor.getTime(i+1) - _tenor.getTime(i)) * nominal;
        }
        _coupons[_coupons.length - 1] += nominal;
    }

    public DefaultableCouponBondForward(TimeDiscretization tenor, double[] coupons, double strike) {
        this(tenor, coupons, strike, 0, 0, Perspective.DEBTOR);
    }

    public DefaultableCouponBondForward(double[] tenor, double[] coupons, double strike) {
        this(new TimeDiscretizationFromArray(tenor), coupons, strike);
    }

    public DefaultableCouponBondForward(double[] tenor, double[] couponRates, double nominal, double strike) {
        this(tenor, couponRates, nominal, strike, 0, 0, Perspective.DEBTOR);
    }

    double _strike;
    TimeDiscretization _tenor;
    double[] _coupons;
    int _debtorIndex;
    int _creditorIndex;


    @Override
    public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
        if (model.getModel() instanceof DefaultableLIBORMarketModel defModel) {
            return getValue(evaluationTime, defModel, model.getProcess());
        } else if (model.getModel() instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, model.getProcess());
        }

        // Case No Defaultable model given:
        return getValue(evaluationTime, (LIBORMarketModel) model.getModel(), model.getProcess());
    }

    @Override
    RandomVariable getValue(double evaluationTime, MultiLIBORVectorModel model, MonteCarloProcess process) throws CalculationException {
        if (evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime larger than 0 is not supported");
        }
        final int terminalIndex = _tenor.getNumberOfTimes();
        RandomVariable payoff = model.getRandomVariableForConstant(_strike);

        TriFunction<MonteCarloProcess, Double, Integer, RandomVariable> numeraire = (p, t, m) -> {
            try {
                if(m<0) {
                    return model.getNumeraire(p, t);
                } else {
                    return model.getDefaultableNumeraire(p, t, m);
                }
            } catch (CalculationException e) {
                throw new RuntimeException(e);
            }
        };
        TriFunction<MonteCarloProcess, Double, Integer, RandomVariable> survivalProb = (p, t, m) -> {
            try {
                if(m<0) {
                    return model.getRandomVariableForConstant(1.0);
                } else {
                    return model.getSurvivalProbability(p, t, m);
                }
            } catch (CalculationException e) {
                throw new RuntimeException(e);
            }
        };
        payoff = payoff.mult(survivalProb.apply(process, _tenor.getFirstTime(), _creditorIndex));
        payoff = payoff.div(numeraire.apply(process, _tenor.getFirstTime(), _debtorIndex));

        for (int i = 1; i < terminalIndex; i++) {
            RandomVariable coupon = model.getRandomVariableForConstant(_coupons[i - 1]);
            switch(getValuationPerspective()) {
                case CREDITOR:
                    coupon = coupon.mult(survivalProb.apply(process, _tenor.getTime(i), _creditorIndex));
                    coupon = coupon.div(numeraire.apply(process, _tenor.getTime(i), _debtorIndex));
                    break;
                case MARKET_IMPLIED:
                case DEBTOR:
                    coupon = coupon.mult(survivalProb.apply(process, _tenor.getTime(0), _creditorIndex));
                    coupon = coupon.div(numeraire.apply(process, _tenor.getTime(i), _debtorIndex));
                    break;
            }
            payoff = payoff.sub(coupon);
        }
        return payoff;
    }

    @Override
    RandomVariable getValue(double evaluationTime, DefaultableLIBORMarketModel model, MonteCarloProcess process) throws CalculationException {
        if (evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime larger than 0 is not supported");
        }

        final int terminalIndex = _tenor.getNumberOfTimes();
        RandomVariable payoff = model.getRandomVariableForConstant(_strike).div(model.getDefaultableNumeraire(process, _tenor.getTime(0)));
        for (int i = 1; i < terminalIndex; i++) {
            RandomVariable coupon = model.getRandomVariableForConstant(_coupons[i - 1]);
            coupon = coupon.div(model.getDefaultableNumeraire(process, _tenor.getTime(i)));
            payoff = payoff.sub(coupon);
        }
        return payoff;
    }

    public RandomVariable getValue(double evaluationTime, LIBORMarketModel model, MonteCarloProcess process) throws CalculationException {
        if (model instanceof DefaultableLIBORMarketModel defModel) {
            return getValue(evaluationTime, defModel, process);
        } else if (model instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, process);
        }

        // Case No Defaultable model given:
        if (evaluationTime > 0) {
            throw new IllegalArgumentException("For now evaluationTime larger than 0 is not supported");
        }

        final int terminalIndex = _tenor.getNumberOfTimes();
        RandomVariable payoff = model.getRandomVariableForConstant(_strike).div(model.getNumeraire(process, _tenor.getTime(0)));
        for (int i = 1; i < terminalIndex; i++) {
            RandomVariable coupon = model.getRandomVariableForConstant(_coupons[i - 1]);
            coupon = coupon.div(model.getNumeraire(process, _tenor.getTime(i)));
            payoff = payoff.sub(coupon);
        }
        return payoff;
    }
}
