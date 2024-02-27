package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationModel;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.conditionalexpectation.RegressionBasisFunctionsProvider;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.ConditionalExpectationEstimator;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

import java.util.ArrayList;
import java.util.Arrays;

public class CancellableLoan extends LoanProduct {

    private final TimeDiscretization _tenor;

    private final double[] _coupons;

    private final double _maturityTime;

    private final int _creditorIndex;

    private final int _debtorIndex;

    private final double _redemption;

    private final Perspective _stoppingPerspective;

    private transient RandomVariable _loanValue;

    private transient RandomVariable _optionValue;

    private transient MonteCarloProcess _mcProcessCache;


    public CancellableLoan(TimeDiscretization tenor, double[] coupons, double maturityTime, double redemption, int creditorIndex, int debtorIndex, Perspective valuationPerspective, Perspective stoppingPerspective) {
        super(valuationPerspective);
        _tenor = tenor;
        _coupons = coupons;
        _maturityTime = maturityTime;
        _creditorIndex = creditorIndex;
        _debtorIndex = debtorIndex;
        _redemption = redemption;
        _stoppingPerspective = stoppingPerspective;
    }


    public CancellableLoan(double[] tenor, double[] coupons, double maturityTime, double redemption, int creditorIndex, int debtorIndex) {
        this(new TimeDiscretizationFromArray(tenor), coupons, maturityTime, redemption, creditorIndex, debtorIndex, Perspective.MARKET_IMPLIED, Perspective.MARKET_IMPLIED);
    }

    @Override
    public RandomVariable getValue(double evaluationTime, MultiLIBORVectorModel market, MonteCarloProcess process) throws CalculationException {

        return getLoanValue(market, process).sub(getOptionValue(market,process));
    }

    public RandomVariable getLoanValue(MultiLIBORVectorModel market, MonteCarloProcess process) throws CalculationException {
        if(_loanValue != null && _mcProcessCache == process) {
            return _loanValue;
        }

        _mcProcessCache = process;
        _loanValue = market.getRandomVariableForConstant(0.0);
        int i=0;
        for(double time: _tenor) {
            RandomVariable cashflow = market.getRandomVariableForConstant(_coupons[i++]);
            switch (getValuationPerspective()) {
                case CREDITOR:
                    cashflow = cashflow.mult(market.getSurvivalProbability(process, time, _creditorIndex));
                case DEBTOR:
                case MARKET_IMPLIED:
                cashflow = cashflow.div(market.getDefaultableNumeraire(process, time, _debtorIndex));
                    break;
            }
            _loanValue = _loanValue.add(cashflow);
        }
        return _loanValue;
    }


    public RandomVariable getOptionValue(MultiLIBORVectorModel market, MonteCarloProcess process) throws CalculationException {
        if(_optionValue != null && _mcProcessCache == process)
            return _optionValue;

        int matIndex = _tenor.getTimeIndex(_maturityTime);
        if(matIndex < 0) {
            matIndex = - matIndex - 1;
        }
        matIndex--;
        RandomVariable newLoan = new RandomVariableFromDoubleArray(0.0);
        RandomVariable stoppingCondition = new RandomVariableFromDoubleArray(0.0);
        _optionValue = new RandomVariableFromDoubleArray(0.0);
        for(int j=_tenor.getNumberOfTimes() - 1; j > matIndex; j--) {
            if(j != _tenor.getNumberOfTimes() - 1) {
                ConditionalExpectationEstimator estimator = new MonteCarloConditionalExpectationRegression(getBasisFunctions(_tenor.getTime(j), market, process));

                // Get conditional expectation for stopping condition
                RandomVariable expectationStoppingCond = estimator.getConditionalExpectation(stoppingCondition);
                // Subtract redemption (redemption needs to be discounted)
                RandomVariable redemption = market.getRandomVariableForConstant(_redemption);
                switch (_stoppingPerspective) {
                    case CREDITOR:
                        redemption = redemption.mult(market.getSurvivalProbability(process, _tenor.getTime(j), _creditorIndex));
                    case DEBTOR:
                    case MARKET_IMPLIED:
                        redemption = redemption.div(market.getDefaultableNumeraire(process, _tenor.getTime(j), _debtorIndex));
                        break;
                }
                expectationStoppingCond = expectationStoppingCond.sub(redemption);

                // Do the same for payoff expectation
                RandomVariable expectationLoan = null;
                if (getValuationPerspective() != _stoppingPerspective) {
                    expectationLoan = estimator.getConditionalExpectation(newLoan);
                    redemption = market.getRandomVariableForConstant(_redemption);
                    switch (_stoppingPerspective) {
                        case CREDITOR:
                            redemption = redemption.mult(market.getSurvivalProbability(process, _tenor.getTime(j), _creditorIndex));
                        case DEBTOR:
                        case MARKET_IMPLIED:
                            redemption = redemption.div(market.getDefaultableNumeraire(process, _tenor.getTime(j), _debtorIndex));
                            break;
                    }
                    expectationLoan = expectationLoan.sub(redemption);
                }
                else expectationLoan = expectationStoppingCond;
                expectationLoan = expectationLoan.floor(0.0);


                // Check if stopping condition is met (i.e. K - E[sum(c_i)](omega) > 0)
                _optionValue = _optionValue.apply(
                        (payoffBefore, payoffCondExp, stopCondition)-> {
                            if(stopCondition > 0.0) return payoffCondExp;
                            else return payoffBefore;},
                        expectationLoan, expectationStoppingCond);
            }
            RandomVariable cashflow = market.getRandomVariableForConstant(_coupons[j]);
            switch (getValuationPerspective()) {
                case CREDITOR:
                    cashflow = cashflow.mult(market.getSurvivalProbability(process, _tenor.getTime(j), _creditorIndex));
                case DEBTOR:
                case MARKET_IMPLIED:
                    cashflow = cashflow.div(market.getDefaultableNumeraire(process, _tenor.getTime(j), _debtorIndex));
                    break;
            }
            newLoan = newLoan.add(cashflow);

            if(getValuationPerspective() == _stoppingPerspective) {
                stoppingCondition = newLoan;
            } else {
                RandomVariable cashflow2 = market.getRandomVariableForConstant(_coupons[j]);
                switch (_stoppingPerspective) {
                    case CREDITOR:
                        cashflow2 = cashflow2.mult(market.getSurvivalProbability(process, _tenor.getTime(j), _creditorIndex));
                    case DEBTOR:
                    case MARKET_IMPLIED:
                        cashflow2 = cashflow2.div(market.getDefaultableNumeraire(process, _tenor.getTime(j), _debtorIndex));
                        break;
                }
                stoppingCondition = stoppingCondition.add(cashflow2);
            }
        }
        return _optionValue;
    }

    @Override
    RandomVariable getValue(double evaluationTime, DefaultableLIBORMarketModel defaultableModel, MonteCarloProcess process) throws CalculationException {
        throw new UnsupportedOperationException("Not yet implemented for single models");
    }

    /**
     * Return the basis functions for the regression suitable for this product.
     *
     * @param evaluationTime The condition time.
     * @param model The model
     * @return The basis functions for the regression suitable for this product.
     * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
     */
    public RandomVariable[] getBasisFunctions(final double evaluationTime, final MultiLIBORVectorModel model, final MonteCarloProcess process) throws CalculationException {

        final ArrayList<RandomVariable> basisFunctions = new ArrayList<>();

        // Constant
        final RandomVariable basisFunction = new RandomVariableFromDoubleArray(1.0);//.getRandomVariableForConstant(1.0);
        basisFunctions.add(basisFunction);

        int liborIndex = model.getLiborPeriodIndex(evaluationTime);
        if(liborIndex < 0) {
            liborIndex = - liborIndex - 1;
        }
        if(liborIndex == model.getNumberOfLiborPeriods()) {
            liborIndex--;
        }

        int timeIndex = process.getTimeIndex(evaluationTime);
        if(timeIndex < 0) {
            timeIndex = - timeIndex - 1;
        }
        if(timeIndex == process.getTimeDiscretization().getNumberOfTimes()) {
            timeIndex--;
        }
        // forward rate to the next period of all models
        final RandomVariable rateShort = model.getNonDefaultableLIBOR(process, timeIndex, liborIndex);
        final RandomVariable discountShort = rateShort.mult(model.getLiborPeriod(liborIndex + 1) - model.getLiborPeriod(liborIndex)).add(1.0).invert();
        basisFunctions.add(discountShort);
        basisFunctions.add(discountShort.pow(2.0));
        for(int i=0; i < model.getNumberOfDefaultableModels(); i++) {
            final RandomVariable defRate = model.getDefaultableLIBOR(process, timeIndex, liborIndex, i);
            final RandomVariable defDiscount = defRate.mult(model.getLiborPeriod(liborIndex + 1) - model.getLiborPeriod(liborIndex)).add(1.0).invert();
            basisFunctions.add(defDiscount);
            basisFunctions.add(defDiscount.pow(2.0));
        }

        // forward rate to the end of the product
        final RandomVariable rateLong = model.getNonDefaultableModel().getForwardRate(model.getNonDefaultableProcess(process), evaluationTime, evaluationTime, _tenor.getLastTime());
        final RandomVariable discountLong = rateLong.mult(_tenor.getLastTime()-evaluationTime).add(1.0).invert();
        basisFunctions.add(discountLong);
        basisFunctions.add(discountLong.pow(2.0));
        for(int i=0; i < model.getNumberOfDefaultableModels(); i++) {
            final RandomVariable defRateLong = model.getNonDefaultableModel().getForwardRate(model.getNonDefaultableProcess(process), evaluationTime, evaluationTime, _tenor.getLastTime());
            final RandomVariable defDiscountLong = defRateLong.mult(_tenor.getLastTime()-evaluationTime).add(1.0).invert();
            basisFunctions.add(defDiscountLong);
            basisFunctions.add(defDiscountLong.pow(2.0));
        }

        // Numeraire
        final RandomVariable numeraire = model.getNumeraire(process, evaluationTime).invert();
        basisFunctions.add(numeraire);
        for(int i=0; i < model.getNumberOfDefaultableModels(); i++) {
            final RandomVariable defNumeraire = model.getDefaultableNumeraire(process, evaluationTime, i).invert();
            basisFunctions.add(numeraire);
        }

        return basisFunctions.toArray(new RandomVariable[0]);
    }

}
