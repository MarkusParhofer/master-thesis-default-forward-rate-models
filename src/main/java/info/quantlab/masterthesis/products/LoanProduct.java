package info.quantlab.masterthesis.products;

import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

public abstract class LoanProduct extends AbstractDefaultableTermStructureProduct {
    public enum Perspective {
        CREDITOR,
        DEBTOR,
        MARKET_IMPLIED
    }

    private final Perspective _valuationPerspective;

    public Perspective getValuationPerspective() {
        return _valuationPerspective;
    }

    LoanProduct(Perspective valuationPerspective) {
        _valuationPerspective = valuationPerspective;
    }

    @Override
    public RandomVariable getValue(double evaluationTime, TermStructureMonteCarloSimulationModel model) throws CalculationException {
        if(model.getModel() instanceof MultiLIBORVectorModel multiModel) {
            return getValue(evaluationTime, multiModel, model.getProcess());
        }
        if(model.getModel() instanceof DefaultableLIBORMarketModel defModel) {
            return getValue(evaluationTime, defModel, model.getProcess());
        }
        if(model.getModel() instanceof LIBORMarketModel nonDefModel)
            return getValue(evaluationTime, nonDefModel, model.getProcess());

        throw new UnsupportedOperationException("Model must be of type MultiLIBORVectorModel, DefaultableLIBORMarketModel or LIBORMarketModel");
    }

    /**
     * Gets the value of the product as a RandomVariable
     * @param evaluationTime Time of valuation
     * @param market Model reflecting the models of all relevant parties: (e.g. debtor and creditor).
     * @param process Process simulating the model.
     * @return The value of the option as RandomVariable, not yet as expectation.
     * @throws CalculationException if the calculation fails (mostly due to the model).
     */
    abstract RandomVariable getValue(double evaluationTime, MultiLIBORVectorModel market, MonteCarloProcess process) throws CalculationException;

    /**
     * Gets the value of the product as a RandomVariable
     * @param evaluationTime Time of valuation
     * @param defaultableModel Model reflecting the non-defaultable model and a defaultable party.
     * @param process Process simulating the model.
     * @return The value of the option as RandomVariable, not yet as expectation.
     * @throws CalculationException if the calculation fails (mostly due to the model).
     */
    abstract RandomVariable getValue(double evaluationTime, DefaultableLIBORMarketModel defaultableModel, MonteCarloProcess process) throws CalculationException;

    /**
     * Gets the value of the product as a RandomVariable
     * @param evaluationTime Time of valuation
     * @param model Model reflecting the non-defaultable model.
     * @param process Process simulating the model.
     * @return The value of the option as RandomVariable, not yet as expectation.
     * @throws CalculationException if the calculation fails (mostly due to the model).
     */
    RandomVariable getValue(double evaluationTime, LIBORMarketModel model, MonteCarloProcess process) throws CalculationException, UnsupportedOperationException {
        throw new UnsupportedOperationException("Loan products should use a defaultable model for valuation!");
    }
}
