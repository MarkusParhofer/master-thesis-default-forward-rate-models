package info.quantlab.masterthesis.products;


import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

/**
 * This class implements the valuation of a zero coupon bond.
 *
 * @author Markus Parhofer
 * @version 1.0
 */
public class DefaultableZeroCouponBond extends AbstractDefaultableTermStructureProduct {

	private double maturity;

	/**
	 * @param maturity The maturity given as double.
	 */
	public DefaultableZeroCouponBond(final double maturity) {
		super();
		this.maturity = maturity;
	}

	/**
	 * This method returns the value random variable of the product within the specified model, evaluated at a given evalutationTime.
	 * Note: For a lattice this is often the value conditional to evalutationTime, for a Monte-Carlo simulation this is the (sum of) value discounted to evaluation time.
	 * Cashflows prior evaluationTime are not considered.
	 *
	 * @param evaluationTime The time on which this products value should be observed.
	 * @param model The model used to price the product.
	 * @return The random variable representing the value of the product discounted to evaluation time
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	@Override
	public RandomVariable getValue(final double evaluationTime, final TermStructureMonteCarloSimulationModel model) throws CalculationException {

		DefaultableLIBORMarketModel defModel = (DefaultableLIBORMarketModel)model.getModel();
		MonteCarloProcess process = model.getProcess();
		// Get random variables
		final RandomVariable	numeraire				= defModel.getDefaultableNumeraire(process, maturity);
		final RandomVariable	monteCarloProbabilities	= model.getMonteCarloWeights(maturity);

		// Calculate numeraire relative value
		RandomVariable values = model.getRandomVariableForConstant(1.0);
		values = values.div(numeraire).mult(monteCarloProbabilities);

		// Convert back to values
		final RandomVariable	numeraireAtEvaluationTime				= model.getNumeraire(evaluationTime);
		final RandomVariable	monteCarloProbabilitiesAtEvaluationTime	= model.getMonteCarloWeights(evaluationTime);
		values = values.mult(numeraireAtEvaluationTime).div(monteCarloProbabilitiesAtEvaluationTime);

		// Return values
		return values;
	}
}
