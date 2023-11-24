package info.quantlab.masterthesis.process;

import java.util.Map;

import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.util.TriFunction;

public class MilsteinSchemeAnalyticDerivative extends AbstractMilsteinScheme implements MonteCarloProcess {
	
	private TriFunction<Integer, Double, RandomVariable[], RandomVariable[]> m_AnalyticDerivative;

	/**
	 * Constructs a Milstein Scheme where the Factor Loadings have analytic derivatives.
	 * 
	 * @see
	 * {@link https://en.wikipedia.org/wiki/Milstein_method}.
	 * 
	 * @param stochasticDriver The Stochastic driver i.e. the Brownian Motion.
	 * @param timeDiscretization Time Discretization on which the process shall be generated
	 * @param model Model for which the Process is generated.
	 * @param analyticDerivative The derivative of the factorLoadings used for the Milstein adjustment.<br>
	 * The function must be of the form: 
	 * <b><code>apply(int componentIndex, double time, RandomVariable[] realizations)</code></b>
	 */
	public MilsteinSchemeAnalyticDerivative(IndependentIncrements stochasticDriver, TimeDiscretization timeDiscretization, ProcessModel model, TriFunction<Integer, Double, RandomVariable[], RandomVariable[]> analyticDerivative) {
		super(stochasticDriver, timeDiscretization, model);
		m_AnalyticDerivative = analyticDerivative;
	}


	@Override
	protected RandomVariable[] getDifferentialOfFactorLoading(int componentIndex, int timeIndex, RandomVariable[] factorLoadings) {
		return m_AnalyticDerivative.apply(componentIndex, getTime(timeIndex), getProcessValues(timeIndex));
	}
	
	@Override
	public MilsteinSchemeAnalyticDerivative getCloneWithModifiedModel(ProcessModel model) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MilsteinSchemeAnalyticDerivative getCloneWithModifiedData(Map<String, Object> dataModified) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public MilsteinSchemeAnalyticDerivative getCloneWithModifiedSeed(int seed) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MilsteinSchemeAnalyticDerivative clone() {
		// TODO Auto-generated method stub
		return null;
	}

}
