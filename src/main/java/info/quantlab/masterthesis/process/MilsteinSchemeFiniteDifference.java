package info.quantlab.masterthesis.process;

import java.util.Arrays;
import java.util.Map;

import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

public class MilsteinSchemeFiniteDifference extends AbstractMilsteinScheme implements MonteCarloProcess {

	// Used locally for multi-threadded calculation.
	public enum FiniteDifferenceMethod {
		FORWARD,
		BACKWARD,
		CENTRAL
	}
	
	private double m_DeltaXForFiniteDifference;
	private FiniteDifferenceMethod m_Method;	
	
	public MilsteinSchemeFiniteDifference(IndependentIncrements stochasticDriver, TimeDiscretization timeDiscretization, ProcessModel model, FiniteDifferenceMethod finiteDifferenceMethod, double hForfiniteDifference) {
		super(stochasticDriver, timeDiscretization, model);
		m_Method = finiteDifferenceMethod;
		m_DeltaXForFiniteDifference = hForfiniteDifference;
	}

	
	@Override
	protected RandomVariable[] getDifferentialOfFactorLoading(int componentIndex, int timeIndex, RandomVariable[] factorLoading) {
		// Copy array otherwise we modify the process values!
		RandomVariable[] adjustedProcess = Arrays.copyOf(getProcessValues(timeIndex), getNumberOfComponents()); 
		
		RandomVariable[] factorLoadingBackward = null;
		RandomVariable[] factorLoadingForward = null;
		
		// Calculate backward adjusted FL:
		if(m_Method == FiniteDifferenceMethod.BACKWARD || m_Method == FiniteDifferenceMethod.CENTRAL) {
			RandomVariable originalValue = adjustedProcess[componentIndex];
			adjustedProcess[componentIndex] = adjustedProcess[componentIndex].sub(m_DeltaXForFiniteDifference);
			factorLoadingBackward = getFactorLoading(timeIndex, componentIndex, adjustedProcess);
			if(m_Method == FiniteDifferenceMethod.BACKWARD) {
				for(int k=0; k < factorLoadingBackward.length; k++)
					factorLoadingBackward[k] = factorLoading[k].sub(factorLoadingBackward[k]).div(m_DeltaXForFiniteDifference);
				
				return factorLoadingBackward;
			}
			// Adjust back to normal state
			adjustedProcess[componentIndex] = originalValue;	
		}
		
		// Calculate forward adjusted FL:
		adjustedProcess[componentIndex] = adjustedProcess[componentIndex].add(m_DeltaXForFiniteDifference);
		factorLoadingForward = getFactorLoading(timeIndex, componentIndex, adjustedProcess);
		if(m_Method == FiniteDifferenceMethod.FORWARD) {
			for(int k=0; k < factorLoadingForward.length; k++)
				factorLoadingForward[k] = factorLoadingForward[k].sub(factorLoading[k]).div(m_DeltaXForFiniteDifference);
			
			return factorLoadingForward;
		}
		
		// Calculate central difference
		for(int k=0; k < factorLoadingForward.length; k++)
			factorLoadingForward[k] = factorLoadingForward[k].sub(factorLoadingBackward[k]).div(2.0d * m_DeltaXForFiniteDifference);
		
		return factorLoadingForward;
	}
	

	@Override
	public MilsteinSchemeFiniteDifference getCloneWithModifiedModel(ProcessModel model) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MilsteinSchemeFiniteDifference getCloneWithModifiedData(Map<String, Object> dataModified) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MilsteinSchemeFiniteDifference getCloneWithModifiedSeed(int seed) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public MilsteinSchemeFiniteDifference clone() {
		// TODO Auto-generated method stub
		return null;
	}
	
}
