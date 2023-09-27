package info.quantlab.masterthesis.defaultablecovariancemodels;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

public class SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread extends AbstractDefaultableLIBORCovarianceModel {

	final double[][] _freeParameters;
	
	public SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread(LIBORMarketModel underlyingUndefaultableModel, double[][] freeParameters) {
		super(underlyingUndefaultableModel);
		if(freeParameters.length != underlyingUndefaultableModel.getNumberOfComponents())
			throw new IllegalArgumentException("freeParameters must have the same number of components as underlyingUndefaultableModel i.e. " + underlyingUndefaultableModel.getNumberOfComponents());
		
		_freeParameters = freeParameters;
	}

	@Override
	public RandomVariable[] getFactorLoading(int timeIndex, int component, RandomVariable[] defaultableRealization, RandomVariable[] undefaultableRealization) {
		final int undefaultableFactors = getUnderlyingUndefaultableModel().getNumberOfFactors();
		
		final double periodLength = getLiborPeriodDiscretization().getTimeStep(component);
		
		RandomVariable[] undefaultableFactorLoading = getUnderlyingUndefaultableModel().getCovarianceModel().getFactorLoading(timeIndex, component, undefaultableRealization);
		RandomVariable[] factorLoading = new RandomVariable[getNumberOfFactors()];
		
		// Calculate from underlying undefaultable Model
		RandomVariable relationFactor = defaultableRealization[component].mult(periodLength).add(1.0).discount(undefaultableRealization[component], periodLength);
		for(int k = 0; k < undefaultableFactors; k++) {
			factorLoading[k] = relationFactor.mult(undefaultableFactorLoading[k]);
		}
		
		// Calculate from free parameters
		relationFactor = defaultableRealization[component].sub(undefaultableRealization[component]);
		for(int k = undefaultableFactors; k < getNumberOfFactors(); k++) {
			factorLoading[k] = relationFactor.mult(_freeParameters[component][k - undefaultableFactors]);
		}
		return factorLoading;
	}

	@Override
	public TimeDiscretization getTimeDiscretization() {
		return getUnderlyingUndefaultableModel().getCovarianceModel().getTimeDiscretization();
	}

	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return getUnderlyingUndefaultableModel().getLiborPeriodDiscretization();
	}

	@Override
	public DefaultableLIBORCovarianceModel getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
		LIBORMarketModel underlyingUndefaultableModel = getUnderlyingUndefaultableModel();
		double[][] freeParameters = _freeParameters;
		if(dataModified != null) {
			underlyingUndefaultableModel = (LIBORMarketModel)dataModified.getOrDefault("underlyingUndefaultableModel", underlyingUndefaultableModel);
			freeParameters = (double[][])dataModified.getOrDefault("freeParameters", freeParameters);
		}
		return new SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread(underlyingUndefaultableModel, freeParameters);
	}

	@Override
	public int getNumberOfFactors() {
		return getUnderlyingUndefaultableModel().getNumberOfFactors() + _freeParameters[0].length;
	}

	@Override
	public DefaultableLIBORCovarianceModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel) {
		return new SimpleDefaultableLIBORCovarianceWithGuaranteedPositiveSpread(newUndefaultableModel, _freeParameters);
	}

}
