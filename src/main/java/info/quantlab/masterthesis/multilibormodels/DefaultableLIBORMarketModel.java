package info.quantlab.masterthesis.multilibormodels;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

public interface DefaultableLIBORMarketModel extends LIBORMarketModel {

	public LIBORMarketModel getUndefaultableLIBORModel();
	
	public RandomVariable[] getDriftOfUndefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);
	
	public RandomVariable[] getDriftOfDefaultableModel(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex,	RandomVariable[] realizationPredictor);
	
	@Override
	public DefaultableLIBORCovarianceModel getCovarianceModel();
	
	@Override
	public DefaultableLIBORMarketModel getCloneWithModifiedCovarianceModel(LIBORCovarianceModel newCovarianceModel);
	
	public DefaultableLIBORMarketModel getCloneWithModifiedUndefaultableModel(LIBORMarketModel newUndefaultableModel);
	
	/**
	 * This method will only get the defaultable LIBOR Rate. 
	 */
	@Override
	public RandomVariable getLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;

	public RandomVariable getUndefaultableLIBOR(MonteCarloProcess process, int timeIndex, int liborIndex) throws CalculationException;

}
