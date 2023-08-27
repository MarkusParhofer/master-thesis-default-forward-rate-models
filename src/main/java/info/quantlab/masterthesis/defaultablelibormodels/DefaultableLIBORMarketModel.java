package info.quantlab.masterthesis.defaultablelibormodels;

import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.stochastic.RandomVariable;
/**
 * Interface for a Defaultable LIBOR Market Model.
 * Each Model must have an underlying LIBORMarketModel which takes the role of the undefaultable Market Situation.
 * 
 * @author Markus Parhofer
 */
public interface DefaultableLIBORMarketModel extends LIBORMarketModel {

	RandomVariable getSpreadAtGivenTime(final double evalTime, final int liborIndex, final double time);
	
	RandomVariable getSpreadAtGivenTimeIndex(final double evalTime, final int liborIndex, final int timeIndex);
	
	LIBORMarketModel getUnderlyingNonDefaultableModel();
	
}
