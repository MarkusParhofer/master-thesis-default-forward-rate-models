/**
 * 
 */
package info.quantlab.masterthesis.functional;

import org.apache.commons.lang3.Validate;

import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * This class is used for some functional methods. They are all helpers for the master thesis.
 * 
 * @author Markus Parhofer
 * @version 1.0
 */
public class Functional {

	/**
	 * No construction of this class allowed. All static methods.
	 */
	private Functional() {}
	
	/**
	 * Creates a new instance of IndependentIncrements that is basically a copy of increments,
	 * but uses only the range factors that is specified.
	 * @param increments The original instance of independent increments
	 * @param firstFactor The first factor to consider in the new instance
	 * @param lastFactor The last factor to consider in the new instance
	 * @return A new instance with reduced factors
	 */
	public static IndependentIncrements getRangeOfFactors(IndependentIncrements increments, int firstFactor, int lastFactor) {
		Validate.isTrue(0 <= firstFactor && firstFactor < lastFactor && lastFactor < increments.getNumberOfFactors(), 
				"The following does not hold: 0 <= firstFactor (=%d) < lastFactor (=%d) < increments.numberOfFactors (=%d)",
				firstFactor, lastFactor, increments.getNumberOfFactors());

		if(firstFactor == 0 && increments.getNumberOfFactors() == lastFactor + 1) {
			return increments;			
		}
		return new IndependentIncrements() {
			private final int numberOfFactors = lastFactor - firstFactor + 1;
			
			@Override
			public RandomVariable getIncrement(int timeIndex, int factor) {
				Validate.isTrue(factor <= lastFactor - firstFactor);
				return increments.getIncrement(timeIndex, factor + firstFactor);
			}

			@Override
			public TimeDiscretization getTimeDiscretization() {
				return increments.getTimeDiscretization();
			}

			@Override
			public int getNumberOfFactors() {
				return numberOfFactors;
			}

			@Override
			public int getNumberOfPaths() {
				return increments.getNumberOfPaths();
			}

			@Override
			public RandomVariable getRandomVariableForConstant(double value) {
				return increments.getRandomVariableForConstant(value);
			}

			@Override
			public IndependentIncrements getCloneWithModifiedSeed(int seed) {
				IndependentIncrements originClone = increments.getCloneWithModifiedSeed(seed);
				return getRangeOfFactors(originClone, firstFactor, lastFactor);
			}

			@Override
			public IndependentIncrements getCloneWithModifiedTimeDiscretization(TimeDiscretization newTimeDiscretization) {
				IndependentIncrements originClone = increments.getCloneWithModifiedTimeDiscretization(newTimeDiscretization);
				return getRangeOfFactors(originClone, firstFactor, lastFactor);
			}
			
		};
	}
	
	/**
	 * Creates a new instance of IndependentIncrements that is basically a copy of increments,
	 * but uses only the first factors.
	 * @param increments The original instance of independent increments
	 * @param numberOfFactors The number factors to consider in the new instance
	 * @return A new instance with reduced factors
	 */
	public static IndependentIncrements getFirstFactors(IndependentIncrements increments, int numberOfFactors) {
		return getRangeOfFactors(increments, 0, numberOfFactors - 1);
	}

}
