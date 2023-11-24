package info.quantlab.masterthesis.functional;

import org.apache.commons.lang3.Validate;

import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

public class FunctionsOnIndependentIncrements {

	private FunctionsOnIndependentIncrements() {}

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

	/**
	 * Creates a new instance of IndependentIncrements that is basically a copy of increments,
	 * but uses only the selected factors. The order of the selection is important and will be used. A duplicate in the array will cause an exception, 
	 * as this would destroy the independence of the increments.
	 * @param increments The original instance of independent increments.
	 * @param selectedFactorIndices The indices of the factors to consider in the new instance. The factors will be ordered as given by this array.
	 * @return A new instance with reordered and possibly reduced factors.
	 */
	public static IndependentIncrements getSelectedFactors(IndependentIncrements increments, int[] selectedFactorIndices) {
		final int originNumberOfFactors = increments.getNumberOfFactors();
		for(int k=0; k < selectedFactorIndices.length; k++) {
			if(selectedFactorIndices[k] >= originNumberOfFactors || selectedFactorIndices[k] < 0)
				throw new IllegalArgumentException("selectedFactorIndices includes a factor index that was not in the original independent increments!");
			for(int j=0; j < k; j++) {
				if(selectedFactorIndices[k] == selectedFactorIndices[j])
					throw new IllegalArgumentException("selectedFactorIndices contains duplicates!");
			}
		}
		
		return new IndependentIncrements() {
			private final int[] factorSel = selectedFactorIndices.clone();
			
			@Override
			public TimeDiscretization getTimeDiscretization() {
				return increments.getTimeDiscretization();
			}
			
			@Override
			public RandomVariable getRandomVariableForConstant(double value) {
				return increments.getRandomVariableForConstant(value);
			}
			
			@Override
			public int getNumberOfPaths() {
				return increments.getNumberOfPaths();
			}
			
			@Override
			public int getNumberOfFactors() {
				return factorSel.length;
			}
			
			@Override
			public RandomVariable getIncrement(int timeIndex, int factor) {
				return increments.getIncrement(timeIndex, factorSel[factor]);
			}
			
			@Override
			public IndependentIncrements getCloneWithModifiedSeed(int seed) {
				IndependentIncrements originClone = increments.getCloneWithModifiedSeed(seed);
				return getSelectedFactors(originClone, factorSel);
			}

			@Override
			public IndependentIncrements getCloneWithModifiedTimeDiscretization(TimeDiscretization newTimeDiscretization) {
				IndependentIncrements originClone = increments.getCloneWithModifiedTimeDiscretization(newTimeDiscretization);
				return getSelectedFactors(originClone, factorSel);
			}
		};
	}
	
	
}
