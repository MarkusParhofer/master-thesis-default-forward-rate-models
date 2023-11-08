package sandbox;

import java.util.Arrays;

import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.BrownianMotionFromRandomNumberGenerator;
import net.finmath.montecarlo.IndependentIncrements;
import net.finmath.montecarlo.IndependentIncrementsFromICDF;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class Playground {

	private static long startTime;
	
	public static void tic() {
		startTime = System.nanoTime();
	}
	
	public static long toc() {
		long result = System.nanoTime();
		result -= startTime;
		System.out.println("Elapsed Time is: " + result/1000 + " Microseconds");
		return result;
	}
	
	public static void main(String[] args) {
		TimeDiscretization tenor = new TimeDiscretizationFromArray(0.0, 20, 0.5);
		IndependentIncrements stochasticDriver = new BrownianMotionFromMersenneRandomNumbers(tenor, 1000, 10000, 1935);
		
		final RandomVariable[] incrementFifteen = stochasticDriver.getIncrement(15);
		final RandomVariable[] incrementSixteen = stochasticDriver.getIncrement(16);
		tic();
		RandomVariable[] firstArray = Arrays.copyOf(incrementFifteen, 2000);
		for(int i=0; i < incrementSixteen.length; i++) {
			firstArray[incrementFifteen.length + i] = incrementSixteen[i];
		}
		toc();
		
		tic();
		RandomVariable[] secondArray = new RandomVariable[incrementFifteen.length + incrementSixteen.length];
		final int separator = incrementFifteen.length;
		Arrays.parallelSetAll(secondArray, i -> i < separator ? incrementFifteen[i] : incrementSixteen[i - separator]);
		toc();
		
		tic();
		RandomVariable[] thirdArray = new RandomVariable[incrementFifteen.length + incrementSixteen.length];
		final int separator2 = incrementFifteen.length;
		Arrays.setAll(secondArray, i -> i < separator2 ? incrementFifteen[i] : incrementSixteen[i - separator2]);
		toc();
		
		for(int i=0; i< incrementFifteen.length + incrementSixteen.length; i++) {
			if(firstArray[i] != secondArray[i]) {
				System.out.println("Problem at index " + i);
			}
		}
	}

}
