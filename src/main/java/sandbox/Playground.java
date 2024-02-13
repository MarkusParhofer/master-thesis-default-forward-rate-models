package sandbox;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

import info.quantlab.debug.Time;
import info.quantlab.easyplot.EasyPlot2D;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.*;
import net.finmath.opencl.montecarlo.RandomVariableOpenCLFactory;
import net.finmath.plots.Named;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class Playground extends Time {

	public static boolean useGPU = false;
	
	public static void main(String[] args) {
		testForLoopWithOuterVar();
	}

	public static void testForLoopWithOuterVar() {
		int outerVar = 0;
		int anotherVar = 0;
		tic();
		for (; outerVar < 100000; outerVar++, anotherVar++) {
			System.out.println("Outer Variable is " + outerVar);
		}
		long firstRes = toc();
		System.out.println("\nOut of loop now!");
		System.out.println("Outer Variable is " + outerVar);
		System.out.println("Another Variable is " + anotherVar);
		System.out.println();
		tic();
		anotherVar = 100000;
		outerVar = 0;
		for (; outerVar < 100000; outerVar++) {
			System.out.println("Outer Variable is " + outerVar);
		}
		System.out.println("\nOut of loop now!");
		System.out.println("Outer Variable is " + outerVar);
		System.out.println("Another Variable is " + anotherVar);
		printTimeResult(firstRes);
		toc();
	}

	public static void testGPURV() {
		int numberOfPathsToPlot = 1;
		RandomVariableFactory factory = null;
		if(useGPU)
			factory = new RandomVariableOpenCLFactory();
		else
			factory = new RandomVariableFromArrayFactory();
		TimeDiscretization td = new TimeDiscretizationFromArray(0.0, 1000, 0.01);
		tic();
		BrownianMotion bm = new BrownianMotionFromMersenneRandomNumbers(td, 1, 250000, 10, factory);
		double sigma = 0.2, mu = 0.1, S0 = 12.0;
		RandomVariable[] path = new RandomVariable[1001];
		path[0] = factory.createRandomVariable(S0);
		for(int timeIndex = 1; timeIndex < 1001; timeIndex++) {
			RandomVariable diffusion = bm.getIncrement(timeIndex - 1, 0).mult(sigma);
			RandomVariable drift = factory.createRandomVariable((mu - sigma * sigma) * (td.getTime(timeIndex) - td.getTime(timeIndex - 1)));
			path[timeIndex] = path[timeIndex - 1].mult(drift.add(diffusion).exp());
		}
		toc();
		DoubleUnaryOperator[] fArray = new DoubleUnaryOperator[numberOfPathsToPlot];

		final DoubleToIntFunction getTimeIndex = operand -> {
			int timeIndex = td.getTimeIndex(operand);
			return timeIndex < 0 ? - timeIndex - 2 : timeIndex;
		};

		for(int i=0; i < numberOfPathsToPlot; i++) {
			final int pathIndex = i;
			fArray[i] = (operand) -> {
				return path[getTimeIndex.applyAsInt(operand)].getAverage();
			};
		}
		List<DoubleUnaryOperator> fList = Arrays.asList(fArray);
		List<Named<DoubleUnaryOperator>> fListNamedNonDefModel =
				fList.stream().map(operator -> new Named<>(
						"Path " + fList.indexOf(operator), operator)).collect(Collectors.toList());
		EasyPlot2D plotPaths = new EasyPlot2D(td.getTime(0), td.getLastTime(), 51, fListNamedNonDefModel);
		plotPaths.setTitle("Sample Paths");
		plotPaths.show();
	}
}
