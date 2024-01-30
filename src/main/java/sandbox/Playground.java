package sandbox;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

import info.quantlab.easyplot.EasyPlot2D;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.*;
import net.finmath.opencl.montecarlo.RandomVariableOpenCLFactory;
import net.finmath.plots.Named;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class Playground {

	public static boolean useGPU = false;
	private static long startTime;
	
	public static void tic() {
		startTime = System.nanoTime();
	}
	
	public static long toc() {
		long result = System.nanoTime();
		result -= startTime;
		System.out.println("Elapsed Time is: " + result/1000000 + " Milliseconds");
		return result;
	}
	
	public static void main(String[] args) {
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
