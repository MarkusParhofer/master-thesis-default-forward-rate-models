package info.quantlab.debug;

public class Time {


	private static long startTime = 0;
	
	public static void tic() {
		startTime = System.nanoTime();
	}
	
	public static long toc(boolean printResult) {
		long result = System.nanoTime();
		result -= startTime;
		if(printResult)
			printTimeResult(result);
		return result;
	}
	
	public static long toc() {
		return toc(true);
	}
	
	public static void printTimeResult(long result) {
		System.out.println("\n\nElapsed Time is " + (result / 1000) + " Microseconds\n");
	}
	
}
