package info.quantlab.debug;

import net.finmath.stochastic.RandomVariable;

public class Debug {

	public static void ln() {
		System.out.println();
	}
	
	public static void log(String s) {
		System.out.print(s);
	}
	
	public static void logf(String s, Object ... args) {
		System.out.printf("[DEBUG] " + s, args);
		System.out.println();
	}

	public static void logln(String s) {
		System.out.println("[DEBUG] " + s);
	}
	
	public static void logVar(String varName, Object var) {
		System.out.printf("[DEBUG] %20s: " + var + "%n", varName);
	}
	
	public static int findNaNPath(RandomVariable rv) {
		if(!Double.isNaN(rv.getAverage()))
			return -1;
		
		if(rv.isDeterministic() && Double.isNaN(rv.doubleValue())) {
			logln("Deterministic Value is NaN!");
			return 0;
		}
			
		int path = 0;
		while(!Double.isNaN(rv.get(path)))
			path++;
		return path;
	}
	
	public static int getNumberOfNaNs(RandomVariable rv) {
		if(rv.isDeterministic() && Double.isNaN(rv.doubleValue())) {
			logln("Deterministic Value is NaN!");
			return 0;
		}
		RandomVariable isNan = rv.isNaN();
		return (int)(isNan.getAverage() * isNan.size());
	}
}
