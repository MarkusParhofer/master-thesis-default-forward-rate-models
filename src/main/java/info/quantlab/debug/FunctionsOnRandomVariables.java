package info.quantlab.debug;

import net.finmath.stochastic.RandomVariable;

public class FunctionsOnRandomVariables {

	private FunctionsOnRandomVariables() {
		// TODO Auto-generated constructor stub
	}
	
	public static int findFirstPathWithValue(RandomVariable rv, double value) {
		if(rv.isDeterministic() && rv.doubleValue() == value) {
			// Debug.logln("Value is deterministic value of RV!");
			return 0;
		} else if(rv.isDeterministic()) {
			// Debug.logln("Value is deterministic value of RV!");
			return -1;
		}
			
		int path = 0;
		while((path < rv.size()) && !(rv.get(path) == value))
			path++;
		return path == rv.size() ? -1 : path - 1;
	}

	public static int findNaNPath(RandomVariable rv) {
		if(!Double.isNaN(rv.getAverage()))
			return -1;
		
		if(rv.isDeterministic() && Double.isNaN(rv.doubleValue())) {
			Debug.logln("Deterministic Value is NaN!");
			return 0;
		}
			
		int path = 0;
		while(!Double.isNaN(rv.get(path)))
			path++;
		return path;
	}
	
	public static int getNumberOfNaNs(RandomVariable rv) {
		if(rv.isDeterministic() && Double.isNaN(rv.doubleValue())) {
			Debug.logln("Deterministic Value is NaN!");
			return 0;
		}
		RandomVariable isNan = rv.isNaN();
		return (int)(isNan.getAverage() * isNan.size());
	}
}
