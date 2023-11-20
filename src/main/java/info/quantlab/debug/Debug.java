package info.quantlab.debug;


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
	
}
