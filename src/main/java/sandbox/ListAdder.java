package sandbox;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class ListAdder {

	List<String> myList;
	
	public ListAdder(String[] stringArray) {
		myList = Arrays.stream(stringArray).map(str -> { return str; }).collect(Collectors.toList());
	}

	public void add(String strToAdd) {
		Stream<String> stringStream = Stream.concat(myList.stream(), Stream.of(strToAdd));
		myList = stringStream.collect(Collectors.toList());
	}
	
	public void print() {
		for(int i = 0; i < myList.size(); i++) {
			if(i==0) {
				System.out.print("{ ");
			} else
			{
				System.out.print(", ");
			}
			System.out.print(myList.get(i));
		}
		System.out.println(" }");
	}
	
	public static void main(String[] args) {
		System.out.println("Command Line Arguments:");
		for(int i = 0; i < args.length; i++) {
			System.out.print(i);
			System.out.print(": ");
			System.out.println(args[i]);
		}
		System.out.println("===============================================================================");
		
		String myArray[] = new String[10];
		myArray[0] = "DE";
		myArray[1] = "EN";
		myArray[2] = "FR";
		myArray[3] = "IT";
		myArray[4] = "BL";
		myArray[5] = "ES";
		myArray[6] = "SWE";
		myArray[7] = "FN";
		myArray[8] = "NO";
		myArray[9] = "BA";
		
		ListAdder tester = new ListAdder(myArray);
		
		tester.print();
		tester.add("Hamburch");
		tester.print();
	}

}
