package DnaDesign.Batch;

import java.util.Scanner;

public class StringSplitter {
	/**
	 * Used for using NUPACK as a designer.
	 */
	public static void main(String[] args){
		Scanner in = new Scanner(System.in);
		String line = in.nextLine();
		int numWrap = new Integer(in.nextLine().split("\\s+")[0]);
		for(int k = 0; k < line.length();){
			for(int y = 0; y < numWrap; y++, k++){
				System.out.print(line.charAt(k));
			}
			if (k < line.length()){
				System.out.println();
			}
		}
	}
}
