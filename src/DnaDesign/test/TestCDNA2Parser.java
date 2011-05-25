package DnaDesign.test;

import java.util.Scanner;

import DnaDesign.parser.CDNA2PublicParser;

public class TestCDNA2Parser {
	public static void main(String[] args){
		Scanner in = new Scanner(System.in);
		while(in.hasNextLine()){
			System.out.println(CDNA2PublicParser.parse(in.nextLine()));
		}
	}
}
