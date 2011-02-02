package DnaDesign.test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Scanner;

import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;

public class FoldingImplTest4 {
	public static void main(String[] args) throws FileNotFoundException{
		File outFile = new File("output.txt");
		System.out.println(outFile.getAbsolutePath());
		System.setOut(new PrintStream(new FileOutputStream(outFile)));
		Scanner in = new Scanner(new File("C:\\Users\\Benjamin\\PROGRAMMING\\Personal\\DNAWork\\Comparison\\NupackMfoldCircDesign.txt"));
		in.nextLine();
		FoldingImpl fl = new FoldingImpl();
		DomainStructureData dsd = new DomainStructureData();
		int[][] domain = new int[2][];
		int[][] domain_markings = new int[2][];

		DomainSequence ds1 = new DomainSequence();
		ds1.setDomains(0, null);
		DomainSequence ds2 = new DomainSequence();
		ds2.setDomains(1, null);
		while(in.hasNextLine()){
			String line2 = in.nextLine();
			String[] line = line2.split(" ");
			int seqLength = line[2].length();
			for(int j = 0; j < 2; j++){
				domain[j] = new int[seqLength];
				domain_markings[j] = new int[seqLength];
				for(int k = 0; k < seqLength; k++){
					domain[j][k] = DomainDesignerImpl.decodeConstraintChar(line[j+2].charAt(k));
				}
			}
			double result = fl.mfeHybridDeltaG_viaUnafold(ds1, ds2, domain, domain_markings);
			System.out.println(line2+" "+(-result));
		}
		System.out.flush();
	}
}
