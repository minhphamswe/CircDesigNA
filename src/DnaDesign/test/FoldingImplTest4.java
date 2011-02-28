package DnaDesign.test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Scanner;

import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.FoldingImpl;

public class FoldingImplTest4 {
	public static void main(String[] args) throws FileNotFoundException{
		File outFile = new File("output.txt");
		CircDesigNAConfig config = new CircDesigNAConfig();
		System.out.println(outFile.getAbsolutePath());
		System.setOut(new PrintStream(new FileOutputStream(outFile)));
		//Error in this file. third column represents pairing of seq1 with seq1! Embarrasing.
		Scanner in = new Scanner(new File("C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\Circ_DesigNAPaper\\AssemblaRepo\\circdesignapaper_w_figures\\scoresComparison_SS.txt"));
		in.nextLine();
		FoldingImpl fl = new FoldingImpl(config);
		DomainStructureData dsd = new DomainStructureData(config);
		int[][] domain = new int[2][];
		int[][] domain_markings = new int[2][];

		DomainSequence ds1 = new DomainSequence();
		ds1.setDomains(0, null);
		DomainSequence ds2 = new DomainSequence();
		ds2.setDomains(1, null);
		while(in.hasNextLine()){
			String line2 = in.nextLine();
			String[] line = line2.split(" ");
			int seqLength = line[1].length();
			for(int j = 0; j < 1; j++){
				domain[j] = new int[seqLength];
				domain_markings[j] = new int[seqLength];
				for(int k = 0; k < seqLength; k++){
					domain[j][k] = config.monomer.decodeConstraintChar(line[j+1].charAt(k));
				}
			}
			//double result = fl.mfeHybridDeltaG_viaUnafold(ds1, ds2, domain, domain_markings);
			double result = fl.foldSingleStranded_viaUnafold(ds1, domain, domain_markings);
			System.out.println(line2+" "+(-result));
		}
		System.out.flush();
	}
}
