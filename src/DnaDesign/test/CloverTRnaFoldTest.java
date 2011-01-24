package DnaDesign.test;

import DnaDesign.DomainDesigner;
import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.NAFolding;
import DnaDesign.impl.FoldingImpl;

public class CloverTRnaFoldTest {
	//The following sequence should fold into a clover.
	//AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT
	//Correct duplexes:
	//1-4x
	public static void main(String[] args){
		NAFolding na = new FoldingImpl();
		FoldingImpl.DEBUG_selfCrosstalkMethod = true;
		String seq = "AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT";
		int[][] domain = new int[1][seq.length()];
		int[][] domainMark= new int[1][seq.length()];
		for(int k = 0; k < seq.length(); k++){
			domain[0][k] = DomainDesigner.decodeConstraintChar(seq.charAt(k));
		}
		DomainStructureData dsd = new DomainStructureData();
		DomainSequence ds = new DomainSequence();
		ds.setDomains(0, null);
		System.out.println(na.mfeSSDeltaG(ds, domain, domainMark));
	}
}
