package DnaDesign.test;

import DnaDesign.DomainDesigner;
import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.FoldingImpl;

public class FoldingImplTest1 {
	public static void main(String[] args){
		String seq = "CAACACAGGAACGTGATGTGGTGTGATG";
		FoldingImpl fl = new FoldingImpl();
		DomainSequence seqS = new DomainSequence();
		DomainStructureData dsd = new DomainStructureData();
		seqS.setDomains(0, dsd);
		int[][] domain = new int[1][];
		domain[0] = new int[seq.length()];
		int[][] domainMark= new int[1][seq.length()];
		for(int k = 0; k < seq.length(); k++){
			domain[0][k] = DomainDesigner.decodeConstraintChar(seq.charAt(k));
		}
		System.out.println(fl.foldSingleStranded(seqS, domain, domainMark));
	}
}
