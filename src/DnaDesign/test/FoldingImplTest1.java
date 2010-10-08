package DnaDesign.test;

import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.T;

import java.util.Arrays;

import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
public class FoldingImplTest1 {
	public static void main(String[] args){
		for(int i = 0; i < 1; i++){
			FoldingImpl fl = new FoldingImpl();
			DomainSequence seqS = new DomainSequence();
			DomainStructureData dsd = new DomainStructureData();
			seqS.setDomains(0, dsd);
			int[][] domain = new int[1][];
			String seq = "GTGATAGA";
			int seqLength;
			if (seq==null){
				seqLength = (int)(Math.random()*10+10);
			} else {
				seqLength = seq.length();
			}
			domain[0] = new int[seqLength];
			int[][] domainMark= new int[1][seqLength];
			if (seq!=null){
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = DomainDesignerImpl.decodeConstraintChar(seq.charAt(k));
				}
			} else {
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = randomChoice(A,C,G,T);
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			final double viaMatrix = fl.foldSingleStranded_viaMatrix(seqS, domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			final double viaUnafold = fl.foldSingleStranded_viaUnafold(seqS, domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			System.out.println(seqLength+" "+-viaMatrix+" "+-viaUnafold);
		}
	}

	private static int randomChoice(int a, int b, int c, int d) {
		switch((int)(Math.random()*4)){
		case 0:return a;
		case 1:return b;
		case 2:return c;
		case 3:return d;
		}
		throw new RuntimeException("Assertion error: randomChoice");
	}
}
