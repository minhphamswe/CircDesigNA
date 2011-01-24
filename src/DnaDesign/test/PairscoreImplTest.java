package DnaDesign.test;

import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.T;

import java.util.Arrays;

import DnaDesign.DnaDefinition;
import DnaDesign.DomainDesigner;
import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.impl.FoldingImpl;
public class PairscoreImplTest {
	public static void main(String[] args){
		for(int i = 0; i < 1000; i++){
			FoldingImpl fl = new FoldingImpl();
			DomainSequence seqS = new DomainSequence();
			DomainSequence seq2S = new DomainSequence();
			DomainStructureData dsd = new DomainStructureData();
			seqS.setDomains(0, null);
			seq2S.setDomains(1, null);
			String Seq1 = "ATGCATGC";//"GTGTTCTTGA";
			String Seq2 = "TACGTACG";//"GAGAGGTGGAA";
			int seqLength1 = (int)(Math.random()*10+10);
			int seqLength2 = (int)(Math.random()*10+10);
			if (Seq1!=null){
				seqLength1 = Seq1.length();
			}
			if (Seq2!=null){
				seqLength2 = Seq2.length();
			}
			int[][] domain = new int[2][];
			int[][] domainMark= new int[2][];
			int[][] domainMark2= new int[2][];
			domain[0] = new int[seqLength1];
			domain[1] = new int[seqLength2];
			domainMark[0] = new int[seqLength1];
			domainMark[1] = new int[seqLength2];
			domainMark2[0] = new int[seqLength1];
			domainMark2[1] = new int[seqLength2];
			for(int k = 0; k < seqLength1; k++){
				if (Seq1==null){
					domain[0][k] = randomChoice(A,C,G,T);
				} else {
					domain[0][k] = DomainDesigner.decodeConstraintChar(Seq1.charAt(k));
				}
			}
			for(int k = 0; k < seqLength2; k++){
				if (Seq2==null){
					domain[1][k] = randomChoice(A,C,G,T);
				} else {
					domain[1][k] = DomainDesigner.decodeConstraintChar(Seq2.charAt(k));
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			for(int k = 0; k < domainMark2.length; k++)Arrays.fill(domainMark2[k],0);
			final double viaMatrix = fl.mfeHybridDeltaG_viaMatrix(seqS, seq2S, domain, domainMark);
			final double viaUnafold = fl.mfeHybridDeltaG_viaUnafold(seqS, seq2S, domain, domainMark2);
			//System.out.println(Arrays.deepToString(domainMark));
			if ((viaUnafold-viaMatrix)>2){
				System.out.println(Arrays.deepToString(domainMark));
				System.out.println(Arrays.deepToString(domainMark2));
				for(int[] domainS : domain){
					for(int k : domainS){
						System.out.print(DnaDefinition.displayBase(k));
					}
					System.out.println();
				}
			}
			System.out.println(seqLength1+" "+seqLength2+" "+-viaMatrix+" "+-viaUnafold);
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
