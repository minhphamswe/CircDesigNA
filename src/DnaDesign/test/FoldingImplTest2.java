package DnaDesign.test;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;
import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
public class FoldingImplTest2 {
	/**
	 * Tests base pairing observables calculations
	 */
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		for(int i = 0; i < 1; i++){
			FoldingImpl fl = new FoldingImpl(config);
			DomainSequence seqS = new DomainSequence();
			DomainStructureData dsd = new DomainStructureData(config);
			seqS.setDomains(0, null);
			int[][] domain = new int[1][];
			String seq = null;// "AAAAATTTTTTGGGGGGCCCCCCCC";
			int seqLength;
			if (seq==null){
				seqLength = (int)(Math.random()*40+10);
			} else {
				seqLength = seq.length();
			}
			domain[0] = new int[seqLength];
			int[][] domainMark= new int[1][seqLength];
			if (seq!=null){
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
				}
			} else {
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = randomChoice(A,C,G,T);
				}
			}

			long now;
			double[][] basePairs;
			if (true){
				int N = seqLength + seqLength;
				basePairs = new double[N][N+1];
				now = System.nanoTime();
				fl.pairPrHybrid(basePairs, seqS, seqS, domain);
			} else {
				int N = seqLength;
				basePairs = new double[N][N+1];
				now = System.nanoTime();
				fl.pairPrSS(basePairs, seqS, domain);
			}
			long time = System.nanoTime()-now;
			System.out.println(seqLength+" "+(time/1e9));
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
