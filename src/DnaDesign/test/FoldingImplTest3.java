package DnaDesign.test;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;

import java.util.Arrays;

import DnaDesign.AbstractDomainDesignTarget;
import DnaDesign.DomainSequence;
import DnaDesign.DomainDefinitions;
import DnaDesign.AbstractDomainDesignTarget.HairpinClosingTarget;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
import DnaDesign.impl.DomainDesignerImpl.HairpinOpening;
public class FoldingImplTest3 {
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		for(int i = 0; i < 1; i++){
			FoldingImpl fl = new FoldingImpl(config);
			DomainDesignerImpl impl = new DomainDesignerImpl(fl,config);
			DomainDefinitions dsd = new DomainDefinitions(config);
			HairpinClosingTarget hairpin = new AbstractDomainDesignTarget(dsd,config).new HairpinClosingTarget(1,0,0|DomainSequence.DNA_COMPLEMENT_FLAG,2,true,null);
			HairpinOpening hairpinOpening = impl.new HairpinOpening(hairpin, null);
			DomainSequence seqS = new DomainSequence();
			seqS.setDomains(0, null);
			String[] seq = new String[3];
			int[][] domain = new int[seq.length][];
			int[][] domainMark= new int[seq.length][];
			seq[0] = "GTGATAGACAC";
			seq[1] = "CTAT";
			seq[2] = "ATAGCATAG";
			for(int j = 0; j < seq.length; j++){
				int seqLength;
				if (seq[j]==null){
					seqLength = (int)(Math.random()*10+10);
				} else {
					seqLength = seq[j].length();
				}
				domain[j] = new int[seqLength];
				domainMark[j] = new int[seqLength];
				if (seq!=null){
					for(int k = 0; k < seqLength; k++){
						domain[j][k] = config.monomer.decodeConstraintChar(seq[j].charAt(k));
					}
				} else {
					for(int k = 0; k < seqLength; k++){
						domain[j][k] = randomChoice(A,C,G,T);
					}
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			final double viaMatrix = hairpinOpening.evalScoreSub(domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			System.out.println(viaMatrix);
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			
			//final double viaUnafold = fl.foldSingleStranded_viaUnafold(seqS, domain, domainMark);
			//System.out.println(Arrays.deepToString(domainMark));
			//System.out.println(seqLength+" "+-viaMatrix+" "+-viaUnafold);
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
