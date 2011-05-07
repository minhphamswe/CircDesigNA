package DnaDesign.test;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;

import java.util.Arrays;

import DnaDesign.DomainSequence;
import DnaDesign.DomainDefinitions;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
public class FoldingImplTest1 {
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		config.setMode(CircDesigNAConfig.RNA_MODE);
		for(int i = 0; i < 1; i++){
			FoldingImpl fl = new FoldingImpl(config);
			DomainSequence seqS = new DomainSequence();
			DomainDefinitions dsd = new DomainDefinitions(config);
			seqS.setDomains(0, null);
			int[][] domain = new int[1][];
			String seq = "GTGATAGACAC";
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
					domain[0][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
				}
			} else {
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = randomChoice(1,config.monomer.getNumMonomers()-1);
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			final double viaMatrix = fl.foldSingleStranded_viaMatrix(seqS, domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			System.out.println(viaMatrix);
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			//final double viaUnafold = fl.foldSingleStranded_viaUnafold(seqS, domain, domainMark);
			//System.out.println(Arrays.deepToString(domainMark));
			//System.out.println(seqLength+" "+-viaMatrix+" "+-viaUnafold);
		}
	}

	private static int randomChoice(int i, int j) {
		return ((int) (Math.random()*(j-i)))+i;
	}
}
