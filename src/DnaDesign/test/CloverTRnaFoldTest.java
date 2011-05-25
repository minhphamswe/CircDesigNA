package DnaDesign.test;

import circdesigna.energy.CircDesigNAMCSFolder;
import circdesigna.energy.NAFolding;
import edu.utexas.cssb.circdesigna.DomainDefinitions;
import edu.utexas.cssb.circdesigna.DomainDesigner;
import edu.utexas.cssb.circdesigna.DomainSequence;
import DnaDesign.Config.CircDesigNAConfig;

public class CloverTRnaFoldTest {
	//The following sequence should fold into a clover.
	//AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT
	//Correct duplexes:
	//1-4x
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		NAFolding na = new CircDesigNAMCSFolder(config);
		String seq = "AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT";
		int[][] domain = new int[1][seq.length()];
		int[][] domainMark= new int[1][seq.length()];
		for(int k = 0; k < seq.length(); k++){
			domain[0][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
		}
		DomainDefinitions dsd = new DomainDefinitions(config);
		DomainSequence ds = new DomainSequence();
		ds.setDomains(0, null);
		System.out.println(na.mfe(ds, domain, domainMark));
	}
}
