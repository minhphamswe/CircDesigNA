package DnaDesign.test;

import circdesigna.energy.CircDesigNAMCSFolder;
import circdesigna.energy.NAFolding;
import edu.utexas.cssb.circdesigna.DomainSequence;
import DnaDesign.Config.CircDesigNAConfig;

public class SelfSimilarityScore {
	public static void main(String[] args){
		CircDesigNAConfig cfg = new CircDesigNAConfig();
		NAFolding fli = new CircDesigNAMCSFolder(cfg);
		while(true){
			int len = (int) (Math.random()*8000+10);
			int[][] domain = new int[1][len];
			int[][] nullMark = new int[1][len];
			for(int k = 0; k < domain.length; k++){
				for(int y = 0; y < domain[k].length; y++){
					domain[k][y] = (int) (Math.random()*4 + 1);
				}
			}
			DomainSequence ds1 = new DomainSequence();
			DomainSequence ds2 = new DomainSequence();
			ds1.setDomains(0,null);
			ds2.setDomains(0 | DomainSequence.DNA_COMPLEMENT_FLAG,null);
			System.out.println(len+" "+fli.mfeNoDiag(ds1, ds2, domain, nullMark));
		}
	}
}
