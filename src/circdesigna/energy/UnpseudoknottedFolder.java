package circdesigna.energy;

import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

/**
 * Uses an N^3 DP algorithm to compute the MFE among all unpseudoknotted folded structures
 * of one or two sequences.
 */
public class UnpseudoknottedFolder extends CircDesigNASystemElement implements NAFolding{
	public UnpseudoknottedFolder(CircDesigNAConfig System) {
		super(System);
	}

	public double mfe(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] domain_markings) {
		
		return 0;
	}

	public double mfe(DomainSequence domainSequence, int[][] domain,
			int[][] domain_markings) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double mfeNoDiag(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double mfeStraight(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings, int markLeft, int markRight, int joffset) {
		// TODO Auto-generated method stub
		return 0;
	}
}
