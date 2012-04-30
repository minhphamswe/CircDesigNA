package circdesigna.energy;

import circdesigna.DomainSequence;

public interface ConstraintsNAFolding extends NAFolding{
	public double mfe(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing);
	public double mfe(DomainSequence seq, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing);
}
