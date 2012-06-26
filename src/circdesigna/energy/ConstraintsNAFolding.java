package circdesigna.energy;

import circdesigna.GeneralizedInteractiveRegion;

public interface ConstraintsNAFolding extends NAFolding{
	public double mfe(GeneralizedInteractiveRegion seq1, GeneralizedInteractiveRegion seq2, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing);
	public double mfe(GeneralizedInteractiveRegion seq, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing);
}
