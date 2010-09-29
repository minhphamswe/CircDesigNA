package DnaDesign;

import java.util.List;

/**
 * Describes objects which support these 3 scoring functions:
 * Pairscore: the likelyhood that seq1 is hybridized with seq2 at equilibrium
 * FoldsingleStranded: the amount of secondary structure exhibited by seq1.
 * AffectedSequenceInvalidScore: Some grab-bag for miscellaneous sequence constraints. GGGG and ATATAT may be covered here.
 * 
 * getLongestHelixLength: returns the length of the longest helix existing in the last secondary structure calculated. 
 * getNumBasesPaired: returns the number of paired bases in the last secondary structure calculated.
 * The above two methods are only allowed to be called after first calling Pairscore or foldSingleStranded.
 */
public interface NAFolding {
	public double pairscore(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] problemAreas);
	public double affectedSequenceInvalidScore(int i, List<DomainSequence> seqs, int[][] domain, int[][] domain_markings);
	public double foldSingleStranded(DomainSequence domainSequence, int[][] domain, int[][] domain_markings);
	public int getLongestHelixLength();
	public int getNumBasesPaired();
}
