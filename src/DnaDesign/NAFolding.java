package DnaDesign;

import java.util.List;

/**
 * Describes objects which support these 3 scoring functions:
 * Pairscore: the likelyhood that seq1 is hybridized with seq2 at equilibrium
 * FoldsingleStranded: the amount of secondary structure exhibited by seq1.
 * AffectedSequenceInvalidScore: Some grab-bag for miscellaneous sequence constraints. GGGG and ATATAT may be covered here.
 */
public interface NAFolding {
	public double pairscore(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] problemAreas);
	public double affectedSequenceInvalidScore(int i, List<DomainSequence> seqs, int[][] domain, int[][] domain_markings);
	public double foldSingleStranded(DomainSequence domainSequence, int[][] domain, int[][] domain_markings);
}
