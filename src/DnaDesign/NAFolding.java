package DnaDesign;

import java.util.List;

/**
 * Describes objects which support many useful nucleic acid folding routines
 * 
 * mfeHybridDeltaG: Locate the minimum free energy structure formed by the interaction of 2 strands ("hybridization energy"), and
 * return its delta G.
 * 
 * mfeSSDeltaG: Locate the minimum free energy conformation of a single nucleic acid strand ("secondary structure energy"), and
 * return its delta G.
 * 
 * Derived values from the above two routines:
 * getLongestHelixLength: returns the length of the longest helix existing in the last MFE calculated. 
 * getNumBasesPaired: returns the number of paired bases in the last MFE calculated.
 * The above two methods are only allowed to be called after first calling Pairscore or foldSingleStranded.
 * 
 * AffectedSequenceInvalidScore: Some grab-bag for miscellaneous sequence constraints. Only the sequences in the provided list
 * are analyzed. For example, penalties for unwanted sequences (GGGG and ATATAT) may be covered here.
 * This score is not required to be in any unit, it should be used as a nonrigorous design penalty. However, it should be
 * used, because many of the routines in this library use parameters derived from highly sequences with good G/C/A/T mixes. An
 * inability to evaluate whether a design is good is a sign that another solution should be sought after.
 * 
 * pairPrHybrid: Fills in an appropriately sized double array of base pairing probabilities. After running,
 * the array[i][j] contains the probability that i of strand #1 is paired with j of strand #2. For convenience, the entire matrix will be
 * filled out and zeroed, up to the length of the two sequences (so the input matrix can be larger)
 * 
 * pairPrSS: See pairPrHybrid, but used for secondary structure formation. array[i][j] contains the probability
 * that i is paired with j
 */
public interface NAFolding {
	public double mfeHybridDeltaG(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] problemAreas);
	public double affectedSequenceInvalidScore(int i, List<DomainSequence> seqs, int[][] domain, int[][] domain_markings);
	public double mfeSSDeltaG(DomainSequence domainSequence, int[][] domain, int[][] domain_markings);
	public int getLongestHelixLength();
	public int getNumBasesPaired();
	public void pairPrHybrid(double[][] pairsOut, DomainSequence seq1, DomainSequence seq2, int[][] domain);
	public void pairPrSS(double[][] pairsOut, DomainSequence seq1, int[][] domain);
}