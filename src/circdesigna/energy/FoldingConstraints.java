package circdesigna.energy;

public class FoldingConstraints {
	public FoldingConstraints(int N) {
		preventPairing = new boolean[N][N];
	}
	/**
	 * Symmetric boolean matrix, if index i,j is true, no structure should be considered which 
	 * pairs bases i and j.
	 */
	public final boolean[][] preventPairing;
}
