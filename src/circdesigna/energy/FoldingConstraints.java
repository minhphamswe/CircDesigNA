package circdesigna.energy;

public class FoldingConstraints {
	public FoldingConstraints(int N) {
		if (N == 0){
			return;
		}
		
		preventPairing = new boolean[N][N];
	}
	/**
	 * Symmetric boolean matrix, if index i,j is true, no structure should be considered which 
	 * pairs bases i and j.
	 */
	public boolean[][] preventPairing;
	public boolean preventPairing(int i, int j) {
		if (preventPairing == null){
			return false;
		}
		return preventPairing[i][j];
	}
}
