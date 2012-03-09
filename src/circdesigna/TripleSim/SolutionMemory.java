package circdesigna.TripleSim;

import java.util.ArrayList;

public class SolutionMemory {
	public static class Solution {
		public Solution(int n) {
			values = new double[n];
			f = new double[n];
		}
		public double[] values;
		public double[] f;
		public double time;
	}
	private ArrayList<Solution> solutions = new ArrayList();
	private Solution filling;
	private int N;
	private int numBuffers;
	public SolutionMemory(int N, int numBuffers){
		this.N = N;
		this.numBuffers= numBuffers;
	}
	public Solution getBuffer(){
		if (filling==null){
			filling = new Solution(N);
		}
		return filling;
	}
	private double latticeSpacing;
	public int size(double h) {
		if (h!=latticeSpacing){
			return 1;
		}
		return solutions.size();
	}
	public void dedicateBuffer(double step){
		solutions.add(filling);
		filling = null;
		if (solutions.size() > numBuffers){
			filling = solutions.remove(0);
		}
		if (step != latticeSpacing){
			latticeSpacing = step;
			//Remove all points not on the new lattice.
			reposition();
		}
	}
	public Solution getSolution(int i) {
		return solutions.get(solutions.size() - 1 + i);
	}
	/**
	 * Removes all memory of the past history, except for the current solution
	 */
	public void reposition() {
		while(solutions.size() >= 2){
			solutions.remove(0);
		}
	}
}
