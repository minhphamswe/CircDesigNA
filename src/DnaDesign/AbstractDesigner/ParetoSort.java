package DnaDesign.AbstractDesigner;

import java.util.ArrayList;

public abstract class ParetoSort <T extends PopulationDesignMember<T>>{
	private boolean[][] dominatedMat;
	private ArrayList<Integer> front = new ArrayList();
	
	/**
	 * Forms at most fnum pareto fronts. Adds penalty * (fi-1) to each item, where fi is its pareto ranking, fi <= fnum.
	 * 
	 */
	public void adjustFitness(FitnessPopulationDesignMember<T>[] in, int s, int e, double penalty, int fnum){
		if (dominatedMat==null || dominatedMat.length < in.length){
			dominatedMat = new boolean[in.length][in.length];
		}
		for(int y = 0; y < dominatedMat.length; y++){
			for(int x = 0; x < dominatedMat.length; x++){
				if (y==x){
					continue;
				}
				boolean isDominated = isDominatedBy(in[y].myKey,in[x].myKey);
				dominatedMat[y][x] = isDominated;
			}
		}
		for(int fi = 0; fi < fnum; fi++){
			front.clear();
			for(int y = 0; y < dominatedMat.length; y++){
				boolean undominated = true;
				for(int x = 0; x < dominatedMat.length; x++){
					if (dominatedMat[y][x]){
						undominated = false;
						break;
					}
				}
				if (!undominated){
					in[y].myScore += penalty;
					front.add(y);
				}
			}
			//Remove the front
			for(Integer y : front){
				for(int x = 0; x < dominatedMat.length; x++){
					dominatedMat[x][y] = false;
				}
			}
		}
	}

	public abstract boolean isDominatedBy(T t, T t2);
}
