package circdesigna.TripleSim;

import java.util.ArrayList;
import java.util.List;

import circdesigna.DomainPolymerGraph;

public class ConnectedComplex {
	public static ArrayList<ArrayList<int[]>> getConnectedComponents(DomainPolymerGraph q) {
		ArrayList<ArrayList<int[]>> connected = new ArrayList<ArrayList<int[]>>();
		List<Integer> rots = q.getStrandRotations();
		rots.add(0,0);
		int[][] strands = new int[rots.size()][];
		rots.add(q.length());
		for(int i = 0; i < strands.length; i++){
			strands[i] = new int[]{rots.get(i),rots.get(i+1)};
		}
		for(int[] e : strands){
			ArrayList<int[]> placed = null;
			for(int inBin = 0; inBin < connected.size(); inBin++){
				ArrayList<int[]> bin = connected.get(inBin);
				for(int[] other : bin){
					if (placed==null){
						if (connected(q,other,e)){
							bin.add(e);
							placed = bin;
							break;
						}
					} else {
						if (connected(q,other,e)){
							placed.addAll(bin);
							connected.remove(inBin--);
							break;
						}
					}
				}
			}
			if (placed == null){
				ArrayList<int[]> neu = new ArrayList();
				neu.add(e);
				connected.add(neu);
			}
		}
		return connected;
	}
	private static boolean connected(DomainPolymerGraph q, int[] input, int[] output) {
		for(int k = input[0]; k < input[1]; k++){
			int target = q.getDomainPair(k);
			if (target >= output[0] && target < output[1]){
				return true;
			}
		}
		return false;
	}

}
