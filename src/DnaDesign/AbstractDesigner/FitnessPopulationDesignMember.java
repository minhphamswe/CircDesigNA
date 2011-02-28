package DnaDesign.AbstractDesigner;

import java.util.Comparator;

public class FitnessPopulationDesignMember <T extends PopulationDesignMember<T>> implements Comparable <FitnessPopulationDesignMember<T>>{

	public T myKey;
	public double myScore;
	public int compareTo(FitnessPopulationDesignMember<T> o) {
		double diff = myScore - o.myScore;
		if (diff==0){
			return 0;
		} else if (diff < 0){
			return -1;
		} else {
			return 1;
		}
	}

}
