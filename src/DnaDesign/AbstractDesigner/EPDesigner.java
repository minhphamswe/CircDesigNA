package DnaDesign.AbstractDesigner;

import java.util.TreeMap;

import DnaDesign.DomainDesigner;

public class EPDesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	public EPDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}
	
	public void runBlockIteration_ (DomainDesigner runner, double endThreshold) {
		TreeMap<Double, T> sorted = new TreeMap();
		for(int i = 0; i < populationSize; i++){
			T toMutate = population_mutable[i];
			SingleDesigner.mutateAndTestAndBackup(toMutate);
			sorted.put(SingleDesigner.getOverallScore(toMutate),toMutate);
		};
		setBestChild(sorted.get(sorted.firstKey()));
	}
}
