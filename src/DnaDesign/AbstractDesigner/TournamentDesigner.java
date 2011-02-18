package DnaDesign.AbstractDesigner;

import java.util.TreeMap;

import DnaDesign.DomainDesigner.ScorePenalty;
import DnaDesign.impl.DomainDesignPMemberImpl;

public abstract class TournamentDesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	public TournamentDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}

	public void tournamentSelect(int numReproduce){
		//Seed the fittest
		TreeMap<Double, T> populationView = new TreeMap();
		for(T q : population_mutable){
			double score = SingleDesigner.getOverallScore(q);
			populationView.put(-score, q); //sort descending
		}
		T fittest = populationView.get(populationView.lastKey());
		for(int k = 0; k < numReproduce && populationView.size()>=2; k++){
			T bottom = populationView.remove(populationView.firstKey());
			T replace = populationView.remove(populationView.lastKey());
			if (replace==bottom){
				throw new RuntimeException("Assertion error");
			}
			bottom.seedFromOther(replace);

			if (ASSERT_SCOREFUNC_ISOLATION){
				DomainDesignPMemberImpl q = (DomainDesignPMemberImpl)(Object)replace;
				//Check: before mutation, all penalties should have change in score of 0.
				for(ScorePenalty s : q.penalties){
					if (s.evalScore(q.domain,q.domain_markings)!=0){
						throw new RuntimeException("FAIL");
					}
				}	
				q = (DomainDesignPMemberImpl)(Object)bottom;
				//Check: before mutation, all penalties should have change in score of 0.
				for(ScorePenalty s : q.penalties){
					if (s.evalScore(q.domain,q.domain_markings)!=0){
						throw new RuntimeException("FAIL");
					}
				}
			}
		}
		setBestChild(fittest);
	}
}
