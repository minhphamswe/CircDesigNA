package DnaDesign.impl;

import static DnaDesign.DomainDesigner.int_urn;
import static DnaDesign.DomainSequence.DNAMARKER_DONTMUTATE;

import java.util.Arrays;

import DnaDesign.DesignerCode;
import DnaDesign.DomainDesigner;
import DnaDesign.DomainStructureData;
import DnaDesign.AbstractDesigner.SingleMemberDesigner;
import DnaDesign.DomainDesigner.ScorePenalty;

/**
 * Implementation of BlockDesigner
 */
public class DomainDesignBlockDesignerImpl extends SingleMemberDesigner<DomainDesignPMemberImpl>{
	public DomainDesignBlockDesignerImpl(int[] mutableDomains, DesignerCode[] mutators, DomainStructureData dsd, DomainDesigner domainDesigner, DomainDesignPMemberImpl backupPMember){
		//Temporary backup member for performing reversions
		this.defaultBackupCache = backupPMember;
		this.dsd = dsd;
		this.mutators = mutators;
		this.dd = domainDesigner;
		this.mutableDomains = mutableDomains;
		mutation_shared = new Mutation();
	}
	private DesignerCode[] mutators;
	private DomainDesignPMemberImpl defaultBackupCache;
	private Mutation mutation_shared;
	int[] mutableDomains;
	private DomainDesigner dd;
	private DomainStructureData dsd;

	public double getOverallScore(DomainDesignPMemberImpl q) {
		double current_score = 0;
		for(ScorePenalty s : q.penalties){
			//Sanity check
			if (s.old_score!=s.cur_score){
				throw new RuntimeException("Get overallScore called in the middle of mutation");
			}
			current_score += s.old_score; 
		}
		return current_score;
	}
	
	private class Mutation {
		public int[] mut_domains = new int[dd.MAX_MUTATION_DOMAINS];
		boolean[] domain_mutated;
		int num_domains_mut;
		boolean revert_mutation, newPointReached;
		
		public void Mutate(DomainDesignPMemberImpl q, DomainDesignPMemberImpl backup, boolean fullBackup){
			domain_mutated = new boolean[q.domain.length];
			num_domains_mut = int_urn(1,Math.min(dd.MAX_MUTATION_DOMAINS,mutableDomains.length));
			Arrays.fill(domain_mutated,false);
			
			//Select random domains to mutate
			for(int k = 0; k < num_domains_mut; k++){
				int mut_domain_search = int_urn(0, mutableDomains.length-1);
				int mut_domain = 0;
				for(int o = 0; o < mutableDomains.length; o++){
					mut_domain = mutableDomains[(mut_domain_search+o)%mutableDomains.length];
					if (domain_mutated[mut_domain]){
						continue;
					}
					break;
				}
				if(domain_mutated[mut_domain]){
					//num_domains_mut--;
					k--;
					continue;
				}
				
				/*
				if (timesSameCount > 0){
					max_mutations = Math.max(MAX_MUTATIONS_LIMIT/timesSameCount,min_mutations);
				}
				*/

				mut_domains[k] = mut_domain;
				domain_mutated[mut_domain] = true;
			}

			if (fullBackup){
				backup.seedFromOther(q);
			}
			
			for(int i = 0; i < num_domains_mut; i++){
				int min_mutations = 1;
				int max_mutations = dd.MAX_MUTATIONS_LIMIT;
				int mut_domain = mut_domains[i];
				if (!fullBackup){
					//Backup
					System.arraycopy(q.domain[mut_domain],0,backup.domain[mut_domain],0,q.domain[mut_domain].length);
					System.arraycopy(q.domain_markings[mut_domain],0,backup.domain_markings[mut_domain],0,q.domain[mut_domain].length);
				}
				//Mutate
				dd.mutateUntilValid(mut_domain, q.domain, q.domain_markings, mutators[mut_domain], min_mutations, max_mutations);
			}
		}
		
		public void Evaluate(DomainDesignPMemberImpl q, boolean ShortcircuitOnRegression){
			//Reset markers.
			for(int k = 0; k < num_domains_mut; k++){
				int mut_domain = mut_domains[k];
				Arrays.fill(q.domain_markings[mut_domain], DNAMARKER_DONTMUTATE);
			}
			double deltaScoreSum = 0;
			
			int priority;
			priorityLoop: for(priority = 0; priority <= 2; priority++){
				for(int k = 0; k < num_domains_mut; k++){
					int mut_domain = mut_domains[k];
					for(int sd : q.scoredElements[mut_domain]){
						ScorePenalty s = q.penalties.get(sd);
						if (s.getPriority()==priority){
							s.evalScore(q.domain,q.domain_markings); //STATE CHANGE
						}
					}
				}

				//Decide whether we improved.

				double deltaScore = 0;
				for(ScorePenalty s : q.penalties){
					if (s.getPriority()==priority){
						deltaScore += s.getCDelta();
					}
				}

				revert_mutation = false;
				newPointReached = false;

				if (deltaScore > 0){
					revert_mutation = true;
				}
				
				/*
				if (bestWorstPenaltyScore[priority-1]<worstPenaltyScore){
					revert_mutation = true;
				} else if (bestWorstPenaltyScore[priority-1]==worstPenaltyScore){
					//Very often, there is a "wall" at a problem spot - allow the rest of the system
					//to also optimize in the background.
					if (deltaScore > 0){
						revert_mutation = true;
					}
				} else {
					newPointReached = true;
				}
				 */

				deltaScoreSum += deltaScore;
				if (ShortcircuitOnRegression){
					if (revert_mutation){
						//Short circuit! Get out of there.
						priority ++;

						break priorityLoop;
					} else {
						//keep the mutations
					}
				}
			}

			if (deltaScoreSum > 0){
				revert_mutation = true;
			}
			if (deltaScoreSum < 0){
				newPointReached = true;
			}
			
			if (!revert_mutation && priority < 3){
				throw new RuntimeException("Not all penalties ran!");
			}
		}

		public void Revert(DomainDesignPMemberImpl q) {
			for(int k = 0; k < num_domains_mut; k++){
				int mut_domain = mut_domains[k];
				//Revert ALL scores.
				for(int sd : q.scoredElements[mut_domain]){
					ScorePenalty s = q.penalties.get(sd);
					s.revert();
				}
				//Have to go back to old sequences..
				System.arraycopy(defaultBackupCache.domain[mut_domain],0,q.domain[mut_domain],0,q.domain[mut_domain].length);
				System.arraycopy(defaultBackupCache.domain_markings[mut_domain],0,q.domain_markings[mut_domain],0,q.domain[mut_domain].length);
			}	
		}
	}
	
	public boolean mutateAndTestAndBackup(DomainDesignPMemberImpl q) {
		//num_domains_mut = int_urn(1,Math.min(Math.min(worstPenalty==null?1:worstPenalty.getNumDomainsInvolved()-1,MAX_MUTATION_DOMAINS),mutableDomains.length-1));
		mutation_shared.Mutate(q, defaultBackupCache, false);
		mutation_shared.Evaluate(q, true);

		//Short circuiting constraint evaluation: If score was improved, keep, otherwise, break.
		//Penalties with higher priority are presumably harder to compute, thus the shortcircuiting advantage

		//Having run some or all of the affected penalty calculations,
		//Revert the states or Dedicate them as need be here.

		if (mutation_shared.revert_mutation){
			//Revert
			/*
			if (ENABLE_MARKINGS){
				for(int[] row : domain_markings){
					for(int q : row){
						System.out.printf("%4d",q);
					}
					System.out.println();
				}
			}
			*/
			mutation_shared.Revert(q);

			//Which priority did we short circuit on?
			/*
			if (priority==2){
				//If we're having trouble with crosstalk (which is really crazy n^2 stuff), 
				//ignore the advice of the self fold server.
				for(int k = 0; k < num_domains_mut; k++){
					mut_domain = mut_domains[k];
					//System.out.println(Arrays.toString(domain_markings[mut_domain]));
					Arrays.fill(q.domain_markings[mut_domain], 1);
				}
			}
			*/
			/*
			for(k = 0; k < num_domains_mut; k++){
				mut_domain = mut_domains[k];
				for(ScorePenalty s : scoredElements[mut_domain]){
					s.evalScore(domain, domain_markings);
					if (s.cur_score!=s.old_score){
						throw new RuntimeException("Assertion failure.");
					}
				}
			}
			*/
			return false;
		} else {
			//Dedicate
			for(ScorePenalty s : q.penalties){
				double l = s.cur_score;
			 	s.dedicate();
			 	
				//Check dedication! Slow!
			 	/*
				if (s.evalScore(q.domain, q.domain_markings)!=0){
					System.out.println(s.getClass());
					throw new RuntimeException("FAIL!");
				}
				*/
			}
			/*
			System.out.println("Current matrix:[");
			for(float[] row : DIR.currentMatrix){
				for(float val : row){
					System.out.printf("%.3f, ",val);
				}
				System.out.println();
			}
			System.out.println("]");
			*/

			/*
			if (OUTPUT_RUNTIME_STATS){
				System.out.println(num_mut_attempts+" , "+best_score);
			}
			//System.out.println(Arrays.toString(worstPenalty.getSeqs()));
			//current_score += deltaScoreSum;
			best_score = current_score;
			displayDomains(domain, false); //Updates public domains

			if (DebugPenaltiesNow){
				printDebugPenalties();
				System.out.println();
				displayDomains(domain, true);
				if (ENABLE_MARKINGS && false){
					for(int[] row : domain_markings){
						System.out.println(Arrays.toString(row));
					}
				}
			}

			double bestWorstPenaltyScoreSum = 0;
			for(double q : bestWorstPenaltyScore){
				//Disallow negative contributions. The overall EndThreshold can still be >0, however.
				bestWorstPenaltyScoreSum += Math.max(q,0);
			}
			if (bestWorstPenaltyScoreSum <= EndThreshold){
				break;
			}
			*/
			//System.out.printf("Components: CX%.3f IN%.3f JX%.3f INT%.3f\n",crosstalkSum,interactionSum,junctionPenalties,selfCrossTalk);

			return true;
		}
	}

	public boolean mutateAndTest(DomainDesignPMemberImpl q, DomainDesignPMemberImpl into) {
		mutation_shared.Mutate(q, into, true);
		mutation_shared.Evaluate(q, false);
		for(ScorePenalty s : q.penalties){
			double l = s.cur_score;
		 	s.dedicate();
		}
		return !mutation_shared.revert_mutation;
	}
}
