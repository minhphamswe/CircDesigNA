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
	
	public DomainDesignBlockDesignerImpl(int num_domain, int[] domain_length, int[] mutableDomains, DesignerCode[] mutators, DomainStructureData dsd, DomainDesigner domainDesigner){
		//Make a back buffer for domain reverts
		old_domains = new int[num_domain][];
		old_domains_markings = new int[num_domain][];
		for(int k = 0; k < num_domain; k++){
			old_domains[k] = new int[domain_length[k]];
			old_domains_markings[k] = new int[domain_length[k]];
		}
		this.dsd = dsd;
		this.mutators = mutators;
		this.dd = domainDesigner;
		this.num_domain = num_domain;
		this.mutableDomains = mutableDomains;
		this.domain_length = domain_length;
	}
	private DesignerCode[] mutators;
	private int[] domain_length;
	private int num_domain;
	private int[][] old_domains, old_domains_markings;
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
	public boolean mutateAndTestAndBackup(DomainDesignPMemberImpl q) {
		//Mutate, in some way.
		int num_mut_attempts = 0, mut_domain, num_domains_mut, min_mutations, max_mutations;
		double worstPenaltyScore = -1, deltaScore = 0, deltaScoreSum = 0;
		double[] bestWorstPenaltyScore = new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
		boolean revert_mutation = false, newPointReached = false;
		boolean[] domain_mutated = new boolean[num_domain];
		int[] mut_domains = new int[dd.MAX_MUTATION_DOMAINS];

		double compareWithDelta = 0;
		int timesSameBeforeBump = 10;
		int timesSameCount = 0;
		double bumpAmount = 10;
		
		//num_domains_mut = int_urn(1,Math.min(Math.min(worstPenalty==null?1:worstPenalty.getNumDomainsInvolved()-1,MAX_MUTATION_DOMAINS),mutableDomains.length-1));
		num_domains_mut = int_urn(1,Math.min(dd.MAX_MUTATION_DOMAINS,mutableDomains.length));
		Arrays.fill(domain_mutated,false);
		for(int k = 0; k < num_domains_mut; k++){
			min_mutations = 1;
			/*
			if (worstPenalty==null){
				mut_domain = int_urn(0, mutableDomains.length-1);
				mut_domain = mutableDomains[mut_domain];
			} else {
				
					mut_domain = int_urn(0, mutableDomains.length-1);
					mut_domain = mutableDomains[mut_domain];
				}while (!worstPenalty.affectedBy(mut_domain));
			}
			*/
			max_mutations = dd.MAX_MUTATIONS_LIMIT;
			

			int mut_domain_search = int_urn(0, mutableDomains.length-1);
			mut_domain = 0;
			for(int o = 0; o < mutableDomains.length; o++){
				mut_domain = mutableDomains[(mut_domain_search+o)%mutableDomains.length];
				if (domain_mutated[mut_domain]){
					continue;
				}
				break;
			}
			if(domain_mutated[mut_domain]){
				break;
			}
			
			/*
			if (timesSameCount > 0){
				max_mutations = Math.max(MAX_MUTATIONS_LIMIT/timesSameCount,min_mutations);
			}
			*/

			mut_domains[k] = mut_domain;
			//Backup
			System.arraycopy(q.domain[mut_domain],0,old_domains[mut_domain],0,domain_length[mut_domain]);
			System.arraycopy(q.domain_markings[mut_domain],0,old_domains_markings[mut_domain],0,domain_length[mut_domain]);
			//Mutate
			dd.mutateUntilValid(mut_domain, q.domain, domain_length, q.domain_markings, mutators[mut_domain], min_mutations, max_mutations);
			domain_mutated[mut_domain] = true;
		}

		//Short circuiting constraint evaluation: If score was improved, keep, otherwise, break.
		//Penalties with higher priority are presumably harder to compute, thus the shortcircuiting advantage

		//Reset markers.
		for(int k = 0; k < num_domains_mut; k++){
			mut_domain = mut_domains[k];
			Arrays.fill(q.domain_markings[mut_domain], DNAMARKER_DONTMUTATE);
		}
		deltaScoreSum = 0;
		int priority;
		
		priorityLoop: for(priority = 0; priority <= 2; priority++){
			for(int k = 0; k < num_domains_mut; k++){
				mut_domain = mut_domains[k];
				for(int sd : q.scoredElements[mut_domain]){
					ScorePenalty s = q.penalties.get(sd);
					if (s.getPriority()==priority){
						s.evalScore(q.domain,q.domain_markings); //STATE CHANGE
					}
				}
			}

			//Decide whether we improved.

			deltaScore = 0;
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
			if (revert_mutation){
				//Short circuit! Get out of there.
				priority ++;

				break priorityLoop;
			} else {
				//Keep the mutations
				bestWorstPenaltyScore[priority] = Math.min(bestWorstPenaltyScore[priority],worstPenaltyScore);
			}
		}

		if (deltaScoreSum > 0){
			revert_mutation = true;
		}
		if (deltaScoreSum < 0){
			newPointReached = true;
		}
		//Having run some or all of the affected penalty calculations,
		//Revert the states or Dedicate them as need be here.

		if (revert_mutation){
			timesSameCount++;
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
			for(int k = 0; k < num_domains_mut; k++){
				mut_domain = mut_domains[k];
				//Revert ALL scores.
				for(int sd : q.scoredElements[mut_domain]){
					ScorePenalty s = q.penalties.get(sd);
					s.revert();
				}
				//Have to go back to old sequences..
				System.arraycopy(old_domains_markings[mut_domain],0,q.domain_markings[mut_domain], 0, q.domain[mut_domain].length);
				System.arraycopy(old_domains[mut_domain], 0, q.domain[mut_domain], 0, q.domain[mut_domain].length);
			}	

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
			if (priority != 3){
				throw new RuntimeException("Assertion failure: self fold layer (layer 2) never ran.");
			}
			
			timesSameCount = 0;
			//Dedicate
			for(ScorePenalty s : q.penalties){
			 	s.dedicate();
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
}
