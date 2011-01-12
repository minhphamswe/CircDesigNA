package DnaDesign.impl;

import java.util.ArrayList;
import java.util.List;

import DnaDesign.AbstractDesigner.PopulationDesignMember;
import DnaDesign.DomainDesigner.ScorePenalty;

/**
 * Concrete population member.
 */
public class DomainDesignPMemberImpl extends PopulationDesignMember<DomainDesignPMemberImpl>{
	public DomainDesignPMemberImpl designerCopyConstructor() {
		return new DomainDesignPMemberImpl(null, null, null, null);
	}
	public void seedFromOther(DomainDesignPMemberImpl pdm) {
		if (penalties==null){
			penalties = new ArrayList(); 
			for(ScorePenalty q : pdm.penalties){
				penalties.add((ScorePenalty)q.clone());
			}
		}
		scoredElements = copy2dInt(pdm.scoredElements,scoredElements);
		domain = copy2dInt(pdm.domain,domain);
		domain_markings = copy2dInt(pdm.domain_markings,domain_markings);
	}
	private static int[][] copy2dInt(int[][] se, int[][] toRet) {
		if (toRet==null || toRet.length < se.length){
			toRet = new int[se.length][];
			for(int k = 0; k < se.length; k++){
				int[] toCopy = se[k];
				toRet[k] = new int[toCopy.length];
			}
		}
		for(int k = 0; k < se.length; k++){
			int[] toCopy = se[k];
			System.arraycopy(toCopy, 0, toRet[k],0,toRet[k].length);
		}
		return toRet;
	}
	public List<ScorePenalty> penalties; 
	public int[][] scoredElements; 
	public int[][] domain; 
	public int[][] domain_markings;
	public DomainDesignPMemberImpl(List<ScorePenalty> penalties, int[][] scoredElements, int[][] domain, int[][] domain_markings){
		this.penalties = penalties;
		this.scoredElements = scoredElements;
		this.domain = domain;
		this.domain_markings = domain_markings;
	}
}
