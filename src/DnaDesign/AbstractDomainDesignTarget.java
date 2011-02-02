package DnaDesign;

import java.util.ArrayList;

/**
 * Represents a design target, in as many ways as possible.
 * 
 * Specifically, this class provides convenient calculations for the score penalty creator
 * to use.
 */
public class AbstractDomainDesignTarget {
	public AbstractDomainDesignTarget(DomainStructureData dsd){
		this.dsd = dsd;
		
	}
	public ArrayList<DomainSequence> wholeStrands = new ArrayList();
	public ArrayList<DomainSequence> generalizedSingleStranded = new ArrayList();
	public ArrayList<DomainSequence> singleDomains = new ArrayList();
	public ArrayList<DomainSequence> pairsOfDomains = new ArrayList();
	
	public static class HairpinClosingTarget {
		public DomainSequence[] stemOnly;
		/**
		 * Use these sequences for printing out this guy.
		 */
		public DomainSequence[] stemAndOpening;
		/**
		 * Outside (true) means deltaDeltaG = (d0d1*d2d3)-(d1*d2)
		 * Inside means deltaDeltaG = (d0d1*d2d3)-(d0*d3)
		 * 
		 * When outside, the "unwanted structure" region is 0...(stemAndOpening.length-stemOnly.length)
		 * When inside, the unwanted structure region is (stemAndOpening.length-stemOnly.length)...stemAndOpening.length
		 */
		public boolean outside;
		public HairpinClosingTarget(int domain0, int domain1, int domain2, int domain3, boolean outside, AbstractComplex dsg){
			DomainSequence sA1 = new DomainSequence();
			DomainSequence sA2 = new DomainSequence();
			DomainSequence s1 = new DomainSequence();
			DomainSequence s2 = new DomainSequence();
			sA1.setDomains(domain0, domain1, dsg);
			sA2.setDomains(domain2, domain3, dsg);
			this.outside = outside;
			if (outside){
				s1.setDomains(domain1, dsg);
				s2.setDomains(domain2, dsg);
			} else {
				s1.setDomains(domain0, dsg);
				s2.setDomains(domain3, dsg);
			}
			stemAndOpening = new DomainSequence[]{sA1,sA2};
			stemOnly = new DomainSequence[]{s1,s2};
		}
		public String toString(int[][] domain){
			StringBuffer toRet = new StringBuffer();
			for(int i = 0; i < stemAndOpening.length; i++){
				toRet.append("Hairpin "+i+" ");
				for(int j = 0; j < stemAndOpening[i].length(domain); j++){
					toRet.append(DnaDefinition.displayBase(stemAndOpening[i].base(j, domain)));
				}
				toRet.append("\n");
			}
			return toRet.toString();
		}
		public boolean equals(Object other){
			if (!(other instanceof HairpinClosingTarget)){
				return false;
			}
			HairpinClosingTarget oth = (HairpinClosingTarget) other;
			for(int i = 0; i < stemAndOpening.length; i++){
				if (!oth.stemAndOpening[i].equals(stemAndOpening[i])){
					return false;
				}
			}
			return true;
		}
	}
	public ArrayList<HairpinClosingTarget> hairpinClosings = new ArrayList();
	
	public static final int overlap_length = 8;
	public ArrayList<DomainSequence> singleDomainsWithOverlap = new ArrayList();
	
	public void clear(){
		generalizedSingleStranded.clear();
		pairsOfDomains.clear();
		hairpinClosings.clear();
		wholeStrands.clear();
		singleDomainsWithOverlap.clear();
		singleDomains.clear();
	}
	
	private DomainStructureData dsd;
	
	public void addTargetStructure(String name, String inputStrand) {
		DomainPolymerGraph dpg = new DomainPolymerGraph(dsd);
		DomainPolymerGraph.readStructure(name, inputStrand, dpg);
		
		DomainDesigner_SharedUtils.utilSingleStrandedFinder(dpg, generalizedSingleStranded);
		DomainDesigner_SharedUtils.utilHairpinClosingFinder(dpg, hairpinClosings);
		DomainDesigner_SharedUtils.utilPairsOfDomainsFinder(dpg, pairsOfDomains);
		DomainDesigner_SharedUtils.utilSingleDomainsFinder(dpg, singleDomains);
		DomainDesigner_SharedUtils.utilSingleDomainsWithOverlap(dpg, singleDomainsWithOverlap, overlap_length);
		
		//Whole strands
		for(String subStrand : inputStrand.split("}")){
			DomainSequence ds = new DomainSequence();
			ds.setDomains(subStrand,dsd,dpg);
			wholeStrands.add(ds);
		}
	}
}