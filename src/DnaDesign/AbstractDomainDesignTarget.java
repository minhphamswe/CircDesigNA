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
	public ArrayList<DomainSequence> pairsOfDomains = new ArrayList();
	public ArrayList<DomainSequence[]> hairpinClosings = new ArrayList();
	
	public static final int overlap_length = 8;
	public ArrayList<DomainSequence> singleDomainsWithOverlap = new ArrayList();
	
	public void clear(){
		generalizedSingleStranded.clear();
		pairsOfDomains.clear();
		hairpinClosings.clear();
		wholeStrands.clear();
		singleDomainsWithOverlap.clear();
	}
	public void removeDuplicates(){
		//TODO: remove substrings
		
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(wholeStrands);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(generalizedSingleStranded);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(pairsOfDomains);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(singleDomainsWithOverlap);
		
		//TODO: remove duplicates from hairpinclosingfinder
	}
	
	private DomainStructureData dsd;
	
	public void addTargetStructure(String name, String inputStrand) {
		DomainPolymerGraph dpg = new DomainPolymerGraph(dsd);
		DomainPolymerGraph.readStructure(name, inputStrand, dpg);
		
		DomainDesigner_SharedUtils.utilSingleStrandedFinder(dpg, generalizedSingleStranded);
		DomainDesigner_SharedUtils.utilHairpinClosingFinder(dpg, hairpinClosings);
		DomainDesigner_SharedUtils.utilPairsOfDomainsFinder(dpg, pairsOfDomains);
		DomainDesigner_SharedUtils.utilSingleDomainsWithOverlap(dpg, singleDomainsWithOverlap, overlap_length);
		
		//Whole strands
		for(String subStrand : inputStrand.split("}")){
			DomainSequence ds = new DomainSequence();
			ds.setDomains(subStrand,dsd,dpg);
			wholeStrands.add(ds);
		}
		
		removeDuplicates();
	}
}