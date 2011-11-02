package circdesigna.TripleSim;

import static circdesigna.TripleSim.TripleSim.getConnectedComponents;

import java.util.ListIterator;

import circdesigna.DomainDefinitions;
import circdesigna.DomainPolymerGraph;
import circdesigna.AbstractComplex.Annotation;
import circdesigna.TripleSim.ReactionGraph3X.BimolecularNode;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;


public class Attachment extends CircDesigNASystemElement{
	public Attachment(CircDesigNAConfig config){
		super(config);
	}
	public void attach(Graph g, BimolecularNode couple) {
		DomainPolymerGraph A = couple.associate[0].structure;
		DomainPolymerGraph B = couple.associate[1].structure;

		String Astr = A.getStructureString();
		String Bstr = B.getStructureString();

		AbstractComplexSet<DomainPolymerGraph> docking = new AbstractComplexSet();
		{
			DomainPolymerGraph C = new DomainPolymerGraph(A.getDomainDefs());
			DomainPolymerGraph.readStructure("AB "+A.getStructureString()+B.getStructureString(), C);
			docking.add(C.getCanonicalForm());
		}
		addAllInsertions(docking, Astr, Bstr, A.getDomainDefs());
		addAllInsertions(docking, Bstr, Astr, A.getDomainDefs());
		
		//Associate a single duplex.
		AbstractComplexSet<DomainPolymerGraph> associations = Associate.getAssociations(docking);

		//Remove disonnected structures.
		ListIterator<DomainPolymerGraph> filter = associations.listIterator();
		while(filter.hasNext()){
			DomainPolymerGraph dpg = filter.next();
			if (getConnectedComponents(dpg).size()!=1){
				//System.out.println("Disconnected: "+dpg.getStructureString());
				filter.remove();
			}
		}
		
		for(DomainPolymerGraph d : associations){
			GraphNode map = g.addSpecies(d);
			GraphEdge reaction = g.addReaction(couple,map);
			
			int longestDuplex = 0;
			String bestDescription = "";
			//Within a set of annotations, we are assured that the rotations are clean.
			for(Annotation q : d.getAnnotationLevel("Associate.*")){
				String myDescription = "";
				int curDuplex = 0;
				for(String e : q.annotations){
					String[] split = e.split("\\s+");
					if (split[0].equals("Associate")){
						myDescription += e.split("\\s+",2)[1];
					}
					if (split[0].equals("Associated")){
						curDuplex = new Integer(split[1]);
					}
				}
				if (curDuplex > longestDuplex){
					longestDuplex = curDuplex;
					bestDescription = myDescription;
				}
			}
			if (longestDuplex == 0){ //Occurs if some of the domains have 0 length.
				continue;
			}			

			Associate.reactionKineticsAssociation(Std, reaction, longestDuplex);
			reaction.type = "Attachment "+bestDescription+" END "+d.getStructureString();
			reaction.reverse.type = "Detachment "+bestDescription+" START "+d.getStructureString();
			
			g.cleanup(map);
			g.cleanup(reaction);
		}
	}


	/**
	 * Inserts A into B
	 */
	private static void addAllInsertions(AbstractComplexSet<DomainPolymerGraph> rotations,String a, String b, DomainDefinitions defs) {
		for(int k = 0; k < b.length()-1; k++){
			if (b.substring(k,k+2).equals("}[")){
				DomainPolymerGraph C = new DomainPolymerGraph(defs);
				DomainPolymerGraph.readStructure("AB " + b.substring(0,k+1)+a+b.substring(k+1), C);
				rotations.add(C.getCanonicalForm());
			}
		}
	}
}
