package circdesigna.TripleSim;

import static circdesigna.TripleSim.TripleSim.getDuplexes;
import static circdesigna.TripleSim.TripleSim.splitIntoComponents;

import java.util.ArrayList;

import circdesigna.DomainPolymerGraph;
import circdesigna.AbstractComplex.Annotation;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

public class Disassociate extends CircDesigNASystemElement{
	public Disassociate(CircDesigNAConfig config){
		super(config);
	}
	public void disassociate(Graph g, GraphNode A) {
		AbstractComplexSet<DomainPolymerGraph> input = new AbstractComplexSet<DomainPolymerGraph>();
		input.add(A.structure);
		
		AbstractComplexSet<DomainPolymerGraph> disassociate = getDisassociations(input);

		//May cause splitting.
		for(DomainPolymerGraph q : disassociate){
			ArrayList<DomainPolymerGraph> components = splitIntoComponents(q);

			int longestDuplex = 0;
			String bestDescription = "";
			//Within a set of annotations, we are assured that the rotations are clean.
			for(Annotation u : q.getAnnotationLevel("Disassociate.*")){
				String myDescription = "";
				int curDuplex = 0;
				for(String s : u.annotations){
					String[] split = s.split("\\s+");
					if (split[0].equals("Disassociate")){
						myDescription += s.split("\\s+",2)[1];
					}
					if (split[0].equals("Disassociated")){
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
			
			if (longestDuplex >= 10){
				continue; 
			}
			
			GraphNode target = g.addSpecies(q);
			
			GraphEdge reaction = g.addReaction(A, target);
			Associate.reactionKineticsAssociation(Std, reaction.reverse, longestDuplex);
			
			reaction.type = "Detach "+bestDescription+" START "+A.structure.getStructureString();
			reaction.reverse.type = "Attach "+bestDescription+" END "+A.structure.getStructureString();
			
			g.cleanup(target);
			g.cleanup(reaction);
		}
	}

	private AbstractComplexSet<DomainPolymerGraph> getDisassociations(AbstractComplexSet<DomainPolymerGraph> input) {
		AbstractComplexSet<DomainPolymerGraph> toRet = new AbstractComplexSet<DomainPolymerGraph>();
		
		for(DomainPolymerGraph q : input){
			for(int[] duplex : getDuplexes(q)){
				DomainPolymerGraph C = new DomainPolymerGraph(q.getDomainDefs());
				DomainPolymerGraph.readStructure("AB "+q.getStructureString(), C);
				C.annotate(q.getAnnotationTree());
				int numNucs = 0;
				StringBuffer dis = new StringBuffer();
				
				for(int y = duplex[0]; y<= duplex[1]; y++){
					dis.append(y+" "+C.getDomainPair(y));
					numNucs += q.getDomainDefs().getDomainLength(q.getDomain(y));
					C.setDomainPair(y, -1);
				}
				
				C.annotate("Disassociate "+dis);
				C.annotate("Disassociated "+numNucs+" nucleotides");
				toRet.add(C);
			}
		}
		
		return toRet;
	}


}
