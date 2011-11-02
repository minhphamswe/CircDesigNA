package circdesigna.TripleSim;

import static circdesigna.CircDesigNA_SharedUtils.isComplements;
import circdesigna.DomainDefinitions;
import circdesigna.DomainPolymerGraph;
import circdesigna.AbstractComplex.Annotation;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

public class Associate extends CircDesigNASystemElement{
	public Associate(CircDesigNAConfig config){
		super(config);
	}
	public void associate(Graph g, GraphNode a) {
		if(a.isBiMolecular()){
			throw new RuntimeException("Association reactions occur intramolecularly.");
		}
		
		AbstractComplexSet<DomainPolymerGraph> input = new AbstractComplexSet<DomainPolymerGraph>();
		input.add(a.structure);
		//Associate a single duplex.
		AbstractComplexSet<DomainPolymerGraph> associations = Associate.getAssociations(input);

		for(DomainPolymerGraph d : associations){
			GraphNode map = g.addSpecies(d);
			GraphEdge reaction = g.addReaction(a,map);
			
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

			reaction.k = Std.kinetics.ckfi * (longestDuplex > 6 ? Math.sqrt(longestDuplex / 6.): 1);
			double dg = Std.getDeltaGPerStackPair()*longestDuplex;
			reaction.reverse.k = Math.exp(dg / (Std.kinetics.R * Std.kinetics.T)) * reaction.k;
			reaction.type = "Associate "+bestDescription+" END "+d.getStructureString();
			reaction.reverse.type = "Disassociate "+bestDescription+" START "+d.getStructureString();
			
			g.cleanup(map);
			g.cleanup(reaction);
		}
	}
	
	public static AbstractComplexSet<DomainPolymerGraph> getAssociations(AbstractComplexSet<DomainPolymerGraph> input) {
		AbstractComplexSet<DomainPolymerGraph> toRet = new AbstractComplexSet<DomainPolymerGraph>();
		for(DomainPolymerGraph q : input){
			DomainDefinitions dd = q.getDomainDefs();
			int[] levels = q.getDomainLevels();
			for(int x = 0; x < levels.length; x++){
				for(int y = x+1; y < levels.length; y++){
					if (Math.abs(y - x)==1){
						continue; //Cannot pair a domain with its neighbor, that produces an empty loop.
					}
					if (x < 0 || y < 0 || levels[x] < 0 || levels[y] < 0)
						continue;
					if (levels[x] != levels[y])
						continue;
					if (q.getDomainPair(x) >= 0 || q.getDomainPair(y) >= 0){
						continue;
					}
					if (!isComplements(q.getDomain(x),q.getDomain(y),q.getDomainDefs())){
						continue;
					}
					DomainPolymerGraph Cclone = new DomainPolymerGraph(q.getDomainDefs());
					DomainPolymerGraph.readStructure("AC "+q.getStructureString(), Cclone);
					int len = 0;
					Cclone.annotate(q.getAnnotationTree());
					if (!Cclone.setDomainPair(x, y)){
						throw new RuntimeException("Could not pair domains "+x+" "+y);
					}	
					len += dd.getDomainLength(q.getDomain(x));
					//Take the whole duplex.
					for(int direction : new int[]{-1,1}){
						for(int i = 1; true; i++){
							int nx = x + i * direction;
							int ny = y - i * direction;
							if (Math.abs(ny - nx)==1){
								break; //Cannot pair a domain with its neighbor, that produces an empty loop.
							}
							if (nx < 0 || ny < 0 || levels[nx] < 0 || levels[ny] < 0)
								break;
							if (levels[nx] != levels[ny])
								break;
							if (q.getDomainPair(nx) >= 0 || q.getDomainPair(ny) >= 0){
								break;
							}
							if (!isComplements(q.getDomain(nx),q.getDomain(ny),q.getDomainDefs())){
								break;
							}
							len += dd.getDomainLength(q.getDomain(nx));
							if (!Cclone.setDomainPair(nx, ny)){
								throw new RuntimeException("Could not pair domains "+nx+" "+ny);
							}
						}
					}
					Cclone.annotate("Associate "+x+" "+y);
					Cclone.annotate("Associated "+len+" nucleotides");

					toRet.add(Cclone);
				}
			}
		}
		return toRet;
	}
	public static void reactionKineticsAssociation(CircDesigNAConfig Std, GraphEdge reaction, int longestDuplex) {
		if (reaction.reverse.towards.isBiMolecular()){
			reaction.k = Std.kinetics.ckf * (longestDuplex > 6 ? Math.sqrt(longestDuplex / 6.): 1);
		} else {
			reaction.k = Std.kinetics.ckfi * (longestDuplex > 6 ? Math.sqrt(longestDuplex / 6.): 1);
		}
		double dg = Std.getDeltaGPerStackPair()*(longestDuplex-1) + 1.96;
		reaction.reverse.k = Math.exp(dg / (Std.kinetics.R * Std.kinetics.T)) * reaction.k;
	}
}
