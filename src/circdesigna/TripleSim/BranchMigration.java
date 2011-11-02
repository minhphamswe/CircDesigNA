package circdesigna.TripleSim;

import static circdesigna.CircDesigNA_SharedUtils.isComplements;
import static circdesigna.DomainSequence.NA_COMPLEMENT_FLAGINV;
import static circdesigna.TripleSim.TripleSim.getDuplexes;

import java.util.ArrayList;
import java.util.Arrays;

import circdesigna.DomainPolymerGraph;
import circdesigna.AbstractComplex.Annotation;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;


public class BranchMigration extends CircDesigNASystemElement{
	public BranchMigration(CircDesigNAConfig config){
		super(config);
	}
	
	public void branchMigrate(Graph g, GraphNode node) {
		if (node.isBiMolecular()){
			throw new RuntimeException("Branch migration can only occur on single molecules");
		}
		AbstractComplexSet<DomainPolymerGraph> set = new AbstractComplexSet<DomainPolymerGraph>();
		
		int[][] duplexes = getDuplexes(node.structure);
		for(int y = 0; y < duplexes.length; y++){
			for(int x = y+1; x < duplexes.length; x++){
				strandDisplace(set, node.structure, duplexes[y], duplexes[x]);
			}
		}

		for(DomainPolymerGraph q : set){
			GraphNode map = g.addSpecies(q);
			GraphEdge reaction = g.addReaction(node,map);
			//All branch migrations with the same start and end products are isoenergetic.
			String myDescription = "";
			int curDuplex = 1;
			for(Annotation a : q.getAnnotationLevel("Displace.*")){
				for(String e : a.annotations){
					String[] split = e.split("\\s+");
					if (split[0].equals("Displace")){
						myDescription += e.split("\\s+",2)[1];
					}
					if (split[0].equals("Displaced")){
						curDuplex = new Integer(split[1]);
					}
				}
				break;
			}
			
			reaction.k = Std.kinetics.ckb/(curDuplex*curDuplex);
			//Reverse is left undefined.
			
			reaction.type = "Branch Migration "+myDescription + " END " + q.getStructureString();
			
			g.cleanup(map);
			g.cleanup(reaction);
		}	
	}
	
	
	private static DomainPolymerGraph clone(DomainPolymerGraph q) {
		DomainPolymerGraph neu = new DomainPolymerGraph(q.getDomainDefs());
		DomainPolymerGraph.readStructure("AC "+q.getStructureString(), neu);
		neu.annotate(q.getAnnotationTree());
		return neu;
	}
	
	private static void strandDisplace(AbstractComplexSet<DomainPolymerGraph> outputs, final DomainPolymerGraph q, int[] d1_, int[] d2_) {
		//The duplexes must be distinct.
		boolean match = true;
		for(int i = 0; i < 4; i++){
			match &= (d1_[i] == d2_[(i+2)%4]);
		}
		if (match || Arrays.equals(d1_, d2_)){
			return;
		}
		
		class modArrayList extends ArrayList<int[]>{
			int lastBrokenBond = -1;
			public boolean add(int[] toAdd){
				if (q.getDomainPair(toAdd[0])>=0){
					return false;
				}
				int pairOfBroken = q.getDomainPair(toAdd[1]);
				if (pairOfBroken>=0){
					if (lastBrokenBond >= 0){
						if (Math.abs(lastBrokenBond - pairOfBroken) != 1){
							return false;
						}
					}
					lastBrokenBond = pairOfBroken;
				} else {
					//throw new RuntimeException("Branch migration getting something for free!");
					return false; //Don't allow associations after branch migration; the kinetics are different if this is the case.
				}
				return super.add(toAdd);
			}
		}
		if (canBranchMigrate(d1_[1],d2_[3]+1, q)){
			int[] d1 = d1_;
			int[] d2 = d2_;
			//'1' CW
			ArrayList<int[]> displaced = new modArrayList();
			DomainPolymerGraph c = clone(q); 
			while (d1[1] >= 0 && isComplements(q.getDomain(d2[3] + 1), q.getDomain(d1[1]), q.getDomainDefs())){
				if (!displaced.add(new int[]{d2[3] + 1, d1[1]})) break;
				d1 = new int[]{d1[0], d1[1] - 1, d1[2] + 1, d1[3]};
				d2 = new int[]{d2[0] - 1, d2[1], d2[2], d2[3] + 1};
			}
			strandDisplaceBases(c,displaced,outputs);
		}
		if (canBranchMigrate(d1_[2]-1, d2_[0], q)){
			int[] d1 = d1_;
			int[] d2 = d2_;
			//'1' CW
			ArrayList<int[]> displaced = new modArrayList();
			DomainPolymerGraph c = clone(q); 
			//'1' CCW
			while (d1[2] > 0 && isComplements(q.getDomain(d1[2] - 1), q.getDomain(d2[0]), q.getDomainDefs())){
				if (!displaced.add(new int[]{d1[2] - 1, d2[0]})) break;
				d1 = new int[]{d1[0], d1[1] + 1, d1[2] - 1, d1[3]};
				d2 = new int[]{d2[0] + 1, d2[1], d2[2], d2[3] - 1};
			}
			strandDisplaceBases(c,displaced,outputs);
		}
		if (canBranchMigrate(d2_[3], d1_[1]+1, q)){
			int[] d1 = d1_;
			int[] d2 = d2_;
			//'1' CW
			ArrayList<int[]> displaced = new modArrayList();
			DomainPolymerGraph c = clone(q); 
			//'2' CW
			while (d2[3] >= 0 && isComplements(q.getDomain(d1[1] + 1), q.getDomain(d2[3]), q.getDomainDefs())){
				if (!displaced.add(new int[]{d1[1] + 1, d2[3]})) break;
				d1 = new int[]{d1[0], d1[1] + 1, d1[2] - 1, d1[3]};
				d2 = new int[]{d2[0] + 1, d2[1], d2[2], d2[3] - 1};
			}
			strandDisplaceBases(c,displaced,outputs);
		}

		if (canBranchMigrate(d2_[0]-1, d1_[2], q)){
			int[] d1 = d1_;
			int[] d2 = d2_; 
			//'2' CCW
			ArrayList<int[]> displaced = new modArrayList();
			DomainPolymerGraph c = clone(q);
			while (d2[0] > 0 && d2[0] > 0 && isComplements(q.getDomain(d2[0] - 1), q.getDomain(d1[2]), q.getDomainDefs())){
				if (!displaced.add(new int[]{d2[0] - 1, d1[2]})) break;
				d1 = new int[]{d1[0], d1[1] - 1, d1[2] + 1, d1[3]};
				d2 = new int[]{d2[0] - 1, d2[1], d2[2], d2[3] + 1};
			}
			strandDisplaceBases(c,displaced,outputs);
		}
	}

	private static boolean canBranchMigrate(int i, int j, DomainPolymerGraph c) {
		if (i < 0 || j < 0){
			return false;
		}
		return ((c.getDomainPair(i)<0) ^ (c.getDomainPair(j)<0)) && spatiallyClose(i, j, c, 0);
	}
	private static boolean spatiallyClose(int i, int j, DomainPolymerGraph c, int depth){
		if (depth >= 2){
			//Not close enough, i and j are separated by more than 2 helixes. 
			return false;
		}
		//True if the 5' most base of domain i is close to the 3' most base of domain j.
		int targetIp1 = c.getDomainPair(i+1);
		if (targetIp1 < 0){
			return false; //Either i ends a strand (can flop on itself), or there is a floppy domain between i and j.
		}
		if (targetIp1+1 == j){
			return true;
		}
		return spatiallyClose(targetIp1,j,c, depth+1);
	}
	private static void strandDisplaceBases(DomainPolymerGraph c, ArrayList<int[]> disp, AbstractComplexSet<DomainPolymerGraph> outputs) {
		String displace = "Displace ";
		int numNucs = 0;
		for(int[] row : disp){
			if (!c.setDomainPair(row[0], row[1])){
				break; //No harm done if fails, just stop the branch migration.
			}
			displace += row[0]+" "+row[1]+" ";
			numNucs += c.getDomainDefs().domainLengths[c.getDomain(row[0]) & NA_COMPLEMENT_FLAGINV];
		}
		if (numNucs == 0){
			return; 
		}
		String displaced = "Displaced "+numNucs+" nucleotides";
		c.annotate(displace);
		c.annotate(displaced);
		outputs.add(c);
	}
}
