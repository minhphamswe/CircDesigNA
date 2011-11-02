package circdesigna.TripleSim;

import static circdesigna.TripleSim.TripleSim.splitIntoComponents;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;

import circdesigna.DomainPolymerGraph;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

/**
 * TripleSim Reaction Graph
 * @author Benjamin
 *
 */
public class ReactionGraph3X {
	public static class GraphEdge {
		public GraphNode towards;
		public double k = 0; //Unspecified = impossible.
		public String type;
		public GraphEdge reverse;
		public String toString(){
			return type+" k = "+k;
		}
	}
	public static class BimolecularNode extends GraphNode{
		public BimolecularNode(String structureString, GraphNode ... associate) {
			super(structureString);
			this.associate = associate;
		}
		public BimolecularNode(GraphNode ... associate) {
			super();
			this.associate = associate;
		}
		public GraphNode[] associate;
		public boolean isMonoMolecular(){
			return false;
		}
		public static String toString(GraphNode ... associate){
			StringBuffer sb = new StringBuffer();
			sb.append(associate[0]);
			sb.append(" + ");
			sb.append(associate[1]);
			return sb.toString(); 
		}
		public String toString(){
			return toString(associate);
		}
		public boolean equals(Object o){
			if (!(o instanceof BimolecularNode)){
				return false;
			}
			BimolecularNode other = (BimolecularNode) o;
			if (other.associate.length != associate.length){
				return false;
			}
			for(int i = 0; i < associate.length; i++){
				if (!other.associate[i].equals(associate[i])){
					return false;
				}
			}
			return true;
		}
	}
	public static class GraphNode implements Comparable<GraphNode>{
		private GraphNode(){
			structureString = null;
			structure = null;
		}
		public GraphNode(String structureString){
			this.structureString = structureString;
			structure = null;
		}
		public GraphNode(DomainPolymerGraph s1){
			structure = s1;
			structureString = structure.getStructureString();
		}
		public boolean visited = false;
		public boolean reachable = false;
		public int index = -1;
		public boolean stable = false;
		public double priority = 0;
		public double[] genData = null; //Used for graph algorithms that want to tack information on.
		public final DomainPolymerGraph structure;
		public final String structureString;
		//There SHOULDN'T be multiedges, check for this.
		public List<GraphEdge> neighbors = new LinkedList();
		public double initialConc;

		public boolean isMonoMolecular(){
			return true;
		}
		public final boolean isBiMolecular(){
			return !isMonoMolecular();
		}
		
		public int compareTo(GraphNode o) {
			/*
			int strandDifference = structure.getStrandRotations().size() - o.structure.getStrandRotations().size();
			if (strandDifference != 0){
				return strandDifference;
			}
			*/
			//If you have a higher priority than me, than I wait on you.
			double sign = Math.signum(o.priority-priority); 
			if (Double.isNaN(sign)){
				throw new RuntimeException("Undefined priorities for "+this+" "+o);
			}
			return (int)sign;
		}
		private int sumLengths(GraphNode[] associate) {
			int len = 0;
			for(GraphNode d : associate){
				len += d.structure.length();
			}
			return len;
		}
		public boolean equals(Object o){
			if (!(o instanceof GraphNode)){
				return false;
			}
			GraphNode other = (GraphNode) o;
			return other.structure == structure && other.structureString.equals(structureString);
		}
		public String toString(){
			return structureString;
		}
	}

	public static class Graph extends CircDesigNASystemElement{
		public Graph(CircDesigNAConfig config){
			super(config);
			unvisited = new PriorityQueue<GraphNode>();
		}
		
		public Graph(CircDesigNAConfig config, Comparator<GraphNode> comparator) {
			super(config);
			unvisited = new PriorityQueue<GraphNode>(64, comparator);
		}

		public List<PulseEvents> events = new ArrayList<PulseEvents>();
		public HashMap<String,GraphNode> allSingles = new HashMap();
		public HashMap<String,BimolecularNode> allDockings = new HashMap();
		public PriorityQueue<GraphNode> unvisited;
		public List<GraphEdge> edges = new ArrayList<GraphEdge>();
		public List<GraphNode> allVisited = new ArrayList<GraphNode>();

		public void visit(GraphNode a) {
			a.visited = true;
			allVisited.add(a);
		}
		public List<GraphNode> getNodes(){
			ArrayList<GraphNode> toRet = new ArrayList(allSingles.size()+allDockings.size());
			toRet.addAll(allSingles.values());
			toRet.addAll(allDockings.values());
			return toRet;
		}
		public int size() {
			return allSingles.size();
		}
		public int countOneWayReactions() {
			int toRet = 0;
			for(GraphEdge p : edges){
				if (p.k > 0){
					toRet++;
				}
				if (p.reverse.k > 0){
					toRet ++;
				}
			}
			return toRet;
		}
		public void cleanup(GraphNode target) {
			if(!Std.saveReactionDescriptions()){
				if (target.structure!=null){
					target.structure.getAnnotationTree().clear();
				}
			}
		}
		public void cleanup(GraphEdge reaction) {
			/*
			if(!Std.saveReactionDescriptions()){
				reaction.type = "FORWARD";
				reaction.reverse.type = "BACKWARD";
			}
			*/
		}
		public GraphNode addSpecies(DomainPolymerGraph neu){
			ArrayList<DomainPolymerGraph> components = splitIntoComponents(neu);

			if (components.size()==2){
				GraphNode B = addSpecies(components.get(0));
				GraphNode C = addSpecies(components.get(1));
				return getDocking(B, C);
			}
			if (components.size() > 2){
				throw new RuntimeException("An operation split a molecule into more than two parts. Not handled.");
			}

			//Otherwise, use the entire input structure.
			neu = neu.getCanonicalForm();
			GraphNode neuSingle = new GraphNode(neu);
			GraphNode already = allSingles.get(neuSingle.toString());
			if (already!=null){
				return already;
			}
			neuSingle.index = allSingles.size();
			allSingles.put(neuSingle.toString(),neuSingle);
			//We want the full n^2, because an instance of a molecule can react with itself.
			for(GraphNode d : allSingles.values()){
				BimolecularNode couple = new BimolecularNode(d,neuSingle);
				couple.index = -1-allDockings.size();
				allDockings.put(couple.toString(),couple);
				//couple.visited = true; //Meh, unnecessary. unvisited only contains singles anyway.
			}
			
			//Stability!
			neuSingle.stable = true;
			List<Integer> strandEnds = neu.getStrandRotations();
			if (!strandEnds.isEmpty()){ //Single strand structures are always stable.
				int[] strands = new int[strandEnds.size()+1];
				int[] strandIndexes = new int[neu.length()];
				{
					int strandIndex = 0;
					for(int i = 0; i < neu.length(); i++){
						if (!strandEnds.isEmpty() && strandEnds.get(0) <= i){
							strandEnds.remove(0);
							strandIndex++;
						}
						strandIndexes[i] = strandIndex;
					}
				}
				for(int i = 0; i < neu.length(); i++){
					int j = neu.getDomainPair(i); 
					if (j >= 0){
						strands[strandIndexes[i]]+=neu.getDomainDefs().getDomainLength(neu.getDomain(i));
						strands[strandIndexes[j]]+=neu.getDomainDefs().getDomainLength(neu.getDomain(j));
					}
				}
				for(int i = 0; i < strands.length; i++){
					if (strands[i] < 6){
						neuSingle.stable = false;
						//System.out.println(neuSingle.structureString+" is not stable.");
						break;
					}
				}
			}
			
			unvisited.add(neuSingle);
			return neuSingle;
		}

		public void recreateUnvisited() {
			unvisited.clear();
			
			for(GraphNode p : allSingles.values()){
				if (!p.visited){
					unvisited.add(p);
				}
			}
		}
		
		public BimolecularNode getDocking(GraphNode A, GraphNode B){
			if (A.isBiMolecular() || B.isBiMolecular()){
				throw new RuntimeException("No docking bimolecular molecules (physically unrealistic).");
			}
			BimolecularNode bimolecularNode = allDockings.get(BimolecularNode.toString(A,B));
			if (bimolecularNode==null){
				bimolecularNode = allDockings.get(BimolecularNode.toString(B,A));
				if (bimolecularNode==null){
					throw new RuntimeException("Invalid state: a pair of molecule nodes with no bimolecular node");
				}
			}
			return bimolecularNode;
		}

		public GraphEdge addReaction(GraphNode A, GraphNode B) {
			if (A==B){
				throw new RuntimeException("No reaction (so far) should change a molecule into itself. Error.");
			}
			for(GraphEdge d : A.neighbors){
				if (d.towards == B){
					return d;
				}
			}
			for(GraphEdge d : B.neighbors){
				if (d.towards == A){
					return d.reverse;
				}
			}
			GraphEdge fw = new GraphEdge();
			GraphEdge back = new GraphEdge();
			fw.towards = B;
			back.towards = A;
			fw.reverse = back;
			back.reverse = fw;
			A.neighbors.add(fw);
			B.neighbors.add(back);
			/*
			//Subset of "recreate unvisited".
			if (!A.visited){
				unvisited.remove(A);
				unvisited.add(A);
			}
			if (!B.visited){
				unvisited.remove(B);
				unvisited.add(B);
			}
			*/
			edges.add(fw);
			return fw;
		}
	}
}
	
