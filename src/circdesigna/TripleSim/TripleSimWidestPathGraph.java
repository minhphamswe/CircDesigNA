package circdesigna.TripleSim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import circdesigna.TripleSim.ReactionGraph3X.BimolecularNode;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;

public class TripleSimWidestPathGraph {

	private Graph g;

	public TripleSimWidestPathGraph(Graph g) {
		this.g = g;
	}

	public void write(String string) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(string));
			
			{
				ArrayList<String> pairs = new ArrayList();
				for(BimolecularNode pair : g.allDockings.values()){
					if (pair.neighbors.isEmpty()){
						continue;
					}
					pairs.add(pair.index+" "+pair.associate[0].index+" "+pair.associate[1].index);
				}
				printLines(pairs,out);
			}
			
			{
				ArrayList<String> singles = new ArrayList();
				for(GraphNode node : g.allSingles.values()){
					singles.add(node.index+" "+node.priority+" "+node.structureString);
				}
				printLines(singles,out);
			}

			{
				ArrayList<String> edges = new ArrayList();
				for(GraphEdge e : g.edges){
					edges.add(e.reverse.towards.index+" "+e.towards.index+" "+e.k+" "+e.type);
					e = e.reverse;
					edges.add(e.reverse.towards.index+" "+e.towards.index+" "+e.k+" "+e.type);
				}
				printLines(edges,out);
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void printLines(ArrayList<String> pairs, PrintWriter out) {
		Collections.sort(pairs);
		for(String p : pairs){
			out.println(p);
		}
		out.println("END");
	}

}
