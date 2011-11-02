package circdesigna.TripleSim;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import circdesigna.DomainPolymerGraph;
import circdesigna.TripleSim.ReactionGraph3X.BimolecularNode;
import circdesigna.TripleSim.ReactionGraph3X.Graph;
import circdesigna.TripleSim.ReactionGraph3X.GraphEdge;
import circdesigna.TripleSim.ReactionGraph3X.GraphNode;
import circdesigna.TripleSim.VeryGoodLA.LUFactorization;
import circdesigna.TripleSim.VeryGoodLA.Matrix;
import circdesigna.TripleSim.VeryGoodLA.MatrixOp;

public class TripleSim {
	/**
	 * Gets the maximal duplexes of q (that is, no duplex is a subduplex of any other.)
	 */
	public static int[][] getDuplexes(DomainPolymerGraph q) {
		ArrayList<int[]> duplexes = new ArrayList();

		for(int i = 0; i < q.length(); i++){
			int[] duplex = new int[4];
			duplex[0] = i;
			duplex[3] = q.getDomainPair(duplex[0]);
			if (duplex[3] < 0){
				continue;
			}
			int j = i+1;
			for(; j < q.length() && q.getDomain(j) >= 0 && q.getDomainPair(j) >= 0; j++){
				if(q.getDomainPair(j) != q.getDomainPair(j-1)-1){
					break;
				}
			}
			duplex[1] = j-1;
			duplex[2] = q.getDomainPair(duplex[1]);
			duplexes.add(duplex);
			i = duplex[1];
		}

		int[][] toRet = new int[duplexes.size()][];
		duplexes.toArray(toRet);
		return toRet;
	}
	public static ArrayList<DomainPolymerGraph> splitIntoComponents(
			DomainPolymerGraph q) {
		ArrayList<DomainPolymerGraph> toRet = new ArrayList();
		
		ArrayList<ArrayList<int[]>> connected = getConnectedComponents(q);
		for(ArrayList<int[]> bin : connected){
			StringBuffer substring = new StringBuffer("AC ");
			Collections.sort(bin, new Comparator<int[]>(){
				public int compare(int[] o1, int[] o2) {
					return o1[0] - o2[0];
				}
			});
			for(int[] e : bin){
				substring.append(q.getStructureString(e[0], e[1]));
			}
			DomainPolymerGraph product = new DomainPolymerGraph(q.getDomainDefs());
			DomainPolymerGraph.readStructure(substring.toString(), product);
			product.annotate(q.getAnnotationTree());
			
			toRet.add(product);
		}
		return toRet;
	}
	public static ArrayList<ArrayList<int[]>> getConnectedComponents(DomainPolymerGraph q) {
		ArrayList<ArrayList<int[]>> connected = new ArrayList<ArrayList<int[]>>();
		int[][] strands = q.getStrands();
		for(int[] e : strands){
			ArrayList<int[]> placed = null;
			for(int inBin = 0; inBin < connected.size(); inBin++){
				ArrayList<int[]> bin = connected.get(inBin);
				for(int[] other : bin){
					if (placed==null){
						if (connected(q,other,e)){
							bin.add(e);
							placed = bin;
							break;
						}
					} else {
						if (connected(q,other,e)){
							placed.addAll(bin);
							connected.remove(inBin--);
							break;
						}
					}
				}
			}
			if (placed == null){
				ArrayList<int[]> neu = new ArrayList();
				neu.add(e);
				connected.add(neu);
			}
		}
		return connected;
	}
	private static boolean connected(DomainPolymerGraph q, int[] input, int[] output) {
		for(int k = input[0]; k < input[1]; k++){
			int target = q.getDomainPair(k);
			if (target >= output[0] && target < output[1]){
				return true;
			}
		}
		return false;
	}
	private static double maxAbs(double[] slope) {
		double m = 0;
		for(double q : slope){
			m = Math.max(m,Math.abs(q));
		}
		return m;
	}
	private double[] lastEndConcentrations;
	public void printLastConcentrations(Graph g){
		GraphNode[] nums = new GraphNode[lastEndConcentrations.length];
		for(GraphNode q : g.allSingles.values()){
			if (q.index < nums.length){
				nums[q.index] = q;
			}
		}
		Arrays.sort(nums, new Comparator<GraphNode>(){
			public int compare(GraphNode o1, GraphNode o2) {
				return -(int)Math.signum(lastEndConcentrations[o1.index] - lastEndConcentrations[o2.index]);
			}
		});
		for(int i = 0; i < nums.length && i < 50; i++){
			System.out.printf("%s %.3e",nums[i].structureString,lastEndConcentrations[nums[i].index]);
			System.out.println();
		}
	}
	/**
	 * Uses the Gragg-Burlish-Stoer multidimensional arbitrary order ODE numerical solver.
	 * 
	 * Trashes the "index" field of all nodes in g for internal use.
	 * 
	 * Stops the simulation after the iteration when at least one unvisited node reaches priority > ignorePriority, or tf time is reached.
	 * 
	 * Returns the actual time that was simulated up to.
	 */
	public double updatePriorities(Graph g, double epsilon, double tf, PrintWriter out, double ignorePriority) {
		double t0 = 0; 
		
		//Reset priorities of all nodes, and set up initial concentration vector.
		int n = g.allSingles.size();
		double[] y = new double[n];
		for(GraphNode u : g.allSingles.values()){
			u.priority = y[u.index] = u.initialConc;

		}
		if (out!=null){
			out.printf("%-10s ","Step");
			out.printf("%-10s ","Time");
			for(int i = 0; i < y.length; i++){
				out.printf("%-10d ",i);
			}
			out.println();
		}
		for(PulseEvents e : g.events){
			e.reset();
		}
				
		//Uses an implicit (stiff) method to integrate the system.
		double step = 1e-4;
		while(t0 < tf){
			boolean hasSolution = false;
			double[] yhat = new double[y.length];

			//System.out.print(t0+" ");
			big:while(!hasSolution){
				boolean invalidState = false;
				double acceptableError = epsilon / (tf / step);
				
				boolean diverged = !EulerTrapezoidal(g, step, yhat, new double[][]{y}, acceptableError);
				if (diverged){
					invalidState = true;
					//throw new RuntimeException("Newton Raphson Non-Convergence with stepsize "+step);
				} else {
					for(double q : yhat){
						if (q < 0){
							invalidState = true;
							//System.out.print("O");
						}
					}
				}

				if (invalidState){
					step /= 2;
				} else {
					step *= 2;
					hasSolution = true;
				}
			}
			
			//Move to next timestep and update priorities
			y = yhat;
			t0 += step;
			for(PulseEvents e : g.events){
				double[] shock = e.handlePulses(y, t0, step);
				t0 = shock[0];
				step = shock[1];
			}
			
			for(GraphNode u : g.allSingles.values()){
				u.priority = Math.max(u.priority,y[u.index]);
				if (Double.isNaN(u.priority)){
					throw new RuntimeException("Undefined priorities.");
				}
				
				if (!u.visited && u.priority > ignorePriority){
					tf = t0; //End simulation.
				}
			}
			
			if (out!=null){
				out.printf("%-10.3e ",step);
				out.printf("%-10.3e ",t0);
				for(double q : y){
					out.printf("%-10.3e ",q);
				}
				out.println();
				out.flush();
			}
		}
		//System.out.println();
		
		//Update data structures to respond to changes in priority.
		g.recreateUnvisited();
		
		lastEndConcentrations = y;
		
		return t0;
	}
	/**
	 * Solves the system g(yhat) = yhat - y - h*f(yhat) = 0 for yhat.
	 * Return true if solution was successful.
	 * 
	 * 
	 */
	private boolean EulerTrapezoidal(Graph g, double h, double[] yhat_, double[][] y_, double epsilon){
		//Evaluate slope at each of the y_
		Matrix[] fys = new Matrix[y_.length];
		for(int k = 0; k < fys.length; k++){
			double[] slope = new double[yhat_.length];
			fys[k] = new Matrix(slope);
			F(slope,g,y_[k]);
		}
		
		//Initial guess: yhat = y + h f(y).
		Matrix y = new Matrix(y_[0]); //BACKED!
		Matrix yhat = new Matrix(yhat_); //BACKED!
		//System.arraycopy(y_[0], 0, yhat_, 0, y_[0].length);
		MatrixOp.add(h, fys[0], y, yhat);
		
		Matrix b = new Matrix(yhat.m,1);
		
		Matrix A = new Matrix(yhat.m,yhat.m);
		
		Matrix I_A = Matrix.eye(A);
		
		//Evaluate jacobian of g at y = A, and its LU factorization.
		DF(g,y_[0],A);
		MatrixOp.add(-h, A, I_A, A); 
		LUFactorization LU = LUFactorization.LUFactorize(A,false);
		
		double[] slope = new double[yhat_.length];
		Matrix fyhat = new Matrix(slope);
		
		Matrix x = new Matrix(yhat.m,1);
		Matrix yhat_new = new Matrix(yhat.m,1);
		Matrix delta = new Matrix(yhat.m,1);
		double lastError = Double.MAX_VALUE;
		for(int i = 0; ; i++){
			//Evaluate f(yhat)
			F(slope, g, yhat_);

			//Evaluate g(yhat) = b = yhat - y - h / 2 f(yhat) - h / 2 f(y)
			MatrixOp.add(-1, y, yhat, b);
			MatrixOp.add(-h, fyhat, b, b);
			//MatrixOp.add(-h/2, fys[0], b, b);
			
			/*
			//Evaluate jacobian of g at yhat = A, and its LU factorization.
			DF(g,yhat_,A);
			MatrixOp.add(-h/2, A, I_A, A); 
			LUFactorization LU = LUFactorization.LUFactorize(A,false);
			*/
			
			//No dependence on h beyond this point
			
			//yhat_new = yhat - (x after solving Ax = b)
			LU.solve(b, x);
			MatrixOp.add(-1, x, yhat, yhat_new);
			
			//Compute delta, yhat_new - yhat,  
			MatrixOp.add(-1, yhat_new, yhat, delta);
			//Compute one-norm.
			double error = 0;
			for(int u = 0; u < delta.m; u++){
				error += Math.abs(delta.get(u, 0));
			}

			//Move to yhat_new.
			MatrixOp.add(0, yhat, yhat_new, yhat);
			
			//Break if delta suff small
			
			if (error <= epsilon){
				break;
			}
			//Return false if error is not decaying
			if (!(error < lastError)){
				return false;
			}
			lastError = error;
			//System.out.println("Error after step "+i+": "+ MatrixOp.twoNorm(delta));
		}
		return true;
	}
	/**
	 * Performs richardson extrapolation on estimates of y(t0 + h) using subdivision of the interval (t0, t0+h) and the T function.
	 * 
	 * Uses y(t0) = y, and writes y(t0+h) over yhat.
	 * @param epsilon 
	 */
	private boolean GBSCore(Graph g, double[] yhat, double[] y, double h, double epsilon) {
		int rows = 4;
		int columns = rows;
		int t = 2; 
		double[] k = new double[]{0,2,4,6,8,10,12,14,16};
		double[][][] tableu = new double[rows][columns][y.length];
		Arrays.fill(yhat,0);
		tableu[rows-1][columns-1] = yhat; 
		boolean toRet = true;
		int L = 1;
		for(int r = rows-columns; r < rows; r++){
			toRet &= T(g, tableu[r][0], y, h, L);
			L *= 2;
			h /= 2;
			for(int c = 1; c <= r; c++){
				axpy(tableu[r][c], tableu[r][c], -1/(Math.pow(t,k[c])-1), tableu[r-1][c-1]);
				axpy(tableu[r][c], tableu[r][c], Math.pow(t,k[c])/(Math.pow(t,k[c])-1), tableu[r][c-1]);
			}
			if (r > 0){
				//Compute error from last two elements of the r-1th column. If they differ less than epsilon (1 norm), then convergence is achieved.
				double error = 0;
				for(int u = 0; u < y.length; u++){
					error += Math.abs(tableu[r][r-1][u] - tableu[r-1][r-1][u]);
				}
				if (error < epsilon){
					System.arraycopy(tableu[r][r],0,tableu[rows-1][columns-1],0,y.length);
					return toRet;
				}	
			}
		}
		return false;
		//return toRet;
	}
	/**
	 * Result is written into y1.
	 * Returns false if some value of y1 goes negative.
	 */
	private boolean T(Graph g, double[] y1, double[] y0, double h, int L){
		h/=2;
		L*=2;
		
		double[] min2 = new double[y0.length];
		double[] min1 = new double[y0.length];
		double[] min0 = new double[y0.length];
		System.arraycopy(y0,0,min2,0,y0.length);
		boolean toRet = true;
		toRet &= axpfy(min1, min2, h, g, min2);
		//POST CONDITION: min2 and min1 mean what they should.
		for(int k = 1; k < L; k++){
			Arrays.fill(min0,0);
			toRet &= axpfy(min0, min2, 2*h, g, min1);
			
			//Cyclic shift the three buffers
			double[] temp = min2;
			min2 = min1;
			min1 = min0;
			min0 = temp;
			
			//POST CONDITION: min2 and min1 mean what they should.
		}

		//S*2 = y1 = min1 + min2 + h * f(min1)
		toRet &= axpfy(y1, min2, h, g, min1);
		toRet &= axpy(y1, min1, 1, y1);
		//Divide the whole thing by two.
		div(y1,y1,2);
		return toRet;
	}
	private void div(double[] y1, double[] y2, int L) {
		for(int k = 0; k < y1.length; k++){
			y1[k] = y2[k] / L;
		}
	}
	/**
	 * Computes a += x + h f(y), returns false if some value of a goes negative.
	 */
	private boolean axpfy(double[] a, double[] x, double h, Graph g, double[] y) {
		double[] slope = new double[x.length];
		F(slope, g, y);
		
		return axpy(a,x,h,slope);
	}
	private void DF(Graph g, double[] y, Matrix A) {
		for(int i = 0; i < A.m; i++){
			for(int j = 0; j < A.n; j++){
				A.set(i, j, 0);
			}
		}
		
		for(GraphEdge rxn : g.edges){
			gedmv(A,rxn,y);
			gedmv(A,rxn.reverse,y);
		}
	}
	private void F(double[] slope, Graph g, double[] y) {
		Arrays.fill(slope,0);
		for(GraphEdge rxn : g.edges){
			gemv(slope,rxn,y);
			gemv(slope,rxn.reverse,y);
		}
	}
	private void gedmv(Matrix a, GraphEdge rxn, double[] y){
		if (rxn.reverse.towards instanceof BimolecularNode){
			gedmv(a,((BimolecularNode)rxn.reverse.towards).associate,rxn,y);
		} else {
			gedmv(a,new GraphNode[]{rxn.reverse.towards},rxn,y);
		}	
	}
	private void gedmv(Matrix a, GraphNode[] reactants, GraphEdge rxn, double[] y) {
		for(GraphNode reactant : reactants){
			//Rate of reaction with respect to me:
			double rate = rxn.k;
			int iOccur = 0;
			for(GraphNode deriv : reactants){
				if (deriv != reactant){
					rate *= y[deriv.index];
				} else {
					iOccur ++;
				}
			}
			if (iOccur == 2){
				rate *= y[reactant.index] * 2; //Derivative of x^2 is 2x
			}
			if (iOccur > 2){
				throw new RuntimeException("Totally not handling reactions between more than 2 reactant species!");
			}
			
			for(GraphNode deriv : reactants){
				a.add(deriv.index,reactant.index,-rate);
			}
			if (rxn.towards instanceof BimolecularNode){
				for(GraphNode product : ((BimolecularNode)rxn.towards).associate){
					a.add(product.index,reactant.index,rate);
				}
			} else {
				a.add(rxn.towards.index,reactant.index,rate);
			}	
		}
	}
	/**
	 * Computes a += rxn( y )
	 */
	private void gemv(double[] a, GraphEdge rxn, double[] y) {
		double prop = rxn.k;
		if (rxn.reverse.towards instanceof BimolecularNode){
			for(GraphNode q : ((BimolecularNode)rxn.reverse.towards).associate){
				prop *= y[q.index];
			}
		} else {
			prop *= y[rxn.reverse.towards.index];
		}
		
		xpy(a,prop,rxn.towards);
		xpy(a,-prop,rxn.reverse.towards);
	}
	private void xpy(double[] a, double h, GraphNode towards) {
		if (towards instanceof BimolecularNode){
			for(GraphNode q : ((BimolecularNode)towards).associate){
				a[q.index] += h;
			}
		} else {
			a[towards.index] += h;
		}
	}
	/**
	 * Returns false if some value of a goes negative. 
	 */
	private boolean axpy(double[] a, double[] x, double h, double[] slope) {
		boolean toRet = true;
		/** AXPY **/
		for(int k = 0; k < a.length; k++){
			a[k] = x[k] + h * slope[k];
			if (a[k] < 0){
				toRet = false;
			}
		}	
		return toRet;
	}	
}
