/*
  Part of the CircDesigNA Project - http://cssb.utexas.edu/circdesigna
  
  Copyright (c) 2010-11 Ben Braun
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation, version 2.1.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General
  Public License along with this library; if not, write to the
  Free Software Foundation, Inc., 59 Temple Place, Suite 330,
  Boston, MA  02111-1307  USA
*/
package circdesigna.geometry;

import java.util.Arrays;

/**
 * Given an irregular polygon with side lengths specified (and ordered), this class calculates
 * the radius of a circle, R, and a set of Theta_i where the Theta_i are the exterior angles of the
 * polygon with the specified side lengths such that a circle of radius R circumscribes the polygon.
 * 
 * Useful for rendering DNA molecules.
 * 
 * Uses a newton raphson iteration procedure to approximate the root of a nonlinear equation in one 
 * variable, R. The Theta_i are uniquely defined once R is calculated.
 */
public class CircumscribedPolygonTool {
	public static class CircumscribedPolygon {
		int count = 0;
		public void resetCounter(){
			count = 0;
		}
		public float[] S;
		public float[] theta;
		public float[] interiorTheta;
		public float[] dIThetaDR;
		public float R;
		public float Rhi;
		public float f, fprime;
		public float next() {
			if (count >= theta.length){
				return 0;
			}
			return theta[count++];
		}
	}
	/**
	 * Modifies input to fill in the theta and R fields.
	 */
	public static void solvePolygonProblem(CircumscribedPolygon input){
		if (input == null){
			return;
		}
		//Make a guess at a good value of R:
		//NOTATION:
		//S: a side length
		//THETA: an arc length
		//R: the radius of the circumscribing circle
		int n = input.S.length;
		
		float avgS = 0;
		for(float k : input.S){
			avgS += k;
		}
		avgS /= n;
		input.R = (float) (avgS/Math.sqrt(2*(1-Math.cos(2*Math.PI/n))));
		
		float maxS = 0;
		for(float k : input.S){
			maxS = Math.max(maxS,k);
		}
		//input.R = (float) maxS / 2 + 0.01f;
		
		/*
		float minS = Float.MAX_VALUE;
		float sumS = 0;
		for(float k : input.S){
			sumS += k;
			minS = Math.min(minS,k);
		}
		input.R = minS / 2;
		input.Rhi = (float) (sumS / (Math.PI*2))*2;
		*/

		input.interiorTheta = new float[n];
		input.dIThetaDR = new float[n]; //For newton raphson.
		input.theta = new float[n];
		for(int i = 0; i < 8; i++){
			input.R = Math.max(input.R,maxS / 2);
			iterateNewtonRaphson(input);
			//R at this point DISAGREES with the theta / interiorTheta.
		}
		
		for(int k = 0; k < n; k++){
			input.theta[k] = (input.interiorTheta[k]+input.interiorTheta[(k+1)%n])/2;
		}
		
		//System.out.println(Arrays.toString(input.theta));
	}
	private static void iterateNewtonRaphson(CircumscribedPolygon input) {
		int n = input.S.length;
		/*
		float mid = (input.Rhi+input.R)/2;
		compute(input, mid);
		float fmid = input.f;
		//System.out.println(input.R+" "+fmid+" "+input.Rhi+" ");
		if (fmid < 0){
			input.R = mid;
		} else {
			input.Rhi = mid;
		}
		*/
		
		
		//Newton raphson
		compute(input,input.R);
		//System.out.println(input.R+" "+input.f);
		if (Float.isNaN(input.f)){
			input.R += 1;
		} else {
			input.R -= input.f / input.fprime;
		}
		
		//BST:
		
	}
	private static void compute(CircumscribedPolygon input, float r) {
		input.f = -(float) (2*Math.PI);
		input.fprime = 0;
		int n = input.S.length;
		for(int k = 0; k < n; k++){
			float a = (float) (1-Math.pow(input.S[k]/r,2)/2);
			float dA = -1/2f*input.S[k]*input.S[k];
			dA *= -2/(r*r*r);
			input.interiorTheta[k] = (float) Math.acos(a);
			input.dIThetaDR[k] = (float)(-1/(Math.sqrt(1-a*a)) * dA);
			
			input.f += input.interiorTheta[k];
			input.fprime += input.dIThetaDR[k];
		}
	}
}
