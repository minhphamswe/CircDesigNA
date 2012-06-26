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
package circdesigna;

import circdesigna.SequenceDesigner.SeqDesignerOption;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;


/**
 * A basic set of options for controlling the sequence designer.
 */
public class CircDesigNAOptions extends CircDesigNASystemElement{
	public CircDesigNAOptions(CircDesigNAConfig Std) {
		super(Std);
	}


	/**
	 * Sets up the default set of options.
	 */
	public static CircDesigNAOptions getDefaultOptions(CircDesigNAConfig Std) {
		return new CircDesigNAOptions(Std);
	}
	

	public SeqDesignerOption.Boolean rule_ccend_option = new SeqDesignerOption.Boolean(){
		public String getDescription() {
			return "Force domains to begin and end with either G or C. Other sequence constraints take priority.";
		}
		private boolean rule_ccend = getDefaultState();
		public boolean getState() {
			return rule_ccend;
		}
		public synchronized void toggle() {
			rule_ccend = !rule_ccend;
		}
		public boolean getDefaultState() {
			return false;
		}
		public void setState(boolean state) {
			rule_ccend = state;
		}
	};
	
	public SeqDesignerOption.Double end_score_threshold = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Stopping condition (solution has lower score)";
		}
		public double getDefaultState(){
			return 0;
		}
		private double EndThreshold = getDefaultState(); 
		public double getState() {
			return EndThreshold;
		}
		public synchronized void setState(double newVal) {
			EndThreshold = newVal;
		}
	};
	
	public SeqDesignerOption.Boolean globalSearch = new SeqDesignerOption.Boolean(){
		public String getDescription() {
			return "Global search (turn off to search for nearby solutions)";
		}
		private boolean sort_markings = getDefaultState();
		public boolean getState() {
			return sort_markings;
		}
		public synchronized void toggle() {
			sort_markings = !sort_markings;
		}
		public boolean getDefaultState() {
			return true;
		}
		public void setState(boolean state) {
			sort_markings = state;
		}
	};	
	
	public SeqDesignerOption.Double bimolecularPenalty = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "\u0394G intermolecular (in kcal/mol). Ref: Zuker, 2003";
		}
		public double getDefaultState(){
			if (Std.isDNAMode()){
				return 1.96;
			} else {
				return 1.70; //RNA has less of a penalty
			}
		}
		private double time = getDefaultState(); 
		public double getState() {
			return time;
		}
		public synchronized void setState(double newVal) {
			time = newVal;
		}
	};
	
	public SeqDesignerOption.Integer population_size = new SeqDesignerOption.Integer(){
		public String getDescription() {
			return "Population size when genetic algorithm is used";
		}
		public int getDefaultState(){
			return 30;
		}
		private int population_size = getDefaultState(); 
		public int getState() {
			return population_size;
		}
		public synchronized void setState(int newVal) {
			if (newVal <= 0){
				throw new RuntimeException("Population size is > 0");
			}
			population_size = newVal;
		}
	};
	
	public SeqDesignerOption.Boolean random_design = new SeqDesignerOption.Boolean(){
		public String getDescription() {
			return "Randomize all population members each iteration.";
		}
		private boolean randomize = getDefaultState();
		public boolean getState() {
			return randomize;
		}
		public synchronized void toggle() {
			randomize = !randomize;
		}
		public boolean getDefaultState() {
			return false;
		}
		public void setState(boolean state) {
			randomize = state;
		}
	};

	
	//Make sure to update this please.
	public final SeqDesignerOption[] options = new SeqDesignerOption[]{
			bimolecularPenalty, rule_ccend_option, globalSearch, population_size, end_score_threshold 
	};
	
}
