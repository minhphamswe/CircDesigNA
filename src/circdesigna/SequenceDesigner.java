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


/**
 * A model of a running sequence designer (compiler). It is designed to provide real-time updates
 * about its state.
 * 
 * This API is used by the graphic user interface to display the real-time state of the designer.
 * It also provides interactivity, such as abort, resume, pause, etc.
 */
public interface SequenceDesigner<T extends CircDesigNAOptions> {
	public interface SeqDesignerOption{
		public java.lang.String getDescription();
		public interface Boolean extends SeqDesignerOption {
			public boolean getState();
			public boolean getDefaultState();
			public void toggle();
			public void setState(boolean state);
		}
		public interface Double extends SeqDesignerOption{
			public double getState();
			public double getDefaultState();
			public void setState(double newVal);
		}
		public interface Integer extends SeqDesignerOption{
			public int getState();
			public int getDefaultState();
			public void setState(int newVal);
		}
		public interface Str extends SeqDesignerOption{
			public java.lang.String getState();
			public java.lang.String getDefaultState();
			public void setState(java.lang.String newVal);
		}
	}
	/**
	 * Returns the generic options object for this designer.
	 */
	public T getOptions();
	/**
	 * Returns the score value of the best solution.
	 */
	public double getBestScore();
	/**
	 * Returns some representation of how many steps the algorithm has taken.
	 * Not enforced to be of any specific standard unit.
	 * Just make it monotonic increasing with time, please;
	 */
	public int getCurrentIteration();
	/**
	 * Use System.getProperty("line.separator") to split multiline result into individual lines.
	 * Should return the Best Result.
	 * 
	 * Must be equal to getResult(getAlternativeResult(BEST_RESULT));
	 */
	public String getResult();
	
	public static abstract class AlternativeResult {
		public static final int 
				BEST = 0, 
				FARTHEST1 = BEST+1, 
				ERROR = FARTHEST1+1,
				LOG = ERROR + 1,
				OTHER=LOG + 1;
		public int TYPE;
		//Used as a "subtype"
		public int ID;
		public abstract String getDescription();
		public String toString(){
			return getDescription();
		}
	}
	
	//Default. To get a specific population member, use 
	public AlternativeResult[] getAlternativeResults();
	public String getResult(AlternativeResult alternative);
	/**
	 * Returns the shared score reporter object (used for generating real time visuals)
	 */
	public DesignIntermediateReporter getDir();
	
	public boolean isRunning();
	public boolean isFinished();
	/**
	 * Returns true if the designer ended due to an error.
	 */
	public boolean isEndConditionError();
	/**
	 * Starts if stopped, resumes if paused
	 */
	public void resume();
	/**
	 * Stops is running.
	 */
	public void pause();
	/**
	 * Stops in all cases; cannot be resumed; cleanup occurs immediately.
	 */
	public void abort();
	/**
	 * resumes / starts the runner, and pauses it after one iteration has completed.
	 * The current iteration can be gotten with getIterationCount() 
	 */
	public void runIteration();
	/**
	 * Returns the current score breakdown.
	 */
	public DesignScoreBreakdown getScoreBreakdown();
}
