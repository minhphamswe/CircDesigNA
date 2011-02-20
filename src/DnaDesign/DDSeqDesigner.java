package DnaDesign;


/**
 * A model of a running sequence designer (compiler). It is designed to provide real-time updates
 * about its state.
 * 
 * This API is used by the graphic user interface to display the real-time state of the designer.
 */
public interface DDSeqDesigner<T extends DesignerOptions> {
	public interface SeqDesignerOption{
		public String getDescription();
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
	public int getIterationCount();
	/**
	 * Use System.getProperty("line.separator") to split multiline result into individual lines.
	 */
	public String getResult();
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
}
