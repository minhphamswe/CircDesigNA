package DnaDesign;

import java.util.List;

public interface DDSeqDesigner {
	public interface SeqDesignerOption{
		public String getDescription();
		public boolean getState();
		public void toggle();
	}
	public List<SeqDesignerOption> getOptions();
	/**
	 * Use System.getProperty("line.separator") to format multiline result
	 */
	public String getResult();
	
	public boolean isRunning();
	public boolean isFinished();
	public float statusVal();
	public float scoreVal();
	
	/**
	 * Starts if stopped, resumes if paused
	 */
	public void resume();
	
	public void pause();
	
	public void abort();
}
