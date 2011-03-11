package DnaDesign;

import DnaDesign.DDSeqDesigner.SeqDesignerOption;


/**
 * A basic set of options for controlling the sequence designer.
 */
public class DesignerOptions {
	/**
	 * Sets up the default set of options.
	 */
	public static DesignerOptions getDefaultOptions() {
		return new DesignerOptions();
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
			return true;
		}
		public void setState(boolean state) {
			rule_ccend = state;
		}
	};
	
	public SeqDesignerOption.Double end_score_threshold = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Stop design when a design candidate with a score less than this value is found.";
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
	
	public SeqDesignerOption.Boolean standardUseGA = new SeqDesignerOption.Boolean(){
		public String getDescription() {
			return "Allow regressive mutations (more standard genetic algorithm, doubles memory usage)";
		}
		private boolean sort_markings = getDefaultState();
		public boolean getState() {
			return sort_markings;
		}
		public synchronized void toggle() {
			sort_markings = !sort_markings;
		}
		public boolean getDefaultState() {
			return false;
		}
		public void setState(boolean state) {
			sort_markings = state;
		}
	};	
	
	public SeqDesignerOption.Double resourcePerMember = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Time to spend in local search loop on each population member before reporting the best candidate or reproducing the fittest member. Defines an \"iteration.\" When time limit is exceeded, any running score evaluations are allowed to finish.";
		}
		public double getDefaultState(){
			return .1;
		}
		private double time = getDefaultState(); 
		public double getState() {
			return time;
		}
		public synchronized void setState(double newVal) {
			if (newVal < 0){
				throw new RuntimeException("Error: time < 0");
			}
			time = newVal;
		}
	};

	
	public SeqDesignerOption.Double bimolecularPenalty = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Delta G (in kcal per mol) of intermolecular structure formations. Ref: Zuker, 2003";
		}
		public double getDefaultState(){
			return 1.96;
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
			return "Population size";
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
	
	public SeqDesignerOption.Integer selfSimilarityPenalty = new SeqDesignerOption.Integer(){
		public String getDescription() {
			return "Apply the \"Self Similarity\" penalty only to domains of length greater than or equal to this value. Negative input disables the penalty for all domains.";
		}
		public int getDefaultState(){
			return 20;
		}
		private int minimumForSelfSimilarity = getDefaultState(); 
		public int getState() {
			return minimumForSelfSimilarity;
		}
		public synchronized void setState(int newVal) {
			this.minimumForSelfSimilarity = newVal;
		}
	};

	
	//Make sure to update this please.
	public final SeqDesignerOption[] options = new SeqDesignerOption[]{
			rule_ccend_option, end_score_threshold, standardUseGA, resourcePerMember, population_size, selfSimilarityPenalty, bimolecularPenalty
	};
	
}
