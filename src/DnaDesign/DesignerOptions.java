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
			return "Force domains to begin / end with G or C, where not specified";
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
	};
	
	public SeqDesignerOption.Boolean sort_markings = new SeqDesignerOption.Boolean(){
		public String getDescription() {
			return "Create a priority queue of bases to mutate, and mutate the problematic ones first";
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
	};
	
	public SeqDesignerOption.Double end_score_threshold = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "When a solution with a score less than this value is found, design will stop.";
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
	
	public SeqDesignerOption.Double resourcePerMember = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Amount of computer time (in seconds) to apply random-walk optimizaiton on each member each iteration (in-progress mutations are allowed to complete). Negative toggles Infinite Resource Tournament";
		}
		public double getDefaultState(){
			return .1;
		}
		private double time = getDefaultState(); 
		public double getState() {
			return time;
		}
		public synchronized void setState(double newVal) {
			time = newVal;
		}
	};

	
	public SeqDesignerOption.Double bimolecularPenalty = new SeqDesignerOption.Double(){
		public String getDescription() {
			return "Delta G (in kcal per mol) of intermolecular structure formations";
		}
		public double getDefaultState(){
			return 1.513;
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
			return "Population size for \"Block Designer.\"";
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
			return "Minimum length of domain to apply \"Self Similarity\" penalty to. Negative input disables.";
		}
		public int getDefaultState(){
			return -1;
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
			rule_ccend_option, sort_markings, end_score_threshold, resourcePerMember, population_size, selfSimilarityPenalty
	};
	
}
