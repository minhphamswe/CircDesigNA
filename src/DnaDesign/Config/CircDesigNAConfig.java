package DnaDesign.Config;

import static DnaDesign.AbstractPolymer.DnaDefinition.P;
import static DnaDesign.AbstractPolymer.DnaDefinition.Z;
import edu.utexas.cssb.circdesigna.DesignSequenceConstraints;
import DnaDesign.AbstractPolymer.DnaDefinition;
import DnaDesign.AbstractPolymer.MonomerDefinition;
import DnaDesign.AbstractPolymer.RnaDefinition;
import DnaDesign.impl.CodonCode;
import DnaDesign.impl.SequenceCode;

/**
 * A configuration class, a reference to an instance of this class is passed
 * around to all CircDesigNASystemElement -s
 */
public class CircDesigNAConfig {
	public static final int DNA_MODE = 0, RNA_MODE = DNA_MODE+1;
	private int current_mode = DNA_MODE;

	public CircDesigNAConfig(){
		setMode(current_mode);

		//Codon table information.
		customCodonTable = new CodonCode(this).defaultTable();
	}

	/**
	 * Returns true if the current mode deals with nucleic acids
	 */
	public boolean isNAmode() {
		return true;
	}
	public void setMode(int mode){
		switch(mode){
		case DNA_MODE:
			monomer = new DnaDefinition();
			break;
		case RNA_MODE:
			monomer = new RnaDefinition();
			break;
		default:
			throw new RuntimeException("Invalid mode");
		}
		current_mode = mode;
	}

	public String getParameterName() {
		if (isDNAMode()){
			return "DNA_mfold2.3";
		} else {
			return "RNA_mfold3.0";
		}
	}

	public boolean isDNAMode() {
		return current_mode==DNA_MODE;
	}

	//Products
	public MonomerDefinition monomer;
	public String customCodonTable;

	public DesignSequenceConstraints getDefaultConstraints() {
		DesignSequenceConstraints defaults = new DesignSequenceConstraints(this);

		//defaults.setMaxConstraint(0, H);
		defaults.setMaxConstraint(0, P);
		defaults.setMaxConstraint(0, Z);
		//defaults.setMaxConstraint(0, D);

		return defaults;
	}

	public SequenceCode getDefaultSequenceCode() {
		SequenceCode init = new SequenceCode();
		init.setConstraints(getDefaultConstraints());
		return init;
	}
}
