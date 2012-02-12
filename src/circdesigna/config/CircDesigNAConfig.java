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
package circdesigna.config;

import static circdesigna.abstractpolymer.DnaDefinition.P;
import static circdesigna.abstractpolymer.DnaDefinition.Z;
import circdesigna.BannedPatterns;
import circdesigna.DesignSequenceConstraints;
import circdesigna.abstractpolymer.DnaDefinition;
import circdesigna.abstractpolymer.MonomerDefinition;
import circdesigna.abstractpolymer.RnaDefinition;
import circdesigna.energy.KineticsDefinition;
import circdesigna.impl.CodonCode;
import circdesigna.impl.SequenceCode;

/**
 * A configuration class, a reference to an instance of this class is passed
 * around to all CircDesigNASystemElement -s
 */
public class CircDesigNAConfig {
	public static final int DNA_MODE = 0, RNA_MODE = DNA_MODE+1;
	private int current_mode = DNA_MODE;
	private boolean saveReactionDescriptions;

	public CircDesigNAConfig(){
		setMode(current_mode);

		//Codon table information.
		customCodonTable = new CodonCode(this).defaultTable();
		bannedPatternsList = new BannedPatterns(this).defaultBannedWords();
		kinetics = new KineticsDefinition();
		saveReactionDescriptions = false;
	}

	public boolean saveReactionDescriptions(){
		return saveReactionDescriptions;
	}
	public void setSaveReactionDescriptions(boolean val){
		saveReactionDescriptions = val;
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
	public String bannedPatternsList;
	public KineticsDefinition kinetics;

	public double getDeltaGPerStackPair() {
		if (isDNAMode()){
			return -1.6;
		} else {
			return -2.0;
		}
	}
	
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

