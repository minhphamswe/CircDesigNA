package circdesigna.energy;

/**
 * Defines constants for a model of nucleic acid interactions.
 */
public class KineticsDefinition {
	
	//UNITS ARE /M/s
	
	//In the following three cases, reverse kinetics are computed using e^(-deltaG / (RT)) = kon / koff
	
	//k(f) = ckf * sqrt(n / 6) for intermolecular association of a single duplex 
	public double ckf = 3.5e6;
	
	//kb = ckb / x^2, where x is the number of bases shifted during branch migration.
	public double ckb = 400;
	
	//k(f) = ckfi * sqrt(n / 6) for intramolecular association of a single duplex
	//can be found by multiplying bimolecular constant by local concentration
	public double ckfi = 3.5e6;
	
	//In kcal / mol
	public double R = .0019858775;
		
	//In kelvins
	public double T = 298;
}
