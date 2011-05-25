package circdesigna.energy;

import java.io.PrintWriter;

import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.Config.CircDesigNASystemElement;
import DnaDesign.plugins.RunUnafoldTool.UnafoldFoldEntry;
import DnaDesign.plugins.RunUnafoldTool.UnafoldRunner;
import edu.utexas.cssb.circdesigna.DomainSequence;

/**
 * Implements MFE prediction and folding score functions
 */
public class UnafoldFolder extends CircDesigNASystemElement implements NAFolding{
	
	/**
	 * Constructors, define parameters and / or a configuration.
	 */
	private NAExperimentDatabase eParams;
	public UnafoldFolder(CircDesigNAConfig sys){
		super(sys);
		eParams = new ExperimentalDuplexParamsImpl(sys);
	}
	public UnafoldFolder(NAExperimentDatabase params, CircDesigNAConfig sys){
		super(sys);
		eParams = params;
	}
	

	/**
	 * UNAFOLD extension: send a request out to unafold.
	 */
	public double mfe(DomainSequence ds, DomainSequence ds2, int[][] domain, int[][] domain_markings) {
		UnafoldRunner ufr = new UnafoldRunner(Std.isDNAMode()?"DNA":"RNA");

		int len = ds.length(domain);
		int len2 = ds2.length(domain);
		PrintWriter out = new PrintWriter(ufr.getArgsFile(0));
		{
			StringBuffer create = new StringBuffer();
			for(int k = 0; k < len; k++){
				create.append(Std.monomer.displayBase(base(ds, k, domain)));
			}
			out.println(">A");
			out.println(create.toString());
		}
		out.close();

		out = new PrintWriter(ufr.getArgsFile(1));
		{
			StringBuffer create = new StringBuffer();
			for(int k = 0; k < len2; k++){
				create.append(Std.monomer.displayBase(base(ds2, k, domain)));
			}
			out.println(">B");
			out.println(create.toString());
		}
		out.close();
		
		double val = 0;
		double PERFECTscore = 0;
		try {
			ufr.runHybridizedJob();
			final UnafoldFoldEntry next = ufr.getResults().iterator().next();
			val = next.mfeDG;
			val = Math.min(val,PERFECTscore);
			if (val == PERFECTscore){ //Unafold does not give structures in this case.
				return val;
			}
			//We have structure:
			for(int k = 0; k < len; k++){
				if (next.pairs[k]>0){
					ds.mark(k, domain, domain_markings);
				}
			}
			for(int k = 0; k < len2; k++){
				if (next.pairs[k+len]>0){
					ds2.mark(k, domain, domain_markings);
				}
			}
			
		} catch( Throwable e){
			e.printStackTrace();
		}
		return val;
	}
	/**
	 * Unafold extension: send a request out to unafold
	 */
	public double mfe(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		UnafoldRunner ufr = new UnafoldRunner(Std.isDNAMode()?"DNA":"RNA");

		int len = seq.length(domain);
		PrintWriter out = new PrintWriter(ufr.getArgsFile(0));
		{
			StringBuffer create = new StringBuffer();
			for(int k = 0; k < len; k++){
				create.append(Std.monomer.displayBase(base(seq, k, domain)));
			}
			out.println(">A");
			out.println(create.toString());
		}
		out.close();
		
		double val = 0;
		double PERFECTscore = 0;
		try {
			ufr.runSingleStrandedJob();
			final UnafoldFoldEntry next = ufr.getResults().iterator().next();
			val = next.mfeDG;
			val = Math.min(val,PERFECTscore);
			if (val == PERFECTscore){ //Unafold does not give structures in this case.
				return val;
			}
			//We have structure:
			for(int k = 0; k < len; k++){
				if (next.pairs[k]>0){
					seq.mark(k, domain, domain_markings);
				}
			}
		} catch( Throwable e){
			e.printStackTrace();
		}
		return val;
		
		/*
		StringBuffer create = new StringBuffer();
		int len = seq.length(domain);
		for(int k = 0; k < len; k++){
			create.append(Std.monomer.displayBase(base(seq, k, domain)));
		}
		String str = create.toString();
		try {
			Process p = Runtime.getRuntime().exec(absPathToHybridSSMinMod+str);
			Scanner in = new Scanner(p.getInputStream());
			double val = 0;
			double PERFECTscore = 0;
			try {
				while(in.hasNextLine()){
					val = new Double(in.nextLine());
					val = Math.min(val,PERFECTscore);
					return val;
				}
			} finally {
				if (val == PERFECTscore){ //"Infinite" Free Energy (?)
					return val;
				}
				in.nextLine(); //Read off "dg" line
				for(int k = 0; k < len; k++){
					char[] arr = in.nextLine().toCharArray();
					//System.out.println(new String(arr));
					int regions = 0;
					int z, end;
					for(z = 0; z < arr.length && regions < 4; z++){
						if (arr[z]=='\t'){
							regions++;
						}
					}
					for(end = z+1; end < arr.length; end++){
						if (arr[end]=='\t'){
							break;
						}
					}
					//System.out.println(new String(arr));
					int num = new Integer(new String(arr,z,end-z));
					//System.out.println(num);
					if (num > 0){
						seq.mark(num-1, domain, domain_markings);
					}
					//Thread.sleep(100);
				}
				in.close();
				p.waitFor();
			}
		} catch (Throwable e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
		*/
	}
	public double mfeNoDiag(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings) {
		throw new RuntimeException("Not available with this folding tool.");
	}
	public double mfeStraight(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings, int markLeft, int markRight, int joffset) {
		throw new RuntimeException("Not available with this folding tool.");
	}
	public void pairPr(double[][] pairsOut, DomainSequence seq1, DomainSequence seq2, int[][] domain) {
		throw new RuntimeException("Not available with this folding tool.");
	}
	public void pairPr(double[][] pairsOut, DomainSequence seq1, int[][] domain) {
		throw new RuntimeException("Not available with this folding tool.");
	}
}
