package DnaDesign.Batch;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import DnaDesign.DomainStructureData;
import DnaDesign.Config.CircDesigNAConfig;
import static DnaDesign.Batch.DesignMultipleTimes.*;

public class ExtractScoresFromRun {
	public static void main(String[] args) throws IOException{
		Scanner in = new Scanner(System.in);
		String file = in.nextLine();

		String Domains = readToEnd(in);
		String Molecules = readToEnd(in);
		
		DomainStructureData dsd = new DomainStructureData(new CircDesigNAConfig());
		dsd.readDomainDefs(Domains, dsd);
		
		DesignMultipleTimes.RunEvaluation(new File(file), Molecules, dsd, -1, false);
	}
}
