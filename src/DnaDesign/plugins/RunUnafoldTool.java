package DnaDesign.plugins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Scanner;
import java.util.TreeMap;

/**
 * Utilities for using the UNAFOLD package.
 * @author Benjamin
 */
public class RunUnafoldTool {
	public static class UnafoldFoldEntry {
		public UnafoldFoldEntry(String strA) {
			this.nameFirst = strA;
			this.nameSecond = null;
		}
		public UnafoldFoldEntry(String strA, String strB) {
			this.nameFirst = strA;
			this.nameSecond = strB;
		}
		public double mfeDG;
		public String nameFirst;
		public String nameSecond;
		/**
		 * -1: Unpaired.
		 * All other values: Paired with that base (0 indexed).
		 */
		public int[] pairs;
	}
	
	public static void loadPaths(){
		if (absPathToHybridMinMod==null){
			absPathToHybridSSMinMod = "\""+new File(System.getProperty("UNAFOLDHOME")).getAbsolutePath()+"hybrid-ss-min.exe\"";
			absPathToHybridMinMod = "\""+new File(System.getProperty("UNAFOLDHOME")).getAbsolutePath()+"hybrid-min.exe\"";
			System.out.println("Paths to unafold: "+absPathToHybridSSMinMod);
		}
	}
	private static String absPathToHybridSSMinMod =  null; //"\"C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\unafold\\hybrid-ss-min.exe\" --NA=DNA ";
	private static String absPathToHybridMinMod = null; //"\"C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\unafold\\hybrid-min.exe\" --NA=DNA ";
	
	public static class UnafoldRunner{
		private ArrayList<UnafoldFoldEntry> data = new ArrayList();
		private TreeMap<Integer, File> argsFiles = new TreeMap();
		private String na_arg = "DNA";
		public UnafoldRunner(String na){
			na_arg = na;
		}
		public OutputStream getArgsFile(int i) {
			File toRet = argsFiles.get(i);
			if (toRet==null){
				toRet = new File("unafoldRunnerTmp"+i+".txt");
				argsFiles.put(i,toRet);
			}
			try {
				return new FileOutputStream(toRet);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				throw new RuntimeException("Could not open temporary file for unafold runner.");
			}
		}
		/**
		 * Input files have the format
		 * name sequence
		 * 
		 * for both files.
		 */
		public void runHybridizedJob() throws IOException, InterruptedException{
			loadPaths();
			Process p = Runtime.getRuntime().exec(absPathToHybridMinMod+" --NA="+na_arg+" \""+argsFiles.get(0).getAbsolutePath()+"\" \""+argsFiles.get(1)+"\"");
			runJob(p);
		}
		/**
		 * Input file has the format
		 * name sequence
		 */
		public void runSingleStrandedJob() throws IOException, InterruptedException{
			loadPaths();
			Process p = Runtime.getRuntime().exec(absPathToHybridSSMinMod+" --NA=DNA \""+argsFiles.get(0)+"\"");
			runJob(p);
		}
		public Collection<UnafoldFoldEntry> getResults(){
			return data;
		}
		public void clear(){
			data.clear();
			argsFiles.clear();
		}
		private void runJob(Process p){
			InputStream in2 = p.getInputStream();
			interpretFromStream(in2);
			try {
				p.waitFor();			
				p.getOutputStream().close();
				p.getInputStream().close();
				p.getErrorStream().close();
			} catch (Throwable e){
				System.err.println("I/O Error in Unafold Extension:");
				e.printStackTrace();
			}
		}
		public void interpretFromStream(InputStream in2){
			Scanner in = new Scanner(in2);
			
			while(in.hasNextLine()){
				//Calculating for 0 and 0, t = 37
				String[] line = in.nextLine().split("\\s+|[,]");
				
				String strA = line[2];
				String strB = line[4];
				
				UnafoldFoldEntry newEntry = new UnafoldFoldEntry(strA,strB);
				data.add(newEntry);
				
				line = in.nextLine().split("\\s+");
				int num = new Integer(line[0]);
				//39	dG = -17.5	dH = -100.9	0-0
				double dg = new Double(line[3]);
				newEntry.mfeDG = dg;
				newEntry.pairs = new int[num];
				for(int k = 0; k < num; k++){
					line = in.nextLine().split("\\s+");
					newEntry.pairs[k] = new Integer(line[4]) - 1;
				}
			}
		}
	}
}
