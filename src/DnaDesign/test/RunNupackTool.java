package DnaDesign.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

public class RunNupackTool {
	public static void runNupack(String seqs, String concs, int maximumComplexSize, String prefix, boolean concentrations, File nupackDir) throws IOException{
		System.err.println(nupackDir.getAbsolutePath()+"/"+prefix);
		nupackDir.mkdirs();
		File nupackList = new File("nupackTest/"+prefix+".in");

		PrintWriter nuPackListOut = new PrintWriter(new FileWriter(nupackList));
		nuPackListOut.print(seqs);
		nuPackListOut.println(maximumComplexSize);
		nuPackListOut.close();

		PrintWriter conOut = new PrintWriter(new FileWriter(new File("nupackTest/"+prefix+".con")));
		conOut.print(concs);
		conOut.close();

		/*
			runProcess(System.getProperty("NUPACKHOME")+"/bin/complexes -material dna -mfe "+prefix,
					new String[]{
				"NUPACKHOME="+System.getProperty("NUPACKHOME")},
			nupackDir);
		 */

		runProcess(System.getProperty("NUPACKHOME")+"/bin/complexes -material dna -pairs "+prefix,
				new String[]{
			"NUPACKHOME="+System.getProperty("NUPACKHOME")},
			nupackDir);
		
		if (concentrations){
			runProcess(System.getProperty("NUPACKHOME")+"/bin/concentrations -sort 0 "+prefix,
					new String[]{
				"NUPACKHOME="+System.getProperty("NUPACKHOME")},
				nupackDir);
		}
	}
	
	public static void runProcess(String string, String[] env,	File nupackDir) throws IOException {
		System.err.println(">"+string);

		Process p = Runtime.getRuntime().exec(string,env,nupackDir);
		if(false){
			final Scanner in2 = new Scanner(p.getInputStream());
			new Thread(){
				public void run(){
					try {
						while(in2.hasNextLine()){
							System.err.println(in2.nextLine());
						}
					} catch (Throwable e){

					}
				}
			}.start();
			try {
				p.waitFor();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			in2.close();
		} else {
			try {
				p.waitFor();
				p.getOutputStream().close();
				p.getInputStream().close();
			    p.getErrorStream().close();
			} catch (Throwable e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			p.destroy();
		}
}

}
