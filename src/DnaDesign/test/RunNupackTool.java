package DnaDesign.test;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;

public class RunNupackTool {
	public static int OPCODE_COMPLEXES = 0, OPCODE_MFE=OPCODE_COMPLEXES+1, OPCODE_COMPLEXES_AND_CONC = OPCODE_MFE+1;
	public static void runNupack(String seqs, String concs, int maximumComplexSize, String prefix, int opcode, File nupackDir) throws IOException{
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
		if (opcode==OPCODE_MFE){
			runProcess(System.getProperty("NUPACKHOME")+"/bin/mfe -multi -material dna "+prefix,
					new String[]{
				"NUPACKHOME="+System.getProperty("NUPACKHOME")},
				nupackDir);
		}

		if (opcode==OPCODE_COMPLEXES_AND_CONC || opcode == OPCODE_COMPLEXES){
			runProcess(System.getProperty("NUPACKHOME")+"/bin/complexes -material dna -pairs "+prefix,
					new String[]{
				"NUPACKHOME="+System.getProperty("NUPACKHOME")},
				nupackDir);
		}
		if (opcode==OPCODE_COMPLEXES_AND_CONC){
			runProcess(System.getProperty("NUPACKHOME")+"/bin/concentrations -sort 0 "+prefix,
					new String[]{
				"NUPACKHOME="+System.getProperty("NUPACKHOME")},
				nupackDir);
		}
	}
	
	public static void runProcess(String string, String[] env,	File nupackDir) throws IOException {
		System.err.println(">"+string);

		final Process p = Runtime.getRuntime().exec(string,env,nupackDir);
		if(true){
			new Thread(){
				public void run(){
					InputStream is = new BufferedInputStream(p.getInputStream());
					int ch = -1;
					try {
						while((ch=is.read())!=-1){
							System.err.print((char)ch);
						}
					} catch (IOException e) {
						//e.printStackTrace();
					}
				}
			}.start();
		}
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
