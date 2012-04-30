package circdesigna;

import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import circdesigna.energy.ExperimentalDuplexParams;

public class ZipExtractor {
	public static ZipInputStream getFile(String name){
		ZipInputStream paramZip = null;
		try {
			paramZip = new ZipInputStream(ExperimentalDuplexParams.class.getResourceAsStream("/"+name));
			//System.out.println("Done (1)");
		} catch (Throwable e){
			//Try loading it as a file.
			try {
				paramZip = new ZipInputStream(new FileInputStream("parameters.zip"));
				//System.out.println("Done (2)");
				//System.out.println("Loaded parameters file from disk.");
			} catch (Throwable f){
				throw new RuntimeException("Could not load the "+name+" file. Please include this file in the working directory!");
			}
		}
		return paramZip;
	}

	public static ByteArrayOutputStream readFully(ZipInputStream paramZip) {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		byte[] buf = new byte[1024];
		int read = -1;
		try {
			while((read=paramZip.read(buf))>0){
				baos.write(buf, 0, read);
			}
		} catch (IOException e) {
			throw new RuntimeException("Could not extract entry."+ e.getMessage());
		}
		return baos;
	}
}
