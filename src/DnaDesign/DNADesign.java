package DnaDesign;

import java.awt.Dimension;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.JFrame;

import processing.core.PApplet;
import processing.xml.XMLElement;
import BulletGame$1.BulletGame$1Engine$L2$1$BulletGameMain;
import BulletGame$1.BulletGame$1Engine$L3$1$BulletGameEditorScreen;
import BulletGame$1.BulletGame$1Engine$L4$1$BulletGameBulletSimulation;
import Deployments.LocalOnly.TaiGameLaunch;
import TaiGameCore.ProceGLHybrid;

public class DNADesign extends TaiGameLaunch{
	public DNADesign(){
	}
	/**
	 * Launch-time details, handled here:
	 * 1) Default config
	 * 2) Load config
	 */
	public void launch(){
		File config = new File("config.ini");
		if (!config.exists()){
			defaultConfig(config);
		}
		PApplet hollowBase = new BulletHell_Lite(){
			int targetFrameRate;
			{
				XMLElement configX = new XMLElement(this, "config.ini");
				int width = new Integer(configX.getChild("widthWindowed").getContent());
				int height = new Integer(configX.getChild("heightWindowed").getContent());
				targetFrameRate =  new Integer(configX.getChild("frameRate").getContent());
				setPreferredSize(new Dimension(width,height));
			}
			/*
			public void setup(){
				size(width,height,OPENGL);
				frameRate(targetFrameRate);
			}
			public void draw(){};
			public void printStackTrace(Exception e){
				death(e);
			}
			*/
		};
		JFrame holder = new JFrame(getGameName());
		setMainGUIElement(holder);
		holder.add(hollowBase);
		holder.setVisible(true);
		holder.pack();
		hollowBase.init();
		holder.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		ProceGLHybrid gameActual = new BulletGame$1Engine$GROUND(holder, hollowBase);
	}
	private static void defaultConfig(File target){
		PrintWriter out;
		try {
			out = new PrintWriter(new FileWriter(target));
			out.println("<paramlist stage=\"0\">");
			out.println("<widthWindowed>640</widthWindowed>");
			out.println("<heightWindowed>480</heightWindowed>");
			out.println("<fullScreen>false</fullScreen>");
			out.println("<widthFullScreen>640</widthFullScreen>");
			out.println("<heightFullScreen>480</heightFullScreen>");
			out.println("<frameRate>60</frameRate>");
			out.println("</paramlist>");
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public String getGameName() {
		return "Bullet Hell $1";
	}


	/**
	 * THE "GROUND" LAYER, ENDING THE SERIES OF L1.. L2.. L3 and so on. MUST BE UPDATED TO INCLUDE ALL FIELDS ABOVE IT!
	 */
	public static class BulletGame$1Engine$GROUND extends BulletGame$1Engine$L4$1$BulletGameBulletSimulation{
		public BulletGame$1Engine$GROUND(JFrame holder, PApplet hold) {
			super(holder, hold);
		}
		public static final int MAINSCREEN = 0,GAMESCREEN=1,EDITORSCREEN=2,GAMESCREEN_FROMEDITOR=3,BENCHMARK_EDITOR=4,GAMESCREEN_FROMBENCHMARKS=5;
		public BulletGameScreen SCREEN(int num) {
			switch(num){
			case MAINSCREEN: 
				//This specifies the first screen loaded:
				return new BulletGame$1Engine$L2$1$BulletGameMain.Stage2Screen();	//break; 
			case GAMESCREEN: return new BulletGame$1Engine$L4$1$BulletGameBulletSimulation.GameScreen();		//break;
			case EDITORSCREEN: return new BulletGame$1Engine$L3$1$BulletGameEditorScreen.EditorScreen();	//break;
			case GAMESCREEN_FROMEDITOR: return new BulletGame$1Engine$L4$1$BulletGameBulletSimulation.GameScreenFromEditor();		//break;
			case BENCHMARK_EDITOR: return new BulletGame$1.BulletGame$1Engine$L3$2$BulletGameBenchmarks.BenchmarkScreen();		//break;
			case GAMESCREEN_FROMBENCHMARKS: return new BulletGame$1Engine$L4$1$BulletGameBulletSimulation.GameScreenFromBenchmarks();		//break;
			}
			throw new RuntimeException("UNASSIGNED SCREEN NUMBER"+num);
		}
	}
}
