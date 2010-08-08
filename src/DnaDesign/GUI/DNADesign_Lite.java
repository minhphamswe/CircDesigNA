package DnaDesign.GUI;

import javax.swing.JFrame;

import processing.core.PApplet;
import BulletGame$1.BulletGame$1Engine$L1$1$OpenglTextRenderer;
import Deployments.JogampPapplet;
import TaiGameCore.ProceGLHybrid;

/**
 * This is a Lite PApplet, which has no config file and just behaves as one would expect
 * @author Benjamin
 */
public class DNADesign_Lite extends JogampPapplet{
	public void setup(){
		/*
		size(new Integer(getParameter("width")),
				new Integer(getParameter("height")),OPENGL);
		background(0);
		*/
		frameRate(60);
		//gameActual = makeActual();
	}
	public ProceGLHybrid makeActual	(){
		//Frame is not available, but could be made available.
		return new DnaDesign$GROUND(null, this);
	}
	private ProceGLHybrid gameActual;
	public void draw(){
		//The world is your oyster, or you could refer to the page-system above.
		g.background(255);
	};
	/**
	 * THE "GROUND" LAYER, ENDING THE SERIES OF L1.. L2.. L3 and so on. MUST BE UPDATED TO INCLUDE ALL FIELDS ABOVE IT!
	 */
	public static class DnaDesign$GROUND extends BulletGame$1Engine$L1$1$OpenglTextRenderer{//TODO: extend bottom layer
		public DnaDesign$GROUND(JFrame holder, PApplet hold) {
			super(holder, hold);
		}
		public static final int MAINSCREEN = 0;
		public BulletGameScreen SCREEN(int num) {
			switch(num){
			case MAINSCREEN: 
				//This specifies the first screen loaded:
				return null;	//break; 
			}
			throw new RuntimeException("UNASSIGNED SCREEN NUMBER"+num);
		}
	}
}
