package DnaDesign;

import java.awt.Color;

import Deployments.JogampPapplet;
import TaiGameCore.ProceGLHybrid;

/**
 * This is a Lite PApplet, which has no config file and just behaves as one would expect
 * @author Benjamin
 */
public class DNADesign_Lite extends JogampPapplet{
	public void setup(){
		size(new Integer(getParameter("width")),
				new Integer(getParameter("height")),OPENGL);
		frameRate(60);
		background(0);
		
		gameActual = makeActual();
	}
	public ProceGLHybrid makeActual(){
		return new DnaDesign.DNADesign(null, DNADesign_Lite.this);
	}
	private ProceGLHybrid gameActual;
	int loadState = 0;
	public void draw(){
		if (loadState==0){
			g.fill(0);
			rect(0,0,100,100);
		}
	};
	public void init(){
		setBackground(new Color(7,9,43));
		super.init();
		setBackground(new Color(125,125,255));
	}
}
