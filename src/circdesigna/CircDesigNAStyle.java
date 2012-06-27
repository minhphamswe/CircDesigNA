package circdesigna;

import java.awt.Color;

public class CircDesigNAStyle {
	/**
	 * Colors to supplement java.awt.Color
	 */
	public static class Colors {
		public static Color TEAL = new Color(0,128,128);
		public static Color teal = new Color(0,128,128);
		public static Color NAVY = new Color(0,0,128);
		public static Color navy = new Color(0,0,128);
		public static Color AO = new Color(0,128,0);
		public static Color ao = new Color(0,128,0);
	}
	
	public Color color;

	public static CircDesigNAStyle getDefaultStyle() {
		CircDesigNAStyle toRet = new CircDesigNAStyle();
		toRet.color = Color.black;
		return toRet;
	}
}
