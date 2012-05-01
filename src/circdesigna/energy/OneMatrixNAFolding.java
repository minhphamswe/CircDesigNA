package circdesigna.energy;

import java.awt.Point;
import java.util.ArrayList;

/*
 * Useful for debugging, a NA Folder with a single matrix / traceback that can be easily visualized
 */
public interface OneMatrixNAFolding extends NAFolding {

	ArrayList<Point> getTraceback();
	double[][] getScoreMatrix(int len1, int len2);

}
