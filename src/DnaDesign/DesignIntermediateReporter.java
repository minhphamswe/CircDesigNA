package DnaDesign;

import java.util.ArrayList;
import java.util.Arrays;

public class DesignIntermediateReporter {
	public DesignIntermediateReporter(){
		molecules = new ArrayList();
	}
	public class DesignIntermediateScore{
		private ArrayList<int[]> scoreAdders;
		public DesignIntermediateScore(){
			this.scoreAdders = new ArrayList();
		}
		private void registerScore(String molA, String molB){
			int[] newRow = new int[]{indexOf(molA),indexOf(molB)};
			for(int[] row : scoreAdders){
				if (Arrays.equals(row, newRow)){
					return;
				}
			}
			scoreAdders.add(newRow);
		}
		public void addScore(double score){
			for(int[] row : scoreAdders){
				currentMatrix[row[0]][row[1]]+=score;
				currentMatrix[row[1]][row[0]]+=score;
			}
		}
	}
	private int indexOf(String mol){
		int index = molecules.indexOf(mol); 
		if (index==-1){
			molecules.add(mol);
			return molecules.size()-1;
		}
		return index;
	}
	private ArrayList<String> molecules;
	public float[][] currentMatrix;
	public boolean currentMatrixSynchronized = true;
	public void beginScoreReport(){
		currentMatrixSynchronized = false;
		int numMolecules = molecules.size();
		if (currentMatrix==null || currentMatrix.length != numMolecules){
			currentMatrix = new float[numMolecules][numMolecules];
		}
		for(float[] row : currentMatrix){
			Arrays.fill(row,0);
		}
	}
	public void endScoreReport(){
		currentMatrixSynchronized = true;
	}
	/**
	 * All scores contain at most 2 distinct molecules.
	 */
	public DesignIntermediateScore chooseDesignIntermediateScore(String moleculeA, String moleculeB){
		DesignIntermediateScore dis = new DesignIntermediateScore();
		for(String molA : moleculeA.split(",")){
			for(String molB : moleculeB.split(",")){
				dis.registerScore(molA, molB);
			}
		}
		return dis;
	}
	public String getMolecule(int k) {
		return molecules.get(k);
	}
}
