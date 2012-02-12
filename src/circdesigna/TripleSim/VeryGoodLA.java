package circdesigna.TripleSim;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * A set of nice linear algebra operations, including fast multiply, block QR and block LU and system solving.
 * @author Benjamin
 */
public class VeryGoodLA {
	public static void main(String[] args){
		int[] nTestArr;
		/*
		int[] nTestArr = new int[128];
		for(int k = 0; k < nTestArr.length; k++){
			nTestArr[nTestArr.length-1-k] = k;
		}
		*/
		long nano;
		nTestArr = new int[]{5};
		for(int N : nTestArr){
			Matrix A = new Matrix(N+1,N);
			for(int y = 0; y < A.m; y++){
				for(int x = 0; x < A.n; x++){
					//REALLY BAD example.
					//A.set(y, x, y==x?1 : 0);
					A.set(y,x,Math.random());
				}
			}

			A.print();
			testQR(A);
			//testLU(A);

			/*
			Matrix C = Matrix.zeros(A);
			nano = System.nanoTime();
			MatrixOp.GEMM(1,A,A,C);
			double timeSmartGemm = (System.nanoTime()-nano)/1e9;
			nano = System.nanoTime();
			MatrixOp.GEPB(1,A,A,C);
			double timeDumbGemm = (System.nanoTime()-nano)/1e9;
			System.out.println(N +" "+timeSmartGemm+" "+timeDumbGemm);
			*/
			
			//A.print();
			//C.print();
		}
	}
	
	private static void testLU(Matrix A) {
		long nano = System.nanoTime();
		LUFactorization LU = LUFactorization.LUFactorize(A,true);
		System.out.println((System.nanoTime()-nano)/1e9);

		LU.print();
		//Test vector against equivalent multiplication
		Matrix p = new Matrix(A.n,1);
		for(int y = 0; y < p.m; y++){
			p.set(y,0,Math.random());
		}
		Matrix prodA = Matrix.zeros(A.partition(0,0,A.m,1));
		Matrix prodLU = Matrix.zeros(A.partition(0,0,A.m,1));
		Matrix prodLU2 = Matrix.zeros(A.partition(0,0,A.m,1));
		MatrixOp.GEMV(1, A, p, prodA);
		LU.multU(p,prodLU);
		LU.multL(prodLU,prodLU2);
		LU.multPinverse(prodLU2);
		//prodA.print();
		//prodLU.print();
		//prodLU2.print();
		{
			Matrix error = Matrix.zeros(prodA);
			MatrixOp.subtract(prodA, prodLU2, error);
			System.out.println("Forward error from LU decomposition: "+MatrixOp.twoNorm(error));
		}
		//Now solve the system for LUx=prodA
		Matrix x = Matrix.zeros(prodA);
		LU.solve(prodA,x);
		{
			Matrix error = Matrix.zeros(prodA);
			MatrixOp.subtract(x, p, error);
			System.out.println("Forward error from solving system: "+MatrixOp.twoNorm(error));
		}
	}

	private static void testQR(Matrix A) {
		LinearLeastSquares lls = new LinearLeastSquares(A);
		
		//Test solution:
		Matrix b = new Matrix(A.m,1);
		for(int y = 0; y < b.m; y++){
			b.set(y,0,Math.random());
		}
		b.print();
		
		Matrix cooeficients = new Matrix(A.n, 1);
		lls.project(b, cooeficients);
		
		Matrix testx2 = Matrix.zeros(b);
		//MatrixOp.GEMV(1, QR.R.partition(0, 0, QR.R.m, QR.R.n, true), testx, testx2);
		//testx2.print();
		MatrixOp.GEMV(1, A, cooeficients, testx2);
		//testx2.print();
		MatrixOp.subtract(b,testx2,testx2);
		//testx.print();
		System.out.println(MatrixOp.twoNorm(testx2));
	}
	public static class LinearLeastSquares {
		/**
		 * Given a matrix of n m-dimensional vectors, pre-compute a projection system
		 * that finds the cooeficients a_i such that the linear combination of the n
		 * vectors, with cooeficients a_i, comes as close as possible (in two norm) to
		 * a given vector b.
		 */
		public LinearLeastSquares(Matrix LIVectors){
			this.A = LIVectors;
			Matrix AtA = new Matrix(A.n, A.n);
			Matrix At = A.partition(0,0,A.m,A.n,true);
			MatrixOp.GEMM(1, At, A, AtA);
			QR = QRFactorization.HouseholderQR(AtA);
		}
		private Matrix A;
		private QRFactorization QR;
		public void project(Matrix b, Matrix cooeficients){
			Matrix At = A.partition(0,0,A.m,A.n,true);
			MatrixOp.GEMV(1, At, b, cooeficients);
			TriangularMatrix.forwardSub(QR.R.partition(0, 0, QR.R.m, QR.R.n, true), cooeficients);
			QR.multQ(cooeficients, cooeficients);
		}
	}

	public static class TriangularMatrix {
		/**
		 * Solves the system Lx = b, where b becomes x after method completes ("in place solve")
		 */
		public static void forwardSub(Matrix L, Matrix b) {
			if (L.isUnitDiagonal()){
				for(int i = 0; i < L.m-1; i++){
					for(int j = 0; j < b.n; j++){
						Matrix b2 = b.partition(i+1, j, b.m-(i+1), 1);
						Matrix l21 = L.partition(i+1, i, L.m-(i+1), 1);
						MatrixOp.add(-b.get(i,j),l21,b2,b2);
					}
				}
			} else {
				if (b.n == 1){
					for(int i = 0; i < L.m; i++){
						b.set(i, 0, b.get(i,0) / L.get(i,i));
						if (i + 1 < L.m){
							for(int j = 0; j < b.n; j++){
								Matrix b2 = b.partition(i+1, j, b.m-(i+1), 1);
								Matrix l21 = L.partition(i+1, i, L.m-(i+1), 1);
								MatrixOp.add(-b.get(i,j),l21,b2,b2);
							}
						}
					}
				} else {
					throw new RuntimeException("Derp");
				}
			}
		}

		/**
		 * Solves the system Ux = b, where b becomes x after method completes ("in place solve")
		 */
		public static void backSub(Matrix U, Matrix b) {
			if (U.isUnitDiagonal()){
				throw new RuntimeException("Derp");
			} else {
				if (b.n > 1){
					throw new RuntimeException("Derp");
				}
				for(int i = U.m-1; i >= 0; i--){
					double u12dotb2 = 0;
					if (i+1 < U.n){
						Matrix u12t = U.partition(i, i+1, 1, U.n-(i+1), true);
						Matrix b2 = b.partition(i+1,0,b.m-(i+1),1);
						u12dotb2 = DotProduct.dot(u12t, b2);
					}
					b.set(i,0,(b.get(i,0) - u12dotb2) / U.get(i,i));
				}
			}
		}
		
	}
	
	public static class LUFactorization{
		public static LUFactorization LUFactorize(Matrix A, boolean pivoting) {
			LUFactorization toRet = new LUFactorization();
			/*
			if (A.m != A.n){
				throw new RuntimeException("LU Factorization only valid for matrixes where m == n");
			}
			*/
			Matrix LU = new Matrix(A);
			
			Matrix permute = new Matrix(A.m,1);
			toRet.permute = permute;
			
			//LU factorization
			int k = 0;
			{
				Matrix A22 = LU.partition(0, 0, LU.m, LU.n);
				while(A22.m>1){
					if (pivoting){
						int pi = PermutationMatrix.maxi(A22.partition(0,0,A22.m,1));
						permute.set(0, 0, pi);
						permute = permute.partition(1, 0, permute.m-1,1);
						//Swaps 
						PermutationMatrix.swapRows(LU, k, pi+k);
					}
					double a11 = A22.get(0,0);
					Matrix a21 = A22.partition(1, 0, A22.m-1, 1);
					MatrixOp.scale(a21, 1./a11, a21);
					Matrix nA22 = A22.partition(1,1,A22.m-1,A22.n-1);
					/*
					for(int x = 0; x < nA22.n; x++){
						//Subtract off multiples of a21
						MatrixOp.GEMV(-1, a21, A22.partition(0,1+x,1,1), nA22.partition(0, x, nA22.m, 1));
					}
					*/
					MatrixOp.GEPP(-1, a21, A22.partition(0,1,1,A22.n-1), nA22);
					A22 = nA22;
					k++;
				}
			}
			
			Matrix A_L = LU.partition(0, 0, LU.m, LU.n);
			Matrix A_U = LU.partition(0, 0, LU.m, LU.n);
			A_L.setUnitDiagonal(true);
			A_U.setUpperTriangular(true);
			A_L.setLowerTriangular(true);
			toRet.A_L = A_L;
			toRet.A_U = A_U;
			
			//LU.print();
			
			return toRet;
		}
		public void solve(Matrix prodA, Matrix x) {
			//Copy prodA into x, forward sub is unit lower triangular.
			MatrixOp.add(0, x, prodA, x);
			multP(x);
			TriangularMatrix.forwardSub(A_L, x);
			TriangularMatrix.backSub(A_U, x);
		}
		public void print() {
			System.out.println("L:");
			A_L.print();
			System.out.println("U:");
			A_U.print();
		}
		private Matrix A_L;
		private Matrix A_U;
		private Matrix permute;

		public void multU(Matrix v, Matrix y) {
			MatrixOp.GEMV(1, A_U, v, y);
		}
		public void multL(Matrix v, Matrix y) {
			MatrixOp.GEMV(1, A_L, v, y);
		}
		public Matrix getL() {
			return A_L;
		}
		public Matrix getU() {
			return A_U;
		}
		public void multP(Matrix a){
			for(int k = 0; k < permute.m; k++){
				int pi1 = k;
				int pi2 = (int) permute.get(k, 0)+k;
				PermutationMatrix.swapRows(a, pi1, pi2);
			}
		}
		public void multPinverse(Matrix a){
			for(int k = permute.m-1; k >= 0; k--){
				int pi1 = k;
				int pi2 = (int) permute.get(k, 0)+k;
				PermutationMatrix.swapRows(a, pi1, pi2);
			}
		}
	}

	private static class QRFactorization{
		private Householder[] transforms;
		private Matrix R;
		public static QRFactorization HouseholderQR(Matrix A) {
			QRFactorization qr = new QRFactorization();
			Matrix R = new Matrix(A);
			Householder[] transforms = new Householder[R.n];
			for(int i = 0; i < R.n; i++){
				Matrix colR = R.partition(0,i,R.m,1);
				for(int j = 0; j < i; j++){
					transforms[transforms.length-1-j].evaluate(colR, colR);
				}
				colR = colR.partition(i,0,colR.m-i,1);
				transforms[transforms.length-1-i] = new Householder(colR,1);
				colR = R.partition(0,i,R.m,1);
				transforms[transforms.length-1-i].evaluate(colR, colR);
			}
			//R.print();
			R.setUpperTriangular(true);
			qr.R = R;
			qr.transforms = transforms;
			//Postcondition: m is R, transforms multiplies to form Q.
			return qr;
		}
		/**
		 * Result testx2 += R*testx
		 */
		public void multR(Matrix testx, Matrix testx2) {
			MatrixOp.GEMV(1,R,testx,testx2);
		}
		public void multQ(Matrix testx, Matrix testx2) {
			for(int i = 0; i < R.n; i++){
				//Go in opposite order as taken when creating R
				transforms[i].evaluate(testx,testx2);
				testx = testx2;
			}
		}
	}
	private static class PermutationMatrix {
		public static int maxi(Matrix a) {
			if (a.n!=1){
				throw new RuntimeException("maxi only works on column vectors");
			}
			int bestI = 0;
			double best = a.get(bestI, 0);
			for(int k = 1; k < a.m; k++){
				double test = a.get(k,0);
				if (test > best){
					bestI = k;
					best = test;
				}
			}
			return bestI;
		}
		
		public static void swapRows(Matrix a, int pi1, int pi2){
			if (pi1==pi2){
				return;
			}
			for(int k = 0; k < a.n; k++ ){
				double tmp = a.get(pi1, k);
				a.set(pi1,k,a.get(pi2,k));
				a.set(pi2,k,tmp);
			}
		}
	}
	public static class MatrixOp {
		public static void subtract(Matrix testx, Matrix testx2, Matrix testx3) {
			add(-1,testx2,testx,testx3);
		}
		
		public static void GEMM(double alpha, Matrix A, Matrix B, Matrix C) {
			//Using fig. 8 from Goto paper.
			if (A.n != B.m){
				throw new IllegalArgumentException("Nonconformant multiplication: n = "+A.n+", k = "+B.m);
			}
			int ib = 256; 
			for(int i = 0; i < A.n; i+=ib){
				int panelSize = Math.min(A.n-i,ib);
				Matrix AVPanel = A.partition(0, i, A.m, panelSize);
				Matrix BHPanel = B.partition(i, 0, panelSize, B.n);
				GEPP(alpha,AVPanel,BHPanel,C);
			}
		}
		private static void GEPP(double alpha, Matrix A, Matrix B, Matrix C) {
			int ib = 16;
			//Pack A
			//A = new Matrix(A);
			Matrix BBlkPacked = null;
			for(int i = 0; i < B.n; i+=ib){
				int blockSize = Math.min(B.n-i,ib);
				Matrix BBlk = B.partition(0, i, B.m, blockSize);
				Matrix CBlk = C.partition(0, i, C.m, blockSize);
				//Pack B
				if (BBlkPacked==null || BBlkPacked.n != BBlk.n){
					BBlkPacked = new Matrix(BBlk);
				} else {
					for(int y = 0; y < BBlk.m; y++){
						for(int x = 0; x < BBlk.n; x++){
							BBlkPacked.set(y,x,BBlk.get(y,x));
						}
					}
				}
				GEPB(alpha,A,BBlkPacked,CBlk);
			}
		}	
		private static void GEPB(double alpha, Matrix A, Matrix B, Matrix C) {
			for(int i = 0; i < A.m; i++){
				for(int j = 0; j < B.n; j++){
					//KERNEL.
					double dot = 0;
					for(int p = 0; p < B.m; p++){
						dot += A.get(i, p)*B.get(p, j);
					}
					C.set(i, j, C.get(i,j) + dot*alpha);
				}
			}
		}

		//Calculates y = y + alpha * A v
		//y cannot be equal to v, because of the method by which y is written to during calculation.
		public static void GEMV(int alpha, Matrix A, Matrix v, Matrix y) {
			//todo check v == y
			if (v.n!=1){
				throw new IllegalArgumentException("Matrix-Vector Multiply: Second argument was not a vector");
			}
			if (v.m!=A.n || A.m>y.m){
				throw new IllegalArgumentException("Matrix-Vector Multiply: Incompatible sizes");
			}
			if (true){
				for(int i = 0; i < A.m; i++){
					Matrix rA;
					Matrix vA;
					boolean isTriangular = true;
					if (A.isUpperTriangular()){
						if (i >= A.n){
							break;
						}
						rA = A.partition(i, i+1, 1,A.n-(i+1), true);
						vA = v.partition(i+1,0,v.m-(i+1),1);
					} else if (A.isLowerTriangular()){
						rA = A.partition(i,0,1,i,true);
						vA = v.partition(0,0,i,1);
					} else {
						isTriangular = false;
						rA = A.partition(i, 0, 1, A.n, true);
						vA = v;
					}
					//TODO: unit diagonal
					double diag = 0.0;
					if (isTriangular){
						double mult = A.get(i,i);
						if (A.isUnitDiagonal){
							mult = 1;
						}
						diag = mult * v.get(i, 0);
					}
					y.set(i, 0, y.get(i, 0)+alpha*(DotProduct.dot(rA,vA)+diag));
				}
			} else {
				for(int i = 0; i < A.n; i++){
					add(v.get(i, 0)*alpha,A.partition(0, i, A.m, 1),y,y);
				}
			}
		}
		public static void add(double alpha, Matrix testx, Matrix testx2, Matrix testx3) {
			for(int y = 0; y < testx.m; y++){
				for(int x = 0; x < testx.n; x++){
					testx3.set(y, x, alpha*testx.get(y, x)+testx2.get(y,x));
				}
			}
		}

		public static void scale(Matrix testx, double alpha,Matrix testx3) {
			for(int y = 0; y < testx.m; y++){
				for(int x = 0; x < testx.n; x++){
					testx3.set(y, x, testx.get(y, x)*alpha);
				}
			}
		}

		public static double twoNorm(Matrix testx) {
			return Math.sqrt(DotProduct.dot(testx, testx));
		}

		public static void TRSLV(Matrix A, Matrix b) {
			if (A.n != b.m){
				throw new IllegalArgumentException("Column dimension of A does not match Row dimension of b");
			}
			if (A.n != A.m || !(A.isLowerTriangular() || A.isUpperTriangular)){
				throw new IllegalArgumentException("TRSLV only works on square upper or lower triangular matrices");
			}
			if (A.isLowerTriangular()){
				for(int i = 0; i < A.n; i++){
					if (!A.isUnitDiagonal){
						b.set(i, 0, b.get(i,0) / A.get(i,i));
					}
					Matrix b2 = b.partition(i+1, 0, b.m-1-i, 1);
					if (b2.m > 0){
						MatrixOp.add(-b.get(i, 0), A.partition(i+1, i, A.m-1-i, 1), b2, b2);
					}
				}
			} else {
				for(int i = A.n-1; i >= 0; i--){
					if (!A.isUnitDiagonal){
						b.set(i, 0, b.get(i,0) / A.get(i,i));
					}
					Matrix b0 = b.partition(0, 0, i, 1);
					if (b0.m > 0){
						MatrixOp.add(-b.get(i, 0), A.partition(0, i, i, 1), b0, b0);
					}
				}
			}
		}
	}
	/**
	 * Can transform such that neither dimension is increased.
	 * In the case that out is too large, a suitable partition of out is returned.
	 */
	private interface InSpaceTransform{
		public Matrix evaluate(Matrix in, Matrix out);
	}
	private static class DotProduct{
		public static double dot(Matrix u, Matrix in) {
			if (u.m != in.m || u.n != 1 || in.n != 1){
				throw new RuntimeException("Invalid arguments: dot product operates on equal length vectors");
			}
			if (true){
				double result = 0;
				for(int i = 0; i < in.m; i++){
					result += u.get(i,0)*in.get(i,0);
				}
				return result;
			} else {
				return dotLog(u,in,0,u.m);
			}
		}
		private static double dotLog(Matrix u, Matrix in, int l, int r){
			int m = (l+r)/2;
			if (m==l){
				return u.get(l, 0)*in.get(l, 0);
			}
			return dotLog(u,in,l,m)+dotLog(u,in,m,r);
		}
	}
	private static class NormalVector {
		public static void normalize(Matrix u,Matrix out){
			if (u.n > 1){
				throw new RuntimeException("Normalize divides a vector (column matrix) by its two norm.");
			}
			double norm = 0;
			for(int k = 0; k < u.m; k++){
				norm += Math.pow(u.get(k, 0),2);
			}
			norm = Math.sqrt(norm);
			for(int k = 0; k < u.m; k++){
				out.set(k,0, u.get(k,0)/norm);
			}
		}
		public static void normalize(Matrix u) {
			normalize(u,u);
		}
		
	}
	private static class Householder implements InSpaceTransform{
		private Matrix u;
		private int vectorsToRight;
		public Householder(Matrix partition, int vectorsToTheRight) {
			if (partition.n > 1){
				throw new RuntimeException("Householder takes as input a column vector, and creates u to reduce it to a*e0");
			}
			u = new Matrix(partition.m,1);
			vectorsToRight = vectorsToTheRight;
			double twonormx = 0;
			int signx0 = 1;
			for(int k = 0; k < partition.m; k++){
				double val = partition.get(k, 0);
				twonormx += val*val;
				if (k==0){
					if (val == 0){
						throw new RuntimeException("Ohnoes!");
					}
					signx0 = val > 0?1:-1;
				}
				u.set(k, 0, -val);
			}
			twonormx = Math.sqrt(twonormx);
			//u = -x-sign(x0)*norm(x)
			u.set(0,0,u.get(0, 0)-signx0*twonormx);
			//Normalize u
			NormalVector.normalize(u);
		}
		public Matrix evaluate(Matrix in, Matrix out) {
			if (in.n>1){
				throw new RuntimeException("Householder evaluates on vectors.");
			}
			if (out.m < in.m){
				throw new RuntimeException("Householder: Output vector smaller than input.");
			}
			
			int partTo = in.m-u.m;
			
			Matrix outUpper = out.partition(0,0,partTo,1);
			Matrix inUpper = in.partition(0,0,partTo,1);
			Matrix outLower = out.partition(partTo,0,u.m,1);
			Matrix inLower = in.partition(partTo,0,u.m,1);
			
			MatrixOp.add(0, outUpper, inUpper, outUpper);
			
			//u.print();
			//in.print();
			double uHx = DotProduct.dot(u,inLower);
			//System.out.println(uHx);
			if (true){
				MatrixOp.add(-2*uHx,u,inLower,outLower);
			} else {
				for(int k = 0; k < u.m; k++){
					outLower.set(k, 0, inLower.get(k, 0) - 2 * u.get(k, 0) * uHx);
				}
			}
			//out.print();
			return out;
		}		
	}
	private static class Matrix_ {
		private final double[] contents;
		public Matrix_(double[] contents){
			this.contents = contents;
		}
		public Matrix_(int length){
			contents = new double[length];
		}
		public double get(int i){
			return contents[i];
		}
		public void set(int i, double val){
			contents[i] = val;
		}
	}
	public static class Matrix {
		private final Matrix_ real;
		/**
		 * Constructs a new linear transformation from F^n->F^m.
		 */
		public Matrix(int m, int n) {
			if (m < 0 || n < 0){
				throw new IllegalArgumentException("Cannot have negative matrix dimensions");
			}
			real = new Matrix_(m*n);
			fullStride = this.numCols = this.n = n;
			this.numRows = this.m = m;
			offx = 0; offy = 0;
			this.transpose = false;
		}
		/**
		 * Returns a matrix of zeros of the same size as p
		 */
		public static Matrix zeros(Matrix p) {
			return new Matrix(p.m,p.n);
		}
		/**
		 * Returns an identity matrix of the same size as p.
		 * Throws an exception if p is not square.
		 */
		public static Matrix eye(Matrix p) {
			if (p.m!=p.n){
				throw new IllegalArgumentException("Identity matrices are square. "+p.m+", "+p.n);
			}
			Matrix eye = new Matrix(p.m,p.m);
			for(int i = 0; i < p.m; i++){
				eye.set(i, i, 1);
			}
			return eye;
		}
		/**
		 * Returns a backed partition of this matrix from i,j to i+rows,j+cols, optionally transposed
		 * in the representation. No memory is modified or copied by this call. Partitions of a matrix
		 * have a different set of flags from the backing matrix, so a submatrix can be declared upperdiagonal
		 * while an outer one is not. 
		 * Like all matrices, a partition by default has no flags.
		 */
		public Matrix partition(int i, int j, int rows, int cols, boolean b) {
			if (rows < 0 || cols < 0){
				throw new IllegalArgumentException("Cannot have negative matrix dimensions");
			}
			//Ensure it fits, if it is nonempty.
			int offset = -1;
			if (rows > 0 && cols > 0){
				offset = ind(i,j);
				int maxOff = ind(i+rows-1,j+cols-1);
			}
			if (transpose){
				int tmp = cols;
				cols = rows;
				rows = tmp;
				b = !b;
			}
			return new Matrix(this,offset,rows,cols,b);
		}
		/**
		 * see partition (i,j,r,col,transpose)
		 */
		public Matrix partition(int i, int j, int rows, int cols) {
			return partition(i,j,rows,cols,false);
		}
		/**
		 * Formatted matrix print.
		 */
		public void print() {
			StringWriter oo = new StringWriter();
			PrintWriter out = new PrintWriter(oo);
			out.println("[");
			for(int y = 0; y < m; y++){
				for(int x = 0; x < n; x++){
					String toPrint = null;
					if (isUnitDiagonal){
						if (x==y){
							toPrint = "<1>";
						}
					}
					if (isLowerTriangular){
						if (y < x){
							toPrint = "<0>";
						}
					}
					if (isUpperTriangular){
						if (x < y){
							toPrint = "<0>";
						}
					}
					if (toPrint==null){
						toPrint = get(y,x)+"";
					}
					out.printf("%-16s"+"\t",toPrint+(x+1<n?",":""));
				}
				//if (y+1 < numrows){
					out.println();
				//}
			}
			out.println("]");
			out.close();
			System.out.print(oo.toString());
		}
		/**
		 * Used by partition.
		 */
		private Matrix(Matrix real, int offset, int m, int n, boolean transpose){
			this.real = real.real;
			fullStride = real.fullStride;
			this.transpose = transpose;
			this.m = transpose?n:m;
			this.n = transpose?m:n;
			this.numRows = m;
			this.numCols = n;
			setup(offset);
		}
		private void setup(int offset) {
			this.packed = false;
			if (fullStride!=0){
				//Empty partitions have undefined offx / offy.
				offx = offset%fullStride;
				offy = offset/fullStride;
			}
		}
		public Matrix(double[] vector){
			this.real = new Matrix_(vector);
			numRows = m = vector.length;
			numCols = n = fullStride = 1;
			transpose = false;
			setup(0);
		}
		/**
		 * Deep copy constructor. The new matrix is NOT backed by the old one.
		 */
		public Matrix(Matrix m2) {
			this(m2.m,m2.n);
			for(int y = 0; y < m; y++){
				for(int x = 0; x < n; x++){
					set(y,x,m2.get(y,x));
				}
			}
		}
		private boolean isUpperTriangular;
		private boolean isLowerTriangular;
		private boolean isUnitDiagonal;
		public boolean isLowerTriangular() {
			return isLowerTriangular;
		}
		public boolean isUnitDiagonal(){
			return isUnitDiagonal;
		}
		public void setUnitDiagonal(boolean isUnitDiagonal){
			this.isUnitDiagonal = isUnitDiagonal;
		}
		public void setLowerTriangular(boolean isLowerTriangular) {
			this.isLowerTriangular = isLowerTriangular;
		}
		public boolean isUpperTriangular() {
			return isUpperTriangular;
		}
		public void setUpperTriangular(boolean isUpperTriangular) {
			this.isUpperTriangular = isUpperTriangular;
		}
		private int offx, offy;
		private int fullStride;
		public final int n, m;
		private boolean packed = true;
		private final int numRows, numCols;
		private boolean transpose;
		public double get(int row, int col){
			return real.get(ind(row,col));
		}
		public void set(int row, int col, double val){
			real.set(ind(row,col),val);
		}
		public void add(int row, int col, double val){
			set(row,col,val + get(row,col));
		}
		/*
		public void shiftPartition(int drow, int dcol){
			offy += drow;
			offx += dcol;
		}
		*/
		private int ind(int row, int col) {
			if (transpose){
				int tmp = row;
				row = col;
				col = tmp;
			}
			if (col < 0 || col >= numCols){
				throw new ArrayIndexOutOfBoundsException(col);
			}
			if (row < 0 || row >= numRows){
				throw new ArrayIndexOutOfBoundsException(row);
			}
			row += offy;
			col += offx;
			return row*fullStride + col;
		}
	}
}
