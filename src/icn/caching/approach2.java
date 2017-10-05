package icn.caching;

/* --------------------------------------------------------------------------
 * Approach2
 * tuan-minh.pham@lip6.fr
 *
 */

import java.text.DecimalFormat;
import java.util.Random;

import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Collections;
import java.util.ArrayList;
import java.util.Properties;
import java.util.Scanner;
import java.util.regex.Pattern;

public class approach2 {
	/*
	 * // public static final String INPUT_FILE = "in.txt"; public static final
	 * String PROBABILITY_FILE = "ProbabilityFile"; public static final String
	 * GENERATED_PROBABILITY_FILE = "prob_generated.txt"; public static final
	 * String HEURISTIC_SOLUTION_FILE = "heuristic_solution.txt"; public static
	 * final String HEURISTIC_OUPUT_FILE = "heuristic_output.txt";
	 * 
	 * public static final String OUTPUT_FILE = "OutputFile"; public static
	 * final String PARAM_NUMBER_NODES = "NumberOfNodes"; public static final
	 * String PARAM_NUMBER_CONTENT = "NumberOfContent"; public static final
	 * String PARAM_UPLINK_CAPACITY = "UplinkCapacity"; public static final
	 * String PARAM_UPLINK_CAPACITY_ARR = "UplinkCapacityArr"; public static
	 * final String PARAM_COST_ROOT_DSLAM = "CostRootDslam"; public static final
	 * String PARAM_COST_DSLAM_CONSUMER = "CostDslamConsumer"; public static
	 * final String PARAM_STORAGE_CAPACITY = "StorageCapacity"; public static
	 * final String PARAM_STORAGE_CAPACITY_I = "StorageCapacityI"; public static
	 * final String PARAM_STORAGE_CAPACITY_I2 = "StorageCapacityI2"; public
	 * static final String PARAM_STORAGE_CAPACITY_H = "StorageCapacityH"; public
	 * static final String PARAM_ZIPF_EXP = "ZipfExponent"; public static final
	 * String PARAM_ZIPF_EXP_ARR = "ZipfExponentArr"; public static final String
	 * PARAM_OPT_AVERAGE = "Optimization_Average"; public static final String
	 * PARAM_OPT_STEPS = "Optimization_Steps";
	 * 
	 * public static final String PARAM_TRACE_SOLUTION = "Trace_Solution";
	 * 
	 * // 1: use the model of transporation problem to compute the optimal
	 * routing // cost // 0: use heristics to compute the routing cost public
	 * static final String PARAM_RCOST_HEURISTIC = "RCost_Heuristic"; public
	 * static final String PARAM_EXHAUSIVE_SEARCH = "Exhausive_Search";
	 * 
	 * // CONSTANTS public static final int EXHAUSIVE_SEARCH_YES = 1; public
	 * static final int RCOST_HEURISTIC = 1; public static final int
	 * RCOST_OPTIMAL = 0; public static final int TRACE_SOLUTION_YES = 1; public
	 * static final int TRACE_SOLUTION_NO = 0;
	 * 
	 * static int traceSolution = TRACE_SOLUTION_NO; static int rcostHeuristic =
	 * RCOST_OPTIMAL; static int exhausiveSearch = EXHAUSIVE_SEARCH_YES;
	 * 
	 * // public static Properties props; public static Settings setting;
	 * 
	 * static int n; // index of set-top boxes I = {0,1,...,n-1} static int m;
	 * // index of content objects J = {0,1,...,m-1} static int c0; // uplink
	 * capacity static int w0; // traffic cost between the DSLAM h and a set-top
	 * box static int w1; // traffic cost between the root r and the DSLAM h
	 * 
	 * // sto[i] = number of content items that a local box can store if i <= n
	 * // sto[n+1] = number of content items that the DSLAM can store static
	 * int[] sto; static int stoI; static int stoH;
	 * 
	 * static int h; // index of the DSLAM static int r; // index of the root
	 * static BufferedWriter out; static double dZipfExp = 0.8; static int
	 * opt_avg; // index of the root static int opt_steps; // index of the root
	 * 
	 * // For an exhausive search only static int[][][] ed; static int[][] ey;
	 * // static int[][] earr; static int[][] eyopt; static double ecopt =
	 * Double.MAX_VALUE; static double ecountCom = 0; static double lowerBound =
	 * Double.MAX_VALUE; static double lowerBound_continuous = Double.MAX_VALUE;
	 */
	public static void main(String[] args) throws IOException {

		/*
		 * Uncomment to test revdoor
		 */
		/*
		 * try { test_revdoor(); } catch (IOException e2) {}
		 */
		/*
		 * String inputFilename = args[0]; if (new File(inputFilename).exists())
		 * { Settings.init(inputFilename); // props = new Properties(); //
		 * props.load(new FileInputStream(inputFilename)); }
		 * 
		 * setting = new Settings();
		 * 
		 * n = setting.getInt(PARAM_NUMBER_NODES); m =
		 * setting.getInt(PARAM_NUMBER_CONTENT); c0 =
		 * setting.getInt(PARAM_UPLINK_CAPACITY); w0 =
		 * setting.getInt(PARAM_COST_DSLAM_CONSUMER); w1 =
		 * setting.getInt(PARAM_COST_ROOT_DSLAM);
		 * 
		 * dZipfExp = setting.getDouble(PARAM_ZIPF_EXP); opt_avg =
		 * setting.getInt(PARAM_OPT_AVERAGE); opt_steps =
		 * setting.getInt(PARAM_OPT_STEPS); rcostHeuristic =
		 * setting.getInt(PARAM_RCOST_HEURISTIC); exhausiveSearch =
		 * setting.getInt(PARAM_EXHAUSIVE_SEARCH);
		 * 
		 * traceSolution = setting.getInt(PARAM_TRACE_SOLUTION); // Read
		 * individual storage capacity // sto =
		 * getCsvInt(PARAM_STORAGE_CAPACITY); // Or, all local boxes have the
		 * same storage capacity stoI =
		 * setting.getInt(PARAM_STORAGE_CAPACITY_I); stoH =
		 * setting.getInt(PARAM_STORAGE_CAPACITY_H); sto = new int[n + 1]; for
		 * (int i = 0; i < n; i++) { sto[i] = stoI; } sto[n] = stoH;
		 * 
		 * h = n; r = n + 1; String outFilename =
		 * setting.getSetting(OUTPUT_FILE); try { out = new BufferedWriter(new
		 * FileWriter(outFilename)); } catch (IOException e) { }
		 */
		//int maxE = 10, maxV = 10, maxF = 10, maxP = 1;
		int maxP = 1;

		// INPUT
		// VNF
		int nF; // number of available functions
		int[] rf;// = new int[maxF-1]; //required computing resource for
					// processing one unit of flow rate

		// k(vsti) = 1: the ith VNF of demand (s,t) is deployed at node v
		int[][][][] k;// = new int[maxV-1][maxV-1][maxV-1][maxF-1];//
                
		// G(V,E)
		int nV, nE;
		int[] ii;// = new int[maxE - 1]; // starting node of link e
		int[] jj;// = new int[maxE - 1]; // terminating node of link e
		double[] ce;// = new int[maxE - 1]; // link capacity
		int[] cv;// = new int[maxV - 1]; // node capacity

		// demands
		int nD; int nD2;
		int[][] h;// = new int[maxV - 1][maxV - 1]; // flow rate of demand (s,t)
		int[][][] f;// = new int [maxV - 1][maxV - 1][maxF]; // required VNFs of
					// demand (s,t)

		// VARIABLES
		double[] w;// = new int[maxE-1]; // metric
		double[][] l;// = new int[maxV-1][maxV-1]; // length of the shortest
						// path

		// x(epst): rate of flow p from node s to node t on link e
		double[][][][] x;// = new float[maxE-1][maxP][maxV-1][maxV-1] ;

		double[][] y;// = new float[maxV-1][maxV-1];

		// u(e,t) = 1: link e is on the shortest path to node t
		double[][] u;// = new int[maxE-1][maxV-1];

		String input_file_param = "3data.txt";
		File inFile1;
		inFile1 = new File(input_file_param);
		Scanner scanner = new Scanner(inFile1, "UTF-8");

		// maxP
		String line = scanner.nextLine();
		Scanner lineScan = new Scanner(line);
		maxP = Integer.parseInt(lineScan.next());

		// nDTest: so demand se chay
		line = scanner.nextLine();
		lineScan = new Scanner(line);
		int nDTest = Integer.parseInt(lineScan.next());
		
		// VNF
		line = scanner.nextLine();
		line = scanner.nextLine();
		lineScan = new Scanner(line);
		nF = Integer.parseInt(lineScan.next());
		rf = new int[nF];

		for (int i = 0; i < nF; i++) {
			line = scanner.nextLine();
			lineScan = new Scanner(line);
			int i1 = Integer.parseInt(lineScan.next());
			rf[i1] = Integer.parseInt(lineScan.next());
		}

		// G(V,E): E
		line = scanner.nextLine();
		line = scanner.nextLine();
		lineScan = new Scanner(line);
		nE = Integer.parseInt(lineScan.next());
		ii = new int[nE]; // starting node of link e
		jj = new int[nE]; // terminating node of link e
		ce = new double[nE]; // link capacity

		for (int i = 0; i < nE; i++) {
			line = scanner.nextLine();
			lineScan = new Scanner(line);
			int e1 = Integer.parseInt(lineScan.next());
			int ie1 = Integer.parseInt(lineScan.next());
			int je1 = Integer.parseInt(lineScan.next());
			int ce1 = Integer.parseInt(lineScan.next());
			ii[e1] = ie1;
			jj[e1] = je1;
			ce[e1] = ce1;
		}

		// G(V,E): V
		line = scanner.nextLine();
		line = scanner.nextLine();
		lineScan = new Scanner(line);
		nV = Integer.parseInt(lineScan.next());
		cv = new int[nV];

		for (int i = 0; i < nV; i++) {
			line = scanner.nextLine();
			lineScan = new Scanner(line);
			int v1 = Integer.parseInt(lineScan.next());
			int cv1 = Integer.parseInt(lineScan.next());
			cv[v1] = cv1;
		}

		// demands
		line = scanner.nextLine();
		line = scanner.nextLine();
		lineScan = new Scanner(line);
		nD = Integer.parseInt(lineScan.next());
		nD2 = nD;
		nD = nDTest;
		h = new int[nV][nV]; // flow rate of demand (s,t)
		f = new int[nV][nV][nF + 1]; // required VNFs of demand (s,t)

		for (int i = 0; i < nV; i++)
			for (int j = 0; j < nV; j++) {
				h[i][j] = 0;
				f[i][j][nF] = 0;
			}

		for (int i = 0; i < nD2; i++) {
			line = scanner.nextLine();
			lineScan = new Scanner(line);
			if (i < nD){
				int s1 = Integer.parseInt(lineScan.next());
				int t1 = Integer.parseInt(lineScan.next());
				int h1 = Integer.parseInt(lineScan.next());
				int nf1 = Integer.parseInt(lineScan.next());
				h[s1][t1] = h1;
				f[s1][t1][nF] = nf1;
				for (int j = 0; j < nf1; j++) {
					int f1 = Integer.parseInt(lineScan.next());
					f[s1][t1][j] = f1;
				}				
			}
		}

		// k(vsti) = 1: the ith VNF of demand (s,t) is deployed at node v
		k = new int[nV][nV][nV][nF];//
		for (int v = 0; v < nV; v++)
			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					for (int i = 0; i < nF; i++) {
						k[v][s][t][i] = 0;
					}

		line = scanner.nextLine();
		for (int i = 0; i < nD; i++) {
			line = scanner.nextLine();
			lineScan = new Scanner(line);
			int s1 = Integer.parseInt(lineScan.next());
			int t1 = Integer.parseInt(lineScan.next());
			for (int j = 0; j < f[s1][t1][nF]; j++) {
				line = scanner.nextLine();
				lineScan = new Scanner(line);
				int nDeploymentNode = Integer.parseInt(lineScan.next());
				for (int inode = 0; inode < nDeploymentNode; inode++) {
					int nn = Integer.parseInt(lineScan.next());
					k[nn][s1][t1][j] = 1;
				}
			}
		}

		w = new double[nE]; // metric
		l = new double[nV][nV]; // length of the shortest path

		// x(epst): rate of flow p from node s to node t on link e
		x = new double[nE][maxP][nV][nV];

		y = new double[nV][nV];

		// u(e,t) = 1: link e is on the shortest path to node t
		u = new double[nE][nV];

		double cost = fCostOptMetric(nF, nV, nE, maxP, rf, ii, jj, ce, cv, h,
				f, k, w, l, x, y, u);

		// Write details of output
		/*
		 * try { out.write("====Server=======================\r\n");
		 * out.write("\r\n Number of servers: " + nS);
		 * out.write("\r\n Min/Max resource: " + minCapS + "/" + maxCapS);
		 * out.write("\r\n Computing resource of servers: \r\n");
		 * writeArrayInt1(out, dS); out.write("\r\n Min/Max cost: " + minCost +
		 * "/" + maxCost);
		 * out.write("\r\n Communication cost among servers: \r\n");
		 * writeArrayInt2(out, w);
		 * 
		 * out.write("\r\n====Old State=======================\r\n"); out.write(
		 * "\r\n Number of apps, min/max number of components of one app: " +
		 * nApp + ", " + minC + "/" + maxC);
		 * out.write("\r\n Number of components: " + nC);
		 * out.write("\r\n Links among components: \r\n"); writeArrayInt2(out,
		 * g);
		 * 
		 * out.write("\r\n====New State=======================\r\n");
		 * out.write("\r\n Number of components: " + nC2);
		 * out.write("\r\n Computing resource required of components: \r\n");
		 * writeArrayInt1(out, r);
		 * out.write("\r\n Links among components: \r\n"); writeArrayInt2(out,
		 * g2);
		 * 
		 * out.write("\r\n====Old Deployment=======================\r\n");
		 * writeArrayInt2(out, u);
		 * 
		 * out.write("\r\n====New Deployment=======================\r\n");
		 * writeArrayInt2(out, aDeploy); out.write("\r\n");
		 * out.write("Optimal Cost: " + Double.toString(c)); } catch
		 * (IOException e1) { e1.printStackTrace(); } out.close();
		 */
	}

	static void writeArrayInt1(BufferedWriter out, int[] array)
			throws IOException {
		DecimalFormat df = new DecimalFormat("#");
		for (int i = 0; i < array.length; i++) {
			out.write(df.format(array[i]) + "\t");
		}
		out.write("\r\n");
	}

	static void writeArrayDouble2(BufferedWriter out, double[][] array)
			throws IOException {
		DecimalFormat df = new DecimalFormat("#.####");
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				out.write(df.format(array[i][j]) + "\t");
			}
			out.write("\r\n");
		}
	}

	static void writeArrayInt2(BufferedWriter out, int[][] array)
			throws IOException {
		DecimalFormat df = new DecimalFormat("#");
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				out.write(df.format(array[i][j]) + "\t");
			}
			out.write("\r\n");
		}
	}

	// INPUT
	// VNF
	// int nF; // number of available functions
	// int [] rf;// = new int[maxF-1]; //required computing resource for
	// processing one unit of flow rate

	// G(V,E)
	// int nV, nE;
	// int [] ii;// = new int[maxE - 1]; // starting node of link e
	// int [] jj;// = new int[maxE - 1]; // terminating node of link e
	// int [] ce;// = new int[maxE - 1]; // link capacity
	// int [] cv;// = new int[maxV - 1]; // node capacity

	// demands
	// int nD;
	// int [][] h;// = new int[maxV - 1][maxV - 1]; // flow rate of demand (s,t)
	// int [][][] f;// = new int [maxV - 1][maxV - 1][maxF]; // required VNFs of
	// demand (s,t)

	static double fCostOptMetric(int nF, int nV, int nE, int maxP, int[] rf,
			int[] ii, int[] jj, double[] ce, int[] cv, int[][] h, int[][][] f,
			int[][][][] k, double[] ow, double[][] ol, double[][][][] ox,
			double[][] oy, double[][] ou) throws IOException {
		double cost = -1;

		try {
			// Create the modeler/solver object
			IloCplex model = new IloCplex();

			OutputStream os = new FileOutputStream("CplexLog_discrete.txt");

			// model.setOut(null);
			model.setOut(os);
			IloIntVar[][] vari = new IloIntVar[1][];
			IloSemiContVar[][] varf = new IloSemiContVar[1][];
			IloRange[][] rng = new IloRange[1][];

			// The created ranges and variables are returned as element 0 of
			// arrays
			// var and rng.
			double mEpsilon = Double.MIN_VALUE;
			
			// VARIABLES
			IloIntVar[] mw = new IloIntVar[nE]; // metric
			IloIntVar[][] ml = new IloIntVar[nV][nV]; // length of the shortest
														// path

			// x(epst): rate of flow p from node s to node t on link e
			IloSemiContVar[][][][] mx = new IloSemiContVar[nE][maxP][nV][nV];
			IloIntVar[][][][] mb = new IloIntVar[nE][maxP][nV][nV];
			IloSemiContVar[][] my = new IloSemiContVar[nV][nV];
			double[][][][] ob = new double[nE][maxP][nV][nV];
			
			// u(e,t) = 1: link e is on the shortest path to node t
			IloIntVar[][] mu = new IloIntVar[nE][nV];

			// k(vsti) = 1: the ith VNF of demand (s,t) is deployed at node v
			// IloIntVar[][][][] mk = new IloIntVar[nV][nV][nV][nF];

			// mre[e]: penanty cost on link e
			//IloSemiContVar[] mre = new IloSemiContVar[nE];
			IloSemiContVar robj;
			
			// IloIntVar[] zi = new IloIntVar[nE + nV*nV + nE*nV + nV*nV*nV*nF];
			// IloSemiContVar[] zf = new IloSemiContVar[nE*maxP*nV*nV + nE*nV +
			// nE];
			// vari[0] = zi;
			// varf[0] = zf;

			// //VARIABLES: Set range of decision variables z
			// metric of nE links
			int maxW = 0;
			for (int i = 0; i < nV; i++)
				for (int j = 0; j < nV; j++) {
					maxW = maxW + h[i][j]; // sum of all demands
				}
			for (int i = 0; i < nE; i++)
				// zi[i] = model.intVar(1, maxW);
				mw[i] = model.intVar(1, maxW);
			// metric of 2 links = 1, maxW when all demands go through one link
			// then splited to two links
			// and...

			// length of the shortest path in metric w
			int maxLength = nE * maxW;
			for (int i = 0; i < nV; i++)
				for (int j = 0; j < nV; j++)
					// zi[nE + i*nV + j] = model.intVar(1, nE);
					if (i != j)
						ml[i][j] = model.intVar(1, maxLength);
					else
						ml[i][j] = model.intVar(0, 0);

			// u(e,t) = 1: link e is on the shortest path to node t
			// int [][] u;// = new int[nE][nV];
			for (int i = 0; i < nE; i++)
				for (int j = 0; j < nV; j++)
					// zi[nE + nV*nV + i*nV + j] = model.intVar(0, 1);
					mu[i][j] = model.intVar(0, 1);

			// k(vsti) = 1: the ith VNF of demand (s,t) is deployed at node v
			// int [][][][] k;// = new int[nV][nV][nV][nF];//
			// for (int v = 0; v < nV; v++)
			// for (int s = 0; s < nV; s++)
			// for (int t = 0; t < nV; t++)
			// for (int i = 0; i < nF; i++)
			// // zi[nE + nV*nV + nE*nV + v*nV*nV*nF + s*nV*nF +
			// // t*nF + i] =
			// // model.intVar(0, 1);
			// mk[v][s][t][i] = model.intVar(0, 1);

			// x(epst): rate of flow p from node s to node t on link e
			// float [][][][] x;// = new float[nE][maxP][nV][nV] ;
			float maxRate = 0;
			for (int i = 0; i < nV; i++)
				for (int j = 0; j < nV; j++)
					if (h[i][j] > maxRate)
						maxRate = h[i][j];

			for (int e = 0; e < nE; e++)
				for (int p = 0; p < maxP; p++)
					for (int s = 0; s < nV; s++)
						for (int t = 0; t < nV; t++)
							// zf[e*maxP*nV*nV + p*nV*nV + s*nV + t] =
							// model.semiContVar((double)0, (double)maxRate,
							// IloNumVarType.Float);
							if (h[s][t] > 0) {
								mx[e][p][s][t] = model.semiContVar((double) 0,
										(double) maxRate, IloNumVarType.Float);
								mb[e][p][s][t] = model.intVar(0, 1);
							} else {
								mx[e][p][s][t] = model.semiContVar((double) 0,
										(double) 0, IloNumVarType.Float);
								mb[e][p][s][t] = model.intVar(0, 1);
							}

			// float [][] y;// = new float[nV][nV];
			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					// zf[nE*maxP*nV*nV + s*nV + t] =
					// model.semiContVar((double)0, (double)maxRate,
					// IloNumVarType.Float);
					if (s != t)
						my[s][t] = model.semiContVar((double) 0,
								(double) maxRate, IloNumVarType.Float);
					else
						my[s][t] = model.semiContVar((double) 0, (double) 0,
								IloNumVarType.Float);

			// float [] re;
			// can chinh lai cho chinh xac sau
			/*
			for (int e = 0; e < nE; e++)
				// zf[nE*maxP*nV*nV + nV*nV + e] =
				// model.semiContVar((double)0, (double)5000,
				// IloNumVarType.Float);
				mre[e] = model.semiContVar((double) 0,
						(double) Float.MAX_VALUE, IloNumVarType.Float);
			*/
			robj = model.semiContVar((double) 0, (double) 1,
					IloNumVarType.Float);
			// OBJECTIVE
			IloNumExpr QQ = model.constant(0);
			//for (int e = 0; e < nE; e++) {
			//	QQ = model.sum(QQ, mre[e]);
			//}
			QQ = model.sum(QQ, robj);
			model.addMinimize(QQ);

			// CONSTRAINTS
			// Flow balance
			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if (h[s][t] > 0)
						for (int p = 0; p < maxP; p++)
							for (int v = 0; v < nV; v++)
								if ((s != t) && (v != s) && (v != t)) {
									IloNumExpr outv = model.constant(0);
									IloNumExpr inv = model.constant(0);
									for (int e = 0; e < nE; e++) {
										if (ii[e] == v)
											outv = model.sum(outv,
													mx[e][p][s][t]);
										if (jj[e] == v)
											inv = model
													.sum(inv, mx[e][p][s][t]);
									}
									model.addEq(outv, inv);
								}

			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if (h[s][t] > 0)
						if (s != t) {
							IloNumExpr outst = model.constant(0);
							for (int p = 0; p < maxP; p++)
								for (int e = 0; e < nE; e++) {
									if (ii[e] == s)
										outst = model
												.sum(outst, mx[e][p][s][t]);
								}

							model.addEq(outst, h[s][t]);

						}

			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if ((h[s][t] > 0) && (s != t)) {
						IloNumExpr inst = model.constant(0);
						for (int p = 0; p < maxP; p++)
							for (int e = 0; e < nE; e++) {
								if (jj[e] == t)
									inst = model.sum(inst, mx[e][p][s][t]);
							}

						model.addEq(inst, h[s][t]);
					}

			// Traffic flow of demand (s,t) will divided into exactly maxP flows
			// Traffic of every flow > 0, but can go through the same path
			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if (h[s][t] > 0)
						for (int p = 0; p < maxP; p++) {
							IloNumExpr t1 = model.constant(0);
							for (int e = 0; e < nE; e++) {
								t1 = model.sum(t1, mx[e][p][s][t]);
							}
							t1 = model.sum(t1, -mEpsilon);
							model.addGe(t1, 0);
						}

			// Link capacity
			for (int e = 0; e < nE; e++) {
				IloNumExpr inst = model.constant(0);
				for (int s = 0; s < nV; s++)
					for (int t = 0; t < nV; t++)
						if (h[s][t] > 0)
							for (int p = 0; p < maxP; p++)
								inst = model.sum(inst, mx[e][p][s][t]);
				IloNumExpr t2 = model.constant(0);
				t2 = model.sum(t2, robj);
				t2 = model.prod(t2, ce[e]);
				t2 = model.prod(t2, -1);
				inst = model.sum(inst, t2);
				model.addLe(inst, 0);
			}

			// Load balancing
			for (int e = 0; e < nE; e++)
				for (int t = 0; t < nV; t++)
					if (ii[e] != t) {
						IloNumExpr inst = model.constant(0);

						for (int s = 0; s < nV; s++)
							if (h[s][t] > 0)
								for (int p = 0; p < maxP; p++)
									inst = model.sum(inst, mx[e][p][s][t]);
						inst = model.prod(inst, -1);
						inst = model.sum(inst, my[ii[e]][t]);

						model.addGe(inst, 0);

						IloNumExpr eright = model.constant(0);
						for (int s = 0; s < nV; s++)
							if (h[s][t] > 0)
								eright = model.sum(eright, h[s][t]);
						IloNumExpr t2 = model.constant(0);
						t2 = model.sum(t2, mu[e][t]);
						t2 = model.prod(t2, -1);
						t2 = model.sum(t2, 1);
						eright = model.prod(t2, eright);

						model.addLe(inst, eright);
					}

			// shortest path routing
			for (int e = 0; e < nE; e++)
				for (int t = 0; t < nV; t++)
					for (int s = 0; s < nV; s++)
						if (h[s][t] > 0) {
							IloNumExpr inst = model.constant(0);
							for (int p = 0; p < maxP; p++)
								inst = model.sum(inst, mx[e][p][s][t]);

							IloNumExpr eright = model.constant(0);
							eright = model.sum(eright, mu[e][t]);
							eright = model.prod(eright, h[s][t]);

							model.addLe(inst, eright);
						}

			for (int e = 0; e < nE; e++)
				for (int t = 0; t < nV; t++) {
					IloNumExpr inst = model.constant(0);
					inst = model.sum(inst, ml[ii[e]][t]);
					inst = model.prod(inst, -1);
					inst = model.sum(inst, ml[jj[e]][t]);
					inst = model.sum(inst, mw[e]);

					IloNumExpr t2 = model.constant(0);
					t2 = model.sum(t2, mu[e][t]);
					t2 = model.prod(t2, -1);
					t2 = model.sum(t2, 1);

					model.addGe(inst, t2);

					t2 = model.prod(t2, maxLength);
					model.addLe(inst, t2);
				}

			// VNFs
			// for (int s = 0; s < nV; s++)
			// for (int t = 0; t < nV; t++)
			// for (int i = 0; i < f[s][t][nF]; i++)
			// // the ith VNF of demand (s,t)
			// for (int p = 0; p < maxP; p++) {
			// IloNumExpr t1 = model.constant(0);
			// IloNumExpr t2 = model.constant(0);
			// for (int e = 0; e < nE; e++) {
			// t1 = model.sum(t1, k[ii[e]][s][t][i]);
			// t1 = model.sum(t1, k[jj[e]][s][t][i]);
			// t2 = model.sum(t2, mx[e][p][s][t]);
			// }
			// t1 = model.prod(t1, t2);
			// model.addGe(t1, t2);
			// }

			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if (h[s][t] > 0)
						for (int i = 0; i < f[s][t][nF]; i++)
							// the ith VNF of demand (s,t)
							for (int p = 0; p < maxP; p++) {
								IloNumExpr t1 = model.constant(0);
								for (int e = 0; e < nE; e++) {
									IloNumExpr t2 = model.constant(0);
									t2 = model.sum(t2, k[ii[e]][s][t][i]);
									t2 = model.sum(t2, k[jj[e]][s][t][i]);
									t2 = model.prod(t2, mx[e][p][s][t]);
									t1 = model.sum(t1, t2);
								}
								t1 = model.sum(t1, -1);
								model.addGe(t1, 0);
							}

			// Do not split flow p
			for (int s = 0; s < nV; s++)
				for (int t = 0; t < nV; t++)
					if (h[s][t] > 0)
						for (int p = 0; p < maxP; p++)
							for (int e = 0; e < nE; e++) 
							if ((ii[e] != s) && (jj[e] != t)){
								IloNumExpr t1 = model.constant(0);
								t1 = model.sum(t1, maxRate);
								t1 = model.prod(t1, mb[e][p][s][t]);
								model.addLe(mx[e][p][s][t], t1);
								
								t1 = model.constant(0);
								for (int e2 = 0; e2 < nE; e2++)
									if (jj[e2] == ii[e]) {
										t1 = model.sum(t1, mx[e2][p][s][t]);
									}
								model.addLe(mx[e][p][s][t], t1);
								
								IloNumExpr t2 = model.constant(0);
								t2 = model.sum(t2, mb[e][p][s][t]);
								t2 = model.prod(t2, -1);
								t2 = model.sum(t2, 1);
								t2 = model.prod(t2, -maxRate);	
								t2 = model.sum(t2, t1);
								
								model.addGe(mx[e][p][s][t], t2);
							}

			// Node capacity
			for (int v = 0; v < nV; v++) {
				IloNumExpr t2 = model.constant(0);

				for (int s = 0; s < nV; s++)
					for (int t = 0; t < nV; t++)
						for (int i = 0; i < f[s][t][nF]; i++)
							if (k[v][s][t][i] == 1) {
								IloNumExpr t1 = model.constant(0);
								for (int p = 0; p < maxP; p++)
									for (int e = 0; e < nE; e++)
										if (jj[e] == v)
											t1 = model.sum(t1, mx[e][p][s][t]);
								t1 = model.prod(t1, rf[f[s][t][i]]);
								t2 = model.sum(t2, t1);
							}
				model.addLe(t2, cv[v]);
			}

			/*
			// Penalty cost
			for (int e = 0; e < nE; e++) {
				IloNumExpr ye = model.constant(0);

				for (int s = 0; s < nV; s++)
					for (int t = 0; t < nV; t++)
						if (h[s][t] > 0)
							for (int p = 0; p < maxP; p++)
								ye = model.sum(ye, mx[e][p][s][t]);
				model.addGe(mre[e], ye);

				IloNumExpr t1 = model.constant(0);
				t1 = model.sum(t1, ye);
				t1 = model.prod(t1, 3.0f);
				t1 = model.sum(t1, (-2.0f / 3) * ce[e]);
				model.addGe(mre[e], t1);

				t1 = model.constant(0);
				t1 = model.sum(t1, ye);
				t1 = model.prod(t1, 10.0f);
				t1 = model.sum(t1, (-16.0f / 3) * ce[e]);
				model.addGe(mre[e], t1);

				t1 = model.constant(0);
				t1 = model.sum(t1, ye);
				t1 = model.prod(t1, 70.0f);
				t1 = model.sum(t1, (-178.0f / 3) * ce[e]);
				model.addGe(mre[e], t1);

				t1 = model.constant(0);
				t1 = model.sum(t1, ye);
				t1 = model.prod(t1, 500.0f);
				t1 = model.sum(t1, (-1468.0f / 3) * ce[e]);
				model.addGe(mre[e], t1);

				t1 = model.constant(0);
				t1 = model.sum(t1, ye);
				t1 = model.prod(t1, 5000.0f);
				t1 = model.sum(t1, (-16318.0f / 3) * ce[e]);
				model.addGe(mre[e], t1);

			}
			*/
			// rng[0] = new IloRange[nS + nNew];

			// write model to file
			model.exportModel("lpex1.lp");

			// solve the model and display the solution if one was found
			if (model.solve()) {

				ow = model.getValues(mw);

				for (int i = 0; i < nV; i++)
					for (int j = 0; j < nV; j++)
						ol[i][j] = model.getValue(ml[i][j]);

				for (int i = 0; i < nE; i++)
					for (int j = 0; j < nV; j++)
						ou[i][j] = model.getValue(mu[i][j]);

				for (int e = 0; e < nE; e++)
					for (int p = 0; p < maxP; p++)
						for (int s = 0; s < nV; s++)
							for (int t = 0; t < nV; t++)
								if (h[s][t] > 0)
									ox[e][p][s][t] = model
											.getValue(mx[e][p][s][t]);

				for (int s = 0; s < nV; s++)
					for (int t = 0; t < nV; t++)
						if (h[s][t] > 0)
							oy[s][t] = model.getValue(my[s][t]);

				for (int s = 0; s < nV; s++)
					for (int t = 0; t < nV; t++)
						if (h[s][t] > 0)
							for (int p = 0; p < maxP; p++)
								for (int e = 0; e < nE; e++) 
								if ((ii[e] != s) && (jj[e] != t)){
									ob[e][p][s][t] = model.getValue(mb[e][p][s][t]);
								}
			}
			cost = model.getObjValue();
			model.end();
		} catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
		return cost;
	}
}
