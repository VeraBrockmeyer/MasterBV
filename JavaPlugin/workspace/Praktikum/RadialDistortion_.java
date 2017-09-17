import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * IJ Plugin für die Berechnung der Radialen Verzerrung anhand von Pixel-Koordinaten paaren
 * @author Vera Brockmeyer
 * @author Artjom Schwabski
 *
 */
public class RadialDistortion_ implements PlugInFilter {

	/**
	 * Speicher für das Vorlagenbild zum Abrufen der Pixel Werte zur Zeichnung der entzerrung
	 */
	private ImagePlus distPicture;

	/**
	 * Abstand der Gitter-Kreuzpunkte in pxl (zur Berechung des optimalen Gitters)
	 */
	private static final int distCross = 111;

	/**
	 * Anzahl der Gitter-Kreuzpunkte vom Mittelpunkt bis zum Eckpunkt in x-Richtung (zur Berechung des optimalen Gitters)
	 */
	private static final int nXCross2Corner = 10;

	/**
	 * Anzahl der Gitter-Kreuzpunkte vom Mittelpunkt bis zum Eckpunkt in y-Richtung (zur Berechung des optimalen Gitters)
	 */
	private static final int nYCross2Corner = 6;

	/**
	 * Anzahl der Gitter-Spalten zur Berechnung des optimalen Gitters
	 */
	private static final int nCol = 21;

	/**
	 * Anzahl der Gitter-Reihen zur Berechnung des optimalen Gitters
	 */
	private static final int nRow = 13;

	



	@Override
	/**
	 * Main Methode des Plugins startet das einlesen der Punkte die Berechung der entzerrung und das Zeichnen des entzerrten Bildes
	 */
	public void run(ImageProcessor img) {
		try {
			ArrayList<PointPair> PointPairs = readData();
			computeDrawRadialTransformation(distPicture, PointPairs);
			//drawTargets(distPicture.getProcessor(), "SourceImage",PointPairs);
			//sourcePicture=computeDrawAffineTransformation(sourcePicture, xPointPairs, yPointPairs);
			//sourcePicture=computeDrawPerpectiveTransformation(sourcePicture, xPointPairs, yPointPairs);

		} catch (Exception exc) {
			IJ.error(exc.getMessage() + exc.getClass() + exc.getCause() + exc.getStackTrace());
		}
	}


	/**
	 * Öffnet die Auswahl eines Textfensters und lädt alle Punkt paare des Text Fensters in SimplePair Objekte, die intern gespeichert werden
	 * Zusätzlich kann hier das Optimale Gitter berechnet werden um in UwrapJ die Zielpunkte auswählen zu können
	 */
	private ArrayList<PointPair> readData() {

		int xCenter = distPicture.getProcessor().getWidth() / 2;
		int yCenter = distPicture.getProcessor().getHeight() / 2;
		
		ArrayList<PointPair> pairs = new ArrayList<PointPair>();
		// Textfenster mit Punktpaare Textdtei wählen:
		java.awt.Window[] non_img_windows = WindowManager.getAllNonImageWindows();
		GenericDialog gd = new GenericDialog("Textdatei auswählen:");
		String[] Windows_names = new String[non_img_windows.length];
		for (int i = 0; i < non_img_windows.length; i++) {
			if (non_img_windows[i] instanceof ij.plugin.frame.Editor) {
				Windows_names[i] = ((ij.plugin.frame.Editor) non_img_windows[i]).getTitle();
			}
		}

		String choise = "";
		gd.addChoice("Textdatei Fenster", Windows_names, choise);

		gd.showDialog();

		// abbrechen wenn keine eingabe:
		if (gd.wasCanceled())
			return null;
		choise = gd.getNextChoice();

		// Punktpaare einlesen:

		for (int i = 0; i < non_img_windows.length; i++) {
			if (((java.awt.Frame) non_img_windows[i]).getTitle() == choise) {

				String text = ((ij.plugin.frame.Editor) non_img_windows[i]).getText();

				String[] lines = text.split("\n");

				// speicher für punkt paarungen

				for (int j = 1; j < lines.length; j++) // erste zeile
					// überspringen
				{
					String[] numbers = lines[j].split("\t"); // spalte zeile
					// anhand des
					// tab zeichens
					
					double x_dist = Integer.parseInt(numbers[1].trim());
					double y_dist = Integer.parseInt(numbers[2].trim());
					
					double x_undist = Integer.parseInt(numbers[3].trim());
					double y_undist = Integer.parseInt(numbers[4].trim());
					
					int index = Integer.parseInt(numbers[0].trim());

					PointPair pp = new PointPair(x_dist, y_dist, x_undist, y_undist);
					pp.r=pp.computeRadius2Center(x_dist, y_dist, xCenter, yCenter);

					
					//
					// //berechne x_target und y_target und gebe sie aus - NUR
					// zur einmaligen generierung der Pointpairs:
//					 int colid = (int)(index / nRow); // 0 - 8 reihe
//					 int rowid = (int)(index - colid * nRow) ; // 0- 15 spalte
//					
//					 int x_offset = xCenter - nXCross2Corner * distCross;
//					 //koorinate mittelpunkt gitter - anzahl der gitterpunkte nach links
//					 int y_offset = yCenter - nYCross2Corner * distCross;
//					 //koorinate mittelpunkt gitter - anzahl der gitterpunkte nach oben
//					
//					 pp.x_undist = colid * distCross+ x_offset;
//					 pp.y_undist = rowid * distCross + y_offset;
//					
//					 ImageProcessor debug = distPicture.getProcessor().duplicate();
//					 debug.drawOval((int)pp.x_undist, (int)pp.y_undist, 3,3);
//
					pairs.add(pp);

				}

			}
		}
		System.out.println("Read Data List");
		System.out.println("index \t x_dist \t y_dist \t x_undist \t y_undist \t radius \n");
		for (int j = 0; j < pairs.size(); j++) {
			PointPair p = pairs.get(j);
			System.out.println(""+j+"\t"+p.x_dist+"\t"+p.y_dist+"\t"+p.x_undist+"\t"+p.y_undist+"\t"+p.r);
		}
		return pairs;
	}


	public static ImagePlus computeDrawPerpectiveTransformation(ImagePlus sourcePicture, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs){
		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
	
		double[][]B = new double[xPointPairs.size()*2][9];
		int counter = 0;
		for (int i = 0; i < B.length; i+=2) {
			
			B[i][0] = -xPointPairs.get(counter).undistorted;
			B[i][1] = -yPointPairs.get(counter).undistorted;
			B[i][2] = -1.;
			B[i][3] = 0.;
			B[i][4] = 0.;
			B[i][5] = 0.;
			B[i][6] = xPointPairs.get(counter).undistorted*xPointPairs.get(counter).distorted;
			B[i][7] = yPointPairs.get(counter).undistorted*xPointPairs.get(counter).distorted;
			B[i][8] = xPointPairs.get(counter).distorted;
			B[i+1][0] = 0.;
			B[i+1][1] = 0.;
			B[i+1][2] = 0.;
			B[i+1][3] = -xPointPairs.get(counter).undistorted;
			B[i+1][4] = -yPointPairs.get(counter).undistorted;
			B[i+1][5] = -1.;
			B[i+1][6] = xPointPairs.get(counter).undistorted*yPointPairs.get(counter).distorted;
			B[i+1][7] = yPointPairs.get(counter).undistorted*yPointPairs.get(counter).distorted;
			B[i+1][8] = yPointPairs.get(counter).distorted;
			counter++;
		}
		RealMatrix A = new Array2DRowRealMatrix(B);
		RealMatrix At = A.transpose();
		RealMatrix AtA = At.multiply(A);

		EigenDecomposition ed = new EigenDecomposition(AtA);
		RealMatrix V = ed.getV();
	
		double[][] pArray = new double[3][3];
		
		pArray[0][0]=V.getEntry(0,V.getRowDimension()-1);
		pArray[0][1]=V.getEntry(1,V.getRowDimension()-1);
		pArray[0][2]=V.getEntry(2,V.getRowDimension()-1);
		
		pArray[1][0]=V.getEntry(3,V.getRowDimension()-1);
		pArray[1][1]=V.getEntry(4,V.getRowDimension()-1);
		pArray[1][2]=V.getEntry(5,V.getRowDimension()-1);
		
		pArray[2][0]=V.getEntry(6,V.getRowDimension()-1);
		pArray[2][1]=V.getEntry(7,V.getRowDimension()-1);
		pArray[2][2]=V.getEntry(8,V.getRowDimension()-1);		
	
		RealMatrix pMat = new Array2DRowRealMatrix(pArray);
	
		for (int y = 0; y < targetImg.getHeight(); y++) {
			for (int x = 0; x < targetImg.getWidth(); x++) {
				double[] target  = {x,y,1};
				RealVector t_vec = new ArrayRealVector(target, false);
				RealVector coord_vec = pMat.operate(t_vec);
				double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
				double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
			
				if (x_coord_vorlage>=0 && x_coord_vorlage < targetImg.getWidth() && y_coord_vorlage>=0 && y_coord_vorlage < targetImg.getHeight()) {
					targetImg.putPixel(x, y, (int) Math.round(sourcePicture.getProcessor().getInterpolatedPixel(x_coord_vorlage , y_coord_vorlage)));
				} 
				else
				{
					targetImg.putPixel(x, y, 0);
				}
			}
		}
		
		for (int i = 0; i < xPointPairs.size(); i++) {
			double[] target  = {xPointPairs.get(i).distorted,yPointPairs.get(i).distorted,1};
			RealVector t_vec = new ArrayRealVector(target, false);
			RealVector coord_vec = pMat.operate(t_vec);
			double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
			double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
			xPointPairs.get(i).distorted = x_coord_vorlage;
			yPointPairs.get(i).distorted = y_coord_vorlage;
		}

		drawTargets(targetImg, "Projective Transformation", xPointPairs, yPointPairs);
		return new ImagePlus("Projected Image",targetImg);
	}




	/**
	 * Startet die Berechnung der Radialen entzerrung und zeichnet das Ergebnis
	 * @param sourcePicture Verzerrtes IJ Bild
	 * @param xPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 * @param yPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 */
	public static void computeDrawRadialTransformation(ImagePlus sourcePicture, ArrayList<PointPair>pairs) {
		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);

		double xCenter = sourcePicture.getProcessor().getWidth() / 2;
		double yCenter = sourcePicture.getProcessor().getHeight() / 2;

		 // the model function components are the distances to current estimated center,
		 // they should be as close as possible to the specified radius
		 MultivariateJacobianFunction distancesToCurrentCenter = new MultivariateJacobianFunction() {
		      public Pair<RealVector, RealMatrix> value(final RealVector point) {

		          Vector2D center = new Vector2D(point.getEntry(0), point.getEntry(1));

		          RealVector value = new ArrayRealVector(observedPoints.length);
		          RealMatrix jacobian = new Array2DRowRealMatrix(observedPoints.length, 2);

		          for (int i = 0; i < observedPoints.length; ++i) {
		              Vector2D o = observedPoints[i];
		              double modelI = Vector2D.distance(o, center);
		              value.setEntry(i, modelI);
		              // derivative with respect to p0 = x center
		              jacobian.setEntry(i, 0, (center.getX() - o.getX()) / modelI);
		              // derivative with respect to p1 = y center
		              jacobian.setEntry(i, 1, (center.getX() - o.getX()) / modelI);
		          }

		          return new Pair<RealVector, RealMatrix>(value, jacobian);

		      }
		  };

		  // the target is to have all points at the specified radius from the center
		  double[] prescribedDistances = new double[observedPoints.length];
		  Arrays.fill(prescribedDistances, radius);

		  // least squares problem to solve : modeled radius should be close to target radius
		  LeastSquaresProblem problem = new LeastSquaresBuilder().
		                                start(new double[] { 100.0, 50.0 }).
		                                model(distancesToCurrentCenter).
		                                target(prescribedDistances).
		                                lazyEvaluation(false).
		                                maxEvaluations(1000).
		                                maxIterations(1000).
		                                build();
		  LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);
		  Vector2D fittedCenter = new Vector2D(optimum.getPoint().getEntry(0), optimum.getPoint().getEntry(1));
		  System.out.println("fitted center: " + fittedCenter.getX() + " " + fittedCenter.getY());
		  System.out.println("RMS: "           + optimum.getRMS());
		  System.out.println("evaluations: "   + optimum.getEvaluations());
		  System.out.println("iterations: "    + optimum.getIterations());
	}

	/**
	 * Zeichnet Punkte an die Stellen der Zielkoorinaten in das übergebene Bild und bringt es auf den Bildschirm
	 * @param ip Bild in das die Punkte gezeichnet werden
	 * @param s Name des Bildes
	 * @param xPointPairs x-Koordinaten 
	 * @param yPointPairs y-Koordinaten
	 */
	public static void drawTargets(ImageProcessor ip, String s, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs) {
		// Punkte in radial entzerrtes Bild malen
		ImageProcessor res = ip.duplicate();
		// Zeichne ziel punkte:
		res.setLineWidth(3);
//		for (int i = 0; i < xPointPairs.size(); i++) {
//
//			res.setColor(Color.BLACK);
//			res.drawOval((int) xPointPairs.get(i).undistorted, (int) yPointPairs.get(i).undistorted, 3, 3);
//			res.setColor(Color.GRAY);
//			res.drawOval((int) xPointPairs.get(i).distorted, (int) yPointPairs.get(i).distorted, 3, 3);
//
//		}

		ImagePlus resImg = new ImagePlus(s, res);
		resImg.show();
		
	}


	/**
	 * Berechnet und zeichnet die Affine entzerrung anhand der x und y Punkt paare
	 * Die X und Y Arrays müssen eine zusammengehörige Reihenfolge haben
	 * @param xPointPairs Vorlage und Ziel-Koorinaten der x-Achse
	 * @param yPointPairs Vorlage und Ziel-Koorinaten der y-Achse
	 * @return Affiner Verzerrungsvektor p
	 */
	public static ImagePlus computeDrawAffineTransformation(ImagePlus sourcePicture, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs) {

		/// Verschiebungsvektor und matrix berechnen

		// Minimum von B*p -c;
		double[][] B = new double[2 * xPointPairs.size()][6];// Matrix B =
		// koordinaten
		// vorlage
		// double[] p;//p = geschter verschiebungsmatrix und vektor
		double[] c = new double[2 * xPointPairs.size()];// c = koordinaten
		// transformiert

		// konstruiert man die (2n,6) -Matrix
		for (int i = 0; i < xPointPairs.size(); i++) {
			// i = 0 => 0,1 => i *2, i*2+1
			// i = 1 => 2,3
			// i = 2 => 4,5
			// i = 3 => 6,7

			//
			c[i * 2] = xPointPairs.get(i).undistorted;// x_i'
			c[i * 2 + 1] = yPointPairs.get(i).undistorted; // y_i'

			B[i * 2][0] = xPointPairs.get(i).distorted; // x_i
			B[i * 2][1] = yPointPairs.get(i).distorted; // y_i
			B[i * 2][2] = 1;
			B[i * 2][3] = 0;
			B[i * 2][4] = 0;
			B[i * 2][5] = 0;

			B[i * 2 + 1][0] = 0;
			B[i * 2 + 1][1] = 0;
			B[i * 2 + 1][2] = 0;
			B[i * 2 + 1][3] = xPointPairs.get(i).distorted; // x_i
			B[i * 2 + 1][4] = yPointPairs.get(i).distorted; // y_i
			B[i * 2 + 1][5] = 1;
		}

		// B*p = c => p = (Bt * B)^-1 * Bt * c

		RealVector c_vec = new ArrayRealVector(c, false);
		RealMatrix B_mat = MatrixUtils.createRealMatrix(B);

		// Bt
		RealMatrix B_mat_transp = B_mat.transpose();

		// Bt * B
		RealMatrix Bt_B = B_mat_transp.multiply(B_mat);

		// (Bt * B)^-1
		RealMatrix Bt_B_inv = MatrixUtils.inverse(Bt_B);

		// (Bt * B)^-1 * Bt * c
		RealVector p = Bt_B_inv.operate(B_mat_transp.operate(c_vec));

		double[] verschiebungsmatrix = new double[] { p.getEntry(0), p.getEntry(1), p.getEntry(3), p.getEntry(4) };
		double[] verschiebungsvektor = new double[] { p.getEntry(2), p.getEntry(5) };
		double[][] verzerrungsmatrix = { { verschiebungsmatrix[0], verschiebungsmatrix[1], verschiebungsvektor[0] },
				{ verschiebungsmatrix[2], verschiebungsmatrix[3], verschiebungsvektor[1] }, { 0, 0, 1 } };

		IJ.log("a11: " + verschiebungsmatrix[0] + " a12: " + verschiebungsmatrix[1] + " a21: " + verschiebungsmatrix[2]
				+ " a22: " + verschiebungsmatrix[3] + " t1: " + verschiebungsvektor[0] + " t2: "
				+ verschiebungsvektor[1]);

		// Berechnen Sie aus den Werten die Umkehrabbildung,

		// erzeugen Sie ein neues Bild gleicher Größe wie das Eingabebild
		ImageProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);

		// berechnen Sie es nach dem Target-to-Source-Verfahren,
		// d.h.berechnen Sie zu jedem Punkt im neuen Bild durch die
		// Umkehrabbildung die Koordinaten im
		// Eingangsbild und holen Sie von dort den Wert des Pixels. Eine
		// Interpolation wird nicht verlangt.

		RealMatrix Bi_mat = MatrixUtils.createRealMatrix(verzerrungsmatrix);
		DecompositionSolver Bi_mat_solver = new LUDecomposition(Bi_mat).getSolver();

		// zielbild durchlaufen und werde aus vorlagebild abrufen
		for (int y = 0; y < targetImg.getHeight(); y++) {
			for (int x = 0; x < targetImg.getWidth(); x++) {
				double[] t_h = { x, y, 1 };
				RealVector t_vec = new ArrayRealVector(t_h, false);
				RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); // Transformation
				// berechnen

				int x_coord_vorlage = (int) Math.round(coord_vec.getEntry(0));
				int y_coord_vorlage = (int) Math.round(coord_vec.getEntry(1));

				if (x_coord_vorlage>=0 && x_coord_vorlage < targetImg.getWidth() && y_coord_vorlage>=0 && y_coord_vorlage < targetImg.getHeight()) {
					targetImg.putPixel(x, y, (int) Math.round(
							sourcePicture.getProcessor().getInterpolatedPixel(x_coord_vorlage, y_coord_vorlage)));
				} else // mit schwarz auffüllen
				{
					targetImg.putPixel(x, y, 0);
				}
			}
		}

		drawTargets(targetImg, "Affine", xPointPairs, yPointPairs);
		
		//Tranformation auf Punkt paare anwenden zur weiteren verwendung
		for (int i = 0; i < xPointPairs.size(); i++) {
			double[] t_h = { xPointPairs.get(i).distorted, yPointPairs.get(i).distorted, 1 };
			RealVector t_vec = new ArrayRealVector(t_h, false);
			RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); // Transformation
																				// berechnen
			xPointPairs.get(i).distorted = coord_vec.getEntry(0);
			yPointPairs.get(i).distorted = coord_vec.getEntry(1);
		}

		return new ImagePlus ("Affine",targetImg);

	}


	/**
	 * Berechnet den Abstand zum Gittermittelpunkt
	 * @param x x-Koordinate des Punktes
	 * @param y y-Koordinate des Punktes
	 * @param xCenter x-Koordinate des Mittelpunktes
	 * @param yCenter y-Koordinate des Mittelpunktes
	 * @return
	 */
	public static double computeRadius2Center(double x, double y, double xCenter, double yCenter) {
		return Math.sqrt((x - xCenter) * (x - xCenter) + (y - yCenter) * (y - yCenter));
	}

	/**
	 * 
	 * @param punkt_paare Array mit SimplePair Objekten in denen die Vorlage- und Ziel-Pixelkoordinaten gespeichert sind
	 * @return Koefizienten der Radialen-Verzerrung nach Levenberg-Marquadt
	 */
	public static double[] compute_radial_dist_koeff(List<SimplePair> punkt_paare) {
		RadialDistFunction qf = new RadialDistFunction(punkt_paare);
		LeastSquaresBuilder lsb = new LeastSquaresBuilder();

		// set model function and its jacobian
		lsb.model(qf.retMVF(), qf.retMMF());
		double[] newTarget = qf.realTargetPoints();

		// set target data
		lsb.target(newTarget);
		double[] newStart = { 1.e-15, 1.e-15, 1.e-15 };
		// set initial parameters
		lsb.start(newStart);
		// set upper limit of evaluation time
		lsb.maxEvaluations(9000);
		// set upper limit of iteration time
		lsb.maxIterations(20000);

		LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		try {
			// do LevenbergMarquardt optimization
			LeastSquaresOptimizer.Optimum lsoo = lmo.optimize(lsb.build());

			// get optimized parameters
			final double[] optimalValues = lsoo.getPoint().toArray();
			// output data
			IJ.log("A: " + optimalValues[0]);
			IJ.log("B: " + optimalValues[1]);
			IJ.log("C: " + optimalValues[2]);
			IJ.log("Iteration number: " + lsoo.getIterations());
			IJ.log("Evaluation number: " + lsoo.getEvaluations());

			return optimalValues;

		} catch (Exception e) {
			System.out.println(e.toString());
			return null;
		}

	}

	@Override
	/**
	 * Speichern des Voralge Bildes bei Aufruf des Plugins und erzeugung einer Kopie
	 */
	public int setup(String arg0, ImagePlus arg1) {
		this.distPicture = arg1;
		return DOES_ALL;
	}
}
