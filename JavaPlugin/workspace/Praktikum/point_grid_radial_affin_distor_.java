import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

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
public class point_grid_radial_affin_distor_ implements PlugInFilter {
	
	/**
	 * Speicher für das Vorlagenbild zum Abrufen der Pixel Werte zur Zeichnung der entzerrung
	 */
	private ImagePlus sourcePicture;

	
	/**
	 * Speicher für die x-Koordinate des Gitter-Mittelpunktes
	 */
	private int xCenter = 1084;// 970;//
	
	/**
	 * Speicher für die y-Koordinate des Gitter-Mittelpunktes
	 */
	private int yCenter = 713;// 652;//
	
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
	
	/**
	 * 
	 */
	private ImagePlus debugImg;
	
	/**
	 * 
	 */
	private ImageProcessor debug;
	
	/**
	 * Speicher für Koorinaten Paare der x-Achse
	 */
	private ArrayList<SimplePair> xPointPairs = new ArrayList<SimplePair>();
	
	/**
	 * Speicher für Koorinaten Paare der y-Achse
	 */
	private ArrayList<SimplePair> yPointPairs = new ArrayList<SimplePair>();


	@Override
	/**
	 * Main Methode des Plugins startet das einlesen der Punkte die Berechung der entzerrung und das Zeichnen des entzerrten Bildes
	 */
	public void run(ImageProcessor img) {
		try {
			readData();
			drawTargets(sourcePicture.getProcessor(), "SourceImage", xPointPairs, yPointPairs);
			// computeAffineTransformation();
			computeDrawRadialTransformation(sourcePicture, xPointPairs, yPointPairs);

		} catch (Exception exc) {
			IJ.error(exc.getMessage() + exc.getClass() + exc.getCause() + exc.getStackTrace());
		}
	}


	/**
	 * Öffnet die Auswahl eines Textfensters und lädt alle Punkt paare des Text Fensters in SimplePair Objekte, die intern gespeichert werden
	 * Zusätzlich kann hier das Optimale Gitter berechnet werden um in UwrapJ die Zielpunkte auswählen zu können
	 */
	private void readData() {
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
			return;
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

					SimplePair xPair = new SimplePair(Integer.parseInt(numbers[1].trim()), // x
																							// source
							Integer.parseInt(numbers[3].trim()), // x' target
							Math.sqrt((Integer.parseInt(numbers[1].trim()) - xCenter)
									* (Integer.parseInt(numbers[1].trim()) - xCenter)
									+ (Integer.parseInt(numbers[2].trim()) - yCenter)
											* (Integer.parseInt(numbers[2].trim()) - yCenter)),
							Integer.parseInt(numbers[0].trim()));// Index );

					SimplePair yPair = new SimplePair(Integer.parseInt(numbers[2].trim()), // y
																							// source
							Integer.parseInt(numbers[4].trim()), // y' target
							Math.sqrt((Integer.parseInt(numbers[1].trim()) - xCenter)
									* (Integer.parseInt(numbers[1].trim()) - xCenter)
									+ (Integer.parseInt(numbers[2].trim()) - yCenter)
											* (Integer.parseInt(numbers[2].trim()) - yCenter)),
							Integer.parseInt(numbers[0].trim()));// Index );
					//
					// //berechne x_target und y_target und gebe sie aus - NUR
					// zur einmaligen generierung der Pointpairs:
					// int colid = (int)(xPair.index / nRow); // 0 - 8 reihe
					// int rowid = (int)(yPair.index - colid * nRow) ; // 0- 15
					// spalte
					//
					// int x_offset = xCenter - nXCross2Corner * distCross;
					// //koorinate mittelpunkt gitter - anzahl der gitterpunkte
					// nach links
					// int y_offset = yCenter - nYCross2Corner * distCross;
					// //koorinate mittelpunkt gitter - anzahl der gitterpunkte
					// nach oben
					//
					// xPair.target = colid * distCross+ x_offset;
					// yPair.target = rowid * distCross + y_offset;
					//
					// debug.drawOval((int)xPair.target, (int)yPair.target, 3,
					// 3);

					xPointPairs.add(xPair);
					yPointPairs.add(yPair);

				}

			}
		}

	}
	
	/**
	 * Startet die Berechnung der Radialen entzerrung und zeichnet das Ergebnis
	 * @param sourcePicture Verzerrtes IJ Bild
	 * @param xPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 * @param yPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 */
	public static void computeDrawRadialTransformation(ImagePlus sourcePicture, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs) {
		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);

		double xCenter = sourcePicture.getProcessor().getWidth() / 2;
		double yCenter = sourcePicture.getProcessor().getHeight() / 2;
		
		// radiale entzerrung berechnen
		double[] y_koeff = compute_radial_dist_koeff(yPointPairs);
		double[] x_koeff = compute_radial_dist_koeff(xPointPairs);

		// Pixel Werte für neues Bild berechnen nach dem "target to source"
		// Verfahren
		for (int y_target = 0; y_target < sourcePicture.getHeight(); y_target++) {
			for (int x_target = 0; x_target < sourcePicture.getWidth(); x_target++) {
				// x_target / (1+ a*r^2 + b*r^4 * c*r^6) = x_distorted(source)

				double radius2Center = computeRadius2Center(x_target, y_target, xCenter, yCenter);

				double x_distorted = (1. / (1. + x_koeff[0] * Math.pow(radius2Center, 2.00)
						+ x_koeff[1] * Math.pow(radius2Center, 4.00) + x_koeff[2] * Math.pow(radius2Center, 6.00))
						* x_target);

				double y_distorted = (1. / (1. + y_koeff[0] * Math.pow(radius2Center, 2.00)
						+ y_koeff[1] * Math.pow(radius2Center, 4.00) + y_koeff[2] * Math.pow(radius2Center, 6.00))
						* y_target);

				sourcePicture.getProcessor().setInterpolationMethod(sourcePicture.getProcessor().BILINEAR);
				if (x_distorted < sourcePicture.getWidth() && y_distorted < sourcePicture.getHeight()) {
					targetImg.putPixel(x_target, y_target, (int) Math
							.round(sourcePicture.getProcessor().getInterpolatedPixel(x_distorted, y_distorted)));
				}
			}
		}
		drawTargets(targetImg, "Radial", xPointPairs, yPointPairs);

		// Punkte berechnen nach radialer entzerrung zur weiteren bearbeitung (optional):
		for (int i = 0; i < xPointPairs.size(); i++) {
			double radius2Center = computeRadius2Center(xPointPairs.get(i).source, yPointPairs.get(i).source, xCenter, yCenter);

			xPointPairs.get(i).source = (1. / (1. + x_koeff[0] * Math.pow(radius2Center, 2.00)
					+ x_koeff[1] * Math.pow(radius2Center, 4.00) + x_koeff[2] * Math.pow(radius2Center, 6.00))
					* xPointPairs.get(i).source);

			yPointPairs.get(i).source = (1. / (1. + y_koeff[0] * Math.pow(radius2Center, 2.00)
					+ y_koeff[1] * Math.pow(radius2Center, 4.00) + y_koeff[2] * Math.pow(radius2Center, 6.00))
					* yPointPairs.get(i).source);

		}

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
		res.setColor(Color.GRAY);
		res.setLineWidth(3);
		for (int i = 0; i < xPointPairs.size(); i++) {

			res.drawOval((int) xPointPairs.get(i).target, (int) yPointPairs.get(i).target, 3, 3);
		}

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
	public static RealVector computeDrawAffineTransformation(ImagePlus sourcePicture, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs) {

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
			c[i * 2] = xPointPairs.get(i).target;// x_i'
			c[i * 2 + 1] = yPointPairs.get(i).target; // y_i'

			B[i * 2][0] = xPointPairs.get(i).source; // x_i
			B[i * 2][1] = yPointPairs.get(i).source; // y_i
			B[i * 2][2] = 1;
			B[i * 2][3] = 0;
			B[i * 2][4] = 0;
			B[i * 2][5] = 0;

			B[i * 2 + 1][0] = 0;
			B[i * 2 + 1][1] = 0;
			B[i * 2 + 1][2] = 0;
			B[i * 2 + 1][3] = xPointPairs.get(i).source; // x_i
			B[i * 2 + 1][4] = yPointPairs.get(i).source; // y_i
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
		ImageProcessor targetImag = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);

		// berechnen Sie es nach dem Target-to-Source-Verfahren,
		// d.h.berechnen Sie zu jedem Punkt im neuen Bild durch die
		// Umkehrabbildung die Koordinaten im
		// Eingangsbild und holen Sie von dort den Wert des Pixels. Eine
		// Interpolation wird nicht verlangt.

		RealMatrix Bi_mat = MatrixUtils.createRealMatrix(verzerrungsmatrix);
		DecompositionSolver Bi_mat_solver = new LUDecomposition(Bi_mat).getSolver();

		// zielbild durchlaufen und werde aus vorlagebild abrufen
		for (int y = 0; y < targetImag.getHeight(); y++) {
			for (int x = 0; x < targetImag.getWidth(); x++) {
				double[] t_h = { x, y, 1 };
				RealVector t_vec = new ArrayRealVector(t_h, false);
				RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); // Transformation
																					// berechnen

				int x_coord_vorlage = (int) Math.round(coord_vec.getEntry(0));
				int y_coord_vorlage = (int) Math.round(coord_vec.getEntry(1));

				if (x_coord_vorlage < targetImag.getWidth() && y_coord_vorlage < targetImag.getHeight()) {
					targetImag.putPixel(x, y, (int) Math.round(
							sourcePicture.getProcessor().getInterpolatedPixel(x_coord_vorlage, y_coord_vorlage)));
				} else // mit schwarz auffüllen
				{
					targetImag.putPixel(x, y, 0);
				}
			}
		}

		drawTargets(targetImag, "Affine", xPointPairs, yPointPairs);

		//Tranformation auf Punkt paare anwenden zur weiteren verwendung
		/*for (int i = 0; i < xPointPairs.size(); i++) {
			double[] t_h = { xPointPairs.get(i).source, yPointPairs.get(i).source, 1 };
			RealVector t_vec = new ArrayRealVector(t_h, false);
			RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); // Transformation
																				// berechnen
			xPointPairs.get(i).source = coord_vec.getEntry(0);
			yPointPairs.get(i).source = coord_vec.getEntry(1);
		}*/
		
		return p;

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
		double[] newStart = { .001, .001, .001 };
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
		this.sourcePicture = arg1;
		debug = sourcePicture.getProcessor().duplicate();
		debug.setLineWidth(2);
		debug.setColor(Color.WHITE);
		
		xCenter = sourcePicture.getProcessor().getWidth() / 2;
		yCenter = sourcePicture.getProcessor().getHeight() / 2;
		return DOES_ALL;
	}
}
