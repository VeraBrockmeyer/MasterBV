import java.awt.Color;
import java.util.ArrayList;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
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
 * IJ Plugin fuer die Berechnung der Radialen Verzerrung anhand von Pixel-Koordinaten paaren
 * @author Vera Brockmeyer
 * @author Artjom Schwabski
 *
 */
public class point_grid_radial_affin_distor_ implements PlugInFilter {
	
	/**
	 * Speicher fuer das Vorlagenbild zum Abrufen der Pixel Werte zur Zeichnung der entzerrung
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
	
	
	/**
	 * Anzahl der Gitter-Reihen zur Berechnung des optimalen Gitters
	 */
	private static final int XgridCenter = 1084;
	
	/**
	 * Anzahl der Gitter-Reihen zur Berechnung des optimalen Gitters
	 */
	private static final int YgridCenter = 714;
		
	/**
	 * Speicher fuer Punkt Paare
	 */
	private ArrayList<PointPair> pointPairs = new ArrayList<PointPair>();
	
	
	


	@Override
	/**
	 * Main Methode des Plugins startet das einlesen der Punkte die Berechung der entzerrung und das Zeichnen des entzerrten Bildes
	 */
	public void run(ImageProcessor img) {
		try {
			
			//compute_draw_optimal_Grid(distPicture.getProcessor()); //für den ersten durchlauf zum einzeichnen des optimalen gitters
			
			this.pointPairs = readData();
			
			drawPointPairs(distPicture.getProcessor(), "SourceImage with PointPairs", pointPairs);
			//distPicture=computeDrawProjectiveTransformation(distPicture, pointPairs);
			computeDrawRadialTransformation(distPicture, pointPairs);

		} catch (Exception exc) {
			IJ.error(exc.getMessage() + exc.getClass() + exc.getCause() + exc.getStackTrace());
		}
	}


	/**
	 * oeffnet die Auswahl eines Textfensters und laedt alle Punkt-paare des Text Fensters in PointPair Objekte
	 */
	private ArrayList<PointPair> readData() {

		ArrayList<PointPair> pointPairs = new ArrayList<PointPair>();
		
		//Fuer ein optimal zentrierte Gitterabbildung kann die Bildmitte berechnet werden. Ansonnsten hier die koordinaten des Gittermittelpunktes im Bild angeben
		int xCenter = distPicture.getProcessor().getWidth() / 2;
		int yCenter = distPicture.getProcessor().getHeight() / 2; 
	
		//ImageProcessor ip = distPicture.getProcessor().duplicate(); //erzeuge kopie des Ausgangsbildes zum Anzeigen der eingelesenen Punkt Paare oder des Radius

		// Textfenster mit Punktpaare Textdtei waehlen:
		java.awt.Window[] non_img_windows = WindowManager.getAllNonImageWindows();
		GenericDialog gd = new GenericDialog("Textdatei auswaehlen:");
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

				// speicher fuer punkt paarungen

				for (int j = 1; j < lines.length; j++) // erste zeile
														// ueberspringen
				{
					String[] numbers = lines[j].split("\t"); // spalte zeile
																// anhand des
																// tab zeichens

					//Erzeuge Punkt bezogen auf Mittelpunkt:
					PointPair pp = new PointPair(
							Integer.parseInt(numbers[1].trim()), // x dist
							Integer.parseInt(numbers[2].trim()), // y dist
							Integer.parseInt(numbers[3].trim()), // x' undist
							Integer.parseInt(numbers[4].trim()), // y' undist
							xCenter,
							yCenter,
							Integer.parseInt(numbers[0].trim()));// Index
					
					//draw undistorted cordinates for debug check:
					//ip.drawString((pp.x_dist + "/ " + pp.y_undist), (int)pp.x_dist+xCenter,(int) pp.y_dist+yCenter);
					
					//draw radius for debug check:
					//ip.drawString( "" + pp.r, (int)pp.x_dist+xCenter,(int) pp.y_dist+yCenter); 
					pointPairs.add(pp);

				}

			}
		}
		
		//Ausgabe der Bildpunkte oder des Radius optional fürs debugging:
		//ImagePlus img = new ImagePlus("Debug Points", ip);
		//img.show();
		
		return pointPairs;
			
	}
	
	/**
	 * Zeichnet das optimale Gitter anhand der Klassenvariablen-Parametern in das übergeben Bild
	 * @param ip Imageprocessor des Bildes
	 */
	private void compute_draw_optimal_Grid(ImageProcessor ip)
	{
		int xCenter = distPicture.getProcessor().getWidth() / 2;
		int yCenter = distPicture.getProcessor().getHeight() / 2; 
		
		 // zur einmaligen generierung der ziel koordinaten gitters:
		for(int colid = 0; colid <= nCol; colid++)
		{
			for(int rowid = 0; rowid <= nRow; rowid++)
			{
				int x_offset = XgridCenter - nXCross2Corner * distCross;				 
				int y_offset = YgridCenter - nYCross2Corner * distCross;

				double x_undist = colid * distCross+ x_offset;
				double y_undist = rowid * distCross + y_offset;
				
				ip.setColor(Color.WHITE);
				ip.drawOval((int)x_undist, (int)y_undist, 5, 5);
			}
		}	 
	}
	
	//TODO Ursprung Verschieben!!!!
	public static ImagePlus computeDrawProjectiveTransformation(ImagePlus sourcePicture, ArrayList<PointPair> pointPairs){
		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
		int xCenter = targetImg.getWidth() / 2;
		int yCenter = targetImg.getHeight() / 2; 
		double[][]B = new double[pointPairs.size()*2][9];
		int counter = 0;
		for (int i = 0; i < B.length; i+=2) {
			double x_dist = pointPairs.get(counter).x_dist;
			double y_dist = pointPairs.get(counter).y_dist;
			
			double x_undist = pointPairs.get(counter).x_undist;
			double y_undist = pointPairs.get(counter).y_undist;
			
			B[i][0] = -x_dist;
			B[i][1] = -y_dist;
			B[i][2] = -1.;
			B[i][3] = 0.;
			B[i][4] = 0.;
			B[i][5] = 0.;
			B[i][6] = x_dist*x_undist;
			B[i][7] = y_dist*x_undist;
			B[i][8] = x_undist;

			B[i+1][0] = 0.;
			B[i+1][1] = 0.;
			B[i+1][2] = 0.;
			B[i+1][3] = -x_dist;
			B[i+1][4] = -y_dist;
			B[i+1][5] = -1.;
			B[i+1][6] = x_dist*y_undist;
			B[i+1][7] = y_dist*y_undist;
			B[i+1][8] = y_undist;
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
		RealMatrix pMatInv = MatrixUtils.inverse(pMat);
		for (int y = 0; y < targetImg.getHeight(); y++) {
			for (int x = 0; x < targetImg.getWidth(); x++) {
				
				double[] target  = {x-xCenter,y-yCenter,1};
				RealVector t_vec = new ArrayRealVector(target, false);
				RealVector coord_vec = pMat.operate(t_vec);
				double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
				double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
			
				x_coord_vorlage+=xCenter;
				y_coord_vorlage+=yCenter;

				if (x_coord_vorlage>=0 && x_coord_vorlage < targetImg.getWidth() && y_coord_vorlage>=0 && y_coord_vorlage < targetImg.getHeight()) {
					targetImg.putPixel(x, y, (int) Math.round(sourcePicture.getProcessor().getInterpolatedPixel(x_coord_vorlage , y_coord_vorlage)));
				} 
				else
				{
					targetImg.putPixel(x, y, 0);
				}
			}
		}
		
		for (int i = 0; i < pointPairs.size(); i++) {
			double[] target  = {pointPairs.get(i).x_dist,pointPairs.get(i).y_dist,1};
			RealVector t_vec = new ArrayRealVector(target, false);
			RealVector coord_vec = pMatInv.operate(t_vec);
			double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
			double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
			pointPairs.get(i).x_dist = x_coord_vorlage;
			pointPairs.get(i).y_dist = y_coord_vorlage;
			
		}

		drawPointPairs(targetImg, "Projective Transformation",pointPairs);
		return new ImagePlus("Projected Image",targetImg);
	}
	
	
	
	/**
	 * Startet die Berechnung der Radialen entzerrung und zeichnet das Ergebnis
	 * @param distPicture Verzerrtes IJ Bild
	 * @param pointPairs Punkt-Paare fuer Start und Ziel Koordinaten 
	 * 
	 */
	public static void computeDrawRadialTransformation(ImagePlus distPicture, ArrayList<PointPair> pointPairs) 
	{
		ShortProcessor undistImg = new ShortProcessor(distPicture.getWidth(), distPicture.getHeight(), true);

		//Gittermittelpunkt bestimmen fuer ideale Gitterabbildung kann dieser errechnet werden ansonsten hier die koordinaten eingeben:
		double xCenter = distPicture.getProcessor().getWidth() / 2;
		double yCenter = distPicture.getProcessor().getHeight() / 2;
		
		// radiale entzerrung berechnen
		double[] koeff = compute_radial_dist_koeff(pointPairs);

		// Pixel Werte fuer neues Bild berechnen nach dem "target to source"
		// Verfahren
		for (int yImg = 0; yImg < distPicture.getHeight(); yImg++) {
			for (int xImg = 0; xImg < distPicture.getWidth(); xImg++) {
				
				//Erzeuge Punkt mit Koordinatentransformation zum Gittermittelpunkt:
				PointPair pp = new PointPair(xImg, yImg, xCenter, yCenter);

				//radiale Entzerrung mit den vorgebenen Koeffizienten: 
				//x_distorted * (1+ a*r^2 + b*r^4 * c*r^6) = x_undistorted 
				
				pp.x_dist = pp.x_undist / (1. + koeff[0] * pp.r + koeff[1] * pp.r*pp.r + koeff[2] * pp.r*pp.r*pp.r);
				pp.y_dist = pp.y_undist / (1. + koeff[0] * pp.r + koeff[1] * pp.r*pp.r + koeff[2] * pp.r*pp.r*pp.r);
				
				//Koordinatenruecktransformation in ImageJ Koordinaten:
				pp.x_dist = pp.x_dist + xCenter;
				pp.y_dist = pp.y_dist + yCenter;
				
				//Entzerrtes Bild zeichnen:
				if (pp.x_dist < distPicture.getWidth() && pp.y_dist < distPicture.getHeight()) 
				{
					undistImg.putPixel(xImg, yImg, distPicture.getProcessor().getPixel((int)Math.round(pp.x_dist), (int)Math.round(pp.y_dist)));
				}
				else
				{
					undistImg.putPixel(xImg, yImg,255);	//füllen mit weiß			
				}
			}
		}

		// Punkte berechnen nach radialer entzerrung zur weiteren bearbeitung (optional):
		for (int i = 0; i < pointPairs.size(); i++) 
		{
			PointPair pp = 	pointPairs.get(i);
			
			pp.x_dist = (1. + koeff[0] * pp.r + koeff[1] * pp.r*pp.r + koeff[2] * pp.r*pp.r*pp.r) * pp.x_dist;
			pp.y_dist = (1. + koeff[0] * pp.r + koeff[1] * pp.r*pp.r + koeff[2] * pp.r*pp.r*pp.r) * pp.y_dist;

		}
		
		drawPointPairs(undistImg, "Radial", pointPairs);

	}

	/**
	 * Zeichnet Punkte an die Stellen der Zielkoorinaten in das uebergebene Bild und bringt es auf den Bildschirm
	 * @param ip Bild in das die Punkte gezeichnet werden
	 * @param s Name des Bildes
	 * @param pointPairs  Zu zeichnende Start und Ziel Koordinaten Punkte
	 */
	public static void drawPointPairs(ImageProcessor ip, String s, ArrayList<PointPair> pointPairs) {
		// Punkte in radial entzerrtes Bild malen
		ImageProcessor res = ip.duplicate();
		res.set(255);
		double xCenter = ip.getWidth() / 2;
		double yCenter = ip.getHeight() / 2;

		// Zeichne ziel punkte:
		res.setLineWidth(3);
		for (int i = 0; i < pointPairs.size(); i++) {
			res.setColor(Color.GRAY);
			res.setLineWidth(2);
			res.drawLine((int) (pointPairs.get(i).x_dist+xCenter), (int) (pointPairs.get(i).y_dist+yCenter), (int) (pointPairs.get(i).x_undist+xCenter), (int) (pointPairs.get(i).y_undist+yCenter));
//			res.drawOval((int) (pointPairs.get(i).x_dist+xCenter), (int) (pointPairs.get(i).y_dist+yCenter), 2, 2);
//			res.setColor(Color.WHITE);
//			res.drawOval((int) (pointPairs.get(i).x_undist+xCenter), (int) (pointPairs.get(i).y_undist+yCenter), 2, 2);
		}
		ImageProcessor res2 = ip.duplicate();
		ImagePlus resImg2 = new ImagePlus(s, res2);
		resImg2.show();
		ImagePlus resImg = new ImagePlus(s+"_onlydiffs", res);
		resImg.show();

	}

	
	/**
	 * 
	 * @param punkt_paare Array mit PointPair Objekten in denen die Vorlage- und Ziel-Pixelkoordinaten gespeichert sind
	 * @return Koefizienten der Radialen-Verzerrung nach Levenberg-Marquadt
	 */
	public static double[] compute_radial_dist_koeff(ArrayList<PointPair> punkt_paare) {
		RadialDistFunction qf = new RadialDistFunction(punkt_paare);
		LeastSquaresBuilder lsb = new LeastSquaresBuilder();

		// set model function and its jacobian
		lsb.model(qf.retMVF(), qf.retMMF());
		// set target data
		double[] undistPoints = qf.realTargetPoints();
		lsb.target(undistPoints);

		double[] newStart = { 0.,0.,0. };

		// set initial parameters
		lsb.start(newStart);
		// set upper limit of evaluation time
		lsb.maxEvaluations(1000);
		// set upper limit of iteration time
		lsb.maxIterations(500);

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
			IJ.log("RMS: "+lsoo.getRMS());
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