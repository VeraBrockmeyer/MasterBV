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
			drawTargets(distPicture.getProcessor(), "SourceImage",PointPairs);
			computeDrawRadialTransformation(distPicture, PointPairs);
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
					
					double x_dist = Integer.parseInt(numbers[1].trim())-xCenter;
					double y_dist = Integer.parseInt(numbers[2].trim())-yCenter;
					
					double x_undist = Integer.parseInt(numbers[3].trim())-xCenter;
					double y_undist = Integer.parseInt(numbers[4].trim())-yCenter;
					
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


//	public static ImagePlus computeDrawPerpectiveTransformation(ImagePlus sourcePicture, List<SimplePair> xPointPairs, List<SimplePair> yPointPairs){
//		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
//	
//		double[][]B = new double[xPointPairs.size()*2][9];
//		int counter = 0;
//		for (int i = 0; i < B.length; i+=2) {
//			
//			B[i][0] = -xPointPairs.get(counter).undistorted;
//			B[i][1] = -yPointPairs.get(counter).undistorted;
//			B[i][2] = -1.;
//			B[i][3] = 0.;
//			B[i][4] = 0.;
//			B[i][5] = 0.;
//			B[i][6] = xPointPairs.get(counter).undistorted*xPointPairs.get(counter).distorted;
//			B[i][7] = yPointPairs.get(counter).undistorted*xPointPairs.get(counter).distorted;
//			B[i][8] = xPointPairs.get(counter).distorted;
//			B[i+1][0] = 0.;
//			B[i+1][1] = 0.;
//			B[i+1][2] = 0.;
//			B[i+1][3] = -xPointPairs.get(counter).undistorted;
//			B[i+1][4] = -yPointPairs.get(counter).undistorted;
//			B[i+1][5] = -1.;
//			B[i+1][6] = xPointPairs.get(counter).undistorted*yPointPairs.get(counter).distorted;
//			B[i+1][7] = yPointPairs.get(counter).undistorted*yPointPairs.get(counter).distorted;
//			B[i+1][8] = yPointPairs.get(counter).distorted;
//			counter++;
//		}
//		RealMatrix A = new Array2DRowRealMatrix(B);
//		RealMatrix At = A.transpose();
//		RealMatrix AtA = At.multiply(A);
//
//		EigenDecomposition ed = new EigenDecomposition(AtA);
//		RealMatrix V = ed.getV();
//	
//		double[][] pArray = new double[3][3];
//		
//		pArray[0][0]=V.getEntry(0,V.getRowDimension()-1);
//		pArray[0][1]=V.getEntry(1,V.getRowDimension()-1);
//		pArray[0][2]=V.getEntry(2,V.getRowDimension()-1);
//		
//		pArray[1][0]=V.getEntry(3,V.getRowDimension()-1);
//		pArray[1][1]=V.getEntry(4,V.getRowDimension()-1);
//		pArray[1][2]=V.getEntry(5,V.getRowDimension()-1);
//		
//		pArray[2][0]=V.getEntry(6,V.getRowDimension()-1);
//		pArray[2][1]=V.getEntry(7,V.getRowDimension()-1);
//		pArray[2][2]=V.getEntry(8,V.getRowDimension()-1);		
//	
//		RealMatrix pMat = new Array2DRowRealMatrix(pArray);
//	
//		for (int y = 0; y < targetImg.getHeight(); y++) {
//			for (int x = 0; x < targetImg.getWidth(); x++) {
//				double[] target  = {x,y,1};
//				RealVector t_vec = new ArrayRealVector(target, false);
//				RealVector coord_vec = pMat.operate(t_vec);
//				double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
//				double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
//			
//				if (x_coord_vorlage>=0 && x_coord_vorlage < targetImg.getWidth() && y_coord_vorlage>=0 && y_coord_vorlage < targetImg.getHeight()) {
//					targetImg.putPixel(x, y, (int) Math.round(sourcePicture.getProcessor().getInterpolatedPixel(x_coord_vorlage , y_coord_vorlage)));
//				} 
//				else
//				{
//					targetImg.putPixel(x, y, 0);
//				}
//			}
//		}
//		
//		for (int i = 0; i < xPointPairs.size(); i++) {
//			double[] target  = {xPointPairs.get(i).distorted,yPointPairs.get(i).distorted,1};
//			RealVector t_vec = new ArrayRealVector(target, false);
//			RealVector coord_vec = pMat.operate(t_vec);
//			double x_coord_vorlage = coord_vec.getEntry(0)/coord_vec.getEntry(2);
//			double y_coord_vorlage = coord_vec.getEntry(1)/coord_vec.getEntry(2);
//			xPointPairs.get(i).distorted = x_coord_vorlage;
//			yPointPairs.get(i).distorted = y_coord_vorlage;
//		}
//
//		drawTargets(targetImg, "Projective Transformation", xPointPairs, yPointPairs);
//		return new ImagePlus("Projected Image",targetImg);
//	}
//
//


	/**
	 * Startet die Berechnung der Radialen entzerrung und zeichnet das Ergebnis
	 * @param sourcePicture Verzerrtes IJ Bild
	 * @param xPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 * @param yPointPairs Punkt-Paare für Start und Ziel Koordinaten der x-Achse
	 */
	public static void computeDrawRadialTransformation(ImagePlus sourcePicture, ArrayList<PointPair>pairs) {
		ShortProcessor targetImg = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
		 double[] xCoeffs = new double [3];
		 double[] yCoeffs = new double [3];
		double xCenter = sourcePicture.getProcessor().getWidth() / 2;
		double yCenter = sourcePicture.getProcessor().getHeight() / 2;
		
		

		System.out.println("Starting X Coeff Estimation ....");
		
		
		 // the model function components are the distances to current estimated center,
		 // they should be as close as possible to the specified radius
		 MultivariateJacobianFunction modelX = new MultivariateJacobianFunction() {
		      public Pair<RealVector, RealMatrix> value(final RealVector coeffs) {

		          RealVector value = new ArrayRealVector(pairs.size());
		          RealMatrix jacobian = new Array2DRowRealMatrix(pairs.size(), 3);

		          for (int i = 0; i < pairs.size(); ++i) {
		        	  PointPair pp = pairs.get(i);
		              double val = (1. + coeffs.getEntry(0)* pp.r * pp.r 
		            		  + coeffs.getEntry(1) * pp.r * pp.r * pp.r * pp.r
		            		  + coeffs.getEntry(2)* pp.r *pp.r * pp.r *pp.r * pp.r * pp.r)
		            		  * pp.x_dist;
		              
		              value.setEntry(i, val);
		              
		            //0 = (1 + a*r^2 + b*r^4 + c*r^6) * x - x'
			        	//dx'/da = r^2*x
			        	//dx'/db = r^4*x
			        	//dx'/dc = r^6*x
			        	
			        	
			          double j0 = Math.pow(pp.r, 2.00) * pp.x_dist; 
			          double j1 = Math.pow(pp.r, 4.00) * pp.x_dist; 
			          double j2 = Math.pow(pp.r, 6.00) * pp.x_dist;
		              
		              jacobian.setEntry(i, 0, j0);
		              jacobian.setEntry(i, 1, j1);
		              jacobian.setEntry(i, 2, j2);
		          }

		          return new Pair<RealVector, RealMatrix>(value, jacobian);

		      }
		  };

		  double[] targetsX = new double[pairs.size()];
		  for (int i = 0; i < pairs.size(); i++) {
			targetsX[i]=pairs.get(i).x_undist;
		}

		  // least squares problem to solve : modeled radius should be close to target radius
		  LeastSquaresProblem problemX = new LeastSquaresBuilder().
		                                start(new double[] { 1.e-10, 1.e-10,1.e-10 }).
		                                model(modelX).
		                                target(targetsX).
		                                lazyEvaluation(false).
		                                maxEvaluations(1000).
		                                maxIterations(500).
		                                build();
		  LeastSquaresOptimizer.Optimum optimumX = new LevenbergMarquardtOptimizer().optimize(problemX);
		  xCoeffs[0] = optimumX.getPoint().getEntry(0);
		  xCoeffs[1] = optimumX.getPoint().getEntry(1);
		  xCoeffs[2] = optimumX.getPoint().getEntry(2);
		  System.out.println("xCoeffs: " + xCoeffs[0] + " " +  xCoeffs[1]+" "+  xCoeffs[2]);
		  System.out.println("RMS: "           + optimumX.getRMS());
		  System.out.println("evaluations: "   + optimumX.getEvaluations());
		  System.out.println("iterations: "    + optimumX.getIterations());
		  
		  
		  //-----------------------------------------------------------------------------------------------------------------------------
		  
		  System.out.println("Starting Y Coeff Estimation ....");
			
			
			 // the model function components are the distances to current estimated center,
			 // they should be as close as possible to the specified radius
			 MultivariateJacobianFunction modelY = new MultivariateJacobianFunction() {
			      public Pair<RealVector, RealMatrix> value(final RealVector coeffs) {

			          RealVector value = new ArrayRealVector(pairs.size());
			          RealMatrix jacobian = new Array2DRowRealMatrix(pairs.size(), 3);

			          for (int i = 0; i < pairs.size(); ++i) {
			        	  PointPair pp = pairs.get(i);
			              double val = (1. + coeffs.getEntry(0)* pp.r * pp.r 
			            		  + coeffs.getEntry(1) * pp.r * pp.r * pp.r * pp.r
			            		  + coeffs.getEntry(2)* pp.r *pp.r * pp.r *pp.r * pp.r * pp.r)
			            		  * pp.y_dist;
			              
			              value.setEntry(i, val);
			              
			            //0 = (1 + a*r^2 + b*r^4 + c*r^6) * x - x'
				        	//dx'/da = r^2*x
				        	//dx'/db = r^4*x
				        	//dx'/dc = r^6*x
				        	
				        	
				          double j0 = Math.pow(pp.r, 2.00) * pp.y_dist; 
				          double j1 = Math.pow(pp.r, 4.00) * pp.y_dist; 
				          double j2 = Math.pow(pp.r, 6.00) * pp.y_dist;
			              
			              jacobian.setEntry(i, 0, j0);
			              jacobian.setEntry(i, 1, j1);
			              jacobian.setEntry(i, 2, j2);
			          }

			          return new Pair<RealVector, RealMatrix>(value, jacobian);

			      }
			  };

			  double[] targetsY = new double[pairs.size()];
			  for (int i = 0; i < pairs.size(); i++) {
				targetsY[i]=pairs.get(i).y_undist;
			}

			  // least squares problem to solve : modeled radius should be close to target radius
			  LeastSquaresProblem problemY = new LeastSquaresBuilder().
			                                start(new double[] { 1.e-10, 1.e-10,1.e-10 }).
			                                model(modelY).
			                                target(targetsY).
			                                lazyEvaluation(false).
			                                maxEvaluations(1000).
			                                maxIterations(500).
			                                build();
			  LeastSquaresOptimizer.Optimum optimumY = new LevenbergMarquardtOptimizer().optimize(problemY);
			  yCoeffs[0] = optimumY.getPoint().getEntry(0);
			  yCoeffs[1] = optimumY.getPoint().getEntry(1);
			  yCoeffs[2] = optimumY.getPoint().getEntry(2);
			  System.out.println("xCoeffs: " + yCoeffs[0] + " " +  yCoeffs[1]+" "+  yCoeffs[2]);
			  System.out.println("RMS: "           + optimumY.getRMS());
			  System.out.println("evaluations: "   + optimumY.getEvaluations());
			  System.out.println("iterations: "    + optimumY.getIterations());
			  
			//-----------------------------------------------------------------------------------------------------------------------
			  
			// Pixel Werte für neues Bild berechnen nach dem "target to source"
				// Verfahren
				for (int y = 0; y < sourcePicture.getHeight(); y++) {
					for (int x = 0; x < sourcePicture.getWidth(); x++) {
						// x_target / (1+ a*r^2 + b*r^4 * c*r^6) = x_distorted(source)

						double radius2Center = computeRadius2Center(x, y, xCenter, yCenter);
						
						double x_distorted = (1./(1. + xCoeffs[0] * Math.pow(radius2Center, 2.00)
								+ xCoeffs[1] * Math.pow(radius2Center, 4.00) + xCoeffs[2] * Math.pow(radius2Center, 6.00))
								* (x-xCenter));

						double y_distorted = (1./ (1. + yCoeffs[0] * Math.pow(radius2Center, 2.00)
								+ yCoeffs[1] * Math.pow(radius2Center, 4.00) + yCoeffs[2] * Math.pow(radius2Center, 6.00))
								* (y-yCenter));
						
						x_distorted+=xCenter;
						y_distorted+=yCenter;
						
						sourcePicture.getProcessor().setInterpolationMethod(sourcePicture.getProcessor().BILINEAR);
						if (x_distorted>=0 && x_distorted < sourcePicture.getWidth() && y_distorted>=0 && y_distorted < sourcePicture.getHeight()) {
							targetImg.putPixel(x, y, (int) Math
									.round(sourcePicture.getProcessor().getInterpolatedPixel(x_distorted, y_distorted)));
						}
						else{
							targetImg.putPixel(x, y,255);
						}
					}
				}
				
				// Punkte berechnen nach radialer entzerrung zur weiteren bearbeitung (optional):
				for (int i = 0; i < pairs.size(); i++) {
					PointPair pp = pairs.get(i);

					pp.x_dist = (1. / (1. + xCoeffs[0] * Math.pow(pp.r, 2.00)
							+ xCoeffs[1] * Math.pow(pp.r, 4.00) + xCoeffs[2] * Math.pow(pp.r, 6.00))
							* pp.x_dist);
					

					pp.y_dist = (1. / (1. + yCoeffs[0] * Math.pow(pp.r, 2.00)
							+ yCoeffs[1] * Math.pow(pp.r, 4.00) + yCoeffs[2] * Math.pow(pp.r, 6.00))
							* pp.y_dist);

				}
				
	  drawTargets(targetImg, "Radial", pairs);
	  
	}

	/**
	 * Zeichnet Punkte an die Stellen der Zielkoorinaten in das übergebene Bild und bringt es auf den Bildschirm
	 * @param ip Bild in das die Punkte gezeichnet werden
	 * @param s Name des Bildes
	 * @param xPointPairs x-Koordinaten 
	 * @param yPointPairs y-Koordinaten
	 */
	public static void drawTargets(ImageProcessor ip, String s, List<PointPair> pairs) {
		// Punkte in radial entzerrtes Bild malen
		ImageProcessor res = ip.duplicate();
		int xCenter = ip.getWidth() / 2;
		int yCenter = ip.getHeight() / 2;
		// Zeichne ziel punkte:
		res.setLineWidth(3);
		for (int i = 0; i < pairs.size(); i++) {
			PointPair pp = pairs.get(i);
			res.setColor(Color.BLACK);
			res.drawOval((int) pp.x_undist+xCenter, (int) pp.y_undist-yCenter, 3, 3);
			res.setColor(Color.GRAY);
			res.drawOval((int) pp.x_dist+xCenter, (int) pp.y_dist-yCenter, 3, 3);
		}

		ImagePlus resImg = new ImagePlus(s, res);
		resImg.show();
		
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

	

	@Override
	/**
	 * Speichern des Voralge Bildes bei Aufruf des Plugins und erzeugung einer Kopie
	 */
	public int setup(String arg0, ImagePlus arg1) {
		this.distPicture = arg1;
		return DOES_ALL;
	}
}
