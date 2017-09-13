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

public class point_grid_radial_affin_distor_Vera_2_ implements PlugInFilter
{
	private ImagePlus sourcePicture;
	private int nCol = 17;//19;
	private int nRow = 11;
	private int xCenter = 971;//1084;
	private int yCenter = 651;//713;
	private int distCross =111;
	private int nXCross2Corner=8;
	private int nYCross2Corner=5;
	private ImagePlus debugImg;
	private ImageProcessor debug;

	
	@SuppressWarnings("deprecation")
	@Override
	public void run(ImageProcessor img) 
	{
		try
		{
			//Textfenster mit Punktpaare Textdtei w�hlen:
			java.awt.Window[] non_img_windows = WindowManager.getAllNonImageWindows();
			GenericDialog gd = new GenericDialog("Textdatei ausw�hlen:");
			String[] Windows_names = new String[non_img_windows.length];
			for(int i=0; i < non_img_windows.length; i++)
			{
				if(non_img_windows[i] instanceof ij.plugin.frame.Editor)
				{
					Windows_names[i] = ((ij.plugin.frame.Editor)non_img_windows[i]).getTitle();
				}
				
			}
			
			String choise = "";
		  	gd.addChoice("Textdatei Fenster", Windows_names, choise);
			
		  	gd.showDialog();
			
			//abbrechen wenn keine eingabe:
			if (gd.wasCanceled()) return;
			choise = gd.getNextChoice();
			
			//Punktpaare einlesen:
			List<PointPair> PointPairs = new ArrayList<PointPair>();
			ArrayList<SimplePair> xPointPairs = new ArrayList<SimplePair>();
			ArrayList<SimplePair> yPointPairs = new ArrayList<SimplePair>();
			
		  	for(int i=0; i < non_img_windows.length; i++)
			{
				if(((java.awt.Frame)non_img_windows[i]).getTitle() == choise)
				{
	
					String text = ((ij.plugin.frame.Editor)non_img_windows[i]).getText();
						
					String[] lines = text.split("\n");
					
					 //speicher f�r punkt paarungen
					
					for(int j = 1; j < lines.length; j++) //erste zeile �berspringen
					{
						String[] numbers = lines[j].split("\t"); //spalte zeile anhand des tab zeichens
						
						SimplePair xPair = new SimplePair(	
								Integer.parseInt(numbers[1].trim()), //x source
								Integer.parseInt(numbers[3].trim()), //x' target 
								Math.sqrt( (Integer.parseInt(numbers[1].trim()) - xCenter) * (Integer.parseInt(numbers[1].trim()) - xCenter)
										+ (Integer.parseInt(numbers[2].trim()) - yCenter) * (Integer.parseInt(numbers[2].trim()) - yCenter)),
								Integer.parseInt(numbers[0].trim()));//Index	);
								
						SimplePair yPair = new SimplePair(	
								Integer.parseInt(numbers[2].trim()), //y source
								Integer.parseInt(numbers[4].trim()), //y' target 
								Math.sqrt( (Integer.parseInt(numbers[1].trim()) - xCenter) * (Integer.parseInt(numbers[1].trim()) - xCenter)
										+ (Integer.parseInt(numbers[2].trim()) - yCenter) * (Integer.parseInt(numbers[2].trim()) - yCenter)),
								Integer.parseInt(numbers[0].trim()));//Index	);
						
						
						PointPair pair = new PointPair(xPair.source, yPair.source, xPair.target, yPair.target, yPair.radius, yPair.index);
						
						
						/*
						//berechne x_target und y_target anstatt es aus der datei auszulesen:
						int colid = (int)(xPair.index / nRow); // 0 - 8 reihe
						int rowid = (int)(yPair.index - colid * nRow)  ; // 0- 15 spalte
				
						int x_offset = xCenter - nXCross2Corner * distCross; //koorinate mittelpunkt gitter - anzahl der gitterpunkte nach links
						int y_offset = yCenter - nYCross2Corner * distCross; //koorinate mittelpunkt gitter - anzahl der gitterpunkte nach oben
						
						xPair.target = colid * distCross+ x_offset;
						yPair.target = rowid * distCross + y_offset;
						*/
						
						xPointPairs.add(xPair);
						yPointPairs.add(yPair);
						PointPairs.add(pair);
//						debug.drawString(""+xPair.index, (int)xPair.target, (int)yPair.target);
//						debug.setColor(Color.GRAY);
//						debug.drawString(""+xPair.index, (int)xPair.source, (int)yPair.source);
//						debug.setColor(Color.WHITE);
//						debug.drawOval((int)xPair.target, (int)yPair.target, 3, 3);
					}
					
				}
			}		
		  	    
//   			ImagePlus debugImg = new ImagePlus("Debug",debug);
//  			debugImg.show();
		  	
			//radiale entzerrung berechnen
			//double[] koeff = entzerungskoeefizienten_radial_berechnen(PointPairs);
			double[] y_koeff= entzerungskoeefizienten_radial_berechnen_simple(yPointPairs);
			double[] x_koeff= entzerungskoeefizienten_radial_berechnen_simple(xPointPairs);
						
			//erzeugen Sie ein neues Bild gleicher Gr��e wie das Eingabebild 
			ShortProcessor sp = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
			
			//Pixel Werte f�r neues Bild berechnen nach dem "target to source" Verfahren
			for (int y_target= 0; y_target < sourcePicture.getHeight();y_target++) 
			{
				for (int x_target= 0; x_target < sourcePicture.getWidth(); x_target++)
				{
					// x_target / (1+ a*r^2 + b*r^4 * c*r^6) = x_distorted(source)

					double radius2Center =computeRadius2Center(x_target, y_target);
					
					
//					double x_distorted =  (1./ (1. 
//		        			+koeff[0] * Math.pow(radius2Center, 2.00) 
//		        			+koeff[1] * Math.pow(radius2Center, 4.00)
//		        			+koeff[2] * Math.pow(radius2Center, 6.00)
//		        			) * x_target);
//					
//					double y_distorted = 
//					(1./(1. 
//		        			+koeff[0] * Math.pow(radius2Center, 2.00) 
//		        			+koeff[1] * Math.pow(radius2Center, 4.00)
//		        			+koeff[2] * Math.pow(radius2Center, 6.00)
//		        			) * y_target);							
					
					double x_distorted =  (1./ (1. 
        			+x_koeff[0] * Math.pow(radius2Center, 2.00) 
        			+x_koeff[1] * Math.pow(radius2Center, 4.00)
        			+x_koeff[2] * Math.pow(radius2Center, 6.00)
        			) * x_target);
			
					double y_distorted = (1./(1. 
        			+y_koeff[0] * Math.pow(radius2Center, 2.00) 
        			+y_koeff[1] * Math.pow(radius2Center, 4.00)
        			+y_koeff[2] * Math.pow(radius2Center, 6.00)
        			) * y_target);	
					
					sourcePicture.getProcessor().setInterpolationMethod(sourcePicture.getProcessor().BILINEAR);
					if(x_distorted < sourcePicture.getWidth() && y_distorted < sourcePicture.getHeight())
					{
						sp.putPixel(x_target, y_target, (int)Math.round(sourcePicture.getProcessor().getInterpolatedPixel(x_distorted, y_distorted)));
					}
				}	
			}
			
			//Zeichne ziel punkte:
//			sp.setColor(Color.WHITE);
//			sp.setLineWidth(3);
//			for (PointPair p : PointPairs) 
//			{
//				sp.drawOval((int)p.x_target, (int)p.y_target, 3, 3);
//
//			}
			
			//Punkte in radial entzerrtes Bild malen
			ImageProcessor resRad = sp.duplicate();
			 //Zeichne ziel punkte:
			resRad.setColor(Color.BLACK);
			resRad.setLineWidth(3);
			for (PointPair point : PointPairs) 
			{
				resRad.drawOval((int)point.x_target, (int)point.y_target, 3, 3);

			}
			ImagePlus radRes = new ImagePlus("Radial Result", resRad);
			radRes.show();
			
			//Am Ende m�ssen Sie das neue Bild noch auf den Bildschirm bringen.
			ImagePlus newImg = new ImagePlus();
			newImg = new ImagePlus("result", sp);
			//
			//Affine entzerrung berechnen
			//
			
			//Punkt berechnen nach radialer entzerrung berechnen:
			for (PointPair p : PointPairs) 
			{
				double radius2Center = computeRadius2Center(p.x_source, p.y_source);
				
				p.x_source =  (1./ (1. 
	        			+x_koeff[0] * Math.pow(radius2Center, 2.00) 
	        			+x_koeff[1] * Math.pow(radius2Center, 4.00)
	        			+x_koeff[2] * Math.pow(radius2Center, 6.00)
	        			) * p.x_source);
				
				p.y_source = (1./(1. 
    			+y_koeff[0] * Math.pow(radius2Center, 2.00) 
    			+y_koeff[1] * Math.pow(radius2Center, 4.00)
    			+y_koeff[2] * Math.pow(radius2Center, 6.00)
    			) * p.y_source);	
						
			}
			
			//Bild anhand der PunktPaare affin entzerren
			entzerrungsmatrix_affin_berechnen(PointPairs, newImg);
			
		} 
		catch(Exception exc)
		{
			IJ.error(exc.getMessage() + exc.getClass() + exc.getCause() + exc.getStackTrace());
		}			
	}
	
	private double computeRadius2Center(double x_target, double y_target) {
	return	Math.sqrt((x_target - xCenter) * (x_target - xCenter)
				+ (y_target - yCenter) * (y_target - yCenter));
		
	}

//	@SuppressWarnings("deprecation")
//	private double[] entzerungskoeefizienten_radial_berechnen(List<PointPair> punkt_paare)
//	{
//		RadialDistFunction qf = new RadialDistFunction(punkt_paare);
//		LeastSquaresBuilder lsb = new LeastSquaresBuilder();
//
//		//set model function and its jacobian
//		lsb.model(qf.retMVF(), qf.retMMF());
//		double[] newTarget = qf.realTargetPoints();
//		
//		//set target data
//		lsb.target(newTarget);
//		double[] newStart = {.001,.001,.001};
//		//set initial parameters
//		lsb.start(newStart);
//		//set upper limit of evaluation time
//		lsb.maxEvaluations(9000);
//		//set upper limit of iteration time
//		lsb.maxIterations(20000);
//
//		//construct LevenbergMarquardtOptimizer 
//		/*LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer(9000, 
//																			0.000000001,
//																			0.0000001,
//																			0.0000001,
//																			Precision.SAFE_MIN);*/
//		
//		LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
//		try
//		{
//			//do LevenbergMarquardt optimization
//			LeastSquaresOptimizer.Optimum lsoo = lmo.optimize(lsb.build());
//			
//			//get optimized parameters
//			final double[] optimalValues = lsoo.getPoint().toArray();			
//			//output data
//			IJ.log("A: " + optimalValues[0]);
//			IJ.log("B: " + optimalValues[1]);
//			IJ.log("C: " + optimalValues[2]);
//			IJ.log("Iteration number: "+lsoo.getIterations());
//			IJ.log("Evaluation number: "+lsoo.getEvaluations());
//			
//			return optimalValues;
//			
//		} 
//		catch (Exception e) 
//		{
//			System.out.println(e.toString());
//			return null;
//		}
//		
//		
//	}
//	

	@SuppressWarnings("deprecation")
	private double[] entzerungskoeefizienten_radial_berechnen_simple(List<SimplePair> punkt_paare)
	{
		RadialDistFunction_simple qf = new RadialDistFunction_simple(punkt_paare);
		LeastSquaresBuilder lsb = new LeastSquaresBuilder();

		//set model function and its jacobian
		lsb.model(qf.retMVF(), qf.retMMF());
		double[] newTarget = qf.realTargetPoints();
		
		//set target data
		lsb.target(newTarget);
		double[] newStart = {.001,.001,.001};
		//set initial parameters
		lsb.start(newStart);
		//set upper limit of evaluation time
		lsb.maxEvaluations(9000);
		//set upper limit of iteration time
		lsb.maxIterations(20000);

		//construct LevenbergMarquardtOptimizer 
		/*LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer(9000, 
																			0.000000001,
																			0.0000001,
																			0.0000001,
																			Precision.SAFE_MIN);*/
		
		LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		try
		{
			//do LevenbergMarquardt optimization
			LeastSquaresOptimizer.Optimum lsoo = lmo.optimize(lsb.build());
			
			//get optimized parameters
			final double[] optimalValues = lsoo.getPoint().toArray();			
			//output data
			IJ.log("A: " + optimalValues[0]);
			IJ.log("B: " + optimalValues[1]);
			IJ.log("C: " + optimalValues[2]);
			IJ.log("Iteration number: "+lsoo.getIterations());
			IJ.log("Evaluation number: "+lsoo.getEvaluations());
			
			return optimalValues;
			
		} 
		catch (Exception e) 
		{
			System.out.println(e.toString());
			return null;
		}
		
		
	}
	
	
	//Berechnet anhand der Punktpaare die Affine Verzerrungsmatrix und gibt diese zur�ck
	private double[][] entzerrungsmatrix_affin_berechnen(List<PointPair> punkt_paare, ImagePlus radial_distorted_pic)
	{
		///Verschiebungsvektor und matrix berechnen
	    
		//Minimum von B*p -c;
		double[][] B = new double[2 * punkt_paare.size()][6];//Matrix B = koordinaten vorlage
		//double[] p;//p = geschter verschiebungsmatrix und vektor
		double[] c = new double[2 * punkt_paare.size()];//c = koordinaten transformiert
	    
	    //konstruiert man die (2n,6) -Matrix
		for(int i = 0; i< punkt_paare.size(); i++)
		{
			//i = 0 => 0,1 => i *2, i*2+1
			//i = 1 => 2,3
			//i = 2 => 4,5
			//i = 3 => 6,7
			
			//
			c[i*2] = punkt_paare.get(i).x_target;//x_i'
			c[i*2+1] = punkt_paare.get(i).y_target; //y_i'
			
			B[i*2][0] = punkt_paare.get(i).x_source; //x_i
			B[i*2][1] = punkt_paare.get(i).y_source; //y_i
			B[i*2][2] = 1;
			B[i*2][3] = 0;
			B[i*2][4] = 0;
			B[i*2][5] = 0;
			
			B[i*2+1][0] = 0;
			B[i*2+1][1] = 0;
			B[i*2+1][2] = 0;
			B[i*2+1][3] = punkt_paare.get(i).x_source;  //x_i
			B[i*2+1][4] = punkt_paare.get(i).y_source; //y_i
			B[i*2+1][5] = 1;
		}
	    
		//B*p = c => p = (Bt * B)^-1 * Bt * c
		
		RealVector c_vec = new ArrayRealVector(c, false);
		RealMatrix B_mat = MatrixUtils.createRealMatrix(B);
		
		//Bt
		RealMatrix B_mat_transp = B_mat.transpose();
		
		//Bt * B
		RealMatrix Bt_B = B_mat_transp.multiply(B_mat);
		
		//(Bt * B)^-1
		RealMatrix Bt_B_inv = MatrixUtils.inverse(Bt_B);
		
		//(Bt * B)^-1 * Bt * c
		RealVector p = Bt_B_inv.operate(B_mat_transp.operate(c_vec));
		
		double[] verschiebungsmatrix = new double[]{ p.getEntry(0), p.getEntry(1), p.getEntry(3), p.getEntry(4) };
		double[] verschiebungsvektor = new double[]{ p.getEntry(2), p.getEntry(5)};
		double[][] verzerrungsmatrix = {{verschiebungsmatrix[0], verschiebungsmatrix[1], verschiebungsvektor[0]}, {verschiebungsmatrix[2], verschiebungsmatrix[3], verschiebungsvektor[1]}, {0,0,1} };
		   
		IJ.log("a11: " + verschiebungsmatrix[0] + " a12: " + verschiebungsmatrix[1] + " a21: " + verschiebungsmatrix[2] + " a22: " + verschiebungsmatrix[3] + " t1: " + verschiebungsvektor[0] + " t2: " + verschiebungsvektor[1]);
		
		//Berechnen Sie aus den Werten die Umkehrabbildung, 
		
		//erzeugen Sie ein neues Bild gleicher Gr��e wie das Eingabebild 
		ImageProcessor processor;
		if(radial_distorted_pic.getProcessor() instanceof ShortProcessor)
		{
			processor = new ShortProcessor(sourcePicture.getWidth(), sourcePicture.getHeight(), true);
		}
		else //other types to be implemented
		{
			return null;
		}
		
		
		
		//berechnen Sie es nach dem Target-to-Source-Verfahren, 
		//d.h.berechnen Sie zu jedem Punkt im neuen Bild durch die Umkehrabbildung die Koordinaten im
		//Eingangsbild und holen Sie von dort den Wert des Pixels. Eine Interpolation wird nicht verlangt.
		
		RealMatrix Bi_mat = MatrixUtils.createRealMatrix(verzerrungsmatrix);
	    DecompositionSolver Bi_mat_solver = new LUDecomposition(Bi_mat).getSolver();
	    
	    //zielbild durchlaufen und werde aus vorlagebild abrufen
        for (int y= 0; y < processor.getHeight();y++) 
		{
			for (int x= 0; x < processor.getWidth(); x++)
			{
				double[] t_h = {x,y,1};
				RealVector t_vec = new ArrayRealVector(t_h, false);
				RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); //Transformation berechnen
				
				int x_coord_vorlage =  (int) Math.round(coord_vec.getEntry(0));
				int y_coord_vorlage =  (int) Math.round(coord_vec.getEntry(1));
				
				if(x_coord_vorlage < radial_distorted_pic.getWidth() && y_coord_vorlage < radial_distorted_pic.getHeight())
				{
					processor.putPixel(x, y, (int)Math.round(radial_distorted_pic.getProcessor().getInterpolatedPixel(x_coord_vorlage, y_coord_vorlage)));
				}
				else //mit schwarz auff�llen				
				{
					processor.putPixel(x, y, 0);
				}			
			}
		}
		
      //Zeichne ziel punkte:
        processor.setColor(Color.BLACK);
        processor.setLineWidth(3);
		for (PointPair point : punkt_paare) 
		{
			processor.drawOval((int)point.x_target, (int)point.y_target, 3, 3);

		}
		
		//Am Ende m�ssen Sie das neue Bild noch auf den Bildschirm bringen.
        ImagePlus newImg = new ImagePlus("Affin Result", processor);
		newImg.show();
		
		return verzerrungsmatrix;
		
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) 
	{
		this.sourcePicture = arg1;
		debug = sourcePicture.getProcessor().duplicate();
		debug.setLineWidth(2);
		debug.setColor(Color.WHITE);
		xCenter=debug.getWidth()/2;
		yCenter=debug.getHeight()/2;
		return DOES_ALL;
	}
}