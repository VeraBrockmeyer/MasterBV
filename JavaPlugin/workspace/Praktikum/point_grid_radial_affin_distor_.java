import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class point_grid_radial_affin_distor_ implements PlugInFilter
{
	ImagePlus sourcePicture;
	
	@Override
	public void run(ImageProcessor img) 
	{
		try
		{
			//Textfenster mit Punktpaare Textdtei wählen:
			java.awt.Window[] non_img_windows = WindowManager.getAllNonImageWindows();
			GenericDialog gd = new GenericDialog("Textdatei auswählen:");
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
		  	for(int i=0; i < non_img_windows.length; i++)
			{
				if(((java.awt.Frame)non_img_windows[i]).getTitle() == choise)
				{
	
					String text = ((ij.plugin.frame.Editor)non_img_windows[i]).getText();
						
					String[] lines = text.split("\n");
					
					List<PointPair> PointPairs = new ArrayList<PointPair>(); //speicher für punkt paarungen
					
					for(int j = 1; j < lines.length; j++) //erste zeile überspringen
					{
						String[] numbers = lines[j].split("\t");
						PointPair txt_pointPair = new PointPair(
								Integer.parseInt(numbers[1].trim()), //x source
								Integer.parseInt(numbers[2].trim()), //y source
								Integer.parseInt(numbers[3].trim()), //x' target
								Integer.parseInt(numbers[4].trim()), //y' target
								
								//r target
								Math.sqrt(Math.pow(Integer.parseInt(numbers[1].trim()) - img.getWidth() /2, 2.00) 
										+ Math.pow(Integer.parseInt(numbers[2].trim()) - img.getHeight() /2, 2.00)
										 ));
						
						PointPairs.add(txt_pointPair);
					}
					
					
					//TODO: radiale entzerrung berechnen
					double[] radial_dist_koeff = entzerungskoeefizienten_radial_berechnen(PointPairs);
					
					//Entzerrung auf punkte anwenden
					for(int j=0; j < PointPairs.size(); j++)
					{
						//x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
						PointPairs.get(j).x_target = (1 
			        			+ radial_dist_koeff[0] * Math.pow(PointPairs.get(i).r_source, 2.00) 
			        			+ radial_dist_koeff[1] * Math.pow(PointPairs.get(i).r_source, 4.00)
			        			+ radial_dist_koeff[2] * Math.pow(PointPairs.get(i).r_source, 6.00)
			        			) * PointPairs.get(i).x_source ;
						
						PointPairs.get(j).y_target = (1 
			        			+ radial_dist_koeff[0] * Math.pow(PointPairs.get(i).r_source, 2.00) 
			        			+ radial_dist_koeff[1] * Math.pow(PointPairs.get(i).r_source, 4.00)
			        			+ radial_dist_koeff[2] * Math.pow(PointPairs.get(i).r_source, 6.00)
			        			) * PointPairs.get(i).y_source ;
					}
					
					
					//erzeugen Sie ein neues Bild gleicher Größe wie das Eingabebild 
					ImagePlus newImg;
					if(sourcePicture.getBitDepth() <= 8)
					{
						newImg =  NewImage.createByteImage ("entzertes_Bild", sourcePicture.getWidth(), sourcePicture.getHeight(), 1, 
						NewImage.FILL_WHITE); 
					}
					else
					{
						newImg =  NewImage.createRGBImage ("entzertes_Bild", sourcePicture.getWidth(), sourcePicture.getHeight(), 1, 
						NewImage.FILL_WHITE); 
					}
					
					for (int y= 0; y < newImg.getHeight();y++) 
					{
						for (int x= 0; x < newImg.getWidth(); x++)
						{
							int x_coord_vorlage =  (int) Math.round((1 
				        			+ radial_dist_koeff[0] * Math.pow(PointPairs.get(i).r_source, 2.00) 
				        			+ radial_dist_koeff[1] * Math.pow(PointPairs.get(i).r_source, 4.00)
				        			+ radial_dist_koeff[2] * Math.pow(PointPairs.get(i).r_source, 6.00)
				        			) * x);
							
							int y_coord_vorlage =  (int) Math.round((1 
				        			+ radial_dist_koeff[0] * Math.pow(PointPairs.get(i).r_source, 2.00) 
				        			+ radial_dist_koeff[1] * Math.pow(PointPairs.get(i).r_source, 4.00)
				        			+ radial_dist_koeff[2] * Math.pow(PointPairs.get(i).r_source, 6.00)
				        			) * y);
							
							if(x_coord_vorlage < sourcePicture.getWidth() && y_coord_vorlage < sourcePicture.getHeight())
							{
								newImg.getStack().getProcessor(1).putPixel(x, y, sourcePicture.getStack().getProcessor(1).getPixel(x_coord_vorlage, y_coord_vorlage));
							}
							else //mit schwarz auffüllen				
							{
								newImg.getStack().getProcessor(1).putPixel(x, y, 0);
							}	
						}	
					}
					//Am Ende müssen Sie das neue Bild noch auf den Bildschirm bringen.
					newImg.show();
					
					
					//Affine entzerrung berechnen
					//this.entzerrungsmatrix_affin_berechnen(PointPairs);
					
					break; //schleifendurchlauf suche durch bilder abbrechen
				}
			}
		} 
		catch(Exception exc)
		{
			IJ.error(exc.getMessage() + exc.getClass() + exc.getCause() + exc.getStackTrace());
		}			
	}
	
	@SuppressWarnings("deprecation")
	private double[] entzerungskoeefizienten_radial_berechnen(List<PointPair> punkt_paare)
	{
		RadialDistFunction qf = new RadialDistFunction(punkt_paare);
		LeastSquaresBuilder lsb = new LeastSquaresBuilder();

		//set model function and its jacobian
		lsb.model(qf.retMVF(), qf.retMMF());
		double[] newTarget = qf.realTargetPoints();
		
		//set target data
		lsb.target(newTarget);
		double[] newStart = {10,10,10};
		//set initial parameters
		lsb.start(newStart);
		//set upper limit of evaluation time
		lsb.maxEvaluations(1000);
		//set upper limit of iteration time
		lsb.maxIterations(10000);

		//construct LevenbergMarquardtOptimizer 
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
	

	//Berechnet anhand der Punktpaare die Affine Verzerrungsmatrix und gibt diese zurück
	private double[][] entzerrungsmatrix_affin_berechnen(List<PointPair> punkt_paare)
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
		
		//erzeugen Sie ein neues Bild gleicher Größe wie das Eingabebild 
		ImagePlus newImg;
		if(sourcePicture.getBitDepth() <= 8)
		{
			newImg =  NewImage.createByteImage ("entzertes_Bild", sourcePicture.getWidth(), sourcePicture.getHeight(), 1, 
			NewImage.FILL_WHITE); 
		}
		else
		{
			newImg =  NewImage.createRGBImage ("entzertes_Bild", sourcePicture.getWidth(), sourcePicture.getHeight(), 1, 
			NewImage.FILL_WHITE); 
		}
		
		
		//berechnen Sie es nach dem Target-to-Source-Verfahren, 
		//d.h.berechnen Sie zu jedem Punkt im neuen Bild durch die Umkehrabbildung die Koordinaten im
		//Eingangsbild und holen Sie von dort den Wert des Pixels. Eine Interpolation wird nicht verlangt.
		
		RealMatrix Bi_mat = MatrixUtils.createRealMatrix(verzerrungsmatrix);
	    DecompositionSolver Bi_mat_solver = new LUDecomposition(Bi_mat).getSolver();
	    
	    //zielbild durchlaufen und werde aus vorlagebild abrufen
        for (int y= 0; y < newImg.getHeight();y++) 
		{
			for (int x= 0; x < newImg.getWidth(); x++)
			{
				double[] t_h = {x,y,1};
				RealVector t_vec = new ArrayRealVector(t_h, false);
				RealVector coord_vec = Bi_mat_solver.getInverse().operate(t_vec); //Transformation berechnen
				
				int x_coord_vorlage =  (int) Math.round(coord_vec.getEntry(0));
				int y_coord_vorlage =  (int) Math.round(coord_vec.getEntry(1));
				
				if(x_coord_vorlage < sourcePicture.getWidth() && y_coord_vorlage < sourcePicture.getHeight())
				{
					newImg.getStack().getProcessor(1).putPixel(x, y, sourcePicture.getStack().getProcessor(1).getPixel(x_coord_vorlage, y_coord_vorlage));
				}
				else //mit schwarz auffüllen				
				{
					newImg.getStack().getProcessor(1).putPixel(x, y, 0);
				}			
			}
		}
		
		
		//Am Ende müssen Sie das neue Bild noch auf den Bildschirm bringen.
		newImg.show();
		
		return verzerrungsmatrix;
		
	}

	@Override
	public int setup(String arg0, ImagePlus arg1) 
	{
		this.sourcePicture = arg1;
		return DOES_ALL;
	}
}
