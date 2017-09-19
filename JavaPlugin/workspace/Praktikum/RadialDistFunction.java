import java.awt.Point;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;

import ij.IJ;

/**
 * Klasse für den Optimierer zur Berechnung der Verzerrung durch Vorgabe der Koeffizienten
 * @author 
 *
 */
@SuppressWarnings("deprecation")
public class RadialDistFunction 
{

	/**
	 * Interner Speicher für die Punkt-Paare
	 */
	ArrayList<PointPair> _pointPairs =  new ArrayList<PointPair>();
	
	/**
	 * Konstruktor
	 * @param point_pairs Start und Ziel Koordinaten der Punkte
	 */
	RadialDistFunction(ArrayList<PointPair> point_pairs)
	{
		this._pointPairs = point_pairs;
	}
	
    /**
     * Gibt die Ziel-Punkt-Koordinaten als double Array aus für den LevenbergMarquadt Optimierer
     * @return target	double array mit den ZielKoordinaten 
     */
    public double[] realTargetPoints() 
    {
        double[] target = new double[_pointPairs.size()*2]; //Speicher für ZielKoordinaten
        int counter =0;
        for (int i = 0; i < target.length; i+=2) 
        {
            target[i] = _pointPairs.get(counter).x_undist; //zielwerte x
            target[i+1] = _pointPairs.get(counter).y_undist; //zielwerte y
            counter++;
        }
        
        return target;
    }
    
    /**
     * Ezeugt die Verzerrungs-Funktion für den Optimierer (x' = (1 + a*r^2 + b*r^4 + c*r^6) * x)
     * @return	Verzerrungs-Funktion in der zu den Keffizienten die Ziel-Werte berechnet werden
     */
    public MultivariateVectorFunction retMVF() 
    {
		return new MultivariateVectorFunction() 
		{
			@Override
			public double[] value (double[] coeffs) throws IllegalArgumentException 
			{
				IJ.log("MVF called:");
		        double[] F = new double[_pointPairs.size()*2 ];
		        int counter=0;
		        //x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        for (int i = 0; i < F.length; i+=2) 
		        {
		        	//IJ.log("Point par values: R " + _pointPairs.get(i).radius + " source" + _pointPairs.get(i).source + " target " + _pointPairs.get(i).target );
		        	PointPair pp = _pointPairs.get(counter);
		        	//radiale verzerrung in x-richtung:
		        	F[i] = (1 + coeffs[0] * pp.r + coeffs[1] * pp.r * pp.r + coeffs[2] * pp.r * pp.r * pp.r) * pp.x_dist ;
		        	F[i+1] = (1 + coeffs[0] * pp.r + coeffs[1] * pp.r * pp.r + coeffs[2] * pp.r * pp.r * pp.r) * pp.y_dist ;
   		        	counter++;
   		        	
//   		        	System.out.println(F[i]);
//   		        	System.out.println(F[i+1]);
   		        	IJ.log(String.format("Calculated values: X_undist = %f  X_dist = %f ; Y_undist = %f Y_dist = %f", 
   		        			pp.x_undist,F[i],pp.y_undist,F[i+1])); 			        
		      	 }
			
		        return F;
		    }			
		};
    	
    }
    
    
    
    /**
     * ERzeugt die jacobi-Matrix-Funktion (Ableitung nach den Koeffizienten)
     * @return	Jacobi-Matrix-Funktion berechnet Ableitung nach den Koeffizienten
     */
    public MultivariateMatrixFunction retMMF() 
    {
    	return new MultivariateMatrixFunction() 
    	{

			@Override
			public double[][] value(double[] point) throws IllegalArgumentException 
			{
                return jacobian(point);
			}

			/**
			 * calculate and set jacobian
			 * @param	variables	parameters of model function
			 * @return	jacobian	jacobian of the model function
			 */
		    private double[][] jacobian(double[] variables) 
		    {
		    	IJ.log("MMF called:");
		        double[][] jacobian = new double[_pointPairs.size()*2][3];
		        int counter =0;
		        for (int i = 0; i < jacobian.length; i+=2) 
		        {
		        	PointPair pp = _pointPairs.get(counter);
		        	//x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        	//dx'/da = r^2*x
		        	//dx'/db = r^4*x
		        	//dx'/dc = r^6*x
		        			        	
		            jacobian[i][0] = pp.r * pp.x_dist; 
		            jacobian[i][1] = pp.r * pp.r * pp.x_dist; 
		            jacobian[i][2] = pp.r * pp.r * pp.r * pp.x_dist; 
		            
		            jacobian[i+1][0] = pp.r * pp.x_dist; 
		            jacobian[i+1][1] = pp.r * pp.r * pp.x_dist; 
		            jacobian[i+1][2] = pp.r * pp.r * pp.r * pp.x_dist; 
		            counter++;
//		         	System.out.println(jacobian[i][0] + "\t" + jacobian[i][1] + "\t" + jacobian[i][2] + "\t" );
//		         	System.out.println(jacobian[i+1][0] + "\t" + jacobian[i+1][1] + "\t" + jacobian[i+1][2] + "\t" );
		            IJ.log(String.format("Jacobian values: x_1 %f x_2 %f x_3 %f ; y_1 %f y_2 %f y_3 %f", jacobian[i][0], jacobian[i][1], jacobian[i][2],  jacobian[i+1][0],  jacobian[i+1][1],  jacobian[i+1][2])); 			           		 		

		        }
		        return jacobian;
		    }
			
		};
    }

}