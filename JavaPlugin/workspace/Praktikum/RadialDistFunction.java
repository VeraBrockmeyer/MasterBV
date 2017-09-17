
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
	List<SimplePair> _pointPairs =  new ArrayList<SimplePair>();
	
	/**
	 * Konstruktor
	 * @param point_pairs Start und Ziel Koordinaten der Punkte
	 */
	RadialDistFunction(List<SimplePair> point_pairs)
	{
		this._pointPairs = point_pairs;
	}
	
    /**
     * Gibt die Ziel-Punkt-Koordinaten als double Array aus für den LevenbergMarquadt Optimierer
     * @return target	double array mit den ZielKoordinaten 
     */
    public double[] realTargetPoints() 
    {
        double[] target = new double[_pointPairs.size()]; //Speicher für ZielKoordinaten
        
        for (int i = 0; i < _pointPairs.size(); i++) 
        {
            target[i] = _pointPairs.get(i).undistorted; //zielwerte x
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
			public double[] value (double[] radial_dist_coeff) throws IllegalArgumentException 
			{
				IJ.log("MVF called:");
		        double[] calculated_target_points = new double[_pointPairs.size() ];
		        
		        //x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        for (int i = 0; i < _pointPairs.size(); ++i) 
		        {
		        	IJ.log("Point par values: R " + _pointPairs.get(i).radius + " source" + _pointPairs.get(i).distorted + " target " + _pointPairs.get(i).undistorted );
		        	
		        	//radiale verzerrung in x-richtung:
		        	calculated_target_points[i] =  _pointPairs.get(i).distorted 
		        			/(1 
		        			+ radial_dist_coeff[0] * Math.pow(_pointPairs.get(i).radius, 2.00) 
		        			+ radial_dist_coeff[1] * Math.pow(_pointPairs.get(i).radius, 4.00)
		        			+ radial_dist_coeff[2] * Math.pow(_pointPairs.get(i).radius, 6.00)
		        			);
		        	
//		        
		        	
		        	 IJ.log(String.format("Calculated values: target = %f target_cal = %f ;", _pointPairs.get(i).undistorted,calculated_target_points[i])); 			           		 		
		        }
		        
		        return calculated_target_points;
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
				// TODO Auto-generated method stub
                return jacobian(point);
			}

			/**
			 * calculate and retjacobian
			 * @param	variables	parameters of model function
			 * @return	jacobian	jacobian of the model function
			 */
		    private double[][] jacobian(double[] variables) 
		    {
		    	IJ.log("MMF called:");
		        double[][] jacobian = new double[_pointPairs.size()][3];
		        
		        for (int i = 0; i < _pointPairs.size(); ++i) 
		        {
		        	//x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        	//dx'/da = r^2*x
		        	//dx'/db = r^4*x
		        	//dx'/dc = r^6*x
		        	double r2 = Math.pow(_pointPairs.get(i).radius, 2.00);
		        	double r4 = Math.pow(_pointPairs.get(i).radius, 4.00);
		        	double r6 = Math.pow(_pointPairs.get(i).radius, 6.00);
		        	double term1 = (1.+variables[0]*r2+variables[1]*r4+variables[2]*r6);
		        	double term2 = -1.*(1/(term1*term1));
		        	
		        	
		            jacobian[i][0] = _pointPairs.get(i).distorted*r2*term2 ; 
		            jacobian[i][1] = _pointPairs.get(i).distorted*r4*term2 ; 
		            jacobian[i][2] = _pointPairs.get(i).distorted*r6*term2 ; 
		        	
//		        	
//		            jacobian[i][0] = Math.pow(_pointPairs.get(i).radius, 2.00)  * _pointPairs.get(i).distorted; 
//		            jacobian[i][1] = Math.pow(_pointPairs.get(i).radius, 4.00)  * _pointPairs.get(i).distorted; 
//		            jacobian[i][2] = Math.pow(_pointPairs.get(i).radius, 6.00)  * _pointPairs.get(i).distorted;
		            
//		         
		            //IJ.log(String.format("Jacobian values: x_1 %f x_2 %f x_3 %f ; y_1 %f y_2 %f y_3 %f", jacobian[i][0], jacobian[i][1], jacobian[i][2],  jacobian[i+_pointPairs.size()][0],  jacobian[i+_pointPairs.size()][1],  jacobian[i+_pointPairs.size()][2])); 			           		 		

		        }
		        return jacobian;
		    }
			
		};
    }

}
