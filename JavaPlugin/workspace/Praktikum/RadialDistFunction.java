
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;

import ij.IJ;

@SuppressWarnings("deprecation")
public class RadialDistFunction 
{

	List<PointPair> _pointPairs =  new ArrayList<PointPair>();
	
	//Konstruktor
	RadialDistFunction(List<PointPair> point_pairs)
	{
		this._pointPairs = point_pairs;
	}
	
    /**
     * return target data as double array by target data
     * @return target	double array of target data
     */
    public double[] realTargetPoints() 
    {
        double[] target = new double[_pointPairs.size() * 2]; //Speicher für ZielKoordinaten
        
        for (int i = 0; i < _pointPairs.size(); i++) 
        {
            target[i] = _pointPairs.get(i).x_target; //zielwerte x
            target[i + _pointPairs.size()] = _pointPairs.get(i).y_target; //zielwerte y
        }
        
        return target;
    }
    
    /**
     * Define model function and return values
     * @return	return the values of model function by input data
     */
    public MultivariateVectorFunction retMVF() 
    {
		return new MultivariateVectorFunction() 
		{
			@Override
			public double[] value (double[] radial_verz_koeffizienten) throws IllegalArgumentException 
			{
				IJ.log("MVF called:");
		        double[] calculated_target_points = new double[_pointPairs.size() * 2];
		        
		        //x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        for (int i = 0; i < _pointPairs.size(); ++i) 
		        {
		        	//radiale verzerrung in x-richtung:
		        	calculated_target_points[i] = (1 
		        			+ radial_verz_koeffizienten[0] * Math.pow(_pointPairs.get(i).r_source, 2.00) 
		        			+ radial_verz_koeffizienten[1] * Math.pow(_pointPairs.get(i).r_source, 4.00)
		        			+ radial_verz_koeffizienten[2] * Math.pow(_pointPairs.get(i).r_source, 6.00)
		        			) * _pointPairs.get(i).x_source ;
		        	
		        	//radiale verzerrung in y-richtung
		        	calculated_target_points[i+_pointPairs.size()] = (1 
		        			+ radial_verz_koeffizienten[0] * Math.pow(_pointPairs.get(i).r_source, 2.00) 
		        			+ radial_verz_koeffizienten[1] * Math.pow(_pointPairs.get(i).r_source, 4.00)
		        			+ radial_verz_koeffizienten[2] * Math.pow(_pointPairs.get(i).r_source, 6.00)
		        			) * _pointPairs.get(i).y_source ;
		        	
		        	 IJ.log(String.format("Calculated values: X_t = %f  X_t_cal = %f ; Y_t = %f Y_t_cal = %f", _pointPairs.get(i).x_target,calculated_target_points[i],_pointPairs.get(i).y_target,calculated_target_points[i+_pointPairs.size()])); 			           		 		
		        }
		        
		        return calculated_target_points;
		    }			
		};
    	
    }
    
    
    
    /**
     * Return the jacobian of the model function
     * @return	return the jacobian
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
		        double[][] jacobian = new double[_pointPairs.size() * 2][3];
		        
		        for (int i = 0; i < _pointPairs.size(); ++i) 
		        {
		        	//x' = (1 + a*r^2 + b*r^4 + c*r^6) * x
		        	//dx'/da = a*r^2*x
		        	//dx'/db = b*r^4*x
		        	//dx'/dc = c*r^6*x
		        	
		            jacobian[i][0] = Math.pow(_pointPairs.get(i).r_source, 2.00)  * _pointPairs.get(i).x_source; 
		            jacobian[i][1] = Math.pow(_pointPairs.get(i).r_source, 4.00)  * _pointPairs.get(i).x_source; 
		            jacobian[i][2] = Math.pow(_pointPairs.get(i).r_source, 6.00)  * _pointPairs.get(i).x_source;
		            
		            jacobian[i+_pointPairs.size()][0] = Math.pow(_pointPairs.get(i).r_source, 2.00)  * _pointPairs.get(i).y_source; 
		            jacobian[i+_pointPairs.size()][1] = Math.pow(_pointPairs.get(i).r_source, 4.00)  * _pointPairs.get(i).y_source; 
		            jacobian[i+_pointPairs.size()][2] = Math.pow(_pointPairs.get(i).r_source, 6.00)  * _pointPairs.get(i).y_source;
		        
		            IJ.log(String.format("Jacobian values: x_1 %f x_2 %f x_3 %f ; y_1 %f y_2 %f y_3 %f", jacobian[i][0], jacobian[i][1], jacobian[i][2],  jacobian[i+_pointPairs.size()][0],  jacobian[i+_pointPairs.size()][1],  jacobian[i+_pointPairs.size()][2])); 			           		 		

		        }
		        return jacobian;
		    }
			
		};
    }

}
