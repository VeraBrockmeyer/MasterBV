
public class PointPair 
{

	PointPair()
	{
		
	}
	
	PointPair(double x_source, double y_source, double x_target, double y_target)
	{
		this.x_source = x_source;
		this.y_source = y_source;
		this.x_target = x_target;
		this.y_target = y_target;
	}
	
	PointPair(double x_source, double y_source, double x_target, double y_target, double r_source, double index)
	{
		this.x_source = x_source;
		this.y_source = y_source;
		this.x_target = x_target;
		this.y_target = y_target;
		
		this.r_source = r_source;
		this.index = index;
	}
	
	double x_source = 0.00;
	double y_source = 0.00;
	double x_target = 0.00;
	double y_target = 0.00;
	
	double r_source;
	
	double index;
}
