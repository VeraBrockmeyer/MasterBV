
public class SimplePair 
{

	SimplePair()
	{
		
	}
	
	SimplePair(double source,double target)
	{
		this.source = source;
		this.target = target;
	}
	
	SimplePair(double source,double target, double radius, double index)
	{
		this.source = source;
		this.target = target;
		this.radius = radius;
		this.index = index;
	}
	
	double source = 0.00;
	double target = 0.00;
	
	double radius;
	
	double index;
}
