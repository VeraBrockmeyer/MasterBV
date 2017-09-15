
/**
 * Klasse zum Speichern von Ziel und Start Koordinaten und des Zugehörigen Radius Wertes
 * @author 
 *
 */
public class SimplePair 
{

	/**
	 * Konstruktor
	 */
	SimplePair()
	{
		
	}
	
	/**
	 * Konstruktor
	 * @param source Start-Koordinate
	 * @param target Ziel-Koordinate
	 */
	SimplePair(double source,double target)
	{
		this.source = source;
		this.target = target;
	}
	
	/**
	 * Konstruktor
	 * @param source Start-Koordinate
	 * @param target Ziel-Koordinate
	 * @param radius Abstand zum Gitter-Mittelpunkt
	 * @param index Punkt-Index der Textdatei
	 */
	SimplePair(double source,double target, double radius, double index)
	{
		this.source = source;
		this.target = target;
		this.radius = radius;
		this.index = index;
	}
	
	/**
	 * Start-Koordinate
	 */
	double source = 0.00;
	
	/**
	 * Ziel- Koordinate
	 */
	double target = 0.00;
	
	/**
	 * Abstand zum Gittermittelpunkt
	 */
	double radius;
	
	/**
	 * Index des Punktes in der Textdatei
	 */
	double index;
}
