
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
	 * @param distorted Start-Koordinate
	 * @param undistorted Ziel-Koordinate
	 */
	SimplePair(double distorted,double undistorted)
	{
		this.distorted = distorted;
		this.undistorted = undistorted;
	}
	
	/**
	 * Konstruktor
	 * @param distorted Start-Koordinate
	 * @param undistorted Ziel-Koordinate
	 * @param radius Abstand zum Gitter-Mittelpunkt
	 * @param index Punkt-Index der Textdatei
	 */
	SimplePair(double distorted,double undistorted, double radius, double index)
	{
		this.distorted = distorted;
		this.undistorted = undistorted;
		this.radius = radius;
		this.index = index;
	}
	
	/**
	 * Start-Koordinate
	 */
	double distorted = 0.00;
	
	/**
	 * Ziel- Koordinate
	 */
	double undistorted = 0.00;
	
	/**
	 * Abstand zum Gittermittelpunkt
	 */
	double radius;
	
	/**
	 * Index des Punktes in der Textdatei
	 */
	double index;
}
