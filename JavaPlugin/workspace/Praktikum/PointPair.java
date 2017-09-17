
public class PointPair {
	/**
	 * Start-X Koordinate
	 */
	double x_dist = 0.00;
	
	/**
	 * Start-Y Koordinate
	 */
	double y_dist = 0.00;
	
	/**
	 * Ziel- X Koordinate
	 */
	double x_undist = 0.00;
	/**
	 * Ziel- Y Koordinate
	 */
	double y_undist = 0.00;
	
	/**
	 * Abstand zum Bildmittelpunkt
	 */
	double r;
	
	/**
	 * Konstruktor
	 * @param distorted Start-Koordinate
	 * @param undistorted Ziel-Koordinate
	 */
	PointPair(double x_dist,double y_dist, double x_undist, double y_undist)
	{
		this.x_dist=x_dist;
		this.y_dist=y_dist;
		this.x_undist=x_undist;
		this.y_undist=y_undist;
	}
	
	/**
	 * Konstruktor
	 * @param distorted Start-Koordinate
	 * @param undistorted Ziel-Koordinate
	 * @param radius Abstand zum Gitter-Mittelpunkt
	 * @param index Punkt-Index der Textdatei
	 */
	PointPair(double x_dist,double y_dist, double x_undist, double y_undist, double radius)
	{
		this.x_dist=x_dist;
		this.y_dist=y_dist;
		this.x_undist=x_undist;
		this.y_undist=y_undist;
		this.r = radius;
	}
	
	/**
	 * Berechnet den Abstand zum Gittermittelpunkt
	 * @param x x-Koordinate des Punktes
	 * @param y y-Koordinate des Punktes
	 * @param x0 x-Koordinate des Mittelpunktes
	 * @param y0 y-Koordinate des Mittelpunktes
	 * @return
	 */
	public double computeRadius2Center(double x, double y, double x0, double y0) {
		return Math.sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
	}
	
	
}
