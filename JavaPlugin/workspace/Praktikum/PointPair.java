
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
	 * Index des PointPairs in Liste
	 */
	double index;
	
	
	
	/**
	 * Konstruktor
	 * @param x_dist X-Coord of distorted Point
	 * @param y_dist Y-Coord of distorted Point
	 * @param x_undist X-Coord of undistorted Point
	 * @param y_undist Y-Coord of undistorted Point
	 * @param x0 X-Coord of image center
	 * @param y0 Y-Coord of image center
	 * @param index Index des PointPairs in list
	 */
	PointPair(double x_dist,double y_dist, double x_undist, double y_undist, double x0, double y0, double index)
	{
		this.x_dist=x_dist - x0;
		this.y_dist=y_dist - y0;
		this.x_undist=x_undist - x0;
		this.y_undist=y_undist - y0;
		this.index = index;
		this.r = computeRadius2Center(this.x_dist, this.y_dist, x0, y0);
	}
	
	/**
	 * Berechnet den Abstand zum Gittermittelpunkt
	 * @param x x-Koordinate des Punktes
	 * @param y y-Koordinate des Punktes
	 * @param x0 x-Koordinate des Mittelpunktes
	 * @param y0 y-Koordinate des Mittelpunktes
	 * @return
	 */
	private double computeRadius2Center(double x, double y, double x0, double y0) {
		return  ((x - x0) * (x - x0) + (y - y0) * (y - y0));
	}
	
	
}
