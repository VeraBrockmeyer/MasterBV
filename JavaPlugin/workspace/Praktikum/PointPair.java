
public class PointPair {
	/**
	 * Start-X Koordinate
	 */
	public double x_dist = 0.00;
	
	/**
	 * Start-Y Koordinate
	 */
	public double y_dist = 0.00;
	
	/**
	 * Ziel- X Koordinate
	 */
	public double x_undist = 0.00;
	
	/**
	 * Ziel- Y Koordinate
	 */
	public double y_undist = 0.00;
	
	/**
	 * Abstand zum Bildmittelpunkt r^2
	 */
	public double r;
	

	/**
	 * Index des PointPairs in Liste
	 */
	public double index;
	
	
	
	/**
	 * Konstruktor speichert den Punkt bezogen auf das Koordinatensystem des übergebenen Mittelpunktes und berechnet den Abstand zum Mittelpunkt
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
		this.r = computeRadius2Center(this.x_dist, this.y_dist);
	}
	
	/**
	 * Konstruktor speichert den Punkt mit Bezugskoordinatensystem und berechnet den Abstand zum Mittelpunkt
	 * @param x_undist X-Coord of undistorted Point
	 * @param y_undist Y-Coord of undistorted Point
	 * @param x0 X-Coord of image center
	 * @param y0 Y-Coord of image center
	 */
	PointPair(double x_undist,double y_undist, double x0, double y0)
	{
		this.x_undist=x_undist - x0;
		this.y_undist=y_undist - y0;
		this.r = computeRadius2Center(this.x_undist, this.y_undist);
	}
	

	/**
	 * Berechnet den Abstand zum Gittermittelpunkt für Koordinaten die auf den Gittermittelpunkt bezogen sind (x^2 + y^2)
	 * @param x x-Koordinate des Punktes Mittelpunktbezogen
	 * @param y y-Koordinate des Punktes Mittelpunktbezogen
	 * @return
	 */
	private static double computeRadius2Center(double x, double y) {
		return x * x + y * y;
	}
	
	
}
