import ij.*;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;
import java.awt.*;

public class Standard_ implements PlugInFilter {

	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")){
			IJ.showMessage("About Standard Plugin", 
			"Standard Plugin zur einfachen Implementierung von Pixelwertaenderungen"); 
			return DONE;
		}
		return DOES_8G;
	}

	public void run(ImageProcessor ip) {
		
		int width = ip.getWidth();			
		int height = ip.getHeight();
		byte[] pixels = (byte[]) ip.getPixels();
				
		for (int y=0; y<height; y++) {
			for (int x=0; x<width; x++) {
				pixels[y*width+x] = (byte)(0xff & ProcessPixel(x, y, width, height, 0xff & pixels[y*width+x]) );
			}
		}
	}

	private int ProcessPixel(int x, int y, int breite, int hoehe, int i){
	    int result;
// hier ist die Berechnung einzufügen
        
    result = i;
// ----------------------------------        
        return result;
	}
}
