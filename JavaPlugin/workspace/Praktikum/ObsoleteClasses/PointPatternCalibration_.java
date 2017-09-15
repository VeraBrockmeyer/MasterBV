import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.text.TextPanel;
import ij.text.TextWindow;

public class PointPatternCalibration_ implements PlugInFilter {


	private TextPanel textInput;
	private ArrayList <Integer> sourcepoints, targetpoints;
	private ImageProcessor ip;
	
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")){
			IJ.showMessage("About Standard Plugin", 
					"Standard Plugin zur einfachen Implementierung von Pixelwertaenderungen"); 
			return DONE;
		}

		this.ip = imp.getProcessor();
		return DOES_RGB | DOES_8G;
	}
		public void run(ImageProcessor arg0) {
		if( !setupDialog() )
			return;
		BufferedReader br = new BufferedReader( new StringReader( textInput.getText() ));
		String line;
		sourcepoints = new ArrayList<Integer>();
		targetpoints = new ArrayList<Integer>();
		try {
			br.readLine();
			while( ( line = br.readLine() ) != null ){
				String [] data = line.split( "\\s+" );
				int leadingEmpty = data[ 0 ].length() == 0 ? 1 : 0;
				if( data.length >= 5 + leadingEmpty ) {
					String outStr = "";
					for( int i = 1 + leadingEmpty; i < data.length; i++ ){
						int val = new Integer( data[ i ] );
						if((leadingEmpty == 0 && i<=2) || (leadingEmpty == 1 && i <=3)) sourcepoints.add(val);
						else targetpoints.add(val);
						outStr = outStr + val;
						if( i < data.length - 1 ) {
							outStr = outStr + "  ";
						}
					}
					//IJ.log( outStr );
				}
			}
		} catch (IOException e) {
			return;
		}
		computeAffTransform();

	}

	private void computeAffTransform() {
			double[][]B = new double[sourcepoints.size()][6];
			double[]c = new double[targetpoints.size()];
			
			for (int i = 0; i < B.length; i+=2) {
				B[i][0] = sourcepoints.get(i);
				B[i][1] = sourcepoints.get(i+1);
				B[i][2] = 1.;
				B[i][3] = 0.;
				B[i][4] = 0.;
				B[i][5] = 0.;
				B[i+1][0] = 0.;
				B[i+1][1] = 0.;
				B[i+1][2] = 0.;
				B[i+1][3] = sourcepoints.get(i);
				B[i+1][4] = sourcepoints.get(i+1);
				B[i+1][5] = 1.;
				c[i]= targetpoints.get(i);
				c[i+1]=targetpoints.get(i+1);
			}
			for (int i = 0; i < B.length; i++) {
				for (int j = 0; j < B[0].length; j++) {
					System.out.print(""+B[i][j]+", ");
				}
			System.out.println();
			}
			
			RealMatrix A = new Array2DRowRealMatrix(B);
			RealMatrix At = A.transpose();
			RealMatrix AtAinv = At.multiply(A);
			AtAinv = MatrixUtils.inverse(AtAinv);
			double[]p = At.operate(c);
			p = AtAinv.operate(p);
			double[][] mat = new double [3][3];
			mat[2][2]=1.;
			for (int i = 0; i < 3; i++) {
				mat[0][i] = p[i];
				mat[1][i] = p[i + 3];
			}
			
				
			RealMatrix realMat = new Array2DRowRealMatrix(mat);
			int width = ip.getWidth();
			int height = ip.getHeight();
			ImageProcessor ipOut = ip.duplicate();
			ipOut.set(0.);
			double [] xIn = new double [3];
			xIn[2]=1.;
			double[] xOut = new double [3];
			
			for (int row1 = 0; row1 < height; row1++) {
				for (int col = 0; col < width; col++) {
					xIn[0]=row1;
					xIn[1]=col;
					xOut = realMat.operate(xIn);
					ipOut.putPixel(col, row1, ip.getPixel((int)xOut[1], (int)xOut[0]));
				}
			}
			ImagePlus imp = new ImagePlus("Result", ipOut);
			imp.show();
			
		}

	private boolean setupDialog() {
		int i, tCnt, wListL;
		int [] wNum;
		java.awt.Frame [] wList = WindowManager.getNonImageWindows();
		wListL = wList.length;
		wNum = new int [ wListL ];
		for( i = 0, tCnt = 0; i < wListL; i++ ){
			if( wList[ i ] instanceof TextWindow ) {
				wNum[ tCnt ] = i;
				tCnt++;
			}
		}
		String [] tList = new String [ tCnt ];
		for( i = 0; i < tCnt; i++ ){
			tList [ i ] = wList[ wNum[ i ] ].getTitle();
		}
		GenericDialog gd = new GenericDialog( "Datenfenster" );
		gd.addChoice("Textfenster", tList, tList[0] );
		gd.showDialog();
		if( gd.wasCanceled() ){
			return false;
		}
		textInput = 
				((TextWindow)( wList[ wNum[ gd.getNextChoiceIndex() ] ] )).getTextPanel();

		return true;
	}

	

}