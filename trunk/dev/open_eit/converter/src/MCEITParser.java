import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class MCEITParser implements OeitLegacyParser {

	public EitEntity parse(String fileName) throws IOException{
		EitEntity mceitInfo = new EitEntity();
		
		File file = new File(fileName);
		
		FileInputStream fin = new FileInputStream(file);
		BufferedInputStream bin = new BufferedInputStream(fin);
		DataInputStream din = new DataInputStream(bin);
		
		int length = (int)(file.length()/4);
		
		if (length%256!=0){
			System.out.println("File length strange - cropping file");
			length = length-length%256;
		}
		
		double[][] input = new double[256][length/256];
		
/*		method for reading Big-endian
		for (int i = 0; i<256; i++){
			for (int j = 0; j<length/256; j++){
				input[i][j] = (double)din.readFloat();
			}
		}*/
		
		byte[] bytes = new byte[4];
		
		int counter = 0;

		for (int i = 0; i<256; i++){
			System.out.print(counter + ": ");
			for (int j = 0; j<length/256; j++){            //METHOD 1
			bytes[0] = din.readByte();
			bytes[1] = din.readByte();
			bytes[2] = din.readByte();
			bytes[3] = din.readByte();
			input[i][j] = (double)ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN ).getFloat();
			counter++;
			System.out.print(" " + input[i][j]);
			}
			System.out.print("\n");
		}
		
		//http://stackoverflow.com/questions/6471177/converting-8-bytes-of-little-endian-binary-into-a-double-precision-float
		
/*		for (int i = 0; i<256; i++){
			for (int j = 0; j<length/256; j++){           //METHOD 2 FOR BIG ENDIAN
				int ch1 = din.read();
				int ch2 = din.read();
				int ch3 = din.read();
				int ch4 = din.read();
				input[i][j] = (double)(ch1 + (ch2 << 8 ) + (ch3 << 16) + (ch4 << 24));
			}
		}*/
		
        din.close();
        bin.close();
        fin.close();
        
/*        if (input[0][38]==0){
        	//104 measurements per frame
        	//write method to transfer 104 measurement array to 208 measurement array
        	input = transfer104To208(input);
        }*/
        input = untwist(input);
        
        mceitInfo.setOutput2D(input);
        mceitInfo.setLength(length);
        
		return mceitInfo;
	}
	
	public double[][] transfer104To208(double[][] array104){
		double[][] array208 = null;
		int x,y;
		x = array104.length;
		y = array104[0].length;
		
		if (x!=256 || y!=256){
			System.out.println("eidors_readdata: expecting an input array of size 208*n");
			return array104;
		}
		
		for (int i = 0; i < x; i++){
			// to be written
		}
		
		return array208;
	}
	
	public double[][] untwist(double[][] twistedArray){
		double[][] untwistedArray = twistedArray;
		
		return untwistedArray;
	}	
}
