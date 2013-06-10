import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class MCEITParser implements OeitLegacyParser {

	public EitEntity parse(String fileName) throws IOException{
		EitEntity mceitInfo = new EitEntity();
		
		File file = new File(fileName);
		
		int length = (int)(file.length()/4);
		
		if (length%256!=0){
			System.out.println("File length strange - cropping file");
			length = length-length%256;
		}
		
        /*
         * Read in source file here
         */
		FileInputStream fin = new FileInputStream(file);
		BufferedInputStream bin = new BufferedInputStream(fin);
		DataInputStream din = new DataInputStream(bin);
		
		double[][] input = new double[256][length/256];
		byte[] bytes = new byte[4];
		int counter = 0;

		for (int i = 0; i<256; i++){
			//System.out.print(counter + ": ");
			for (int j = 0; j<length/256; j++){            //METHOD 1
				bytes[0] = din.readByte();
				bytes[1] = din.readByte();
				bytes[2] = din.readByte();
				bytes[3] = din.readByte();
				input[i][j] = (double)ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN ).getFloat();
				counter++;
				//System.out.print(" " + input[i][j]);
			}
			//System.out.print("\n");
		}
		
        din.close();
        bin.close();
        fin.close();
        
        double[][] injectionData = new double[208][length/256];
        double[][] current = new double[14][length/256];
        double[][] voltage = new double[14][length/256];
        double[][] aux = new double[14][length/256];
		
        /*
         * Read in Config file here, using DOM
         */
		getConfigInfo();
		
        //input = untwist(input);
        mceitInfo.setOutput2D(input);
        mceitInfo.setLength(length);
        /*        if (input[0][38]==0){
    	//104 measurements per frame
    	//write method to transfer 104 measurement array to 208 measurement array
    	input = transfer104To208(input);
    }*/
        
		return mceitInfo;
	}
	
	public void getConfigInfo(){
		Document doc = XMLReader.returnDOM("sample1.xml");
		
		NodeList inputs = doc.getElementsByTagName("inputs");
		NodeList in = ((Element)inputs.item(0)).getElementsByTagName("stim");
		
		NodeList outputs = doc.getElementsByTagName("outputs");
		NodeList out = ((Element)outputs.item(0)).getElementsByTagName("meas");
		
		NodeList stimTypes = doc.getElementsByTagName("stim_type");
		NodeList measTypes = doc.getElementsByTagName("meas_type");
		
		String stimType = null;
		String measType = null;
		
		for (int i = 0; i < in.getLength(); i++){
			System.out.println("Electrode_number: " + ((Element)in.item(i)).getAttribute("elec")
					+ " Gain: " + ((Element)in.item(i)).getAttribute("gain"));
			stimType = ((Element)in.item(i)).getAttribute("type");
			System.out.println("Stim_Info: ");
			for (int j = 0; j < stimTypes.getLength(); j++){
				if (stimType.equals(((Element)stimTypes.item(j)).getAttribute("name"))){
					System.out.println("Amplitude:" + ((Element)stimTypes.item(j)).getAttribute("amplitude")
							+ " Frequency:" + ((Element)stimTypes.item(j)).getAttribute("frequency")
							+ " Physical_Property:" + ((Element)stimTypes.item(j)).getAttribute("physical_property")
							+ " Wave:" + ((Element)stimTypes.item(j)).getAttribute("wave"));
				}
			}
			measType = ((Element)out.item(i)).getAttribute("type");
			System.out.println("Output_info:");
			for (int j = 0; j < measTypes.getLength(); j++){
				if (measType.equals(((Element)measTypes.item(j)).getAttribute("name"))){
					System.out.println("Demod_frequency:" + ((Element)measTypes.item(j)).getAttribute("demod_frequency")
							+ " Physical_Property:" + ((Element)measTypes.item(j)).getAttribute("physical_property")
							+ " Offset_gain:" + ((Element)measTypes.item(j)).getAttribute("offset_gain")
							+ " Signal_gain:" + ((Element)measTypes.item(j)).getAttribute("signal_gain"));
				}
				Node fields = ((Element)measTypes.item(j)).getElementsByTagName("fields").item(0);
				NodeList fieldList = ((Element)fields).getElementsByTagName("field");
				
				System.out.println("Fields:");
				for (int k = 0; k < fieldList.getLength(); k++){
					System.out.println("Name:" + ((Element)fieldList.item(k)).getAttribute("name")
							+ " Type:" + ((Element)fieldList.item(k)).getAttribute("type")
							+ " Length:" + ((Element)fieldList.item(k)).getAttribute("length"));
				}
			}
			System.out.println("------------------------------------------------");
		}
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
