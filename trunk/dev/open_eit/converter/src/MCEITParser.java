import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
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
		
		int[][] input = new int[256][length/256];
		byte[] bytes = new byte[4];
		int counter = 0;
		
		

		for (int i = 0; i<256; i++){
			//System.out.print(counter + ": ");
			for (int j = 0; j<length/256; j++){            //METHOD 1
				bytes[0] = din.readByte();
				bytes[1] = din.readByte();
				bytes[2] = din.readByte();
				bytes[3] = din.readByte();
				input[i][j] = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN ).getInt();
				counter++;
				//System.out.print(" " + input[i][j]);
			}
			//System.out.print("\n");
		}
		
        din.close();
        bin.close();
        fin.close();
        
        double[] data = new double[256*(length/256)];
		
/*        for (int i = 0; i< 256*(length/256); i++){
    		if (i%14==0 && i<208*(length/256)){
    			data[i]= input[208+i/16][i%length];
    		}else if (i%15==0 && i<224*(length/256)){
    			data[i] = input[224+i/16][i%length];
    		}else if (i%16==0 && i<240*(length/256)){
    			data[i] = input[240+i/16][i%length];
    		}else{
    			data[i] = input[i/length][i%length];
    		}
        }*/
 
        
        /*
         * Read in Config file here, using DOM
         */
        
		getData();
		
        //input = untwist(input);
        mceitInfo.setLength(length);
        /*        if (input[0][38]==0){
    	//104 measurements per frame
    	//write method to transfer 104 measurement array to 208 measurement array
    	input = transfer104To208(input);
    }*/
        
		return mceitInfo;
	}
	
	public void getData(){
		
		//double[] data = new double[256*(length/256)];
	
		Document doc = XMLReader.returnDOM("sample1.xml");
		NodeList measTasks = doc.getElementsByTagName("meas_task");
		
		NodeList measTypeList = doc.getElementsByTagName("meas_type");
		NodeList electrodeList = doc.getElementsByTagName("electrode");
		
		String stimType = null;
		String measType = null;
		
		for (int i = 0; i<measTasks.getLength(); i++){
			Element measTaskE = (Element) measTasks.item(i);
			
			String start = null, stop = null;
			
			start = measTaskE.getAttribute("start");
			stop = measTaskE.getAttribute("stop");
			
			NodeList meas = measTaskE.getElementsByTagName("meas");
			NodeList stim = measTaskE.getElementsByTagName("stim");
			
			NodeList stimTypeList = doc.getElementsByTagName("stim_type");
			
			String elec1 = null, elec2 = null;
			
			String amplitude=null, frequency=null, physical_property=null, wave=null;
			
			Element electrode1 = (Element) stim.item(0);
			Element electrode2 = (Element) stim.item(1);
			stimType = electrode1.getAttribute("type");
			
			elec1 = electrode1.getAttribute("elec");
			elec2 = electrode2.getAttribute("elec");
			
			for (int k = 0; k<stimTypeList.getLength(); k++){
				Element stimulation = (Element)stimTypeList.item(k);
				if (stimType.equals(stimulation.getAttribute("name"))){
					amplitude = stimulation.getAttribute("amplitude");
					frequency = stimulation.getAttribute("frequency");
					physical_property = stimulation.getAttribute("physical_property");
					wave = stimulation.getAttribute("wave");
				}
			}
			
			for (int j = 0; j<meas.getLength(); j++){
				Element measE = (Element) meas.item(j);
				measType = measE.getAttribute("type");
				NodeList elect = measE.getElementsByTagName("e");
				String firstOrdinal = ((Element)elect.item(0)).getAttribute("ordinal");
				String secondOrdinal = ((Element)elect.item(1)).getAttribute("ordinal");
				
				String name1 = null, type1 = null, 
						position1 = null, name2 = null, 
						type2 = null, position2 = null;
				
				for (int k = 0; k<electrodeList.getLength(); k++){
					Element electrode = (Element)electrodeList.item(k);
					if (firstOrdinal.equals(electrode.getAttribute("ordinal"))){
						name1 = electrode.getAttribute("name");
						type1 = electrode.getAttribute("type");
						position1 = electrode.getAttribute("position");
					}
					if (secondOrdinal.equals(electrode.getAttribute("ordinal"))){
						name2 = electrode.getAttribute("name");
						type2 = electrode.getAttribute("type");
						position2 = electrode.getAttribute("position");
					}
				}
				
				String demod_freq = null, physical_property2 = null, 
						offset_gain = null, signal_gain = null;
				
				for (int k = 0; k<measTypeList.getLength(); k++){
					Element measurement = (Element)measTypeList.item(k);
					if (measType.equals(measurement.getAttribute("name"))){
						demod_freq = measurement.getAttribute("demod_freq");
						physical_property2 = measurement.getAttribute("physical_property");
						offset_gain = measurement.getAttribute("offset_gain");
						signal_gain = measurement.getAttribute("signal_gain");
					}
				}
				System.out.println("<<TASK>>" 
						+ " Start:" + start
						+ " Stop:" + stop
						+" <<<STIMULATION>>>"
						+ "Elects:" + elec1
						+ " " + elec2
						+ " type:" + stimType + "("
						+ "amplitude:" + amplitude
						+ " frequency:" + frequency
						+ " physical_property:" + physical_property
						+ " wave" + wave + ")"
						+ " <<<MEASUREMENT>>>" 
						+ "Ordinals:" + firstOrdinal + "("
						+ "name:" + name1
						+ " type:" + type1
						+ " position:" + position1 + ")"
						+ " " + secondOrdinal + "("
						+ "name:" + name2
						+ " type:" + type2
						+ " position:" + position2 + ")"
						+ " type:" + measType  + "("
						+ "demod_freq:" + demod_freq
						+ " physical_property:" + physical_property2
						+ " offset_gain:" + offset_gain
						+ " signal_gain:" + signal_gain + ")"
						+ "   <measurement here>");
			}
		}
		
		//return data;
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
