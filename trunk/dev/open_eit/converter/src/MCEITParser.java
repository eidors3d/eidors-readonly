import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * MCEITParser
 * Chengbo He June 2013
 * Carleton University
 */

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
         * Read in .get file here
         */
		FileInputStream fin = new FileInputStream(file);
		BufferedInputStream bin = new BufferedInputStream(fin);
		DataInputStream din = new DataInputStream(bin);
		
		int[][] input = new int[256][length/256];
		byte[] bytes = new byte[4];
		
		for (int i = 0; i<256; i++){
			for (int j = 0; j<length/256; j++){
				bytes[0] = din.readByte();
				bytes[1] = din.readByte();
				bytes[2] = din.readByte();
				bytes[3] = din.readByte();
				input[i][j] = ByteBuffer.wrap(bytes).order(ByteOrder.BIG_ENDIAN ).getInt();
			}
		}
		
        din.close();
        bin.close();
        fin.close();
        
        int[] data = new int[length];
        
        /*
         *  Rearranges data and reshapes to a 1D array
         */
        int a = 0, b = 0;
        for (int i = 0; i< length-1; i++){
    		if ((i+3)%16==0 && a<16*(length/256)){
    			data[i]= input[208+a/(length/256)][a%(length/256)];
    			data[i+1] = input[224+a/(length/256)][a%(length/256)];
    			data[i+2] = input[240+a/(length/256)][a%(length/256)];
    			a++;
    			i=i+2;
    		}else{
    			data[i] = input[b/(length/256)][b%(length/256)];
    			b++;
    		}
        }
        
        /*
         * Read in XML file here, using DOM and prints data
         */
		showData(data);
		writeToOeit(data,length,fileName);
		
        //input = untwist(input);
        mceitInfo.setLength(length);
        /*        if (input[0][38]==0){
    	//104 measurements per frame
    	//write method to transfer 104 measurement array to 208 measurement array
    	input = transfer104To208(input);
    }*/
        
		return mceitInfo;
	}
	
	public void showData(int[] data){
		
		//double[] data = new double[256*(length/256)];
	
		//Get NodeLists here
		XMLReader reader = new XMLReader();
		Document doc = reader.returnDOM("sample.xml"); // Returns DOM
		NodeList measTasks = doc.getElementsByTagName("acquisition");
		NodeList measTypeList = doc.getElementsByTagName("meas_type");
		NodeList electrodeList = doc.getElementsByTagName("electrode");
		NodeList stimTypeList = doc.getElementsByTagName("stim_type");
		
		String stimType = null;
		String measType = null;
		
		for (int i = 0; i<measTasks.getLength(); i++){
			Element measTaskE = (Element) measTasks.item(i);
			String start = null, stop = null;
			
			start = measTaskE.getAttribute("start");
			stop = measTaskE.getAttribute("stop");
			
			NodeList meas = measTaskE.getElementsByTagName("meas");
			NodeList stim = measTaskE.getElementsByTagName("stim");
			
			/*
			 * --STIM--
			 */
			Element stimE = (Element)stim.item(0);
			stimType = stimE.getAttribute("type");
			NodeList elect0 = stimE.getElementsByTagName("elec");
			String firstOrdinal0 = ((Element)elect0.item(0)).getAttribute("id");
			String secondOrdinal0 = ((Element)elect0.item(1)).getAttribute("id");
			String firstGain0 = ((Element)elect0.item(0)).getAttribute("gain");
			String secondGain0 = ((Element)elect0.item(1)).getAttribute("gain");
			
			/*
			 * GET STIM ORDINAL INFO
			 */
			String name01 = null, type01 = null, position01 = null, 
					name02 = null, type02 = null, position02 = null;

			for (int k = 0; k<electrodeList.getLength(); k++){
				Element electrode = (Element)electrodeList.item(k);
				if (firstOrdinal0.equals(electrode.getAttribute("id"))){
					name01 = electrode.getAttribute("name");
					type01 = electrode.getAttribute("type");
					position01 = electrode.getAttribute("position");
				}
				if (secondOrdinal0.equals(electrode.getAttribute("id"))){
					name02 = electrode.getAttribute("name");
					type02 = electrode.getAttribute("type");
					position02 = electrode.getAttribute("position");
				}
			}
				
			/*
			 * GET STIM TYPES INFO
			 */
			String amplitude=null, frequency=null, physical_property=null, wave=null;
			
			for (int k = 0; k<stimTypeList.getLength(); k++){
				Element stimulation = (Element)stimTypeList.item(k);
				if (stimType.equals(stimulation.getAttribute("name"))){
					amplitude = stimulation.getAttribute("amplitude");
					frequency = stimulation.getAttribute("frequency");
					physical_property = stimulation.getAttribute("physical_property");
					wave = stimulation.getAttribute("wave");
				}
			}
			
			/*
			 * --MEAS--
			 */
			Element measE = (Element) meas.item(0);
			measType = measE.getAttribute("type");
			NodeList elect = measE.getElementsByTagName("elec");
			String firstOrdinal = ((Element)elect.item(0)).getAttribute("id");
			String secondOrdinal = ((Element)elect.item(1)).getAttribute("id");
			String firstGain = ((Element)elect.item(0)).getAttribute("gain");
			String secondGain = ((Element)elect.item(1)).getAttribute("gain");
			
			/*
			 * GET MEAS ORDINAL INFO
			 */
			String name1 = null, type1 = null, position1 = null, 
					name2 = null, type2 = null, position2 = null;

			for (int k = 0; k<electrodeList.getLength(); k++){
				Element electrode = (Element)electrodeList.item(k);
				if (firstOrdinal.equals(electrode.getAttribute("id"))){
					name1 = electrode.getAttribute("name");
					type1 = electrode.getAttribute("type");
					position1 = electrode.getAttribute("position");
				}
				if (secondOrdinal.equals(electrode.getAttribute("id"))){
					name2 = electrode.getAttribute("name");
					type2 = electrode.getAttribute("type");
					position2 = electrode.getAttribute("position");
				}
			}
			
			String demod_freq = null, physical_property2 = null, 
					meas_gain = null;
			/*
			 * GET MEAS TYPES INFO
			 */
			int j = 0;
			for (int k = 0; k<measTypeList.getLength(); k++){
				Element measurement = (Element)measTypeList.item(k);
				if (measType.equals(measurement.getAttribute("name"))){
					demod_freq = measurement.getAttribute("demod_frequency");
					physical_property2 = measurement.getAttribute("physical_property");
					meas_gain = measurement.getAttribute("gain");
					j++;
				}
			}
			Element measurement = (Element)measTypeList.item(j);
			/*
			 * GET FIELD INFO
			 */
			NodeList fieldList = measurement.getElementsByTagName("field");
			String[] fieldName = new String[fieldList.getLength()];
			String[] fieldType = new String[fieldList.getLength()];
			
			for (int k = 0; k<fieldList.getLength(); k++){
				fieldName[k] = ((Element)fieldList.item(k)).getAttribute("name");
				fieldType[k] = ((Element)fieldList.item(k)).getAttribute("type");
			}
			
			/*
			 * PRINT
			 */
			System.out.println("<<TASK>> " 
					+ "Start:" + start
					+ " Stop:" + stop
					+" <<<STIMULATION>>> "
					+ "Ordinals:" + firstOrdinal0 + "{"
					+ "name:" + name01
					+ " gain:" + firstGain0
					+ " type:" + type01
					+ " position:" + position01 + "}"
					+ " " + secondOrdinal0 + "{"
					+ "name:" + name02
					+ " gain:" + secondGain0
					+ " type:" + type02
					+ " position:" + position02 + "}"
					+ " type:" + stimType + "{"
					+ "amplitude:" + amplitude
					+ " frequency:" + frequency
					+ " physical_property:" + physical_property
					+ " wave:" + wave + "}"
					+ " <<<MEASUREMENT>>> " 
					+ "Ordinals:" + firstOrdinal + "{"
					+ "name:" + name1
					+ " gain:" + firstGain
					+ " type:" + type1
					+ " position:" + position1 + "}"
					+ " " + secondOrdinal + "{"
					+ "name:" + name2
					+ " gain:" + secondGain
					+ " type:" + type2
					+ " position:" + position2 + "}"
					+ " type:" + measType  + "{"
					+ "demod_freq:" + demod_freq
					+ " physical_property:" + physical_property2
					+ " measurement_type_gain:" + meas_gain + "}"
					+ " <<<FIELDS>>> "
					+ "name:" + fieldName[0]
					+ " type:" + fieldType[0]
					+ " <<<VALUE>>> "
					+ "raw:0x" + Integer.toHexString(0x1000000 | data[i]).substring(1)
					+ " interpreted:" + data[i]);
		}
		//return data;
	}
	public void writeToOeit(int data[], int length, String fileName) throws IOException{
		String oeitFileName = fileName.substring(0,fileName.length()-4);	
		
		DataOutputStream os = new DataOutputStream(new FileOutputStream(oeitFileName + ".oeit"));
			
			for (int i = 0; i < length; i++){
				os.writeInt(data[i]);
			}
			
			os.close();
	}
	
	public double[][] transfer104To208(double[][] array104){ //Unnecessary right now
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
	
	public double[][] untwist(double[][] twistedArray){ //Unnecessary right now
		double[][] untwistedArray = twistedArray;
		
		return untwistedArray;
	}
}
