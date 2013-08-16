import java.io.File;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;


public class Oeit {

	public static void main(String[] args){
		TransformerFactory transformerFactory = TransformerFactory.newInstance();
		
		try {
			Transformer transformer = transformerFactory.newTransformer();
			DOMSource source = new DOMSource(returnOeit());
			StreamResult result = new StreamResult(new File("sample.xml"));
			StreamResult resultOnConsole = new StreamResult(System.out);
			transformer.transform(source, result);
			transformer.transform(source, resultOnConsole);
			System.out.println("file saved.");
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TransformerConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TransformerException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static Document returnOeit () throws ParserConfigurationException{
		//int subjectType = 0; // 0: person  1: tank  2: volcano
		
		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
		
		Document doc = docBuilder.newDocument();
		
		//ROOT
		Element root = doc.createElement("oeit");
		doc.appendChild(root);
		root.setAttribute("xmlns", "http://org.open-eit/schemas/oeit.v0.9");
		root.setAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
		
		//SUBJECT BRANCH
		Element subject = doc.createElement("subject");
		subject.setAttribute("id", "subject");
		
		Element patient = doc.createElement("patient");
			//PATIENT BRANCH
			patient.setAttribute("id", "1");
			
			Element first_name = doc.createElement("first_name");
			first_name.appendChild(doc.createTextNode("text"));
			Element last_name = doc.createElement("last_name");
			last_name.appendChild(doc.createTextNode("text"));
			patient.appendChild(first_name);
			patient.appendChild(last_name);
			//PATIENT BRANCH END
		Element age = doc.createElement("age");
		age.appendChild(doc.createTextNode("0.0"));
		Element weight = doc.createElement("weight");
		weight.appendChild(doc.createTextNode("0.0 g"));
		Element height = doc.createElement("height");
		height.appendChild(doc.createTextNode("0.0 m"));
		Element radius = doc.createElement("radius");
		radius.appendChild(doc.createTextNode("0.0 m"));
		Element water_level = doc.createElement("water_level");
		water_level.appendChild(doc.createTextNode("0.0 m"));
		Element water_conductivity = doc.createElement("water_conductivity");
		water_conductivity.appendChild(doc.createTextNode("0.0 S/m"));
		Element manufacturer = doc.createElement("manufacturer");
		manufacturer.appendChild(doc.createTextNode("text"));
		
		subject.appendChild(patient);
		subject.appendChild(age);
		subject.appendChild(weight);
		subject.appendChild(height);
		subject.appendChild(radius);
		subject.appendChild(water_level);
		subject.appendChild(water_conductivity);
		subject.appendChild(manufacturer);
		//SUBJECT BRANCH END
		
		
		//DEVICE_LIST BRANCH
		Element device_list = doc.createElement("device_list");
		Element device = doc.createElement("device");
			//DEVICE BRANCH
			device.setAttribute("id", "0");
			Element model = doc.createElement("model");
			model.appendChild(doc.createTextNode("text"));
			Element serial_number = doc.createElement("serial_number");
			serial_number.appendChild(doc.createTextNode("text"));
			Element firmware = doc.createElement("firmware");
			//FIRMWARE BRANCH
				Element name = doc.createElement("name");
				name.appendChild(doc.createTextNode("text"));
				Element version = doc.createElement("version");
				version.appendChild(doc.createTextNode("text"));
				
				firmware.appendChild(name);
				firmware.appendChild(version);
				firmware.appendChild(manufacturer);
			//FIRMWARE BRANCH END
			Element data_acquisition_software = doc.createElement("data_acquisition_software");
			data_acquisition_software.appendChild(doc.createTextNode("text"));
			
			device.appendChild(manufacturer);
			device.appendChild(model);
			device.appendChild(serial_number);
			device.appendChild(firmware);
			device.appendChild(data_acquisition_software);
			//DEVICE BRANCH END
		device_list.appendChild(device);
		//DEVICE_LIST BRANCH END
		
		//ELECTRODE_TYPE_LIST BRANCH
		Element electrode_type_list = doc.createElement("electrode_type_list");
			Element electrode_type = doc.createElement("electrode_type");
			electrode_type.setAttribute("id", "electrode_type");
			//ELECTRODE_TYPE BRANCH
			Element material = doc.createElement("material");
			material.appendChild(doc.createTextNode("text"));
			Element diameter = doc.createElement("diameter");
			diameter.appendChild(doc.createTextNode("0.0 m"));
			Element thickness = doc.createElement("thickness");
			thickness.appendChild(doc.createTextNode("0.0 m"));
			Element width = doc.createElement("width");
			width.appendChild(doc.createTextNode("0.0 m"));
			Element contact_impedance = doc.createElement("contact_impedence");
			
			electrode_type.appendChild(manufacturer);//defined previously
			electrode_type.appendChild(model);//defined previously
			electrode_type.appendChild(material);
			electrode_type.appendChild(diameter);
			electrode_type.appendChild(thickness);
			electrode_type.appendChild(height);
			electrode_type.appendChild(width);
			electrode_type.appendChild(contact_impedance);
			
		electrode_type_list.appendChild(electrode_type);
		//ELECTRODE_TYPE_LIST BRANCH END
		
		//ELECTRODE_LIST BRANCH
		Element electrode_list = doc.createElement("electrode_list");
		electrode_list.setAttribute("id", "electrode_list");
		electrode_list.setAttribute("coordinate_system","Cartesian");
		electrode_list.setAttribute("position_description","Absolute");
		
			Element electrode = doc.createElement("electrode");
			//ELECTRODE BRANCH
			electrode.setAttribute("id", "elec1");
			electrode.setAttribute("type", "elect_type1");
			electrode.setAttribute("position", "[0,0,0]");
			electrode.setAttribute("orientation", "");//missing
			electrode.setAttribute("position_accuracy", "0.0 m");
			electrode.setAttribute("contact_impedance", "");//missing
			//ELECTRODE BRANCH END
			
		electrode_list.appendChild(electrode);
		//ELECTRODE LIST END
		
		//STIM_TYPE_LIST BRANCH
		Element stim_type_list = doc.createElement("stim_type_list");
		Element stim_type = doc.createElement("stim_type");
		
		stim_type.setAttribute("id","stim1");
		stim_type.setAttribute("amplitude","0.0 A");
		stim_type.setAttribute("frequency","0.0 Hz");
		stim_type.setAttribute("wave","0.0");
		stim_type.setAttribute("physical_property","property");
		
		stim_type_list.appendChild(stim_type);
		//STIM_TYPE_LIST BRANCH END
		
		//MEAS_TYPE_LIST BRANCH
		Element meas_type_list = doc.createElement("meas_type_list");
		Element meas_type = doc.createElement("meas_type");
		
		meas_type.setAttribute("id", "meas_type");
		meas_type.setAttribute("frequency", "0.0 Hz");//missing
		meas_type.setAttribute("physical_property", "");//missing
		meas_type.setAttribute("multiplier", "0.0");
		meas_type.setAttribute("offset", "0.0");
			
			//FIELDS FRANCH
			Element fields = doc.createElement("fields");
			Element field = doc.createElement("field");
			
			field.setAttribute("name","Amplitude");
			field.setAttribute("type","uint64");
			
			fields.appendChild(field);
			meas_type.appendChild(fields);
			//FIELDS BRANCH END
		
		meas_type_list.appendChild(meas_type);
		//MEAS_TYPE_LIST BRANCH END
		
		//FRAME_TYPE_LIST BRANCH
		Element frame_type_list = doc.createElement("frame_type_list");
		Element frame_type = doc.createElement("frame_type");
		
			//ACQUISITION LIST BRANCH
			Element acquisition = doc.createElement("acquisition");
			acquisition.setAttribute("id","acquisition");
				
				//STIM_LIST BRANCH
				Element stim_list = doc.createElement("stim_list");
				Element stim = doc.createElement("stim");
				stim.setAttribute("id","stim");
				Element elec = doc.createElement("elec");
				
				stim.appendChild(elec);
				stim_list.appendChild(stim);
						
				//STIM_LIST BRANCH END
				
				//MEAS_LIST BRANCH
				Element meas_list = doc.createElement("meas_list");
				Element meas = doc.createElement("meas");
				meas.setAttribute("id","stim");
				
				meas.appendChild(elec);
				meas_list.appendChild(stim);
				//MEAS_LIST BRANCH END
				
				acquisition.appendChild(stim_list);
				acquisition.appendChild(meas_list);
			frame_type.appendChild(acquisition);
			//ACQUISITION LIST BRANCH END
			
		frame_type_list.appendChild(frame_type);
		//FRAME_TYPE_LIST BRANCH END

		//STREAM BRANCH
		Element streams = doc.createElement("streams");
		streams.setAttribute("id", "streams");
		streams.setAttribute("device", "");//missing
		
			Element process = doc.createElement("process");
			//PROCESS BRANCH
			process.setAttribute("id","process");
			
			Element repeat = doc.createElement("repeat");
			repeat.setAttribute("number","0");
			
			Element decode = doc.createElement("decode");
			decode.setAttribute("frame", ""); //missing
			
			process.appendChild(repeat);
			process.appendChild(decode);
			//PROCESS BRANCH END
			
			//FILES BRANCH
			Element files = doc.createElement("files");
			Element file = doc.createElement("file");
			
			file.setAttribute("start_frame","0");
			file.setAttribute("end_frame","0");
			file.setAttribute("start_time","0");
			file.setAttribute("end_time","0");
			
			files.appendChild(file);
			//FILES BRANCH END
			
			//LOGS BRANCH
			Element logs = doc.createElement("logs");
			
			Element log = doc.createElement("log");
			log.setAttribute("start_frame","0");
			log.setAttribute("end_frame","0");
			log.setAttribute("start_time","0");
			log.setAttribute("end_time","0");
			
			logs.appendChild(log);
			//LOGS BRANCH END
			
		streams.appendChild(process);
		streams.appendChild(files);
		streams.appendChild(logs);
		//STREAM BRANCH END
		
		root.appendChild(subject);
		root.appendChild(electrode_type_list);
		root.appendChild(electrode_list);
		root.appendChild(device_list);
		root.appendChild(stim_type_list);
		root.appendChild(meas_type_list);
		root.appendChild(frame_type_list);
		root.appendChild(streams);
		
		//TEST AREA
		//TEST AREA END
		
		return doc;
	}

}
