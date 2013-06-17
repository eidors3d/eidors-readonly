import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import java.io.File;

public class XMLReader {
	
	/*
	 * This method reads a config xml file and returns a document object model (DOM)
	 */
	public Document returnDOM(String filename){
		
		try{
			File fXmlFile = new File(filename);
			
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(fXmlFile);
			
			doc.getDocumentElement().normalize(); //optional
			
			return doc;
		}catch (Exception e){
			System.out.println(filename + " not found.");
			return null;
		}
	}
}
