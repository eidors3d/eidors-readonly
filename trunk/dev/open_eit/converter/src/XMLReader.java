import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;

public class XMLReader {

	public static Document readXML(String filename){
		
		try{
			File fXmlFile = new File(filename);
			
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(fXmlFile);
			
			doc.getDocumentElement().normalize(); //optional
			
			return doc;
		}catch (Exception e){
			return null;
		}

	}
	
	public static void main(String argv[]) {
		readXML("sample1.xml");
	}
}
