package org.openeit.oeit.v1.validation;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Source;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

public class XmlHelper {
	public static void validate(String xml, InputStream xsd)
			throws ParserConfigurationException, SAXException, IOException {
		Source source = new StreamSource(xsd);
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setXIncludeAware(true);
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document document = builder.parse(xml);

		Schema schema = SchemaFactory.newInstance(
				XMLConstants.W3C_XML_SCHEMA_NS_URI).newSchema(source);
		schema.newValidator().validate(new DOMSource(document));
	}

	public static void validateSet(Map<String, String> xmlset, InputStream xsd)
			throws ParserConfigurationException, SAXException, IOException {

		for (String xml : xmlset.values())
			validate(xml, xsd);
	}
}
