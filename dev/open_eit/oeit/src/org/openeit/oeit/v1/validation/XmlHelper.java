package org.openeit.oeit.v1.validation;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import javax.xml.XMLConstants;
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
	public static boolean isXmlValid(String xml, InputStream xsd) {
		try {
			Source source = new StreamSource(xsd);
			Document document = DocumentBuilderFactory.newInstance()
					.newDocumentBuilder().parse(xml);
			Schema schema = SchemaFactory.newInstance(
					XMLConstants.W3C_XML_SCHEMA_NS_URI).newSchema(source);
			schema.newValidator().validate(new DOMSource(document));
		} catch (SAXException | IOException | ParserConfigurationException e1) {
			return false;
		}

		return true;
	}

	public static boolean isXmlSetValid(Map<String, String> xmlset,
			InputStream xsd) {
		boolean result = !xmlset.isEmpty();

		for (String xml : xmlset.values())
			result &= isXmlValid(xml, xsd);

		return result;
	}
}
