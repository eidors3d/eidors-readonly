import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * @author Chengbo He
 */

/*
 * search code for "to do" for unfinished code.
 */

public class TestCases {

	public static final double VALID_THRESHOLD = 0.8;
	public static final double INVALID_THRESHOLD = 0.9;
	public static final double SKIP_THRESHOLD = 1.0;
	
	enum TestDataType {
		VALID, INVALID, SKIP
	}

	static public Map<String, Map<String, TestDataType>> cases = new ConcurrentHashMap<String, Map<String, TestDataType>>();
	static public List<String> rule_order = new ArrayList<String>();
	static public List<String> ruleList = new ArrayList<String>();
	public Document doc;
	static private Document source_doc;
	
	static {
		ConcurrentHashMap<String, TestDataType> map = new ConcurrentHashMap<String, TestDataType>();
		
		for (int i = 1; i<100; i++){
			for (int j = 1; j<100; j++){
				ruleList.add("" + i + "." + j);
			}
		}
		
		for (String rule:ruleList){
			double r = Math.random();
			TestDataType tt = TestDataType.SKIP;
			
			if(r < VALID_THRESHOLD)
				tt = TestDataType.VALID;
			else if(r < INVALID_THRESHOLD)
				tt = TestDataType.INVALID;
			
			map.put(rule, tt);
		}
		
		try {
			source_doc = Oeit.returnOeit();
		} catch (ParserConfigurationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		map.put("1.1", TestDataType.VALID);
		map.put("1.2", TestDataType.INVALID);
		cases.put("TestCase1", map);

		map = new ConcurrentHashMap<String, TestDataType>();

		map.put("1.1", TestDataType.VALID);
		map.put("1.2", TestDataType.SKIP);
		cases.put("TestCase2", map);
		
		rule_order.add("1.1");
		rule_order.add("1.2");
		rule_order.add("4.3");
		rule_order.add("1.3");
	}

	public void build(Map<String, TestDataType> rules)
			throws NoSuchMethodException, SecurityException,
			IllegalAccessException, IllegalArgumentException,
			InvocationTargetException, ParserConfigurationException {
		
		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
		
		doc = docBuilder.newDocument();

		doc.importNode(source_doc.getDocumentElement(), true);
		
		for (String r : rule_order) {
			String methodname = "";

			if (rules.containsKey(r)){

				switch (rules.get(r)){
				case VALID:
					methodname = "ruleAddValidDataFor" + r.replace('.', '_');
					break;
				case INVALID:
					methodname = "ruleAddInvalidDataFor" + r.replace('.', '_');
					break;
				case SKIP:
					methodname = "ruleSkip" + r.replace('.', '_');
					break;
				}
				Method m = this.getClass().getDeclaredMethod(methodname);
				m.invoke(this, null);
			}
		}

		// finalize_build_and_save();
	}

	//Rule 1.1
	public void ruleAddValidDataFor1_1() {
		
	}
	public void ruleAddInvalidDataFor1_1() {
		System.out.println("ruleAddInvalidDataFor1_1");
	}
	public void ruleSkip1_1() {
		System.out.println("ruleSkip1_1");
	}

	//Rule 1.2
	public void ruleAddValidDataFor1_2() {
		System.out.println("ruleAddValidDataFor1_2");
	}
	public void ruleAddInvalidDataFor1_2() {
		System.out.println("ruleAddInvalidDataFor1_2");
	}
	public void ruleSkip1_2() {
		System.out.println("ruleSkip1_2");
	}
	
	//Rule 3.1
	public void ruleAddValidDataFor3_1() {
		System.out.println("ruleAddValidDataFor3_1");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.setAttribute("id", "valid");
	}
	public void ruleAddInvalidDataFor3_1() {
		System.out.println("ruleAddInvalidDataFor3_1");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.setAttribute("id", "invalid");
	}
	public void ruleSkip3_1() {
		System.out.println("ruleSkip3_2");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.removeAttribute("id");
	}
	
	//Rule 3.2
	public void ruleAddValidDataFor3_2() {
		System.out.println("ruleAddValidDataFor3_1");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor3_2() {
		System.out.println("ruleAddInvalidDataFor3_2");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.setAttribute("someattribute", "invalid");
	}
	public void ruleSkip3_2() {
		System.out.println("ruleSkip3_2");
		
		System.out.println("Cannot skip rule 3.2");
	}
	
	//Rule 3.3
	public void ruleAddValidDataFor3_3() {
		System.out.println("ruleAddValidDataFor3_3");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor3_3() {
		System.out.println("ruleAddInvalidDataFor3_3");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		Element invalid = doc.createElement("stim_list");
		acquisition.appendChild(invalid);
	}
	public void ruleSkip3_3() {
		System.out.println("ruleSkip3_3");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.removeChild(acquisition.getFirstChild());
	}
	
	//Rule 3.4
	public void ruleAddValidDataFor3_4() {
		System.out.println("ruleAddValidDataFor3_4");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor3_4() {
		System.out.println("ruleAddInvalidDataFor3_4");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		Element invalid = doc.createElement("meas_list");
		acquisition.appendChild(invalid);
	}
	public void ruleSkip3_4() {
		System.out.println("ruleSkip3_4");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.removeChild(acquisition.getElementsByTagName("meas_list").item(0));
	}
	
	//Rule 3.5
	public void ruleAddValidDataFor3_5() {
		System.out.println("ruleAddValidDataFor3_5");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.removeChild(acquisition.getElementsByTagName("user_data").item(0));
	}
	public void ruleAddInvalidDataFor3_5() {
		System.out.println("ruleAddInvalidDataFor3_5");
		
		//Do nothing
	}
	public void ruleSkip3_5() {
		System.out.println("ruleSkip3_5");
		
		//Do nothing
	}
	
	//Rule 3.6
	public void ruleAddValidDataFor3_6() {
		System.out.println("ruleAddValidDataFor3_6");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor3_6() {
		System.out.println("ruleAddInvalidDataFor3_6");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		acquisition.appendChild(invalid);
	}
	public void ruleSkip3_6() {
		System.out.println("ruleSkip3_6");
		
		System.out.println("rule 3_6 cannot be skipped.");
	}
	
	//Rule 3.7
	public void ruleAddValidDataFor3_7() {
		System.out.println("ruleAddValidDataFor3_6");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor3_7() {
		System.out.println("ruleAddInvalidDataFor3_7");
		
		NodeList temp = doc.getElementsByTagName("acquisition");
		Element acquisition = (Element)temp.item(0);
		acquisition.appendChild(doc.createTextNode("some text"));
	}
	public void ruleSkip3_7() {
		System.out.println("ruleSkip3_7");
		
		System.out.println("rule 3_7 cannot be skipped.");
	}
	
	//Rule 4.1
	public void ruleAddValidDataFor4_1() {
		System.out.println("ruleAddValidDataFor4_1");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor4_1() {
		System.out.println("ruleAddInvalidDataFor4_1");
		
		NodeList temp = doc.getElementsByTagName("age");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip4_1() {
		System.out.println("ruleSkip4_1");
		
		System.out.println("rule 4_1 cannot be skipped");
	}
	
	//Rule 4.2
	public void ruleAddValidDataFor4_2() {
		System.out.println("ruleAddValidDataFor4_2");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor4_2() {
		System.out.println("ruleAddInvalidDataFor4_2");
		
		NodeList temp = doc.getElementsByTagName("age");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip4_2() {
		System.out.println("ruleSkip4_2");
		
		NodeList temp = doc.getElementsByTagName("age");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 4.3
	public void ruleAddValidDataFor4_3() {
		System.out.println("ruleAddValidDataFor4_3");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor4_3() {
		System.out.println("ruleAddInvalidDataFor4_3");
		
		NodeList temp = doc.getElementsByTagName("age");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip4_3() {
		System.out.println("ruleSkip4_3");
		
		System.out.println("rule 4_3 cannot be skipped");
	}
	
	//Rule 5.1
	public void ruleAddValidDataFor5_1() {
		System.out.println("ruleAddValidDataFor5_1");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor5_1() {
		System.out.println("ruleAddInvalidDataFor5_1");
		
		NodeList temp = doc.getElementsByTagName("contact_impedance");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip5_1() {
		System.out.println("ruleSkip5_1");
		
		System.out.println("rule 5_1 cannot be skipped");
	}
	
	//Rule 5.2
	public void ruleAddValidDataFor5_2(){
		System.out.println("ruleAddValidDataFor5_2");
		//Do nothing
	}
	public void ruleAddInvalidDataFor5_2(){
		System.out.println("ruleAddInvalidDataFor5_2");
		
		NodeList temp = doc.getElementsByTagName("contact_impedance");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip5_2(){
		System.out.println("ruleSkip5_2");
		
		NodeList temp = doc.getElementsByTagName("contact_impedance");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 5.3
	public void ruleAddValidDataFor5_3(){
		System.out.println("ruleAddValidDataFor5_3");
		//Do nothing
	}
	public void ruleAddInvalidDataFor5_3(){
		System.out.println("ruleAddInvalidDataFor5_3");
		
		NodeList temp = doc.getElementsByTagName("contact_impedance");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip5_3(){
		System.out.println("ruleSkip5_3");
		
		System.out.println("rule 5_3 cannot be skipped");
	}
	
	//Rule 6.1
	public void ruleAddValidDataFor6_1(){
		System.out.println("ruleAddValidDataFor6_1");
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_2(){
		System.out.println("ruleAddInvalidDataFor6_1");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip6_1(){
		System.out.println("ruleSkip6_1");
		
		System.out.println("rule 6_1 cannot be skipped");
	}
	
	//Rule 6.2
	public void ruleAddValidDataFor6_2(){
		System.out.println("ruleAddValidDataFor6_2");
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_1(){
		System.out.println("ruleAddInvalidDataFor6_2");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("name");
		element.appendChild(invalid);
	}
	public void ruleSkip6_2(){
		System.out.println("ruleSkip6_2");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 6.3
	public void ruleAddValidDataFor6_3(){
		System.out.println("ruleAddValidDataFor6_3");
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_3(){
		System.out.println("ruleAddInvalidDataFor6_3");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("version");
		element.appendChild(invalid);
	}
	public void ruleSkip6_3(){
		System.out.println("ruleSkip6_3");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("version").item(0));
	}
	
	//Rule 6.4
	public void ruleAddValidDataFor6_4(){
		System.out.println("ruleAddValidDataFor6_4");
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_4(){
		System.out.println("ruleAddInvalidDataFor6_4");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("manufacturer");
		element.appendChild(invalid);
	}
	public void ruleSkip6_4(){
		System.out.println("ruleSkip6_4");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("manufacturer").item(0));
	}
	
	//Rule 6.5
	public void ruleAddValidDataFor6_5(){
		System.out.println("ruleAddValidDataFor6_5");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor6_5(){
		System.out.println("ruleAddInvalidDataFor6_5");
		
		//Do nothing
	}
	public void ruleSkip6_5(){
		System.out.println("ruleSkip6_5");
		
		//Do nothing
	}
	
	//Rule 6.6
	public void ruleAddValidDataFor6_6(){
		System.out.println("ruleAddValidDataFor6_6");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_6(){
		System.out.println("ruleAddInvalidDataFor6_6");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip6_6(){
		System.out.println("ruleSkip6_6");
		
		System.out.println("rule 6_6 cannot be skipped");
	}
	
	//Rule 6.7
	public void ruleAddValidDataFor6_7(){
		System.out.println("ruleAddValidDataFor6_7");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor6_7(){
		System.out.println("ruleAddInvalidDataFor6_7");
		
		NodeList temp = doc.getElementsByTagName("data_acquisition_software");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip6_7(){
		System.out.println("ruleSkip6_7");
		
		System.out.println("rule 6_7 cannot be skipped");
	}
	
	//Rule 7.1
	public void ruleAddValidDataFor7_1(){
		System.out.println("ruleAddValidDataFor7_1");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor7_1(){
		System.out.println("ruleAddInvalidDataFor7_1");
		
		NodeList temp = doc.getElementsByTagName("decode");
		Element element = (Element)temp.item(0);
		element.setAttribute("frame","invalid_value");
	}
	public void ruleSkip7_1(){
		System.out.println("ruleSkip7_1");
		
		System.out.println("rule 7_1 cannot be skipped");
	}
	
	//Rule 7.2
	public void ruleAddValidDataFor7_2(){
		System.out.println("ruleAddValidDataFor7_2");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor7_2(){
		System.out.println("ruleAddInvalidDataFor7_2");
		
		NodeList temp = doc.getElementsByTagName("decode");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip7_2(){
		System.out.println("ruleSkip7_2");
		
		System.out.println("rule 7_2 cannot be skipped");
	}
	
	//Rule 7.3
	public void ruleAddValidDataFor7_3(){
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor7_3(){
		
		NodeList temp = doc.getElementsByTagName("decode");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip7_3(){
		
		System.out.println("rule 7_3 cannot be skipped");
	}
	
	//Rule 7.4
	public void ruleAddValidDataFor7_4(){
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor7_4(){
		
		NodeList temp = doc.getElementsByTagName("decode");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip7_4(){
		
		System.out.println("rule 7_4 cannot be skipped");
	}
	
	//Rule 8.1
	public void ruleAddValidDataFor8_1(){
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor8_1(){
		
		NodeList temp = doc.getElementsByTagName("description");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip8_1(){
		
		System.out.println("rule 8_1 cannot be skipped");
	}
	
	//Rule 8.2
	public void ruleAddValidDataFor8_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor8_2(){		
		NodeList temp = doc.getElementsByTagName("description");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
		element.appendChild(doc.createTextNode(")*&(*^^%$$$$123"));
	}
	public void ruleSkip8_2(){
		NodeList temp = doc.getElementsByTagName("description");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 8.3
	public void ruleAddValidDataFor8_3(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor8_3(){		
		NodeList temp = doc.getElementsByTagName("description");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip8_3(){
		System.out.println("rule 8_3 cannot be skipped");
	}
	
	//Rule 9.1
	public void ruleAddValidDataFor9_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_1(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip9_1(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 9.2
	public void ruleAddValidDataFor9_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_2(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip9_2(){
		System.out.println("rule 9_2 cannot be skipped");
	}
	
	//Rule 9.3
	public void ruleAddValidDataFor9_3(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_3(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("manufacturer");
		element.appendChild(invalid);
	}
	public void ruleSkip9_3(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("manufacturer").item(0));
	}
	
	//Rule 9.4
	public void ruleAddValidDataFor9_4(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_4(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("model");
		element.appendChild(invalid);
	}
	public void ruleSkip9_4(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("model").item(0));
	}
	
	//Rule 9.5
	public void ruleAddValidDataFor9_5(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_5(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("serial_number");
		element.appendChild(invalid);
	}
	public void ruleSkip9_5(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("serial_number").item(0));
	}
	
	//Rule 9.6
	public void ruleAddValidDataFor9_6(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_6(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("firmware");
		element.appendChild(invalid);
	}
	public void ruleSkip9_6(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("firmware").item(0));
	}
	
	//Rule 9.7
	public void ruleAddValidDataFor9_7(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_7(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("data_acquisition_software");
		element.appendChild(invalid);
	}
	public void ruleSkip9_7(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("data_acquisition_software").item(0));
	}
	
	//Rule 9.8
	public void ruleAddValidDataFor9_8(){
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor9_8(){		
		System.out.println("rule 9_8 cannot be invalid");
	}
	public void ruleSkip9_8(){
		//Do nothing
	}
	
	//Rule 9.9
	public void ruleAddValidDataFor9_9(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_9(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip9_9(){
		System.out.println("rule 9_9 cannot be skipped");
	}
	
	//Rule 9.10
	public void ruleAddValidDataFor9_10(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor9_10(){		
		NodeList temp = doc.getElementsByTagName("device");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip9_10(){
		System.out.println("rule 9_10 cannot be skipped");
	}
	
	//Rule 10.1
	public void ruleAddValidDataFor10_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor10_1(){		
		//To Do
	}
	public void ruleSkip10_1(){
		System.out.println("rule 10_1 cannot be skipped");
	}
	
	//Rule 10.2
	public void ruleAddValidDataFor10_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor10_2(){		
		NodeList temp = doc.getElementsByTagName("device_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip10_2(){
		System.out.println("rule 10_2 cannot be skipped");
	}
	
	//Rule 10.3
	public void ruleAddValidDataFor10_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor10_3(){		
		NodeList temp = doc.getElementsByTagName("device_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("device").item(0));
	}
	public void ruleSkip10_3(){
		System.out.println("rule 10_3 cannot be skipped");
	}
	
	//Rule 10.4
	public void ruleAddValidDataFor10_4(){
		NodeList temp = doc.getElementsByTagName("device_list");
		Element element = (Element)temp.item(0);
		Element device = doc.createElement("device");
		element.appendChild(device);
	}
	public void ruleAddInvalidDataFor10_4(){		
		NodeList temp = doc.getElementsByTagName("device_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("device").item(0));
	}
	public void ruleSkip10_4(){
		System.out.println("rule 10_4 cannot be skipped");
	}
	
	//Rule 10.5
	public void ruleAddValidDataFor10_5(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor10_5(){		
		NodeList temp = doc.getElementsByTagName("device_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip10_5(){
		System.out.println("rule 10_5 cannot be skipped");
	}
	
	//Rule 11.1
	public void ruleAddValidDataFor11_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_1(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip11_1(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 11.2
	public void ruleAddValidDataFor11_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_2(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("type","invalid");
	}
	public void ruleSkip11_2(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("type");
	}
	
	//Rule 11.3
	public void ruleAddValidDataFor11_3(){
		//to do
	}
	public void ruleAddInvalidDataFor11_3(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("position_cartesion","invalid");
		element.setAttribute("position_cylindrical","invalid");
		element.setAttribute("position_spherical","invalid");
		element.setAttribute("position_gps","invalid");
	}
	public void ruleSkip11_3(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("position_cartesion");
		element.removeAttribute("position_cylindrical");
		element.removeAttribute("position_spherical");
		element.removeAttribute("position_gps");
	}
	
	//Rule 11.4
	public void ruleAddValidDataFor11_4(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_4(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("orientation","invalid");
	}
	public void ruleSkip11_4(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("orientation");
	}
	
	//Rule 11.5
	public void ruleAddValidDataFor11_5(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_5(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("position_accuracy","invalid");
	}
	public void ruleSkip11_5(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("position_accuracy");
	}
	
	//Rule 11.6
	public void ruleAddValidDataFor11_6(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_6(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("contact_impedance","invalid");
	}
	public void ruleSkip11_6(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.removeAttribute("contact_impedance");
	}
	
	//Rule 11.7
	public void ruleAddValidDataFor11_7(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor11_7(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip11_7(){
		System.out.println("rule 11.7 cannot be skipped");
	}
	
	//Rule 11.8
	public void ruleAddValidDataFor11_8(){
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor11_8(){		
		//do nothing
	}
	public void ruleSkip11_8(){
		//do nothing
	}
	
	//Rule 11.9
	public void ruleAddValidDataFor11_9(){
		//do nothing
	}
	public void ruleAddInvalidDataFor11_9(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip11_9(){
		System.out.println("Rule 11.9 cannot be skipped");
	}
	
	//Rule 11.10
	public void ruleAddValidDataFor11_10(){
		//do nothing
	}
	public void ruleAddInvalidDataFor11_10(){		
		NodeList temp = doc.getElementsByTagName("electrode");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip11_10(){
		System.out.println("Rule 11.10 cannot be skipped");
	}
	
	//Rule 12.1
	public void ruleAddValidDataFor12_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_1(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip12_1(){
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 12.2
	public void ruleAddValidDataFor12_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_2(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("coordinate_system","invalid");
	}
	public void ruleSkip12_2(){
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("coordinate_system");
	}
	
	//Rule 12.3
	public void ruleAddValidDataFor12_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_3(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("position_description","invalid");
	}
	public void ruleSkip12_3(){
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("position_description");
	}
	
	//Rule 12.4
	public void ruleAddValidDataFor12_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_4(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip12_4(){
		System.out.println("rule 12.4 cannot be skipped");
	}
	
	//Rule 12.5
	public void ruleAddValidDataFor12_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_5(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("electrode").item(0));
	}
	public void ruleSkip12_5(){
		System.out.println("rule 12.5 cannot be skipped");
	}
	
	//Rule 12.6
	public void ruleAddValidDataFor12_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_6(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip12_6(){
		System.out.println("rule 12.6 cannot be skipped");
	}
	
	//Rule 12.7
	public void ruleAddValidDataFor12_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor12_7(){		
		NodeList temp = doc.getElementsByTagName("electrode_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip12_7(){
		System.out.println("rule 12.7 cannot be skipped");
	}
	
	//Rule 13.1
	public void ruleAddValidDataFor13_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_1(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip13_1(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 13.2
	public void ruleAddValidDataFor13_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_2(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip13_2(){
		System.out.println("rule 13.2 cannot be skipped");
	}
	
	//Rule 13.3
	public void ruleAddValidDataFor13_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_3(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("manufacturer");
		element.appendChild(invalid);
	}
	public void ruleSkip13_3(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("manufacturer").item(0));
	}
	
	//Rule 13.4
	public void ruleAddValidDataFor13_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_4(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("model");
		element.appendChild(invalid);
	}
	public void ruleSkip13_4(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("model").item(0));
	}
	
	//Rule 13.5
	public void ruleAddValidDataFor13_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_5(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("material");
		element.appendChild(invalid);
	}
	public void ruleSkip13_5(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("material").item(0));
	}
	
	//Rule 13.6
	public void ruleAddValidDataFor13_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_6(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("diameter");
		element.appendChild(invalid);
	}
	public void ruleSkip13_6(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("diameter").item(0));
	}
	
	//Rule 13.7
	public void ruleAddValidDataFor13_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_7(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("thickness");
		element.appendChild(invalid);
	}
	public void ruleSkip13_7(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("thickness").item(0));
	}
	
	//Rule 13.8
	public void ruleAddValidDataFor13_8(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_8(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("height");
		element.appendChild(invalid);
	}
	public void ruleSkip13_8(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("height").item(0));
	}
	
	//Rule 13.9
	public void ruleAddValidDataFor13_9(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_9(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("width");
		element.appendChild(invalid);
	}
	public void ruleSkip13_9(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("width").item(0));
	}
	
	//Rule 13.10
	public void ruleAddValidDataFor13_10(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_10(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("contact_impedance");
		element.appendChild(invalid);
	}
	public void ruleSkip13_10(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("contact_impedance").item(0));
	}
	
	//Rule 13.11
	public void ruleAddValidDataFor13_11(){
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor13_11(){		
		//Do nothing
	}
	public void ruleSkip13_11(){
		//Do nothing
	}
	
	//Rule 13.12
	public void ruleAddValidDataFor13_12(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_12(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip13_12(){
		System.out.println("rule 13.12 cannot be skipped");
	}
	
	//Rule 13.13
	public void ruleAddValidDataFor13_13(){
		//do nothing
	}
	public void ruleAddInvalidDataFor13_13(){		
		NodeList temp = doc.getElementsByTagName("electrode_type");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip13_13(){
		System.out.println("rule 13.13 cannot be skipped");
	}
	
	//Rule 14.1
	public void ruleAddValidDataFor14_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor14_1(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip14_1(){
		System.out.println("rule 14.1 cannot be skipped");
	}
	
	//Rule 14.2
	public void ruleAddValidDataFor14_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor14_2(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("electrode_type").item(0));
	}
	public void ruleSkip14_2(){
		System.out.println("rule 14.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 14.3
	public void ruleAddValidDataFor14_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor14_3(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip14_3(){
		System.out.println("rule 14.3 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 14.4
	public void ruleAddValidDataFor14_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor14_4(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip14_4(){
		System.out.println("rule 14.4 cannot be skipped");
	}
	
	//Rule 15.1
	public void ruleAddValidDataFor15_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_1(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("name", ":invalid:");
	}
	public void ruleSkip15_1(){
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
	}
	
	//Rule 15.2
	public void ruleAddValidDataFor15_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_2(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("type", "invalid");
	}
	public void ruleSkip15_2(){
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("type");
	}
	
	//Rule 15.3
	public void ruleAddValidDataFor15_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_3(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("time", "invalid");
	}
	public void ruleSkip15_3(){
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.removeAttribute("time");
	}
	
	//Rule 15.4
	public void ruleAddValidDataFor15_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_4(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("description");
		element.appendChild(invalid);
	}
	public void ruleSkip15_4(){
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("description").item(0));
	}
	
	//Rule 15.5
	public void ruleAddValidDataFor15_5(){
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor15_5(){		
		//do nothing
	}
	public void ruleSkip15_5(){
		//do nothing
	}
	
	//Rule 15.6
	public void ruleAddValidDataFor15_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_6(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip15_6(){
		System.out.println("rule 15.6 cannot be skipped");
	}
	
	//Rule 15.7
	public void ruleAddValidDataFor15_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor15_7(){		
		NodeList temp = doc.getElementsByTagName("electrode_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip15_7(){
		System.out.println("rule 15.7 cannot be skipped");
	}
	
	//Rule 16.1
	public void ruleAddValidDataFor16_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor16_1(){		
		NodeList temp = doc.getElementsByTagName("event_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip16_1(){
		System.out.println("rule 16.1 cannot be skipped");
	}
	
	//Rule 16.2
	public void ruleAddValidDataFor16_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor16_2(){		
		NodeList temp = doc.getElementsByTagName("event_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("log").item(0));
	}
	public void ruleSkip16_2(){
		System.out.println("rule 14.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 16.3
	public void ruleAddValidDataFor16_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor16_3(){		
		NodeList temp = doc.getElementsByTagName("event_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip16_3(){
		System.out.println("rule 16.3 cannot be skipped");
	}
	
	//Rule 16.4
	public void ruleAddValidDataFor16_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor16_4(){		
		NodeList temp = doc.getElementsByTagName("event_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip16_4(){
		System.out.println("rule 16.4 cannot be skipped");
	}
	
	//Rule 17.1
	public void ruleAddValidDataFor17_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor17_1(){		
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
		element.setAttribute("name", ":invalid:");
	}
	public void ruleSkip17_1(){
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
	}
	
	//Rule 17.2
	public void ruleAddValidDataFor17_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor17_2(){		
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.removeAttribute("type");
		element.setAttribute("type", "invalid");
	}
	public void ruleSkip17_2(){
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.removeAttribute("type");
	}
	
	//Rule 17.3
	public void ruleAddValidDataFor17_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor17_3(){		
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("description");
		element.appendChild(invalid);
	}
	public void ruleSkip17_3(){
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("description").item(0));
	}
	
	//Rule 17.4
	public void ruleAddValidDataFor17_4(){
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor17_4(){		
		//do nothing
	}
	public void ruleSkip17_4(){
		//do nothing
	}
	
	//Rule 17.5
	public void ruleAddValidDataFor17_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor17_5(){		
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip17_5(){
		System.out.println("rule 17.5 cannot be skipped");
	}
	
	//Rule 17.6
	public void ruleAddValidDataFor17_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor17_6(){		
		NodeList temp = doc.getElementsByTagName("event");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip17_6(){
		System.out.println("rule 17.6 cannot be skipped");
	}
	
	//Rule 18.1
	public void ruleAddValidDataFor18_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor18_1(){		
		NodeList temp = doc.getElementsByTagName("fields");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip18_1(){
		System.out.println("rule 16.1 cannot be skipped");
	}
	
	//Rule 18.2
	public void ruleAddValidDataFor18_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor18_2(){		
		NodeList temp = doc.getElementsByTagName("fields");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("field").item(0));
	}
	public void ruleSkip18_2(){
		System.out.println("rule 14.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 18.3
	public void ruleAddValidDataFor18_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor18_3(){		
		NodeList temp = doc.getElementsByTagName("fields");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip18_3(){
		System.out.println("rule 18.3 cannot be skipped");
	}
	
	//Rule 18.4
	public void ruleAddValidDataFor18_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor18_4(){		
		NodeList temp = doc.getElementsByTagName("fields");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip18_4(){
		System.out.println("rule 18.3 cannot be skipped");
	}
	
	//Rule 19.1
	public void ruleAddValidDataFor19_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_2(){
		
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip19_1(){
		System.out.println("rule 6_1 cannot be skipped");
	}
	
	//Rule 19.2
	public void ruleAddValidDataFor19_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_1(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("name");
		element.appendChild(invalid);
	}
	public void ruleSkip19_2(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 19.3
	public void ruleAddValidDataFor19_3(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_3(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("version");
		element.appendChild(invalid);
	}
	public void ruleSkip19_3(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("version").item(0));
	}
	
	//Rule 19.4
	public void ruleAddValidDataFor19_4(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_4(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("manufacturer");
		element.appendChild(invalid);
	}
	public void ruleSkip19_4(){
		
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("manufacturer").item(0));
	}
	
	//Rule 19.5
	public void ruleAddValidDataFor19_5(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor19_5(){
		//Do nothing
	}
	public void ruleSkip19_5(){
		//Do nothing
	}
	
	//Rule 19.6
	public void ruleAddValidDataFor19_6(){
		System.out.println("ruleAddValidDataFor6_6");
		
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_6(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip19_6(){
		System.out.println("rule 6_6 cannot be skipped");
	}
	
	//Rule 19.7
	public void ruleAddValidDataFor19_7(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor19_7(){
		NodeList temp = doc.getElementsByTagName("firmware");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip19_7(){
		System.out.println("rule 6.7 cannot be skipped");
	}
	
	//Rule 20.1
	public void ruleAddValidDataFor20_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_1(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
		element.setAttribute("name","invalid");
	}
	public void ruleSkip20_1(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
	}
	
	//Rule 20.2
	public void ruleAddValidDataFor20_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_2(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_frame");
		element.setAttribute("start_frame","invalid");
	}
	public void ruleSkip20_2(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_frame");
	}
	
	//Rule 20.3
	public void ruleAddValidDataFor20_3(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_3(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_frame");
		element.setAttribute("end_frame","invalid");
	}
	public void ruleSkip20_3(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_frame");
	}
	
	//Rule 20.4
	public void ruleAddValidDataFor20_4(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_4(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_time");
		element.setAttribute("start_time","invalid");
	}
	public void ruleSkip20_4(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_time");
	}
	
	//Rule 20.5
	public void ruleAddValidDataFor20_5(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_5(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_time");
		element.setAttribute("end_time","invalid");
	}
	public void ruleSkip20_5(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_time");
	}
	
	//Rule 20.6
	public void ruleAddValidDataFor20_6(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_6(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip20_6(){
		System.out.println("rule 20.6 cannot be skipped");
	}
	
	//Rule 20.7
	public void ruleAddValidDataFor20_7(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor20_7(){
		//Do nothing
	}
	public void ruleSkip20_7(){
		//Do nothing
	}
	
	//Rule 20.8
	public void ruleAddValidDataFor20_8(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_8(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip20_8(){
		System.out.println("rule 20.8 cannot be skipped");
	}
	
	//Rule 20.9
	public void ruleAddValidDataFor20_9(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor20_9(){
		NodeList temp = doc.getElementsByTagName("file");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip20_9(){
		System.out.println("rule 20.9 cannot be skipped");
	}
	
	//Rule 21.1
	public void ruleAddValidDataFor21_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor21_1(){		
		NodeList temp = doc.getElementsByTagName("files");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip21_1(){
		System.out.println("rule 16.1 cannot be skipped");
	}
	
	//Rule 21.2
	public void ruleAddValidDataFor21_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor21_2(){		
		NodeList temp = doc.getElementsByTagName("files");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("file").item(0));
	}
	public void ruleSkip21_2(){
		System.out.println("rule 21.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 21.3
	public void ruleAddValidDataFor21_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor21_3(){		
		NodeList temp = doc.getElementsByTagName("files");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip21_3(){
		System.out.println("rule 21.3 cannot be skipped");
	}
	
	//Rule 21.4
	public void ruleAddValidDataFor21_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor21_4(){		
		NodeList temp = doc.getElementsByTagName("files");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip21_4(){
		System.out.println("rule 21.4 cannot be skipped");
	}
	
	//Rule 22.1
	public void ruleAddValidDataFor22_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor22_1(){		
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip22_1(){
		System.out.println("rule 22.1 cannot be skipped");
	}
	
	//Rule 22.2
	public void ruleAddValidDataFor22_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor22_2(){		
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip22_2(){
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 22.3
	public void ruleAddValidDataFor22_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor22_3(){		
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip22_3(){
		System.out.println("rule 22.3 cannot be skipped");
	}
	
	//Rule 23.1
	public void ruleAddValidDataFor23_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor23_1(){		
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
		element.setAttribute("id", "invalid");
	}
	public void ruleSkip23_1(){
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 23.2
	public void ruleAddValidDataFor23_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor23_2(){	
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("acquisition").item(0));
	}
	public void ruleSkip23_2(){
		System.out.println("rule 23.2 cannot be skipped");
	}
	
	//Rule 23.3
	public void ruleAddValidDataFor23_3(){
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor23_3(){	
		//do nothing
	}
	public void ruleSkip23_3(){
		//do nothing
	}
	
	//Rule 23.4
	public void ruleAddValidDataFor23_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor23_4(){	
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip23_4(){
		System.out.println("rule 23.4 cannot be skipped");
	}
	
	//Rule 23.5
	public void ruleAddValidDataFor23_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor23_5(){	
		NodeList temp = doc.getElementsByTagName("frame_type");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip23_5(){
		System.out.println("rule 23.5 cannot be skipped");
	}
	
	//Rule 24.1
	public void ruleAddValidDataFor24_1(){
		//to do
	}
	public void ruleAddInvalidDataFor24_1(){	
		//to do
	}
	public void ruleSkip24_1(){
		//to do
	}
	
	//Rule 24.2
	public void ruleAddValidDataFor24_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor24_2(){	
		NodeList temp = doc.getElementsByTagName("frame_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip24_2(){
		System.out.println("rule 24.2 cannot be skipped");
	}
	
	//Rule 24.3
	public void ruleAddValidDataFor24_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor24_3(){	
		NodeList temp = doc.getElementsByTagName("frame_type_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("frame_type").item(0));
	}
	public void ruleSkip24_3(){
		System.out.println("rule 24.3 cannot be skipped");
	}
	
	//Rule 24.4
	public void ruleAddValidDataFor24_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor24_4(){	
		NodeList temp = doc.getElementsByTagName("frame_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip24_4(){
		System.out.println("rule 24.4 cannot be skipped");
	}
	
	//Rule 24.5
	public void ruleAddValidDataFor24_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor24_5(){	
		NodeList temp = doc.getElementsByTagName("frame_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip24_5(){
		System.out.println("rule 24.4 cannot be skipped");
	}
	
	//Rule 25.1
	public void ruleAddValidDataFor25_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor25_1(){		
		NodeList temp = doc.getElementsByTagName("height");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip25_1(){
		System.out.println("rule 22.1 cannot be skipped");
	}
	
	//Rule 25.2
	public void ruleAddValidDataFor25_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor25_2(){		
		NodeList temp = doc.getElementsByTagName("height");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("invalid"));
	}
	public void ruleSkip25_2(){
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 25.3
	public void ruleAddValidDataFor25_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor25_3(){		
		NodeList temp = doc.getElementsByTagName("height");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip25_3(){
		System.out.println("rule 22.3 cannot be skipped");
	}
	
	//Rule 26.1
	public void ruleAddValidDataFor26_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor26_1(){		
		NodeList temp = doc.getElementsByTagName("last_name");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip26_1(){
		System.out.println("rule 22.1 cannot be skipped");
	}
	
	//Rule 26.2
	public void ruleAddValidDataFor26_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor26_2(){		
		NodeList temp = doc.getElementsByTagName("last_name");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip26_2(){
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 26.3
	public void ruleAddValidDataFor26_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor26_3(){		
		NodeList temp = doc.getElementsByTagName("last_name");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip26_3(){
		System.out.println("rule 22.3 cannot be skipped");
	}
	
	//Rule 27.1
	public void ruleAddValidDataFor27_1(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_1(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
		element.setAttribute("name","invalid");
	}
	public void ruleSkip27_1(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("name");
	}
	
	//Rule 27.2
	public void ruleAddValidDataFor27_2(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_2(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_frame");
		element.setAttribute("start_frame","invalid");
	}
	public void ruleSkip27_2(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_frame");
	}
	
	//Rule 27.3
	public void ruleAddValidDataFor27_3(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_3(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_frame");
		element.setAttribute("end_frame","invalid");
	}
	public void ruleSkip27_3(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_frame");
	}
	
	//Rule 27.4
	public void ruleAddValidDataFor27_4(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_4(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_time");
		element.setAttribute("start_time","invalid");
	}
	public void ruleSkip27_4(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("start_time");
	}
	
	//Rule 27.5
	public void ruleAddValidDataFor27_5(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_5(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_time");
		element.setAttribute("end_time","invalid");
	}
	public void ruleSkip27_5(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.removeAttribute("end_time");
	}
	
	//Rule 27.6
	public void ruleAddValidDataFor27_6(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_6(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip27_6(){
		System.out.println("rule 27.6 cannot be skipped");
	}
	
	//Rule 27.7
	public void ruleAddValidDataFor27_7(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor27_7(){
		//Do nothing
	}
	public void ruleSkip27_7(){
		//Do nothing
	}
	
	//Rule 27.8
	public void ruleAddValidDataFor27_8(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_8(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip27_8(){
		System.out.println("rule 27.8 cannot be skipped");
	}
	
	//Rule 27.9
	public void ruleAddValidDataFor27_9(){
		//Do nothing
	}
	public void ruleAddInvalidDataFor27_9(){
		NodeList temp = doc.getElementsByTagName("log");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip27_9(){
		System.out.println("rule 27.9 cannot be skipped");
	}
	
	//Rule 28.1
	public void ruleAddValidDataFor28_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor28_1(){		
		NodeList temp = doc.getElementsByTagName("logs");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip28_1(){
		System.out.println("rule 28.1 cannot be skipped");
	}
	
	//Rule 28.2
	public void ruleAddValidDataFor28_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor28_2(){		
		NodeList temp = doc.getElementsByTagName("logs");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("log").item(0));
	}
	public void ruleSkip28_2(){
		System.out.println("rule 28.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 28.3
	public void ruleAddValidDataFor28_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor28_3(){		
		NodeList temp = doc.getElementsByTagName("logs");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip28_3(){
		System.out.println("rule 28.3 cannot be skipped");
	}
	
	//Rule 28.4
	public void ruleAddValidDataFor28_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor28_4(){		
		NodeList temp = doc.getElementsByTagName("logs");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip28_4(){
		System.out.println("rule 28.4 cannot be skipped");
	}
	
	//Rule 29.1
	public void ruleAddValidDataFor29_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor29_1(){		
		NodeList temp = doc.getElementsByTagName("manufacturer");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip29_1(){
		System.out.println("rule 22.1 cannot be skipped");
	}
	
	//Rule 29.2
	public void ruleAddValidDataFor29_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor29_2(){		
		NodeList temp = doc.getElementsByTagName("manufacturer");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip29_2(){
		NodeList temp = doc.getElementsByTagName("first_name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 29.3
	public void ruleAddValidDataFor29_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor29_3(){		
		NodeList temp = doc.getElementsByTagName("manufacturer");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip29_3(){
		System.out.println("rule 29.3 cannot be skipped");
	}
	
	//Rule 30.1
	public void ruleAddValidDataFor30_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor30_1(){		
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
		element.setAttribute("id", "invalid");
	}
	public void ruleSkip30_1(){
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 30.2
	public void ruleAddValidDataFor30_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor30_2(){	
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("elec").item(0));
	}
	public void ruleSkip30_2(){
		System.out.println("rule 30.2 cannot be skipped");
	}
	
	//Rule 30.3
	public void ruleAddValidDataFor30_3(){
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor30_3(){	
		//do nothing
	}
	public void ruleSkip30_3(){
		//do nothing
	}
	
	//Rule 30.4
	public void ruleAddValidDataFor30_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor30_4(){	
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip30_4(){
		System.out.println("rule 30.4 cannot be skipped");
	}
	
	//Rule 30.5
	public void ruleAddValidDataFor30_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor30_5(){	
		NodeList temp = doc.getElementsByTagName("meas");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip30_5(){
		System.out.println("rule 30.5 cannot be skipped");
	}
	
	//Rule 31.1
	public void ruleAddValidDataFor31_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor31_1(){		
		NodeList temp = doc.getElementsByTagName("meas_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip31_1(){
		System.out.println("rule 31.1 cannot be skipped");
	}
	
	//Rule 31.2
	public void ruleAddValidDataFor31_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor31_2(){		
		NodeList temp = doc.getElementsByTagName("meas_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("meas").item(0));
	}
	public void ruleSkip31_2(){
		System.out.println("rule 31.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 31.3
	public void ruleAddValidDataFor31_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor31_3(){		
		NodeList temp = doc.getElementsByTagName("meas_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip31_3(){
		System.out.println("rule 31.3 cannot be skipped");
	}
	
	//Rule 31.4
	public void ruleAddValidDataFo31_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor31_4(){		
		NodeList temp = doc.getElementsByTagName("meas_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip31_4(){
		System.out.println("rule 31.4 cannot be skipped");
	}
	
	//Rule 32.1
	public void ruleAddValidDataFo32_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_1(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
		element.setAttribute("id","invalid");
	}
	public void ruleSkip32_1(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 32.2
	public void ruleAddValidDataFo32_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_2(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("frequency");
		element.setAttribute("frequency","invalid");
	}
	public void ruleSkip32_2(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("frequency");
	}
	
	//Rule 32.3
	public void ruleAddValidDataFo32_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_3(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("physical_property");
		element.setAttribute("physical_property","invalid");
	}
	public void ruleSkip32_3(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("physical_property");
	}
	
	//Rule 32.4
	public void ruleAddValidDataFo32_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_4(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("multiplier");
		element.setAttribute("multiplier","invalid");
	}
	public void ruleSkip32_4(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("multiplier");
	}
	
	//Rule 32.5
	public void ruleAddValidDataFo32_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_5(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("offset");
		element.setAttribute("offset","invalid");
	}
	public void ruleSkip32_5(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("offset");
	}
	
	//Rule 32.6
	public void ruleAddValidDataFo32_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_6(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip32_6(){
		System.out.println("Rule 32.6 cannot be skipped");
	}
	
	//Rule 32.7
	public void ruleAddValidDataFo32_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_7(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("fields");
		element.appendChild(invalid);
	}
	public void ruleSkip32_7(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("fields").item(0));
	}
	
	//Rule 32.8
	public void ruleAddValidDataFor32_8(){
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor32_8(){	
		//do nothing
	}
	public void ruleSkip32_8(){
		//do nothing
	}
	
	//Rule 32.9
	public void ruleAddValidDataFo32_9(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_9(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip32_9(){
		System.out.println("Rule 32.9 cannot be skipped");
	}
	
	//Rule 32.10
	public void ruleAddValidDataFo32_10(){
		//do nothing
	}
	public void ruleAddInvalidDataFor32_10(){		
		NodeList temp = doc.getElementsByTagName("meas_type");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip32_10(){
		System.out.println("Rule 32.10 cannot be skipped");
	}
	
	//Rule 33.1
	public void ruleAddValidDataFor33_1(){
		//to do
	}
	public void ruleAddInvalidDataFor33_1(){	
		//to do
	}
	public void ruleSkip33_1(){
		//to do
	}
	
	//Rule 33.2
	public void ruleAddValidDataFor33_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor33_2(){	
		NodeList temp = doc.getElementsByTagName("meas_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip33_2(){
		System.out.println("rule 33.2 cannot be skipped");
	}
	
	//Rule 33.3
	public void ruleAddValidDataFor33_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor33_3(){	
		NodeList temp = doc.getElementsByTagName("meas_type_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("frame_type").item(0));
	}
	public void ruleSkip33_3(){
		System.out.println("rule 33.3 cannot be skipped");
	}
	
	//Rule 33.4
	public void ruleAddValidDataFor33_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor33_4(){	
		NodeList temp = doc.getElementsByTagName("meas_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip33_4(){
		System.out.println("rule 24.4 cannot be skipped");
	}
	
	//Rule 33.5
	public void ruleAddValidDataFor33_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor33_5(){	
		NodeList temp = doc.getElementsByTagName("meas_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip33_5(){
		System.out.println("rule 33.4 cannot be skipped");
	}
	
	//Rule 34.1
	public void ruleAddValidDataFor34_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor34_1(){		
		NodeList temp = doc.getElementsByTagName("model");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip34_1(){
		System.out.println("rule 34.1 cannot be skipped");
	}
	
	//Rule 34.2
	public void ruleAddValidDataFor34_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor34_2(){		
		NodeList temp = doc.getElementsByTagName("model");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip34_2(){
		NodeList temp = doc.getElementsByTagName("model");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 34.3
	public void ruleAddValidDataFor34_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor34_3(){		
		NodeList temp = doc.getElementsByTagName("model");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip34_3(){
		System.out.println("rule 34.3 cannot be skipped");
	}
	
	//Rule 35.1
	public void ruleAddValidDataFor35_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor35_1(){		
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip35_1(){
		System.out.println("rule 35.1 cannot be skipped");
	}
	
	//Rule 35.2
	public void ruleAddValidDataFor35_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor35_2(){		
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip35_2(){
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 35.3
	public void ruleAddValidDataFor35_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor35_3(){
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip35_3(){
		System.out.println("rule 35.3 cannot be skipped");
	}
	
	//Rule 36.1
	public void ruleAddValidDataFor36_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_1(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.setAttribute("xmlns","invalid");
	}
	public void ruleSkip36_1(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeAttribute("xmlns");
	}
	
	//Rule 36.2
	public void ruleAddValidDataFor36_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_2(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.setAttribute("xmlns:xi","invalid");
	}
	public void ruleSkip36_2(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeAttribute("xmlns:xi");
	}
	
	//Rule 36.3
	public void ruleAddValidDataFor36_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_3(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip36_3(){
		System.out.println("rule 36.3 cannot be skipped.");
	}
	
	//Rule 36.4
	public void ruleAddValidDataFor36_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_4(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("subject");
		element.appendChild(invalid);
	}
	public void ruleSkip36_4(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("subject").item(0));
	}
	
	//Rule 36.5
	public void ruleAddValidDataFor36_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_5(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("device_list");
		element.appendChild(invalid);
	}
	public void ruleSkip36_5(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("device_list").item(0));
	}
	
	//Rule 36.6
	public void ruleAddValidDataFor36_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_6(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("electrode_type_list");
		element.appendChild(invalid);
	}
	public void ruleSkip36_6(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("electrode_type_list").item(0));
	}
	
	//Rule 36.7
	public void ruleAddValidDataFor36_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_7(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("electrode_list");
		element.appendChild(invalid);
	}
	public void ruleSkip36_7(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("electrode_list").item(0));
	}
	
	//Rule 36.8
	public void ruleAddValidDataFor36_8(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_8(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("stim_type_list");
		element.appendChild(invalid);
	}
	public void ruleSkip36_8(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("stim_type_list").item(0));
	}
	
	//Rule 36.9
	public void ruleAddValidDataFor39_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor39_4(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("meas_type_list");
		element.appendChild(invalid);
	}
	public void ruleSkip39_4(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("meas_type_list").item(0));
	}
	
	//Rule 36.10
	public void ruleAddValidDataFor36_10(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_10(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("frame_type_list");
		element.appendChild(invalid);
	}
	public void ruleSkip36_10(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("frame_type_list").item(0));
	}
	
	//Rule 36.11
	public void ruleAddValidDataFor36_11(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_11(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("streams");
		element.appendChild(invalid);
	}
	public void ruleSkip36_11(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("streams").item(0));
	}
	
	//Rule 36.12
	public void ruleAddValidDataFor36_12(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_12(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip36_12(){
		System.out.println("rule 36.12 cannot be skipped");
	}
	
	//Rule 36.13
	public void ruleAddValidDataFor36_13(){
		//do nothing
	}
	public void ruleAddInvalidDataFor36_13(){
		NodeList temp = doc.getElementsByTagName("oeit");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip36_13(){
		System.out.println("rule 36.13 cannot be skipped");
	}
	
	//Rule 37.1
	public void ruleAddValidDataFor37_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_1(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip37_1(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 37.2
	public void ruleAddValidDataFor37_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_2(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip37_2(){
		System.out.println("rule 37.2 cannot be skipped");
	}
	
	//Rule 37.3
	public void ruleAddValidDataFor37_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_3(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("first_name");
		element.appendChild(invalid);
	}
	public void ruleSkip37_3(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.removeChild(doc.getElementsByTagName("first_name").item(0));
	}
	
	//Rule 37.4
	public void ruleAddValidDataFor37_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_4(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("last_name");
		element.appendChild(invalid);
	}
	public void ruleSkip37_4(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.removeChild(doc.getElementsByTagName("last_name").item(0));
	}
	
	//Rule 37.5
	public void ruleAddValidDataFor37_5(){
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor37_5(){	
		//do nothing
	}
	public void ruleSkip37_5(){
		//do nothing
	}
	
	//Rule 37.6
	public void ruleAddValidDataFo37_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_6(){		
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip37_6(){
		System.out.println("Rule 32.9 cannot be skipped");
	}
	
	//Rule 37.7
	public void ruleAddValidDataFo37_7(){
		//do nothing
	}
	public void ruleAddInvalidDataFor37_7(){		
		NodeList temp = doc.getElementsByTagName("patient");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip37_7(){
		System.out.println("Rule 32.7 cannot be skipped");
	}
	
	//Rule 38.1
	public void ruleAddValidDataFor38_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_1(){
		NodeList temp = doc.getElementsByTagName("process");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip38_1(){
		NodeList temp = doc.getElementsByTagName("process");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 38.2
	public void ruleAddValidDataFor38_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_2(){
		NodeList temp = doc.getElementsByTagName("process");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip38_2(){
		System.out.println("rule 38.2 cannot be skipped");
	}
	
	//Rule 38.3
	public void ruleAddValidDataFor38_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_3(){
		//to do
	}
	public void ruleSkip38_3(){
		//to do
	}
	
	//Rule 38.4
	public void ruleAddValidDataFor38_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_4(){
		//to do
	}
	public void ruleSkip38_4(){
		//to do
	}
	
	//Rule 38.5
	public void ruleAddValidDataFor38_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_5(){
		NodeList temp = doc.getElementsByTagName("process");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip38_5(){
		System.out.println("rule 38.5 cannot be skipped");
	}
	
	//Rule 38.6
	public void ruleAddValidDataFor38_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor38_6(){
		NodeList temp = doc.getElementsByTagName("process");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip38_6(){
		System.out.println("rule 38.6 cannot be skipped");
	}
	
	//Rule 39.1
	public void ruleAddValidDataFor39_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor39_1(){		
		NodeList temp = doc.getElementsByTagName("radius");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip39_1(){
		System.out.println("rule 39.1 cannot be skipped");
	}
	
	//Rule 39.2
	public void ruleAddValidDataFor39_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor39_2(){		
		NodeList temp = doc.getElementsByTagName("radius");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip39_2(){
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 39.3
	public void ruleAddValidDataFor39_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor39_3(){
		NodeList temp = doc.getElementsByTagName("radius");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip39_3(){
		System.out.println("rule 39.3 cannot be skipped");
	}
	
	//Rule 40.1
	public void ruleAddValidDataFor40_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_1(){
		NodeList temp = doc.getElementsByTagName("repeat");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip40_1(){
		NodeList temp = doc.getElementsByTagName("repeat");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 40.2
	public void ruleAddValidDataFor40_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_2(){
		NodeList temp = doc.getElementsByTagName("repeat");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip40_2(){
		System.out.println("rule 40.2 cannot be skipped");
	}
	
	//Rule 40.3
	public void ruleAddValidDataFor40_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_3(){
		//to do
	}
	public void ruleSkip40_3(){
		//to do
	}
	
	//Rule 40.4
	public void ruleAddValidDataFor40_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_4(){
		//to do
	}
	public void ruleSkip40_4(){
		//to do
	}
	
	//Rule 40.5
	public void ruleAddValidDataFor40_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_5(){
		NodeList temp = doc.getElementsByTagName("repeat");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip40_5(){
		System.out.println("rule 40.5 cannot be skipped");
	}
	
	//Rule 40.6
	public void ruleAddValidDataFor40_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor40_6(){
		NodeList temp = doc.getElementsByTagName("repeat");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip40_6(){
		System.out.println("rule 40.6 cannot be skipped");
	}
	
	//Rule 41.1
	public void ruleAddValidDataFor41_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor41_1(){		
		NodeList temp = doc.getElementsByTagName("serial_number");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip41_1(){
		System.out.println("rule 41.1 cannot be skipped");
	}
	
	//Rule 41.2
	public void ruleAddValidDataFor41_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor41_2(){		
		NodeList temp = doc.getElementsByTagName("serial_number");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip41_2(){
		NodeList temp = doc.getElementsByTagName("name");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 41.3
	public void ruleAddValidDataFor41_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor41_3(){
		NodeList temp = doc.getElementsByTagName("serial_number");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip41_3(){
		System.out.println("rule 41.3 cannot be skipped");
	}
	
	//Rule 42.1
	public void ruleAddValidDataFor42_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor42_1(){		
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
		element.setAttribute("id", "invalid");
	}
	public void ruleSkip42_1(){
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 42.2
	public void ruleAddValidDataFor42_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor42_2(){	
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("elec").item(0));
	}
	public void ruleSkip42_2(){
		System.out.println("rule 42.2 cannot be skipped");
	}
	
	//Rule 42.3
	public void ruleAddValidDataFor42_3(){
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
	}
	public void ruleAddInvalidDataFor42_3(){	
		//do nothing
	}
	public void ruleSkip42_3(){
		//do nothing
	}
	
	//Rule 42.4
	public void ruleAddValidDataFor42_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor42_4(){	
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip42_4(){
		System.out.println("rule 42.4 cannot be skipped");
	}
	
	//Rule 42.5
	public void ruleAddValidDataFor42_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor42_5(){	
		NodeList temp = doc.getElementsByTagName("stim");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip42_5(){
		System.out.println("rule 42.5 cannot be skipped");
	}
	
	//Rule 43.1
	public void ruleAddValidDataFor43_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor43_1(){		
		NodeList temp = doc.getElementsByTagName("stim_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip43_1(){
		System.out.println("rule 43.1 cannot be skipped");
	}
	
	//Rule 43.2
	public void ruleAddValidDataFor43_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor43_2(){		
		NodeList temp = doc.getElementsByTagName("stim_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("stim").item(0));
	}
	public void ruleSkip43_2(){
		System.out.println("rule 43.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 43.3
	public void ruleAddValidDataFor43_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor43_3(){		
		NodeList temp = doc.getElementsByTagName("stim_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip43_3(){
		System.out.println("rule 43.3 cannot be skipped");
	}
	
	//Rule 43.4
	public void ruleAddValidDataFor43_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor43_4(){		
		NodeList temp = doc.getElementsByTagName("stim_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip43_4(){
		System.out.println("rule 43.4 cannot be skipped");
	}
	
	//Rule 44.1
	public void ruleAddValidDataFor44_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_1(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip44_1(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 44.2
	public void ruleAddValidDataFor44_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_2(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("amplitude","invalid");
	}
	public void ruleSkip44_2(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("amplitude");
	}
	
	//Rule 44.3
	public void ruleAddValidDataFor44_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_3(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("frequency","invalid");
	}
	public void ruleSkip44_3(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("frequency");
	}
	
	//Rule 44.4
	public void ruleAddValidDataFor44_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_4(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("wave","invalid");
	}
	public void ruleSkip44_4(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("wave");
	}
	
	//Rule 44.5
	public void ruleAddValidDataFor44_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_5(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("physical_property","invalid");
	}
	public void ruleSkip44_5(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.removeAttribute("physical_property");
	}
	
	//Rule 44.6
	public void ruleAddValidDataFor44_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_6(){		
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip44_6(){
		System.out.println("rule 44.6 cannot be skipped");
	}
	
	//Rule 44.7
	public void ruleAddValidDataFor44_7(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor44_7(){
		//Do nothing
	}
	public void ruleSkip44_7(){
		//Do nothing
	}
	
	//Rule 44.8
	public void ruleAddValidDataFor44_8(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_8(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip44_8(){
		System.out.println("rule 44.8 cannot be skipped");
	}
	
	//Rule 44.9
	public void ruleAddValidDataFor44_9(){
		//do nothing
	}
	public void ruleAddInvalidDataFor44_9(){
		NodeList temp = doc.getElementsByTagName("stim_type");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip44_9(){
		System.out.println("rule 44.9 cannot be skipped");
	}
	
	//Rule 45.1
	public void ruleAddValidDataFor45_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor45_1(){
		NodeList temp = doc.getElementsByTagName("stim_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip45_1(){
		System.out.println("rule 45.1 cannot be skipped");
	}
	
	//Rule 45.2
	public void ruleAddValidDataFor45_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor45_2(){
		NodeList temp = doc.getElementsByTagName("stim_type_list");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip45_2(){
		System.out.println("rule 45.2 cannot be skipped");
	}
	
	//Rule 45.3
	public void ruleAddValidDataFor45_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor45_3(){
		NodeList temp = doc.getElementsByTagName("stim_type_list");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("stim_type").item(0));
	}
	public void ruleSkip45_3(){
		System.out.println("rule 45.3 cannot be skipped");
	}
	
	//Rule 45.4
	public void ruleAddValidDataFor45_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor45_4(){
		NodeList temp = doc.getElementsByTagName("stim_type_list");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("invalid");
		element.appendChild(invalid);
	}
	public void ruleSkip45_4(){
		System.out.println("rule 45.4 cannot be skipped");
	}
	
	//Rule 45.5
	public void ruleAddValidDataFor45_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor45_5(){
		NodeList temp = doc.getElementsByTagName("stim_type_list");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip45_5(){
		System.out.println("rule 45.5 cannot be skipped");
	}
	
	//Rule 46.1
	public void ruleAddValidDataFor46_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_1(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip46_1(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 46.2
	public void ruleAddValidDataFor46_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_2(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip46_2(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 46.3
	public void ruleAddValidDataFor46_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_3(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip46_3(){
		System.out.println("Rule 46.3 cannot be skipped");
	}
	
	//Rule 46.4
	public void ruleAddValidDataFor46_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_4(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("process");
		element.appendChild(invalid);
	}
	public void ruleSkip46_4(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("process").item(0));
	}
	
	//Rule 46.5
	public void ruleAddValidDataFor46_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_5(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("process");
		element.appendChild(invalid);
	}
	public void ruleSkip46_5(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("files").item(0));
	}
	
	//Rule 46.6
	public void ruleAddValidDataFor46_6(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_6(){
		NodeList temp = doc.getElementsByTagName("files");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("process");
		element.appendChild(invalid);
	}
	public void ruleSkip46_6(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("files").item(0));
	}
	
	//Rule 46.7
	public void ruleAddValidDataFor46_7(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		Element user_data = doc.createElement("user_data");
		element.appendChild(user_data);
		
	}
	public void ruleAddInvalidDataFor46_7(){
		//Do nothing
	}
	public void ruleSkip46_7(){
		//Do nothing
	}
	
	//Rule 46.8
	public void ruleAddValidDataFor46_8(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_8(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip46_8(){
		System.out.println("rule 44.8 cannot be skipped");
	}
	
	//Rule 46.9
	public void ruleAddValidDataFor46_9(){
		//do nothing
	}
	public void ruleAddInvalidDataFor46_9(){
		NodeList temp = doc.getElementsByTagName("stream");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip46_9(){
		System.out.println("rule 46.9 cannot be skipped");
	}
	
	//Rule 47.1
	public void ruleAddValidDataFor47_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor47_1(){		
		NodeList temp = doc.getElementsByTagName("streams");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute", "some_value");
	}
	public void ruleSkip47_1(){
		System.out.println("rule 47.1 cannot be skipped");
	}
	
	//Rule 47.2
	public void ruleAddValidDataFor47_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor47_2(){		
		NodeList temp = doc.getElementsByTagName("streams");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getElementsByTagName("stream").item(0));
	}
	public void ruleSkip47_2(){
		System.out.println("rule 47.2 cannot be skipped (already skipped in invalid)");
	}
	
	//Rule 47.3
	public void ruleAddValidDataFor47_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor47_3(){		
		NodeList temp = doc.getElementsByTagName("streams");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("stream");
		element.appendChild(invalid);
	}
	public void ruleSkip47_3(){
		System.out.println("rule 47.3 cannot be skipped");
	}
	
	//Rule 47.4
	public void ruleAddValidDataFor47_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor47_4(){		
		NodeList temp = doc.getElementsByTagName("streams");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("stream"));
	}
	public void ruleSkip47_4(){
		System.out.println("rule 47.4 cannot be skipped");
	}
	
	//Rule 48.1
	public void ruleAddValidDataFor48_1(){
		//to do
	}
	public void ruleAddInvalidDataFor48_1(){		
		//to do
	}
	public void ruleSkip48_1(){
		//to do
	}
	
	//Rule 48.2
	public void ruleAddValidDataFor48_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor48_2(){		
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.setAttribute("id","invalid");
	}
	public void ruleSkip48_2(){
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.removeAttribute("id");
	}
	
	//Rule 48.3
	public void ruleAddValidDataFor48_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor48_3(){		
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.setAttribute("type","invalid");
	}
	public void ruleSkip48_3(){
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.removeAttribute("type");
	}
	
	//Rule 48.4
	public void ruleAddValidDataFor48_4(){
		//do nothing
	}
	public void ruleAddInvalidDataFor48_4(){		
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip48_4(){
		System.out.println("rule 48.4 cannot be skipped");
	}
	
	//Rule 48.5
	public void ruleAddValidDataFor48_5(){
		//do nothing
	}
	public void ruleAddInvalidDataFor48_5(){		
		NodeList temp = doc.getElementsByTagName("subject");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("some_text"));
	}
	public void ruleSkip48_5(){
		System.out.println("rule 48.5 cannot be skipped");
	}
	
	//Rule 48.6
	public void ruleAddValidDataFor48_6(){
		//do do
	}
	public void ruleAddInvalidDataFor48_6(){		
		//to do
	}
	public void ruleSkip48_6(){
		//to do
	}
	
	//Rule 49.1
	public void ruleAddValidDataFor49_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor49_1(){		
		NodeList temp = doc.getElementsByTagName("thickness");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip49_1(){
		System.out.println("rule 49.1 cannot be skipped");
	}
	
	//Rule 49.2
	public void ruleAddValidDataFor49_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor49_2(){		
		NodeList temp = doc.getElementsByTagName("thickness");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode("invalid"));
	}
	public void ruleSkip49_2(){
		System.out.println("rule 49.2 cannot be skipped");
	}
	
	//RULE 50 TO DO
	
	//Rule 51.1
	public void ruleAddValidDataFor51_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor51_1(){		
		NodeList temp = doc.getElementsByTagName("version");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip51_1(){
		System.out.println("rule 51.1 cannot be skipped");
	}
	
	//Rule 51.2
	public void ruleAddValidDataFor51_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor51_2(){		
		NodeList temp = doc.getElementsByTagName("version");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip51_2(){
		NodeList temp = doc.getElementsByTagName("version");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 51.3
	public void ruleAddValidDataFor51_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor51_3(){		
		NodeList temp = doc.getElementsByTagName("version");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip51_3(){
		System.out.println("rule 51.3 cannot be skipped");
	}
	
	//Rule 52.1
	public void ruleAddValidDataFor52_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor52_1(){		
		NodeList temp = doc.getElementsByTagName("water_level");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip52_1(){
		System.out.println("rule 52.1 cannot be skipped");
	}
	
	//Rule 52.2
	public void ruleAddValidDataFor52_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor52_2(){		
		NodeList temp = doc.getElementsByTagName("water_level");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip52_2(){
		NodeList temp = doc.getElementsByTagName("water_level");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 52.3
	public void ruleAddValidDataFor52_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor52_3(){		
		NodeList temp = doc.getElementsByTagName("water_level");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip52_3(){
		System.out.println("rule 51.3 cannot be skipped");
	}
	
	//Rule 53.1
	public void ruleAddValidDataFor53_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor53_1(){		
		NodeList temp = doc.getElementsByTagName("water_conductivity");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip53_1(){
		System.out.println("rule 52.1 cannot be skipped");
	}
	
	//Rule 53.2
	public void ruleAddValidDataFor53_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor53_2(){		
		NodeList temp = doc.getElementsByTagName("water_conductivity");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip53_2(){
		NodeList temp = doc.getElementsByTagName("water_conductivity");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 53.3
	public void ruleAddValidDataFor53_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor53_3(){		
		NodeList temp = doc.getElementsByTagName("water_conductivity");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip53_3(){
		System.out.println("rule 53.3 cannot be skipped");
	}
	
	//Rule 54.1
	public void ruleAddValidDataFor54_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor54_1(){		
		NodeList temp = doc.getElementsByTagName("width");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip54_1(){
		System.out.println("rule 54.1 cannot be skipped");
	}
	
	//Rule 54.2
	public void ruleAddValidDataFor54_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor54_2(){		
		NodeList temp = doc.getElementsByTagName("width");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip54_2(){
		NodeList temp = doc.getElementsByTagName("width");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 54.3
	public void ruleAddValidDataFor54_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor54_3(){		
		NodeList temp = doc.getElementsByTagName("width");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip54_3(){
		System.out.println("rule 54.3 cannot be skipped");
	}
	
	//Rule 55.1
	public void ruleAddValidDataFor55_1(){
		//do nothing
	}
	public void ruleAddInvalidDataFor55_1(){		
		NodeList temp = doc.getElementsByTagName("weight");
		Element element = (Element)temp.item(0);
		element.setAttribute("some_attribute","some_value");
	}
	public void ruleSkip55_1(){
		System.out.println("rule 55.1 cannot be skipped");
	}
	
	//Rule 54.2
	public void ruleAddValidDataFor55_2(){
		//do nothing
	}
	public void ruleAddInvalidDataFor55_2(){		
		NodeList temp = doc.getElementsByTagName("weight");
		Element element = (Element)temp.item(0);
		element.appendChild(doc.createTextNode(":invalid:"));
	}
	public void ruleSkip55_2(){
		NodeList temp = doc.getElementsByTagName("weight");
		Element element = (Element)temp.item(0);
		element.removeChild(element.getFirstChild());
	}
	
	//Rule 54.3
	public void ruleAddValidDataFor55_3(){
		//do nothing
	}
	public void ruleAddInvalidDataFor55_3(){		
		NodeList temp = doc.getElementsByTagName("weight");
		Element element = (Element)temp.item(0);
		Element invalid = doc.createElement("some_element");
		element.appendChild(invalid);
	}
	public void ruleSkip55_3(){
		System.out.println("rule 55.3 cannot be skipped");
	}
}
