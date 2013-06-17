package org.openeit.oeit.v1.validation;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.Assume.assumeNotNull;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import org.openeit.oeit.Oeit;

@RunWith(Parameterized.class)
public class OeitFileTest {

	private static Map<FilePattern, String> __filePatterns = null;
	private ZipFileHelper _zfh = null;

	@Parameters
	public static Collection<Object[]> generateData() {
		List<Object[]> p = new ArrayList<Object[]>();

		for (String s : Oeit.getInputFiles()) {
			p.add(new Object[] { s });
		}

		return p;
	}

	public OeitFileTest(String file) {
		_zfh = new ZipFileHelper(file);
	}

	@BeforeClass
	public static void setUpClass() {
		Map<FilePattern, String> map = new ConcurrentHashMap<FilePattern, String>();
		map.put(FilePattern.VERSION, "^" + Pattern.quote("version.txt") + "$");
		map.put(FilePattern.HEADER_DIR, "^" + Pattern.quote("header") + "$");
		map.put(FilePattern.DEVICE, "^" + Pattern.quote("header/device.xml")
				+ "$");
		map.put(FilePattern.ELECTRODE,
				"^" + Pattern.quote("header/electrode.xml") + "$");
		map.put(FilePattern.SUBJECT, "^" + Pattern.quote("header/subject.xml")
				+ "$");
		map.put(FilePattern.PATIENT, "^" + Pattern.quote("header/patient.xml")
				+ "$");
		map.put(FilePattern.RAW_DIR, "^" + Pattern.quote("raw") + "$");
		map.put(FilePattern.RAW, "^" + "raw/.*" + "$");
		map.put(FilePattern.EIT_DIR, "^" + Pattern.quote("eit") + "$");
		map.put(FilePattern.EIT_MANIFEST,
				"^" + Pattern.quote("eit/manifest.xml") + "$");
		map.put(FilePattern.EIT_DATA_DIR, "^" + Pattern.quote("eit/data") + "$");
		map.put(FilePattern.EIT_DATA, "^" + "eit/data/data_\\d{5}\\.sframes"
				+ "$");
		map.put(FilePattern.EIT_DATA_BASE,
				"^" + Pattern.quote("eit/data/data_00000.sframes") + "$");
		map.put(FilePattern.EIT_CONFIG_DIR, "^" + Pattern.quote("eit/config")
				+ "$");
		map.put(FilePattern.EIT_CONFIG, "^" + "eit/config/config_\\d{5}\\.xml"
				+ "$");
		map.put(FilePattern.EIT_CONFIG_BASE,
				"^" + Pattern.quote("eit/config/config_00000.xml") + "$");
		map.put(FilePattern.EIT_LOG_DIR, "^" + Pattern.quote("eit/log") + "$");
		map.put(FilePattern.EIT_LOG, "^" + Pattern.quote("eit/log/log.xml"));
		map.put(FilePattern.OPTIONAL_DATA_DIR,
				"^" + Pattern.quote("optional_data") + "$");
		map.put(FilePattern.OPTIONAL_DATA, "^" + "optional_data/.*" + "$");
		map.put(FilePattern.AUXILIARY_DIR, "^" + "auxiliary\\d{2}" + "$");
		map.put(FilePattern.AUXILIARY_DATA, "^" + "auxiliary\\d{2}/.*" + "$");

		__filePatterns = Collections.unmodifiableMap(map);
	}

	@AfterClass
	public static void tearDownClass() {

	}

	@Before
	public void setUp() throws Exception {
		_zfh.open();
	}

	@After
	public void tearDown() throws Exception {
		_zfh.close();
	}

	@Test
	public void test_file_contains_version_txt() {
		assertThat("version.txt exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.VERSION)), is(true));
	}

	@Test
	public void test_file_contains_device_xml() {
		assertThat("header/device.xml exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.DEVICE)), is(true));
	}

	@Test
	public void test_file_contains_electrode_xml() {
		assertThat("header/electrode.xml exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.ELECTRODE)), is(true));
	}

	@Test
	public void test_file_contains_subject() {
		assertThat("header/subject.xml exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.SUBJECT)), is(true));
	}

	@Test
	public void test_file_contains_manifest() {
		assertThat("eit/manifest.xml exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.EIT_MANIFEST)), is(true));
	}

	@Test
	public void test_file_contains_data() {
		assertThat("eit/data/data_00000.sframes exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.EIT_DATA_BASE)), is(true));
	}

	@Test
	public void test_file_contains_config() {
		assertThat("eit/config/config_00000.xml exists",
				_zfh.doesFileContainPattern(__filePatterns
						.get(FilePattern.EIT_CONFIG_BASE)), is(true));
	}

	@Test
	public void test_file_does_not_contain_unknown_files() {
		assertThat("unknown files exist",
				_zfh.getEntriesNotForPatterns(__filePatterns.values()).size(),
				is(0));
	}

	@Test
	public void test_version_is_1_x_y() throws Exception {
		assertThat("version",
				_zfh.readTextFile(__filePatterns.get(FilePattern.VERSION))
						.matches("^1\\.\\d{1,2}(\\.\\d{1,3})?(\\n)?$"), is(true));
	}

	@Test
	public void test_device_xml_conforms_to_schema() throws Exception {
		assertThat("header/device.xml schema", XmlHelper.isXmlValid(_zfh
				.readTextFile(__filePatterns.get(FilePattern.DEVICE)), this
				.getClass().getResourceAsStream("device.xsd")), is(true));
	}

	@Test
	public void test_electrode_xml_conforms_to_schema() throws Exception {
		assertThat("header/electrode.xml schema", XmlHelper.isXmlValid(
				_zfh.readTextFile(__filePatterns.get(FilePattern.ELECTRODE)),
				this.getClass().getResourceAsStream("electrode.xsd")), is(true));
	}

	@Test
	public void test_subject_xml_conforms_to_schema() throws Exception {
		assertThat("header/subject.xml schema", XmlHelper.isXmlValid(
				_zfh.readTextFile(__filePatterns.get(FilePattern.SUBJECT)),
				this.getClass().getResourceAsStream("subject.xsd")), is(true));
	}

	@Test
	public void test_patient_xml_conforms_to_schema() throws Exception {
		String s = _zfh.readTextFile(__filePatterns.get(FilePattern.SUBJECT));

		assumeNotNull(s);

		assertThat(
				"header/subject.xml schema",
				XmlHelper.isXmlValid(s,
						this.getClass().getResourceAsStream("patient.xsd")),
				is(true));
	}

	@Test
	public void test_manifest_xml_conforms_to_schema() throws Exception {
		assertThat("eit/manifest.xml schema",
				XmlHelper.isXmlValid(_zfh.readTextFile(__filePatterns
						.get(FilePattern.EIT_MANIFEST)), this.getClass()
						.getResourceAsStream("manifest.xsd")), is(true));
	}

	@Test
	public void test_config_NNNNN_xml_conforms_to_schema() throws Exception {
		Map<String, String> configs = _zfh.readTextFiles(__filePatterns
				.get(FilePattern.EIT_CONFIG));

		for (String s : configs.keySet()) {
			assertThat(s + " schema", XmlHelper.isXmlValid(configs.get(s), this
					.getClass().getResourceAsStream("configuration.xsd")),
					is(true));

		}
	}

}