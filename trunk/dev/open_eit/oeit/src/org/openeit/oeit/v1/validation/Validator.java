package org.openeit.oeit.v1.validation;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import javax.xml.parsers.ParserConfigurationException;

import org.openeit.oeit.IValidator;
import org.openeit.oeit.ValidationLevel;
import org.xml.sax.SAXException;

public class Validator implements IValidator {

	private ConcurrentHashMap<ValidationLevel, List<String>> _messages;
	private ZipFileHelper _zf = null;
	private List<? extends ZipEntry> _entries = null;
	private String _filespec;

	public Validator() {
	}

	private void addError(String s) {
		_messages.get(ValidationLevel.ERROR).add(
				"File = '" + _filespec + "' --> " + s);
	}

	private void addInfo(String s) {
		_messages.get(ValidationLevel.INFO).add(
				"File = '" + _filespec + "' --> " + s);
	}

	private void addWarning(String s) {
		_messages.get(ValidationLevel.WARNING).add(
				"File = '" + _filespec + "' --> " + s);
	}

	private void initialize(String filespec) {
		_messages = new ConcurrentHashMap<ValidationLevel, List<String>>();

		for (ValidationLevel level : ValidationLevel.values())
			_messages.put(level, new LinkedList<String>());

		_filespec = filespec;
		_zf = new ZipFileHelper(filespec);
		_entries = null;

	}

	@Override
	public Map<ValidationLevel, List<String>> validate(String filespec) {

		try {
			initialize(filespec);
			validateOeitFile();
			validateOeitXml();
		} catch (ZipException e) {
			addError("Is not a valid zip file.");
		} catch (Exception e) {
			addError(e.getMessage());
		} finally {
			try {
				_zf.close();
			} catch (Exception e) {
			}
		}

		return _messages;
	}

	private void validateOeitXml() throws IOException,
			ParserConfigurationException, SAXException {
		String text = _zf.readTextFile("oeit.xml");
		XmlHelper.validate(text, this.getClass()
				.getResourceAsStream("oeit.xsd"));
	}

	private void validateOeitFile() throws Exception {
		if (!_filespec.endsWith(".oeit"))
			addWarning("Does not end with .oeit");

		validateOeitFileContainsOeitXml();
		validateOeitFileDoesNotContainRootFiles();
		validateOeitFileContainsNonEmptyEit0000Directory();
		validateOeitFileDoesNotContainNonEitSframeFiles();
		validateOeitFileDoesNotContainNonAuxiliarySframeFiles();
		validateOeitFileContainsInfoDevices();
		validateOeitFileContainsInfoElectrodeTypes();
		validateOeitFileContainsInfoElectrodes();
		validateOeitFileContainsInfoFrameTypes();
		validateOeitFileContainsInfoMeasTypes();
		validateOeitFileContainsInfoStimTypes();
		validateOeitFileContainsInfoStreams();
		validateOeitFileContainsInfoSubject();
	}

	private void validateOeitFileContainsInfoDevices() {
		validateOeitFileContainsSubfile("info/devices.xml", false);
	}

	private void validateOeitFileContainsInfoElectrodeTypes() {
		validateOeitFileContainsSubfile("info/electrode_types.xml", false);
	}

	private void validateOeitFileContainsInfoElectrodes() {
		validateOeitFileContainsSubfile("info/electrodes.xml", false);
	}

	private void validateOeitFileContainsInfoFrameTypes() {
		validateOeitFileContainsSubfile("info/frame_types.xml", false);
	}

	private void validateOeitFileContainsInfoMeasTypes() {
		validateOeitFileContainsSubfile("info/meas_types.xml", false);
	}

	private void validateOeitFileContainsInfoStimTypes() {
		validateOeitFileContainsSubfile("info/stim_types.xml", false);
	}

	private void validateOeitFileContainsInfoStreams() {
		validateOeitFileContainsSubfile("info/streams.xml", false);
	}

	private void validateOeitFileContainsInfoSubject() {
		validateOeitFileContainsSubfile("info/subject.xml", false);
	}

	private void validateOeitFileContainsNonEmptyEit0000Directory() {
		_entries = _zf.getEntriesForPattern("^eit/0000/.*$");

		if (_entries.size() == 0) {
			addError("Contains empty 'eit/0000/' directory");
		}
	}

	private void validateOeitFileContainsOeitXml() {
		validateOeitFileContainsSubfile("oeit.xml", true);
	}

	private void validateOeitFileContainsSubfile(String subfile,
			boolean required) {
		if (!_zf.containsPattern(subfile)) {
			if (required)
				addError("Must contain '" + subfile + "'");
			else
				addWarning("Should contain '" + subfile + "'");
		}
	}

	private void validateOeitFileDoesNotContainNonAuxiliarySframeFiles() {
		_entries = _zf.getEntriesForPattern("^auxiliary/\\d{4}/.*$");
		_entries.removeAll(_zf
				.getEntriesForPattern("^auxiliary/\\d{4}/\\d{5}.sframes$"));

		for (ZipEntry e : _entries) {
			addWarning("Contains auxiliary-non-sframe file '" + e.getName()
					+ "'");
		}
	}

	private void validateOeitFileDoesNotContainNonEitSframeFiles() {
		_entries = _zf.getEntriesForPattern("^eit/\\d{4}/.*$");
		_entries.removeAll(_zf
				.getEntriesForPattern("^eit/\\d{4}/\\d{5}.sframes$"));

		for (ZipEntry e : _entries) {
			addWarning("Contains non-sframe file '" + e.getName() + "'");
		}
	}

	private void validateOeitFileDoesNotContainRootFiles() {
		_entries = _zf.getEntriesForPattern("^[^/]*$");

		if (_entries.size() > 1) {
			for (ZipEntry e : _entries) {
				if (e.getName() != "oeit.xml")
					addWarning("Contains root file '" + e.getName() + "'");
			}
		}
	}
}
