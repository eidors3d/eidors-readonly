package org.openeit.oeit.v1.validation;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class ZipFileHelper implements Closeable {
	private String _name = "";
	private ZipFile _zf = null;
	private List<? extends ZipEntry> _ze = null;

	public ZipFileHelper(String file) {
		_name = file;
	}

	public void open() throws IOException {
		_zf = new ZipFile(_name);
		_ze = Collections.unmodifiableList(Collections.list(_zf.entries()));
	}

	public void close() throws IOException {
		if (null != _zf)
			_zf.close();
	}

	public boolean doesFileContainPattern(String pattern) {
		return !getEntriesForPattern(pattern).isEmpty();
	}

	public List<? extends ZipEntry> getEntriesForPattern(String pattern) {
		List<ZipEntry> entries = new LinkedList<ZipEntry>();

		for (ZipEntry e : _ze) {
			if (e.getName().matches(pattern))
				entries.add(e);
		}

		return Collections.unmodifiableList(entries);
	}

	public List<? extends ZipEntry> getEntriesForPatterns(
			Collection<String> collection) {
		List<ZipEntry> entries = new LinkedList<ZipEntry>();

		for (String p : collection) {
			entries.addAll(getEntriesForPattern(p));
		}

		return Collections.unmodifiableList(entries);
	}

	public List<? extends ZipEntry> getEntriesNotForPatterns(
			Collection<String> collection) {
		List<ZipEntry> entries = new LinkedList<ZipEntry>();
		entries.addAll(_ze);
		entries.removeAll(getEntriesForPatterns(collection));

		return Collections.unmodifiableList(entries);
	}

	public String readTextFile(String pattern) throws IOException {
		for (ZipEntry e : _ze) {
			if (e.getName().matches(pattern)) {
				BufferedReader in = new BufferedReader(new InputStreamReader(
						_zf.getInputStream(e)));
				StringBuilder builder = new StringBuilder();

				for (String line = in.readLine(); line != null; line = in
						.readLine()) {
					builder.append(line + "\n");
				}

				return builder.toString();
			}
		}

		return null;
	}

	public Map<String, String> readTextFiles(String pattern) throws IOException {
		Map<String, String> results = new ConcurrentHashMap<String, String>();

		for (ZipEntry e : getEntriesForPattern(pattern))
			results.put(e.getName(),
					readTextFile("^" + Pattern.quote(e.getName()) + "$"));

		return Collections.unmodifiableMap(results);
	}
}
