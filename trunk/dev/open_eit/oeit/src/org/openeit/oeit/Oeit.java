package org.openeit.oeit;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

public class Oeit {

	private static List<String> __inputFiles = new LinkedList<String>();

	public static List<? extends String> getInputFiles() {
		return Collections.unmodifiableList(__inputFiles);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length == 0)
			displayHelp();
		else {
			handleCommand(args);
		}
	}

	public static void displayHelp() {
		System.out.println("Usage:  java -jar oiet.jar <command>");
		System.out.println("<command> = ");
		System.out.println("\tvalidate\t\t- validate OEIT files ");
		System.out.println("\thelp <command>\t\t- get help for <command>");
	}

	public static void handleCommand(String[] args) {
		switch (args[0]) {

		case "validate":
			validateOeitFile(args);
			break;
		case "help":
			displayHelp();
			break;
		default:
			displayHelp();
			break;
		}
	}

	public static void validateOeitFile(String[] args) {
		if (args.length != 2) {
			System.out
					.println("Usage:  java -jar oeit.jar validate <filespec>");
		}

		if (!resolveInputFiles(args[1])) {
			System.out.println("Invalid file specification for validation!");

			System.out
					.println("Usage:  java -jar oeit.jar validate <filespec>");
		}

		Result result = org.junit.runner.JUnitCore
				.runClasses(org.openeit.oeit.v1.validation.OeitFileTest.class);

		for (Failure f : result.getFailures()) {
			System.out.println(f.getTestHeader() + ": " + f.getMessage());
		}
	}

	private static boolean resolveInputFiles(String name) {
		File f = new File(name);

		if (f.exists()) {
			if (f.isFile()) {
				__inputFiles.add(f.getAbsolutePath());
				return true;
			} else if (f.isDirectory()) {
				for (String s : f.list()) {
					if (new File(s).isFile())
						__inputFiles.add(s);
				}

				return true;
			}
		}

		return false;
	}
}
