all: MceitTest

MceitTest: GetParser.jar 
	java -cp GetParser.jar RunParser

GetParser.jar: 
	jar cf $< src/*
