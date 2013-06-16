all: MceitTest

MceitTest: GetParser.jar 
	java -cp GetParser.jar MceitTest

GetParser.jar: 
	jar cf $< src/*
