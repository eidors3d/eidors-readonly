all: MceitTest
CLASS=src/CarefusionParser.class \
      src/DixtalEncodeParser.class \
      src/DixtalParser.class \
      src/DraegereitParser.class \
      src/DraegergetParser.class \
      src/EitEntity.class \
      src/IIRCParser.class \
      src/ITCParser.class \
      src/Landquart1Parser.class \
      src/Landquart2Parser.class \
      src/MCEITParser.class \
      src/OeitLegacyParser.class \
      src/OeitParserFactory.class \
      src/RunParser.class \
      src/SheffieldParser.class \
      src/UTCCalParser.class \
      src/UTCDataParser.class \
      src/UTCSeqParser.class \
      src/XMLReader.class \
      test/CarefusionTest.class \
      test/DixtalEncodeTest.class \
      test/DixtalTest.class \
      test/DraegereitTest.class \
      test/DraegergetTest.class \
      test/IIRCTest.class \
      test/ITSTest.class \
      test/Landquart1Test.class \
      test/Landquart2Test.class \
      test/MceitTest.class \
      test/SheffieldTest.class \
      test/UTCCalTest.class \
      test/UTCDataTest.class \
      test/UTCSeqTest.class

SRC = $(patsubst %.class, %.java, $(CLASS))

# This is for UNIX
# In ununtu, install apt-get install junit4
CLASSPATH = :./usr/share/java/junit4.jar

$(CLASS):
	javac -cp $(CLASSPATH) $(SRC)

MceitTest: GetParser.jar 
	java -cp GetParser.jar MceitTest

GetParser.jar: $(CLASS)
	jar cf $@ $(CLASS)
