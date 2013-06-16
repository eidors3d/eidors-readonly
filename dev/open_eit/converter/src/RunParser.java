import java.io.FileNotFoundException;
import java.io.IOException;

public class RunParser {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		OeitLegacyParser parserFactory = OeitParserFactory
				.createParser(OeitParserFactory.MCEIT);
		EitEntity entity = parserFactory.parse("sample.get");
	}
}
