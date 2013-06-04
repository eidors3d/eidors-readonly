import java.io.FileNotFoundException;
import java.io.IOException;


public interface OeitLegacyParser {
	public EitEntity parse(String fileName) throws FileNotFoundException, IOException;
}

