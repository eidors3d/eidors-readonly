import java.io.IOException;
import junit.framework.Assert;

import org.junit.Test;

public class DraegergetTest {

	@Test
	public void test() throws IOException {
		OeitLegacyParser parserFactory = OeitParserFactory
				.createParser(OeitParserFactory.DRAEGER_GET);
		EitEntity entity = parserFactory.parse("sample");

		// Assert.assertEquals(0.0, entity.getFrameLength());
		Assert.assertEquals(0.0, entity.getCurrent());
		Assert.assertEquals(0.0, entity.getLength());

		double[][] data = entity.getOutput2D();

		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < entity.getLength() / 256; j++) {
				Assert.assertEquals(0.0, data[i][j]);
			}
		}
	}
}