import java.io.IOException;

import junit.framework.Assert;

import org.junit.Test;

public class UTCDataTest {

	@Test
	public void test() throws IOException {
		OeitLegacyParser parserFactory = OeitParserFactory
				.createParser("UTC_DATA");
		EitEntity entity = parserFactory.parse("sample");

		Assert.assertEquals(0.0, entity.getLength());

		double[] data = entity.getOutput();

		for (int i = 0; i < entity.getLength(); i++) {
			Assert.assertEquals(0.0, data[i]);
		}
	}

}
