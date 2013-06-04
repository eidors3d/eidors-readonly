import java.io.IOException;

import junit.framework.Assert;

import org.junit.Test;

public class IIRCTest {
	@Test
	public void test() throws IOException {
		OeitLegacyParser parserFactory = OeitParserFactory.createParser("IIRC");
		EitEntity entity = parserFactory.parse("sample");

		Assert.assertEquals(0.0, entity.getLength());

		double[] data = entity.getOutput();

		for (int i = 0; i < 256; i++) {
			Assert.assertEquals(0.0, data[i]);
		}
	}
}
