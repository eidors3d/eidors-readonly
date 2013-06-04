import java.io.IOException;

import junit.framework.Assert;

import org.junit.Test;

public class Landquart1Test {

	@Test
	public void test() throws IOException {
		OeitLegacyParser parserFactory = OeitParserFactory.createParser("lq1");
		EitEntity entity = parserFactory.parse("sample");

		Assert.assertEquals(0.0, entity.getLength());

		double[] data = entity.getOutput();

		for (int i = 0; i < entity.getLength(); i++) {
			Assert.assertEquals(0.0, data[i]);
		}
	}

}
