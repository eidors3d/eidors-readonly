import java.io.IOException;

import junit.framework.Assert;

import org.junit.Test;

public class SheffieldTest {
	@Test
	public void test() throws IOException {
		OeitLegacyParser parserFactory = OeitParserFactory
				.createParser("sheffield");
		EitEntity entity = parserFactory.parse("sample");

		Assert.assertEquals(0.0, entity.getVoltage());
		Assert.assertEquals(0.0, entity.getCurrent());
		Assert.assertEquals(0.0, entity.getLength());

		double[][] data = entity.getOutput2D();

		for (int i = 0; i < 0/* unknown */; i++) {
			for (int j = 0; j < entity.getLength(); j++) {
				Assert.assertEquals(0.0, data[i][j]);
			}
		}
	}
}
