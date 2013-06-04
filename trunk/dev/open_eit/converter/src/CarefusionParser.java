import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;


public class CarefusionParser implements OeitLegacyParser{
	
	public EitEntity parse(String fileName) throws IOException {
		EitEntity carefusionInfo = new EitEntity();
		
		InputStream is = new FileInputStream(fileName);
		DataInputStream dis = new DataInputStream(is);
		
		double[] d = new double[180];
		
		for (int i = 1; i<180; i++){
			d[i] = dis.readDouble();
		}
		
		dis.close();
		
		return carefusionInfo;
	}
}
