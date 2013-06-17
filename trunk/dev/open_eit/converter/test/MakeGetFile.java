import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public class MakeGetFile {
	
	public static void main(String[] args) throws IOException {
		int i, j, value;
		DataOutputStream os = new DataOutputStream(new FileOutputStream("sample.get"));
		
		for (int frame = 0; frame<1; frame++){
			for (i = 0; i<16; i++){
				for (j = 0; j<13; j++){
					value = (frame << 24) + ((i & 0xF) << 8) + (j & 0x0F);
					System.out.println(value);
					os.writeInt(value);
				}
			}
			
			for (j = 1; j<4; j++){
				for (i = 0; i<16; i++){
					value = (frame << 24) + (j<<16) + i;
					os.writeInt(value);
				}
			}
		}
		
		os.close();
	}
}
