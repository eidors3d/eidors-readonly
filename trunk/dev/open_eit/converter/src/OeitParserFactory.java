/*   - MCEIT (Goettingen / Viasys) "get" file format 
%        format = "GET" or "MCEIT"
%    - Draeger "eit" file format (post 2008 format for Draeger equipment)
%        format = "draeger-eit"
%    - Draeger "get" file format (older - pre 2007 - format for Draeger equipment)
%        format = "GET" or "draeger-get"
%    - New Carefusion "EIT" file format
%        format = "EIT" or "carefusion"
%    - Sheffield MK I "RAW" file format
%        format = "RAW" or "sheffield"
%    - ITS (International Tomography Systems)
%        format = "ITS" or "p2k"
%    - IIRC (Impedance Imaging Research Center, Korea)
%        format = "txt" or "IIRC"
%    - University of Cape Town formats
%        format = "UCT_SEQ"  UCT sequence file
%           - Output: [ stimulation, meas_select]= eidors_readdata(fname, 'UCT_SEQ')
%        format = "UCT_CAL"  UCT calibration file
%           - Output: [vv, no_cur_caldata_raw ]= eidors_readdata( fname, 'UCT_CAL' )
%                 where no_cur_caldata_raw is data captured with no current
%        format = "UCT_DATA"  UCT data frame file
%           - Output: [vv]= eidors_readdata( fname, 'UCT_DATA' )
%    - Landquart, Switzerland EIT equipment 'LQ1' (pre - 2011)
%    - Landquart, Switzerland EIT equipment 'LQ2' (2013 - ?)
%        files are 'EIT' files, but are not autodetected
%
%    - Dixtal file format, from Dixtal inc, Brazil
%        format = 'DIXTAL_encode' extract encoder from provided Dll
%           - Output: [encodepage] = eidors_readdata( path,'DIXTAL_encode');
%              where path= '/path/to/Criptografa_New.dll' (provided with system)
%        format = 'DIXTAL'
%           - output: [vv] = eidors_readdata( fname, 'DIXTAL', [], encodepage );
%
 */

public class OeitParserFactory {
	public static final String MCEIT = "MCEIT";
	public static final String DRAEGER_EIT = "draeger-eit";
	public static final String DRAEGER_GET = "draeger-get";

	public static OeitLegacyParser createParser(String parserType) {
		OeitLegacyParser factory = null;
		if (MCEIT.equalsIgnoreCase(parserType)) {
			factory = new MCEITParser();
		} else if (DRAEGER_EIT.equals(parserType)) {
			factory = new DraegereitParser();
		} else if (DRAEGER_GET.equals(parserType)) {
			factory = new DraegergetParser();
		} else {
			// new NotImplementedException("Invalid parser type requested.");
		}
		return factory;
	}

}
