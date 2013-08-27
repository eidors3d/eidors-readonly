package org.openeit.oeit;

import java.util.List;
import java.util.Map;

public interface IValidator {
	public Map<ValidationLevel, List<String>> validate(String filespec);
}
