package org.openeit.oeit;

public class ValidatorFactory {
	public static IValidator getValidator(){
		return (IValidator) new org.openeit.oeit.v1.validation.Validator();
	}
}
