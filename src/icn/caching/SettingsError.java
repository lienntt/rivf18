package icn.caching;

/**
 * Settings related error
 *
 */
public class SettingsError extends SimError {

	public SettingsError(String cause) {
		super(cause);
	}
	
	public SettingsError(String cause, Exception e) {
		super(cause,e);
	}
	
	public SettingsError(Exception e) {
		super(e);
	}

}
