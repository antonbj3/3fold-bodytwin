"""Guard explicitly synthetic locators; other physical provenance needs review."""
from stress_port import measured_query as original_measured_query

def measured_query(port, R):
    tokens = ('our_fixture', 'our_own_fixture', 'simulated')
    if any((any((token in str(locator).lower() for token in tokens)) for locator in port.get('measurement_locators', []))):
        return {'status': 'UNKNOWN', 'reason': 'known synthetic measurement locator'}
    return original_measured_query(port, R)
