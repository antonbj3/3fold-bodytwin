"""Reject impossible DC values, invalid uncertainty and observations predating freeze."""
import datetime
import math

def checked_intake(record, original_intake, frozen):
    if record.get('kind') == 'resin_DC_30min':
        value = record.get('value')
        uncertainty = record.get('uncertainty')
        if isinstance(value, bool) or not isinstance(value, (int, float)) or (not math.isfinite(value)) or (not 0 <= value <= 100):
            return {'status': 'REJECTED', 'reason': 'DC must be finite in [0,100] percent'}
        if isinstance(uncertainty, bool) or not isinstance(uncertainty, (int, float)) or (not math.isfinite(uncertainty)) or (uncertainty < 0):
            return {'status': 'REJECTED', 'reason': 'Nonnegative finite DC uncertainty in percentage points required'}
        try:
            measured = datetime.datetime.fromisoformat(record['measured_at'].replace('Z', '+00:00'))
            freeze = datetime.datetime.fromisoformat(frozen['frozen_at'].replace('Z', '+00:00'))
            if measured.tzinfo is None or freeze.tzinfo is None or measured < freeze:
                raise ValueError('Observation predates freeze or lacks timezone')
        except (KeyError, ValueError, TypeError):
            return {'status': 'REJECTED', 'reason': 'Timestamp must be explicit and on or after prediction freeze'}
    return original_intake(record)
