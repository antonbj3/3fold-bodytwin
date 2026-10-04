# HUM robotport interval extension

This catalog extends the current202 with two conditional HUM model answers.
All old records and their source classes remain unchanged. These are not native
21DoF-hand predictions or empirically validated surgical grasp capacities.

Run the commands in ../CATALOG_HUM_20261003.json to validate or compose a204-port
view. They print to stdout and never modify the historical PORT.json.
The adapter invokes the exact ROBOT_PORT_CONNECT engine_r1 node validator,
then separately enforces interval/provenance/status semantics and replays the
two-contact mathematical certificate. Existing verify_r5.py remains untouched.
Its original immutable terminal-handoff assertions are historical, not this
extension's numerical validator.

The scalar-output rejection test is test_point_value_rejected_where_interval_required.
Run: python3 -B -m unittest discover -s HUM_ROBOTPORT_20261003 -p test_hum_robotport.py -v
The optional Contact port v1 check requires field_engine.contact_port_v1 from
the pinned1216032 implementation. Other tests use the Python standard library.
Every source has an explicit file and hash. Full primary sources and original
raw artifacts remain at their recorded read-only source locations in this lane;
no private subject records or licensed anatomy assets are part of this patch.

Force ranges:12–20N per finger for the original model;11.6–13.2N total shear
capacity for the dry flat-scalpel/single-glove mean-friction envelope at10N
per finger. Published geometry is regime data, not an inferred contact-gap law.
Needle-driver mean20.11N is not covered by the stated friction-only closures.
UNKNOWN remains on empirical validity, normal-load transfer, torque, peak
pressure and rate friction. BodyTwin graph admission of the external CAD
consumer namespace requires coordinator review. No graph status is promoted.
