"""Compose existing geometry/guide/protocol ports without a clinical closure."""
from __future__ import annotations
import hashlib
import json
import math
import zipfile
from pathlib import Path
import numpy as np
from .vendor.x8_geometry import voxel_cylinder_bracket
from .vendor.x87_science import guide as reviewed_guide

class ContractError(ValueError):
    pass

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path):
    return json.loads(Path(path).read_text())

def classify(interval, threshold=2.0):
    if interval is None:
        return 'UNKNOWN'
    if interval[0] >= threshold:
        return 'AT_OR_ABOVE_REFERENCE'
    if interval[1] < threshold:
        return 'BELOW_REFERENCE'
    return 'UNRESOLVED_NUMERIC_BRACKET'

def finite_number(x, name, positive=False):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or (not math.isfinite(x)) or (positive and x <= 0):
        raise ContractError(name + ': finite ' + ('positive ' if positive else '') + 'number required')
    return float(x)

def validate_pose(pose, frame):
    if pose.get('frame') != frame:
        raise ContractError('CBCT frame must match source: zyx, mm, spacing, origin and direction')
    try:
        e = np.asarray(pose['entry_zyx_mm'], float)
        a = np.asarray(pose['axis_zyx'], float)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError('explicit entry and unit axis required') from exc
    if e.shape != (3,) or a.shape != (3,) or (not np.isfinite(e).all()) or (not np.isfinite(a).all()):
        raise ContractError('pose must contain finite 3-vectors')
    if abs(float(np.linalg.norm(a)) - 1) > 1e-12:
        raise ContractError('axis is not unit length; no silent normalization')
    return (e, a)

def validate_frame(frame):
    if frame.get('order') != 'zyx' or frame.get('units') != 'mm':
        raise ContractError('source frame requires zyx coordinates in mm')
    sp = np.asarray(frame.get('spacing_mm'), float)
    origin = np.asarray(frame.get('origin_mm'), float)
    if sp.shape != (3,) or not np.isfinite(sp).all() or (sp <= 0).any():
        raise ContractError('positive finite voxel spacing required')
    if origin.shape != (3,) or not np.isfinite(origin).all():
        raise ContractError('finite explicit voxel origin required')
    if frame.get('direction') != [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]:
        raise ContractError('reviewed voxel-box operator requires axis-aligned identity direction; resample/transform with validated geometry before registration')

def guide_components(profile, confidence, radius, extra):
    g = reviewed_guide(profile, 1 - confidence)
    rot = g['radius_term_mm'] * (radius + extra) / 2
    total = max(g['entry_bound_mm'], g['apex_bound_mm']) + rot
    return dict(entry_mm=g['entry_bound_mm'], apex_mm=g['apex_bound_mm'], angle_deg=g['angle_bound_deg'], rotation_mm=rot, combined_displacement_mm=total, confidence=confidence, component_failure_bound=(1 - confidence) / 3, combination='max(entry,apex) + rotation; entry and apex are overlapping global pose errors', statistical_scope='Conditional on sample moments equalling population moments and implant-error transport to the drill; arbitrary dependence allowed by union bound', resolution='PHENOMENOLOGICAL', source=profile['source'])

def compact_bracket(b):
    return dict(interval_mm=[b['lower_mm'], b['upper_mm']], width_mm=b['gap_mm'], active_voxel_boxes=b['active_voxel_boxes'], witness=b['witness'], real_arithmetic_enclosure='convex supporting plane lower / feasible point upper', rigorous_IEEE_enclosure='MISSING', resolution='PER_TOOTH')

class SafetyModule:
    """The bundle root is explicit and can move without changing scientific hashes.

    Inputs are registered in the same CBCT frame. Radius is the research cylinder
    radius, not an inferred commercial implant diameter. No risk distribution or
    physical certificate is inferred from a digital gap.
    """

    def __init__(self, bundle_root=None):
        self.root = Path(bundle_root or Path(__file__).resolve().parents[1]).resolve()
        self.verify_bundle()
        self.sites = {s['site_id']: s for s in load(self.root / 'data/sites.json')}
        self.profiles = {p['id']: p for p in load(self.root / 'data/guide_profiles.json')['profiles']}
        self.protocols = {p['id']: p for p in load(self.root / 'data/protocols.json')['protocols']}
        self._points = {}
        self._geometry = {}

    def verify_bundle(self):
        for r in load(self.root / 'SOURCE_MANIFEST.json'):
            if r.get('local') and sha(self.root / r['local']) != r['sha256']:
                raise ContractError('source hash mismatch: ' + r['local'])
        for p in self.root.glob('PREREG_*.json'):
            if sha(p) != p.with_suffix(p.suffix + '.sha256').read_text().strip():
                raise ContractError('prereg hash mismatch: ' + p.name)

    def default_pose(self, site_id):
        s = self.sites[site_id]
        return dict(entry_zyx_mm=s['pose']['entry_zyx_mm'], axis_zyx=s['pose']['axis_zyx'], frame=s['frame'])

    def register_site(self, *, site_id, points_path, points_sha256, frame, nominal_key, source_locator, source_bindings=None):
        """Register local, complete voxel occupancy explicitly; no copying/data search.

        NPZ keys ending _voxels_zyx contain integer voxel-center indices, shape
        (N,3). The origin is the physical position of index zero, not a corner.
        Retained extra keys are releases, never independent anatomy by default.
        An external reader must establish completeness and provenance beforehand.
        """
        if not isinstance(site_id, str) or not site_id or site_id in self.sites:
            raise ContractError('new nonempty site ID required; existing source cannot be overwritten')
        if not isinstance(source_locator, str) or not source_locator.strip():
            raise ContractError('explicit independent dataset/source locator required')
        validate_frame(frame)
        p = Path(points_path).resolve()
        if not p.is_file() or p.stat().st_size > 50000000:
            raise ContractError('registered point source must exist and be <=50 MB; larger sources need an external bounded adapter')
        if sha(p) != points_sha256:
            raise ContractError('registered source hash mismatch')
        with zipfile.ZipFile(p) as archive:
            if sum((i.file_size for i in archive.infolist())) > 50000000:
                raise ContractError('expanded point archive exceeds registered-source memory contract')
        with np.load(p, allow_pickle=False) as z:
            keys = [k for k in z.files if k.endswith('_voxels_zyx')]
            if nominal_key not in keys:
                raise ContractError('nominal occupancy key absent')
            total = 0
            for k in keys:
                a = z[k]
                if a.ndim != 2 or a.shape[1] != 3 or len(a) == 0 or (not np.issubdtype(a.dtype, np.integer)) or (a < 0).any():
                    raise ContractError('occupancy must be nonempty nonnegative integer voxel center indices')
                total += a.nbytes
            if total > 50000000:
                raise ContractError('expanded occupancy exceeds registered-source memory contract')
        site = dict(site_id=site_id, case=site_id, fdi=None, image_group=None, pose=None, frame=frame, points_local=str(p), points_sha256=points_sha256, nominal_key=nominal_key, annotation_kind='RELEASE_REVISION_UNION' if len(keys) > 1 else 'SINGLE_MASK_NO_REVISION_MEASUREMENT', segment='UNKNOWN', bone_binding=None, external_source_bindings=source_bindings or [dict(locator=source_locator, sha256=points_sha256)], source_locator=source_locator)
        self.sites[site_id] = site
        return dict(site_id=site_id, status='REGISTERED_DIGITAL_OCCUPANCY', annotation_kind=site['annotation_kind'], physical_safety='UNKNOWN', completeness='CALLER_SOURCE_READER_ASSUMPTION')

    def points(self, site):
        path = self.root / site['points_local']
        if sha(path) != site['points_sha256']:
            raise ContractError('canal point source hash mismatch')
        point_key = (site['points_sha256'], json.dumps(site['frame'], sort_keys=True), site.get('nominal_key', 'tf2_voxels_zyx'))
        if point_key not in self._points:
            with np.load(path, allow_pickle=False) as a:
                versions = {k: a[k].astype(float) * np.array(site['frame']['spacing_mm']) + np.array(site['frame']['origin_mm']) for k in a.files if k.endswith('_voxels_zyx')}
            nominal = versions.get(site.get('nominal_key', 'tf2_voxels_zyx'))
            if nominal is None or len(nominal) == 0:
                raise ContractError('missing nonempty TF2 occupancy')
            union = np.unique(np.vstack(list(versions.values())), axis=0)
            self._points[point_key] = (nominal, union, list(versions))
        return self._points[point_key]

    def geometry(self, site, e, a, length, radius, extra):
        (nominal, union, versions) = self.points(site)
        key = (site['points_sha256'], json.dumps(site['frame'], sort_keys=True), site.get('nominal_key', 'tf2_voxels_zyx'), tuple(e), tuple(a), length, radius, extra)
        if key not in self._geometry:
            sp = np.array(site['frame']['spacing_mm'])
            d0 = compact_bracket(voxel_cylinder_bracket(nominal, sp, e, a, length, radius))
            d1 = compact_bracket(voxel_cylinder_bracket(union, sp, e, a, length, radius))
            d2 = compact_bracket(voxel_cylinder_bracket(union, sp, e, a, length + extra, radius)) if extra is not None else None
            for b in [d0, d1, d2]:
                if b and (b['width_mm'] < -1e-10 or b['width_mm'] > 1e-06):
                    raise ContractError('geometry bracket outside frozen numeric tolerance')
            self._geometry[key] = (d0, d1, d2, versions)
        return self._geometry[key]

    def protocol(self, system, length, length_reference, datum):
        p = self.protocols.get(system)
        if length_reference != 'actual_cylinder_length':
            return (None, 'UNKNOWN_LABEL_TO_ACTUAL_LENGTH: manufacturer-labelled length is not accepted as geometric length')
        if p is None:
            return (None, 'UNKNOWN_SYSTEM: no reviewed depth datum')
        if system.startswith('ZimVie'):
            required = dict(fixture_label_mm=13.0, actual_implant_mm=12.6, actual_mark_mm=13.7, platform_depth_mm=1.0, tip_mm=1.2)
            if datum != required or length != 12.6:
                return (None, 'UNKNOWN_DATUM_TRANSFER: only reviewed 13mm label /12.6mm actual /1mm platform example supported')
        return (p, 'KNOWN_PROTOCOL_SCENARIO')

    def query(self, *, site_id, cbct_pose, guide_type, implant_system, length_mm, radius_mm=2.0, confidence=0.95, length_reference='actual_cylinder_length', protocol_datum=None):
        if site_id not in self.sites:
            raise ContractError('unknown source site; register source occupancy and frame before querying')
        site = self.sites[site_id]
        validate_frame(site['frame'])
        (e, a) = validate_pose(cbct_pose, site['frame'])
        length = finite_number(length_mm, 'length_mm', True)
        radius = finite_number(radius_mm, 'radius_mm', True)
        if radius > 2:
            raise ContractError('reviewed guide radius contract r<=2mm')
        if confidence not in (0.9, 0.95):
            raise ContractError('reviewed confidence levels are .90 or .95')
        (p, protocol_status) = self.protocol(implant_system, length, length_reference, protocol_datum)
        extra = p['extra_mm'] if p else None
        (d0, d1, d2, versions) = self.geometry(site, e, a, length, radius, extra)
        profile = self.profiles.get(guide_type)
        g = guide_components(profile, confidence, radius, extra) if profile and extra is not None else None
        (n, u) = (d0['interval_mm'], d1['interval_mm'])
        has_revision = len(versions) > 1
        revision_loss = [max(0.0, n[0] - u[1]), max(0.0, n[1] - u[0])] if has_revision else None
        swept = d2['interval_mm'] if d2 else None
        drill_loss = [max(0.0, u[0] - swept[1]), max(0.0, u[1] - swept[0])] if swept else None
        guide_interval = [max(0.0, swept[0] - g['combined_displacement_mm']), swept[1] + g['combined_displacement_mm']] if g else None
        debts = [dict(term='anatomy_vs_annotation', status='UNKNOWN', replacement='same-site registered independent canal-wall boundary with measurement uncertainty'), dict(term='joint_signed_achieved_drill_pose', status='UNKNOWN', replacement='registered planned/achieved drill poses and signed distance loss in same image'), dict(term='sample_to_population_and_drill_transport', status='CONSTITUTIVE_CLOSURE', replacement='paired guide-specific achieved drill measurements with confidence bounds on moments'), dict(term='commercial_implant_and_drill_shape', status='UNKNOWN', replacement='measured true profile, diameter, datum and depth/stop tolerance'), dict(term='bone_deflection_and_material_response', status='UNKNOWN', replacement='same-site calibrated bone mechanical and drill response'), dict(term='neurosensory_outcome', status='UNKNOWN', replacement='paired observed outcome; digital clearance is not an injury endpoint'), dict(term='rigorous_float_enclosure', status='MISSING', replacement='outward-rounded implementation of complete geometry operator')]
        if not p:
            debts.append(dict(term='protocol_depth', status=protocol_status, replacement='reviewed version and matching depth datum'))
        if not profile:
            debts.append(dict(term='guide_profile', status='UNKNOWN', replacement='reviewed guide-specific moments or local measurements'))
        if not has_revision:
            debts.append(dict(term='annotation_revision_spread', status='UNKNOWN', replacement='registered retained revision or independent reread; single-mask spread cannot be reported as zero'))
        bone_pose = dict(entry_zyx_mm=e.tolist(), axis_zyx=a.tolist(), length_mm=length, radius_mm=radius)
        bound_bone = site['bone_binding'] if site['bone_binding'] and site['bone_binding']['pose'] == bone_pose else None
        conditional_class = classify(guide_interval)
        if conditional_class == 'UNRESOLVED_NUMERIC_BRACKET':
            conditional_class = 'INTERVAL_CROSSES_REFERENCE'
        conditional = dict(interval_mm=guide_interval, margin_to_2mm_interval_mm=[guide_interval[0] - 2, guide_interval[1] - 2] if guide_interval else None, class_vs_2mm=conditional_class, confidence=confidence if g else None, coverage='CONDITIONAL_MODEL_ONLY; empirical joint/anatomical coverage UNKNOWN', resolution='PHENOMENOLOGICAL')
        return dict(schema='implant-safety-query-v1', claim_type='capability', site_id=site_id, query=dict(cbct_pose=cbct_pose, length_mm=length, radius_mm=radius, guide_type=guide_type, implant_system=implant_system, length_reference=length_reference, protocol_datum=protocol_datum, confidence=confidence), reference_rule=dict(clearance_mm=2.0, nominal_interval_mm=n, class_vs_2mm=classify(n), scope='digital reference comparison, no patient approval', resolution='PER_TOOTH'), terms=dict(annotation_revision_loss=dict(interval_mm=revision_loss, status='OBSERVED_RELEASE_REVISIONS_ONLY' if has_revision else 'UNKNOWN_NO_REVISION_MEASUREMENT', versions=versions, resolution='PER_TOOTH'), drill_extra_depth=dict(maximum_or_example_mm=extra, status=protocol_status, source=p, interpretation='Swept full-radius envelope contains hypothetical tool; actual shape and achieved depth unknown', resolution='PHENOMENOLOGICAL'), directed_tool_gap_loss=dict(interval_mm=drill_loss, resolution='PER_TOOTH'), guide_entry=dict(bound_mm=g['entry_mm'] if g else None, resolution='PHENOMENOLOGICAL'), guide_apex=dict(bound_mm=g['apex_mm'] if g else None, resolution='PHENOMENOLOGICAL'), guide_angle=dict(bound_deg=g['angle_deg'] if g else None, rotation_mm=g['rotation_mm'] if g else None, resolution='PHENOMENOLOGICAL'), guide_combination=g, anatomy_error=dict(interval_mm=None, status='UNKNOWN', resolution='PER_POINT')), digital_geometry=dict(nominal=d0, revision_union=d1, tool_envelope=d2, tool_class_vs_2mm=classify(swept), telescoping_rule='d0-(d0-d1)-(d1-d2)=d2; combined interval uses d2 directly, preserving shared endpoints'), combined_conditional=conditional, physical_safety=dict(status='UNKNOWN', interval_mm=None, injury_probability=None, reason='No independent anatomy or paired signed achieved drill error; finite scenario interval cannot close these ports'), bone_observation=dict(status='BOUND_AT_EXACT_SOURCE_POSE' if bound_bone else 'UNKNOWN_CHANGED_POSE_OR_LENGTH', observation=bound_bone), unknowns=debts, source=dict(points_local=site['points_local'], sha256=site['points_sha256'], image_group=site['image_group'], segment=site['segment'] if bone_pose == site['pose'] else 'UNKNOWN_FOR_CHANGED_POSE', bindings=site['external_source_bindings']), resolution='PER_TOOTH', edge_resolution='PER_POINT', time_scale='SIMULTANEOUS', numerical_precision_mm=1e-06, physical_spatial_resolution_mm=0.3, assumptions=['Generic finite solid cylinder, r<=2mm; no commercial shape inference', 'Local snapshot contains complete retained digital canal occupancy', 'No additional Gaussian channel term, signed depth term or independent registration term is inferred from global apex error', 'Guide interval assumes source moments are population moments; no validated population or clinical transfer'])
