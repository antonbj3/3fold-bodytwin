"""Opt-in adapters for conditional contact bands and exact scalar equilibria."""
from .contact_port_v1 import Gap, Unit, Source, Assurance, ContactCandidate, Branch, BranchSet, Witness, canonical, digest, exact

def from_contact_band(band, *, source, condition, participants, frame, normals=None):
    """Preserve represented bounds and projected JOINT sigma, never tighten them.

    `condition` must state the upstream z/event/threshold; it cannot be inferred
    from ContactBand, which does not retain z or its original covariance.
    """
    n = len(band.lower_mm)
    if len(participants) != n or (normals is not None and len(normals) != n):
        raise ValueError('explicit geometry binding must match band')
    if any((len(getattr(band, k)) != n for k in ('upper_mm', 'sigma_tot_mm', 'classes', 'completion_unknown', 'sign_conflict'))):
        raise ValueError('inconsistent band arrays')
    return tuple((ContactCandidate(str(i), participants[i], frame, Gap(float(band.lower_mm[i]), float(band.upper_mm[i]), Unit.MM, source, float(band.sigma_tot_mm[i]), Assurance.CONDITIONAL, condition, bool(band.completion_unknown[i] or band.sign_conflict[i])).to(Unit.M), None if normals is None else normals[i]) for i in range(n)))

def from_univariate(cert, *, source):
    from .univariate_equilibrium import verify
    if not verify(cert):
        raise ValueError('invalid polynomial certificate')
    w = Witness('univariate_equilibrium', canonical(cert), digest(cert), source)
    branches = tuple((Branch('root:' + str(i), canonical({'box': r['box'], 'scope': cert['scope']}), (w,)) for (i, r) in enumerate(cert['roots'])))
    return BranchSet(branches, cert['complete'], cert['scope'], w if cert['complete'] else None, cert.get('reason', ''))

def from_shell_equilibria(cert, *, source, replay):
    """Consumer replay must check the energy, both supports and force readout.

    A root-box check alone is insufficient. Full source input is retained in
    every witness; instability is preserved, not silently filtered out.
    """
    if not replay(cert):
        raise ValueError('invalid reduced-shell consumer certificate')
    if cert.get('model_scope') != 'REDUCED_MODEL_ONLY':
        raise ValueError('unsupported shell scope')
    w = Witness('reduced_shell_contact_set', canonical(cert), digest(cert), source)
    bs = tuple((Branch('shell:' + str(i), canonical(s), (w,)) for (i, s) in enumerate(cert['states'])))
    return BranchSet(bs, cert['complete'], cert['model_scope'], w if cert['complete'] else None, cert['status'])

def from_aligner_port(port, *, source):
    """Lossless legacy transport; standalone model_branches lack coverage replay."""
    if port.get('schema') != 'dental-aligner-branch/v1':
        raise ValueError('unsupported legacy port')
    w = Witness('legacy_aligner_port', canonical(port), digest(port), source)
    bs = tuple((Branch(k, canonical(v), (w,)) for (k, v) in port['model_branches'].items()))
    return BranchSet(bs, False, 'REDUCED_MODEL_ONLY', reason='LEGACY_COVERAGE_NOT_REPLAYED')
