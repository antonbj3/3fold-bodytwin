"""Minimal adapter; replace this function with a field/commercial/AI generator.

An STL wrapper should return triangle_mesh_v1 with labelled outer and intaglio
arrays in task coordinates. Unsupported proof scopes honestly remain UNKNOWN.
"""
from .generators.baselines import parametric

def generate(task):
    return parametric(task)
