"""Minimal offline example for the public GenCAD v2 task interface."""
import numpy as np

def generate(task):
    inner = np.asarray(task['preparation_z']) + 0.08
    outer = inner + task['requirements']['wall_mm'] + 0.05
    return dict(task_id=task['task_id'], status='DESIGN', units='mm', frame=task['frame'], outer_vertices=np.c_[task['xy'], outer].tolist(), inner_vertices=np.c_[task['xy'], inner].tolist(), faces=task['faces'].tolist())
