import igl
import numpy as np
from scipy.io import savemat
import os
from PyRMT import RMTMesh

SamplesDir = "../../samples"

v, f = igl.read_triangle_mesh(os.path.join(SamplesDir, "cat2.off"))
v = np.asfortranarray(v)
f = np.asfortranarray(f)

m = RMTMesh(v, f)
m.make_manifold()
m = m.remesh(6000)
m.clean_up()

bm = m.baryc_map(v)

igl.write_triangle_mesh(os.path.join(SamplesDir, "cat2-remesh.off"), m.vertices, m.triangles)
savemat(os.path.join(SamplesDir, "cat2-remesh.mat"), {"U" : bm})