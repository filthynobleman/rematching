![CMake Windows](https://github.com/filthynobleman/rematching/actions/workflows/cmake-windows.yml/badge.svg)
![CMake Ubuntu](https://github.com/filthynobleman/rematching/actions/workflows/cmake-ubuntu.yml/badge.svg)
![CMake Mac OS](https://github.com/filthynobleman/rematching/actions/workflows/cmake-macos.yml/badge.svg)

# ReMatching
This repository implements the remeshing algorithm presented in *ReMatching: Low-Resolution Representations for Scalable Shape Correspondence*.  
This branch contains the Python bindings for the core C++ implementation.
- DOI: https://arxiv.org/abs/2305.09274
- PDF: https://arxiv.org/abs/2305.09274

## Installation instructions
(Optional: create a virtual environment and activate it.)
`PyRMT` python bindings can simply be installed by cloning the project. Inside the root directory of the project one simply needs to execute
```
python setup.py install
```

Alternatively, one can do:

```bash
pip install git+https://github.com/filthynobleman/rematching.git@python-binding
```

For building the project, you need [CMake](https://cmake.org/) and a C++ compiler compliant with [C++17 standard](https://en.cppreference.com/w/cpp/compiler_support/17).  
The Python bindings also require a working installation of [Python 3](https://www.python.org/) and [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html).  


## Usage
The directory `python/` contains a sample Python application `sample.py` for testing the method and showcasing the interface.
The sample application also requires a working installation of [libigl for Python](https://libigl.github.io/libigl-python-bindings/), [numpy](https://numpy.org/), and [scipy](https://scipy.org/).  


```python
import igl
import numpy as np
from scipy.io import savemat
import os
from PyRMT import RMTMesh
```
The application imports all the required headers.

```python
SamplesDir = "../../samples"

v, f = igl.read_triangle_mesh(os.path.join(SamplesDir, "cat2.off"))
v = np.asfortranarray(v)
f = np.asfortranarray(f)
```
We load the mesh using the `igl.read_triangle_mesh()` function, obtaining a matrix of vertices `v` and a matrix of triangles `f`.  
> :warning: The `np.asfortranarray()` is required for converting the array into Fortran contiguous arrays. This is required for the current implementation of the binding, but could be removed in future releases.

```python
m = RMTMesh(v, f)
```
We initialize a `RMTMesh` structure with the vertices and faces. They are accessible as `m.vertices` and `m.triangles`.  
Notice that the matrices are copied for initializing the mesh, thus changing `v` and `f` does not affect `m`, and vice versa.

```python
m.make_manifold()
m = m.remesh(6000)
m.clean_up()
```
The mesh is preprocessed to ensure working with manifold meshes, remeshed to 6000 vertices, and then postprocessed to remove small disconnected components and degenerate faces.

```python
bm = m.baryc_map(v)
```
We compute the **sparse** matrix `bm` for mapping scalar functions from the low-resolution mesh to the high-resolution mesh.  
Approximately, `v = bm @ m.vertices`.

```python
igl.write_triangle_mesh(os.path.join(SamplesDir, "cat2-remesh.off"), m.vertices, m.triangles)
savemat(os.path.join(SamplesDir, "cat2-remesh.mat"), {"U" : bm})
```
Finally, the remeshed surface and the sparse barycentric mapping matrix are saved with proper file formats.
