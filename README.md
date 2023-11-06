# ReMatching
This repository implements the remeshing algorithm presented in *ReMatching: Low-Resolution Representations for Scalable Shape Correspondence*.  
- DOI: https://arxiv.org/abs/2305.09274
- PDF: https://arxiv.org/abs/2305.09274

## Building instructions
The building process is entirely carried out with CMake. If you have not already cloned the repository recursively, or if you have not updated the submodules, please run
```
git submodule update --init --recursive --remote
```

For building the project, you need [CMake](https://cmake.org/) and a C++ compiler compliant with [C++17 standard](https://en.cppreference.com/w/cpp/compiler_support/17).  

First you need to build the dependencies. The only library that needs to be built is [*CUT*](https://github.com/filthynobleman/cut). Execute the following commands
```
mkdir build-cut
cd build-cut
cmake ../ext/cut -DBUILD_SAMPLES=OFF -DCMAKE_INSTALL_PREFIX="../install"
cmake --build . --config release --parallel
cmake --install .
```
This will install the *CUT* library into the `install` directory of the project's root.  

Go back to the root folder and execute
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="../install"
cmake --build . --config release --parallel
cmake --install .
```
This will produce the application executables and the install files that are used to compile
the MATLAB Mex functions.


## Integration with MATLAB
Inside the project's root folder there is the directory `matlab`. The directory contains
- `compile.m`: a MATLAB script that compiles the Mex functions realizing the remeshing and resampling;
- `rmt_remesh.m`: a MATLAB function that wraps the Mex function for the remeshing algorithm;
- `rmt_resample.m`: a MATLAB function that wraps the Mex function for the resampling algorithm.
- `rmt_wmap.m`: a MATLAB function that wraps the Mex function for computing the mapping that brings scalar functions from the remesh to the full-resolution shape.

The compilation script `compile.m` assumes that the installation of the *CUT* library and this repository are carried out with the default settings, and that the compilation is ran from the `matlab` directory. **If this is not the case, please update the paths accordingly.**


## Applications
The building process should produce two executables called `Remesh` and `BatchRemesh`, which operate a remeshing on, respectively, a single mesh or an entire dataset.

### Remeshing a single shape
The `Remesh` application applies the remeshing algorithm to a single 3D model. To run it, please execute the following command
```
Remesh input_mesh num_samples [-o|--output out_mesh] [-r|--resample] [-e|--evaluate]
```
The semantics of the arguments is the following:
 - `input_mesh` is the path to a **triangular** mesh. Currently, only `OBJ`, `OFF` and `PLY` file formats are supported.
 - `num_samples` is the number of vertices that the output mesh must have.
 - `out_mesh` is the path where the output is saved, by default the basename of the input mesh in the current working directory. The output format is inferred from the path name.
 - `-r` applies a resampling strategy during the preprocessing for a more uniform remeshing.
 - `-e` evaluates the resulting mesh with different metrics.
Together with the output mesh, a file in ASCII Market file format (`.mat`) is produced, which contains the triplets to build a sparse matrix that can transfer scalar functions from the remeshed shape to the original meshes using barycentric interpolation (from now on, referred to as the _weight map_).  

**WARNING:** please, be aware that the program currently does not warn about overwriting. So, if your input mesh is in the current working directory and you don't provide an output path, the input mesh will be silently overwritten. Again, if you run the algorithm multiple times without specifying the output path, only the last output mesh and the last _weight map_ will be saved.  

Alternatively, the program can be fed with a configuration file in `JSON` file format
```
Remesh -f|--file config_file
```
The configuration file must at least contain the string attribute `input_mesh` and the integer numeric attribute `num_samples`. Optionally, you can provide:
 - the string attribute `out_mesh`;
 - the boolean attribute `resample`;
 - the boolean attribute `evaluate`;
 
The semantics of the attribute is the same as the command line.  

The program also supports the help command as
```
Remesh -h|--help
```

### Remeshing a dataset
The `BatchRemesh` application applies the remeshing algorithm to multiple shapes and/or with multiple options. To run it, please execute the following command
```
BatchRemesh config_file
```
The program accepts a configuration file in `JSON` file format. The file must contain at least the following attributes:
 - `input_dir`, which tells to the program where the meshes of the dataset are located;
 - `output_dir`, which tells to the program where to save the remeshed shapes;
 - `meshes`, which tells to the program on which meshes the remeshing must be applied.
The `meshes` attribute can be a string (meaning the algorithm must be applied on a single mesh), a JSON array of strings (listing the meshes on which the remeshing must be applied) or it can be the string `"*"` (meaning the remeshing must be applied to all the meshes in the input directory).  
**TODO:** the configuration file still misses a support for regex, which is planned to be added in future.  

Additionally, either the attribute `num_samples` or `resolution` must be set. If `num_samples` is set, it must be an integer value and all the meshes will be remeshed to have that number of vertices. If `resolution` is set, it must be a floating point value `0 < r <= 1` and all the meshes will be remeshed to have a fraction `r` of their original vertices.  
If both `num_samples` and `resolution` are set, the configuration file must specify the boolean attribute `fixed_size`. If set to true `num_samples` is used and `resolution` is ignored, while if set to false `resolution` is used and `num_samples` is ignored.  
The attributes `num_samples` and `resolution` can also be lists of values. In that case, the remeshing is applied to all the specified meshes multiple times, one for each specified value of `num_samples` or `resolution`. The output meshes will be divided in subdirectories.  

The configuration file can optionally be fed with the attributes `resample` and `evaluate`, whose meaning is the same as the single run execution.  

The program also generates a CSV file `batch.csv` in the output directory containing the statistics of the meshes, the number of output vertices and the time needed to remesh the shape and to perform every step of the algorithm. If the attribute `evaluate` is set to true, the CSV also contains the evaluation metrics for each shape.  

The program also supports the help command as
```
BatchRemesh -h|--help
```
