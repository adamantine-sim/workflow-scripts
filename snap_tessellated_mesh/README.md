# snap\_to\_tessellated\_mesh
`snap_to_tesselated_mesh` snaps a tessellated vtk mesh to the geometries surface given by an iges file.

## Installation
Installing `snap_to_tessellated_mesh` requires:
* deal.II: 9.5 or later
* OpenCASCADE (compatible with deal.II)

You need to compile deal.II with OpenCASCADE support.

To configure `snap_to_tessellated_mesh` use:
```CMake
cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D DEAL_II_DIR=/path/to/dealii \
/path/to/source/dir
```
Then simply use `make`. This will compile `snap_to_tessellated_mesh` and create an executable
called `snap_to_tessellated_mesh.exe` in the build directory

## Run
After compiling `snap_to_tessellated_mesh`, you can create a snapped vtk mesh file via
```bash
./snap_to_tessellated_mesh input.vtk input.iges output.vtk
```
Note that the name of the input and output files is totally arbitrary, `input.vtk` is as
valid as `input_vtk`.

## License
`snap_to_tessellated_mesh` is distributed under the 3-Clause BSD License.

## Questions
If you have any question, find a bug, or have feature request please open an
issue.
