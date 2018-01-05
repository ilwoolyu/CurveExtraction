# CurveExtraction

## Description
Automatic sulcal/gyral curve extraction on the human cortical surfaces.
### Input
surface file (.vtk): triangular 3D mesh
### Output
curve files (.scurve or/and .gcurve): indices of the selected sulcal/gyral points<br />
surface file (.vtk): triangluar 3D mesh with curve information that stores the number of branch.
### Usage
The following command line will generate "output.scurve" and "output.scurve.vtk"<br />
```
CurveExtraction -i input.vtk -o output --sulcus
```
See more information
```
CurveExtraction -h
```

## Requirements
<a href="https://github.com/ilwoolyu/MeshLib">MeshLib (general mesh processing)</a><br />
<a href="https://github.com/ilwoolyu/SlicerExecutionModel">SlicerExecutionModel (CLI)</a>

## References
<ol>
<li>Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, IEEE Transactions on Medical Imaging, in press.</li>
<li>Lyu, I., Kim, S., Styner, M., <a href="http://dx.doi.org/10.1117/12.2078291">Automatic Sulcal Curve Extraction on the Human Cortical Surface</a>, SPIE Medical Imaging 2015, SPIE9413, 94130P-1-94130P-7.</li>
</ol>
