# CurveExtraction

## Description
Automatic sulcal/gyral curve extraction on the human cortical surfaces.
### Input
surface file (.vtk): triangular 3D mesh
### Output
curve files (.scurve or .gcurve): indices of the selected sulcal/gyral points
### Usage
CurveExtraction -i input.vtk -o output --sulcus

## Requirements
<a href="https://github.com/ilwoolyu/MeshLib" target="_blank">MeshLib (general mesh processing)</a><br />
<a href="https://github.com/ilwoolyu/SlicerExecutionModel" target="_blank">SlicerExecutionModel (CLI)</a>

## References
Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction, IEEE Transactions on Medical Imaging, in press.<br />
Lyu, I., Kim, S., Styner, M., Automatic Sulcal Curve Extraction on the Human Cortical Surface, SPIE Medical Imaging 2015, SPIE9413, 94130P-1-94130P-7.
