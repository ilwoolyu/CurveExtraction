# CurveExtraction

## Description
A proper geometric representation of the cortical regions is a fundamental task for cortical shape analysis and landmark extraction. However, a significant challenge has arisen due to the highly variable, convoluted cortical folding patterns. Here, we propose a novel topological graph representation for automatic sulcal curve extraction (TRACE). In practice, the reconstructed surface suffers from noise influences introduced during image acquisition/surface reconstruction. In the presence of noise on the surface, TRACE determines stable sulcal fundic regions by employing the line simplification method that prevents the sulcal folding pattern from being significantly smoothed out. The sulcal curves are then traced over the connected graph in the determined regions by the Dijkstra’s shortest path algorithm.

### Input
* surface file (.vtk): triangular 3D mesh
### Output
* curve files (.scurve or/and .gcurve): indices of the selected sulcal/gyral points<br />
* surface file (.vtk): triangular 3D mesh with curve information that stores the number of branch.
### Usage
The following command line will generate "output.scurve" and "output.scurve.vtk"<br />
```
CurveExtraction -i input.vtk -o output --sulcus
```
To enable multi-thread support (OpenMP)<br />
```
CurveExtraction --nThreads <# of threads>
```
To ignore vtk output<br />
```
CurveExtraction --novtk
```
To allow only one juction point at multiple branch in output<br />
```
CurveExtraction --nojunction
```
The actual geodesic trajectories of sulcal curves can be obtained:<br />
```
CurveExtraction --geodesic
```
See more options:
```
CurveExtraction -h
```
Please refer to our papers [1,2] for technical details (theory, parameters, methodological validation, etc.).
## Requirements
<a href="https://github.com/ilwoolyu/MeshLib">MeshLib (general mesh processing)</a><br />
<a href="https://github.com/Slicer/SlicerExecutionModel">SlicerExecutionModel (CLI)</a>

## References
<ol>
<li>Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, <i>IEEE Transactions on Medical Imaging</i>, 37(7), 1653-1663, 2018.</li>
<li>Lyu, I., Kim, S., Styner, M., <a href="http://dx.doi.org/10.1117/12.2078291">Automatic Sulcal Curve Extraction on the Human Cortical Surface</a>, <i>SPIE Medical Imaging 2015</i>, SPIE9413, 94130P-1-94130P-7.</li>
</ol>
