# TRACE

![image](https://user-images.githubusercontent.com/9325798/47694824-3c39ce00-dbcd-11e8-8bf5-ddd317d773a6.png)

## Description
A proper geometric representation of the cortical regions is a fundamental task for cortical shape analysis and landmark extraction. However, a significant challenge has arisen due to the highly variable, convoluted cortical folding patterns. Here, we propose a novel topological graph representation for automatic sulcal curve extraction (TRACE). In practice, the reconstructed surface suffers from noise influences introduced during image acquisition/surface reconstruction. In the presence of noise on the surface, TRACE determines stable sulcal fundic regions by employing the line simplification method that prevents the sulcal folding pattern from being significantly smoothed out. The sulcal curves are then traced over the connected graph in the determined regions by the Dijkstra’s shortest path algorithm.

## Installation
You can download and compile the source code using <a href="https://cmake.org/">CMake</a>. Or you can pull <a href="https://hub.docker.com/r/ilwoolyu/cmorph/">docker image</a>:
```
$ docker pull ilwoolyu/cmorph:1.0
```
## Usage
### Input
* surface file (.vtk): triangular 3D mesh
### Output
* curve files (.scurve or/and .gcurve): indices of the selected sulcal/gyral points
* barycentric curve files (.scurve.bary or/and .gcurve.bary): barycentric coefficients and edges for geodesic trajectories of sulcal curves - *the most accurate curve information*
* surface file (.vtk): triangular 3D mesh with curve information that stores the number of branch.
### Commands
The following command line will generate "output.scurve" and "output.scurve.vtk"
```
$ CurveExtraction -i input.vtk -o output --sulcus
```
By default, a junction point appears in every branch it belongs to. To allow only one juction point in output:
```
$ CurveExtraction --nojunction
```
To applying spatial or shape operator smoothing:
```
$ CurveExtraction --smoothIter <# of iters> --smoothTensorIter <# of iters>
```
To obtain the complete geodesic trajectories of sulcal curves:
```
$ CurveExtraction --cart (Cartesian coordinates)
$ CurveExtraction --bary (barycentric coordinates)
```
>**Note**: *.scurve only contains curves traced along the extracted sulcal points. The intermediate points may not be fully captured. *These outputs are not recommended for detailed curve representations.*

To ignore vtk output
```
$ CurveExtraction --noVtk
```
>**Note**: vtk outputs are designed for visual quality check.

To enable multi-thread support (OpenMP)
```
$ CurveExtraction --nThreads <# of threads>
```
See more options:
```
$ CurveExtraction -h
```
In Docker, you need a sudo acces. To run, type:
```
$ docker run \
         -v <LOCAL_INPUT_PATH>:/INPUT/ \
         -v <LOCAL_OUTPUT_PATH>:/OUTPUT/ \
         --rm ilwoolyu/cmorph:1.0 \
         CurveExtraction -i /INPUT/input.vtk -o /OUTPUT/output --sulcus
```
Please refer to our papers [[1](#ref1),[2](#ref2)] for technical details (theory, parameters, methodological validation, etc.).

## Requirements for build
<a href="https://github.com/ilwoolyu/MeshLib">MeshLib (general mesh processing)</a><br />
<a href="https://github.com/Slicer/SlicerExecutionModel">SlicerExecutionModel (CLI)</a>

## References
<ol>
<li><a id="ref1"></a>Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, <i>IEEE Transactions on Medical Imaging</i>, 37(7), 1653-1663, 2018.</li>
<li><a id="ref2"></a>Lyu, I., Kim, S., Styner, M., <a href="http://dx.doi.org/10.1117/12.2078291">Automatic Sulcal Curve Extraction on the Human Cortical Surface</a>, <i>SPIE Medical Imaging 2015</i>, SPIE9413, 94130P-1-94130P-7.</li>
</ol>
