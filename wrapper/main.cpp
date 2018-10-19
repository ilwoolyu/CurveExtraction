#include <cstdlib>
#include <omp.h>
#include "main.h"

int main(int argc, char *argv[])
{
	PARSE_ARGS;

	if (argc < 2)
	{
		cout << "Usage: " << argv[0] << " " << "--help" << endl;
		return -1;
	}
	if (!sulc && !gyr)
	{
		cout << "Set extraction on either sulcal or gyral curves! --sulcus / --gyrus" << endl;
		return -1;
	}

	if (nThreads == 0) nThreads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	omp_set_num_threads(nThreads);

	cout << "----- Configuration Summary -----" << endl;
	cout << "- Input surface: " << input << endl;
	cout << "- Output curve: " << output << endl;
	if (!inputPoint.empty()) cout << "- Input point: " << inputPoint << endl;
	cout << "- Output detected points: " << interm << endl;
	cout << "- Sulcal curve extraction: " << sulc << endl;
	cout << "- Gyral curve extraction: " << gyr << endl;
	if (sulc) cout << "- Sulcal seed threshold: " << sseed << endl;
	if (gyr) cout << "- Gyral seed threshold: " << gseed << endl;
	cout << "- Curve simpliciation: " << simp << endl;
	if (iter > 0) cout << "- NN Surface smoothing: " << iter << endl;
	if (iterTensor > 0) cout << "- NN Tensor matrix smoothing: " << iterTensor << endl;
	cout << "- Threshold for line simplification: " << lsThreshold << endl;
	cout << "- Radius for end point detection: " << eprad << endl;
	cout << "- Geodesic distance for curve delineation: " << nhdist << endl;
	cout << "- Pruning distance: " << prune << endl;
	cout << "- Exclusion of junction points: " << junc << endl;
	cout << "- Number of OpenMP threads: " << nThreads << endl;
	cout << "---------------------------------" << endl;

	extraction(input, output, inputPoint, interm, sseed, gseed, sulc, gyr, simp, iter, iterTensor, !junc, nThreads, eprad, nhdist, lsThreshold, prune, noVtk, geodesic, bary);

	return 0;
}
