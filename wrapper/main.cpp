#include <cstdlib>
#include <omp.h>
#include "PARSE_ARGS.h"
#include "main.h"

int main(int argc, char** argv)
{
	PARSE_ARGS(argc, argv);

	if (!sulc && !gyr)
	{
		cout << "Specify at least one extraction type! --sulcus / --gyrus" << endl;
		return -1;
	}

#ifdef _USE_OPENMP
	if (nThreads == 0)
	{
		const char *env = getenv("OMP_NUM_THREADS");
		nThreads = (env != NULL) ? max(atoi(env), 1): 1;
	}
	omp_set_num_threads(nThreads);
#endif

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

	extraction(input, output, inputPoint, interm, sseed, gseed, sulc, gyr, simp, iter, iterTensor, !junc, nThreads, eprad, nhdist, lsThreshold, prune, noVtk, cart, bary);

	return 0;
}
