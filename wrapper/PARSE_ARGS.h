#include <cstring>
#include <vector>
#include "CLI11.hpp"

std::string input;
std::string inputPoint;
std::string output;
bool interm = false;
float gseed = 0.03;
float sseed = -0.05;
float eprad = 4;
float nhdist = 4;
float lsThreshold = 2.5;
float prune = 5.0;
int iter = 0;
int iterTensor = 3;
int nThreads = 0;
bool sulc = false;
bool gyr = false;
bool simp = false;
bool junc = false;
bool noVtk = false;
bool cart = false;
bool bary = false;

void PARSE_ARGS(int argc, char **argv)
{
    
    std::string desc("TRACE: Topological Graph Representation for Automatic Sulcal Curve Extraction v1.1.1\n"
					 "Author: Ilwoo Lyu\n"
					 "Please refer to the following paper for details:\n"
					 "TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction, IEEE Transactions on Medical Imaging, 37(7), 1653-1663, 2018.\n"
					 );

    CLI::App app(desc);

	app.add_option("-i,--inputSurface", input, "provides input surface")->required()->check(CLI::ExistingFile);
	app.add_option("--inputPoints", inputPoint, "provides input sulcal/gyral points prior to the curve delineation (no extension required): *.spoint | *.gpoint")->check(CLI::ExistingFile);
	app.add_option("-o,--outputCurve", output, "provides a prefix of output curves (no extension required): *.scurve | *.gcurve | *.scurve.vtk | *.gcurve.vtk")->required();
	app.add_flag("--interm", interm, "enables output sulcal/gyral points prior to the curve delineation, which follow the output name: *.spoint | *.gpoint");
	app.add_option("--gyralSeed", gseed, "sets initial seed threshold for gyral point detection", true);
	app.add_option("--sulcalSeed", sseed, "sets initial seed threshold for sulcal point detection", true);
	app.add_option("--epRad", eprad, "sets a radius for end point detection", true)->check(CLI::NonNegativeNumber);
	app.add_option("--nhDist", nhdist, "sets a geodesic distance of neighborhoods for curve delineation", true)->check(CLI::NonNegativeNumber);
	app.add_option("--lsTh", lsThreshold, "sets a threshold for line simplification", true)->check(CLI::NonNegativeNumber);
	app.add_option("--prune", prune, "sets a threshold for curve pruning", true)->check(CLI::NonNegativeNumber);
	app.add_option("-s,--smoothIter", iter, "sets NN smoothing: # of iterations", true)->check(CLI::NonNegativeNumber);
	app.add_option("--smoothTensorIter", iterTensor, "sets NN tensor matrix smoothing: # of iterations", true)->check(CLI::NonNegativeNumber);
	app.add_option("--nThreads", nThreads, "sets # of OpenMP threads")->check(CLI::NonNegativeNumber);
	app.add_flag("--sulcus", sulc, "enables sulcal curve extraction");
	app.add_flag("--gyrus", gyr, "enables sulcal curve extraction");
	app.add_flag("--simplification", simp, "enables curve simplification to reduce the number of sulcal points of the delineated curves");
	app.add_flag("--nojunction", junc, "excludes junction points in the output file (scurve or gcurve)");
	app.add_flag("--noVtk", noVtk, "ignores vtk output");
	app.add_flag("--cart", cart, "writes geodesic curve trajectory (Cartesian coorditates)");
	app.add_flag("--bary", bary, "writes geodesic curve trajectory (barycentric coorditates)");

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		exit(app.exit(e));
	}
}

