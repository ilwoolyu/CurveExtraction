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
    
	std::string desc("TRACE: Topological Graph Representation for Automatic Sulcal Curve Extraction "
					 TRACE_VERSION "\n"
					 "Author: Ilwoo Lyu\n"
					 "Please refer to the following paper for details:\n"
					 "TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction, IEEE Transactions on Medical Imaging, 37(7), 1653-1663, 2018.\n"
					 );

	CLI::App app(desc);

	app.add_option("-i,--inputSurface", input, "Specify input surface")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("--inputPoints", inputPoint, "Specify input sulcal/gyral points prior to the curve delineation (no extension required): *.spoint | *.gpoint")->check(CLI::ExistingFile)->group("Inputs");
	app.add_flag("--sulcus", sulc, "Execute sulcal curve extraction");
	app.add_flag("--gyrus", gyr, "Execute gyral curve extraction");

	app.add_option("-o,--outputCurve", output, "Specify a prefix of output curves (no extension required): *.scurve | *.gcurve | *.scurve.vtk | *.gcurve.vtk")->required()->group("Outputs");
	app.add_flag("--cart", cart, "Write geodesic curve trajectory (Cartesian coorditates)")->group("Outputs");
	app.add_flag("--bary", bary, "Write geodesic curve trajectory (barycentric coorditates)")->group("Outputs");
	app.add_flag("--noVtk", noVtk, "Ignore vtk output")->group("Outputs");
	app.add_flag("--interm", interm, "Specify output sulcal/gyral points prior to the curve delineation, which follow the output name: *.spoint | *.gpoint")->group("Outputs");
	app.add_flag("--nojunction", junc, "Exclude junction points in the output file (scurve or gcurve)")->group("Outputs");

	app.add_option("--gyralSeed", gseed, "Specify initial seed threshold for gyral point detection", true)->group("TRACE parameters");
	app.add_option("--sulcalSeed", sseed, "Specify initial seed threshold for sulcal point detection", true)->group("TRACE parameters");
	app.add_option("--epRad", eprad, "Specify a radius for end point detection", true)->check(CLI::NonNegativeNumber)->group("TRACE parameters");
	app.add_option("--nhDist", nhdist, "Specify a geodesic distance of neighborhoods for curve delineation", true)->check(CLI::NonNegativeNumber)->group("TRACE parameters");
	app.add_option("--lsTh", lsThreshold, "Specify a threshold for line simplification", true)->check(CLI::NonNegativeNumber)->group("TRACE parameters");
	app.add_option("--prune", prune, "Specify a threshold for curve pruning", true)->check(CLI::NonNegativeNumber)->group("TRACE parameters");
	app.add_flag("--simplification", simp, "Execute curve simplification to reduce the number of sulcal points of the delineated curves")->group("TRACE parameters");

	app.add_option("-s,--smoothIter", iter, "Specify NN smoothing: # of iterations", true)->check(CLI::NonNegativeNumber)->group("Preprocessing");
	app.add_option("--smoothTensorIter", iterTensor, "Specify NN tensor matrix smoothing: # of iterations", true)->check(CLI::NonNegativeNumber)->group("Preprocessing");

	app.add_option("--nThreads", nThreads, "Specify the number of OpenMP threads")->check(CLI::NonNegativeNumber)->group("Multi-threading");

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		exit(app.exit(e));
	}
}

