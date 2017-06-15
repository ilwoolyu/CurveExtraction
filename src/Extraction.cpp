/*************************************************
*	Extraction.cpp
*
*	Release: Mar 2015
*	Update: Jun 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <cstring>
#include <omp.h>
#include "Extraction.h"
#include "SurfaceUtil.h"

void extraction(string input, string output, string inputPoint, bool interm, float sseed, float gseed, bool sulc, bool gyr, bool simp, int iter, int iterTensor, bool junc, int nThreads)
{
	if (!inputPoint.empty()) nThreads = 1;
	if (sulc) se = new SulcalPoint*[nThreads];
	if (gyr) ge = new GyralPoint*[nThreads];
	
	mesh = new Mesh();
	mesh->openFile(input.c_str());
	if (iter > 0) SurfaceUtil::smoothing(mesh, iter);
	mesh->centering();

	#pragma omp parallel for
	for (int i = 0; i < nThreads; i++)
	{
		if (sulc)
		{
			se[i] = new SulcalPoint(mesh, iterTensor);
			se[i]->setSeed(sseed);
		}

		if (gyr)
		{
			ge[i] = new GyralPoint(mesh, iterTensor);
			ge[i]->setSeed(gseed);
		}
	}
	if (sulc)
	{
		seed_s = se[0]->getSeed();
		direction_s = se[0]->getDirection();
	}
	if (gyr)
	{
		seed_g = ge[0]->getSeed();
		direction_g = ge[0]->getDirection();
	}
	
	int n = mesh->nVertex();

	if (sulc)
	{
		isValley = new bool[n];
		memset(isValley, 0, sizeof(bool) * n);
	}
	if (gyr)
	{
		isRidge = new bool[n];
		memset(isRidge, 0, sizeof(bool) * n);
	}
	
	likelihood = new float[n]; memset(likelihood, 0, sizeof(float) * n);

	cout << "Point Detection.. ";

	if (inputPoint.empty())	// point computation
	{
		cout << endl;
		int done = 0;
		int threadID = 0;
		#pragma omp parallel for private(threadID)
		for (int i = 0; i < n; i++)
		{
			threadID = omp_get_thread_num();
			if (sulc && seed_s[i])
			{
				se[threadID]->localMaxima(i, NULL, NULL, NULL);
				if (se[threadID]->isSPoint(i)) isValley[i] = true;
			}
			if (gyr && seed_g[i])
			{
				ge[threadID]->localMaxima(i, NULL, NULL, NULL);
				if (ge[threadID]->isGPoint(i)) isRidge[i] = true;
			}
			#pragma omp critical
			{
				done++;
				if (done % 1000 == 0)
				{
					cout << "\r" << done << "/" << n;
					fflush(stdout);
				}
			}
		}
		cout << "\r" << n  << "/" << n << endl;

		if (sulc)
		{
			if (interm)
			{
				char spoint[1024]; sprintf(spoint, "%s.spoint", output.c_str());
				FILE *fp = fopen(spoint, "w");
				for (int i = 0; i < mesh->nVertex(); i++)
					fprintf(fp, "%d\n", (int)isValley[i]);
				fclose(fp);
			}
			for (int i = 1; i < nThreads; i++) delete se[i];
			sc = new SulcalCurve(mesh, isValley, direction_s, likelihood);
		}
		
		if (gyr)
		{
			if (interm)
			{
				char gpoint[1024]; sprintf(gpoint, "%s.gpoint", output.c_str());
				FILE *fp = fopen(gpoint, "w");
				for (int i = 0; i < mesh->nVertex(); i++)
					fprintf(fp, "%d\n", (int)isRidge[i]);
				fclose(fp);
			}
			for (int i = 1; i < nThreads; i++) delete ge[i];
			gc = new GyralCurve(mesh, isRidge, direction_g, likelihood);
		}
	}
	else
	{
		if (sulc)
		{
			char spoint[1024]; sprintf(spoint, "%s.spoint", inputPoint.c_str());
			sc = new SulcalCurve(mesh, spoint, direction_s);
			sc->getSeedPoint(isValley);
		}
		
		if (gyr)
		{
			char gpoint[1024]; sprintf(gpoint, "%s.gpoint", inputPoint.c_str());
			gc = new GyralCurve(mesh, gpoint, direction_g);
			gc->getSeedPoint(isRidge);
		}
	}
	cout << "Done" << endl;
	fflush(stdout);
	
	if (sulc)
	{
		sc->run();
		sc->showInfo();
	}
	if (gyr)
	{
		gc->run();
		gc->showInfo();
	}
	fflush(stdout);

	if (simp)
	{
		cout << "Simplifying curves.. ";
		if (sulc) sc->SimplifyCurves();
		if (gyr) gc->SimplifyCurves();
		cout << "Done" << endl;
		fflush(stdout);

		if (sulc) sc->showInfo();
		if (gyr) gc->showInfo();
		fflush(stdout);
	}

	cout << "Writing an output file.. ";
	if (sulc)
	{
		char sout[1024];
		sprintf(sout, "%s.scurve", output.c_str());
		sc->saveSulcalCurves(sout, junc);
		sprintf(sout, "%s.scurve.vtk", output.c_str());
		sc->saveVTK(sout);
	}
	if (gyr)
	{
		char gout[1024];
		sprintf(gout, "%s.gcurve", output.c_str());
		gc->saveGyralCurves(gout, junc);
		sprintf(gout, "%s.gcurve.vtk", output.c_str());
		gc->saveVTK(gout);
	}
	cout << "Done" << endl;
	fflush(stdout);
	
	// free resources
	if (sulc)
	{
		delete se[0];
		delete [] se;
		delete sc;
		delete [] isValley;
	}
	if (gyr)
	{
		delete ge[0];
		delete [] ge;
		delete gc;
		delete [] isRidge;
	}
	delete [] likelihood;
	delete mesh;
}

