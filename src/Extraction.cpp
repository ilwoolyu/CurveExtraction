/*************************************************
*	Extraction.cpp
*
*	Release: Mar 2015
*	Update: Oct 2017
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

void extraction(string input, string output, string inputPoint, bool interm, float sseed, float gseed, bool sulc, bool gyr, bool simp, int iter, int iterTensor, bool junc, int nThreads, float eprad, float nhdist, float lsThreshold, float prune, bool noVTK, bool cart, bool bary)
{
	if (!inputPoint.empty()) nThreads = 1;
	if (sulc) se = new SulcalPoint*[nThreads];
	if (gyr) ge = new GyralPoint*[nThreads];
	
	mesh = new Mesh();
	mesh->openFile(input.c_str());
	if (iter > 0) SurfaceUtil::smoothing(mesh, iter);

	#pragma omp parallel for
	for (int i = 0; i < nThreads; i++)
	{
		if (sulc)
		{
			se[i] = new SulcalPoint(mesh, iterTensor, lsThreshold);
			se[i]->setSeed(sseed);
		}

		if (gyr)
		{
			ge[i] = new GyralPoint(mesh, iterTensor, lsThreshold);
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
	}
	cout << "Done" << endl;
	fflush(stdout);

	if (sulc)
	{
		if (inputPoint.empty())
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
		else
		{
			char spoint[1024]; sprintf(spoint, "%s.spoint", inputPoint.c_str());
			sc = new SulcalCurve(mesh, spoint, direction_s);
			sc->getSeedPoint(isValley);
		}
		sc->setThreshold(nhdist, nhdist, eprad, prune);
		sc->run();
		sc->showInfo();
		if (simp)
		{
			fflush(stdout);
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
		char sout[1024];
		sprintf(sout, "%s.scurve", output.c_str());
		sc->saveSulcalCurves(sout, junc);
		if (!noVTK)
		{
			sprintf(sout, "%s.scurve.vtk", output.c_str());
			sc->saveVTK(sout);
		}
		if (cart)
		{
			sprintf(sout, "%s.scurve.cart", output.c_str());
			sc->saveGeodesicPath(sout, false);
		}
		if (bary)
		{
			sprintf(sout, "%s.scurve.bary", output.c_str());
			sc->saveGeodesicPath(sout, true);
		}
		cout << "Done" << endl;
		delete se[0];
		delete [] se;
		delete sc;
		delete [] isValley;
	}
		
	if (gyr)
	{
		if (inputPoint.empty())
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
		else
		{
			char gpoint[1024]; sprintf(gpoint, "%s.gpoint", inputPoint.c_str());
			gc = new GyralCurve(mesh, gpoint, direction_g);
			gc->getSeedPoint(isRidge);
		}
		gc->setThreshold(nhdist, nhdist, eprad, prune);
		gc->run();
		gc->showInfo();
		if (simp)
		{
			fflush(stdout);
			cout << "Simplifying curves.. ";
			if (sulc) sc->SimplifyCurves();
			if (gyr) gc->SimplifyCurves();
			cout << "Done" << endl;
			fflush(stdout);

			if (sulc) sc->showInfo();
			if (gyr) gc->showInfo();
			fflush(stdout);
		}
		fflush(stdout);
		cout << "Writing an output file.. ";
		fflush(stdout);
		char gout[1024];
		sprintf(gout, "%s.gcurve", output.c_str());
		gc->saveGyralCurves(gout, junc);
		if (!noVTK)
		{
			sprintf(gout, "%s.gcurve.vtk", output.c_str());
			gc->saveVTK(gout);
		}
		if (cart)
		{
			sprintf(gout, "%s.gcurve.cart", output.c_str());
			gc->saveGeodesicPath(gout, false);
		}
		if (bary)
		{
			sprintf(gout, "%s.gcurve.bary", output.c_str());
			gc->saveGeodesicPath(gout, true);
		}
		cout << "Done" << endl;
		fflush(stdout);
		delete ge[0];
		delete [] ge;
		delete gc;
		delete [] isRidge;
	}
		
	delete [] likelihood;
	delete mesh;
}

