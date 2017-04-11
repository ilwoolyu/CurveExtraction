/*************************************************
*	Extraction.h
*
*	Release: Mar 2015
*	Update: Apr 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include "GyralPoint.h"
#include "GyralCurve.h"
#include "SulcalPoint.h"
#include "SulcalCurve.h"

Mesh *mesh;
SulcalPoint **se;
SulcalCurve *sc;
GyralPoint **ge;
GyralCurve *gc;

const bool *seed_g, *seed_s;
const float *direction_g;
const float *direction_s;
bool *isValley;
bool *isRidge;
float *likelihood;

