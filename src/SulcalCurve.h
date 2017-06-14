/*************************************************
*	SulcalCurve.h
*
*	Release: Nov 2013
*	Update: Mar 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <vector>
#include "Mesh.h"
#include "Geodesic/Geodesic.h"

using namespace std;

class SulcalCurve
{
private:
	struct curveList;
	struct curveElem;

	struct curveElem
	{
		const float *v;
		const float *normal;
		const float *direction;	// principal direction
		int id;
		int vid;
		int nJunctions;
		bool isEndPoint;	// is endpoint in the current group?
		bool isJunction;	// is endpoint a junction?
		bool traced;
		bool deleted;
		curveList *header;
		Vector orientation;
		float likelihood;
		bool operator < (const curveElem& elem) const
		{
			return (Vector(v) < Vector(elem.v));
		};
	};
	struct curveList
	{
		vector<curveElem *> item;
		curveList *next;
		curveList *prev;
		float length;
	};

	bool *m_candEndPoint;	// candidate
	bool *m_sulcalPoint;	// grouped points
	int m_nPoints;
	curveList *m_list;
	const Mesh *m_mesh;
	curveElem **m_curveElem;
	Geodesic *m_geodesic;
	float **m_dist;

public:
	SulcalCurve(void);
	SulcalCurve(const char *mesh, const bool *valley);
	SulcalCurve(const Mesh *mesh, const char *valley, const float *curvature = NULL, const float *likelihood = NULL);
	SulcalCurve(const Mesh *mesh, const bool *valley, const float *curvature = NULL, const float *likelihood = NULL);
	~SulcalCurve(void);
	void run(void);
	void grouping(float threshold1 = 2.5f, float threshold2 = 2.5f, float threshold3 = 2.0f);	// th1: delineation, th2: neighbor to be deleted, th3: endpoint radius
	void refineCurves(float threshold = 3.0f);
	void getSeedPoint(bool *isValley);
	void saveSulcalPoint(const char *filename);
	void saveSulcalCurves(const char *filename, bool incJunc = true);
	void saveVTK(const char *filename);
	void showInfo(void);
	void SimplifyCurves(float threshold = 2.5f, float inner = cos(PI / 10.0f));
	int nCurves(void);
	bool testIntegrity(void);

private:
	int delineation(float threshold1, float threshold2, float threshold3);	// th3: endpoint radius
	int detectEndPoints(float threshold, float inner = 0);
	void detectNearestPoints(float threshold = 3.0f);
	void reverseCurveOrder(curveList *list);
	void joinCurves(curveElem *elem1, curveElem *elem2, float inner = cos(PI / 4.0f));
	void joinCurves(float threshold);
	void extendCurves(curveElem *elem, float threshold, float inner1 = cos(PI / 3.0f), float inner2 = cos(PI / 3.0f));
	void deleteCurves(curveList *list);
	void deleteCurves(float threshold);
	void deleteCurveElem(curveElem *elem);
	void deleteNearestPoints(float threshold);
	void deleteNearestPoints(curveList *list, float threshold);
	void updateOrientation(curveList *list, curveElem *_elem = NULL);
	void separateBranch(void);
	void findHeader(curveElem *elem);
	float distCurveElem(curveElem *elem1, curveElem *elem2);
	curveElem * closestPoint(curveElem *elem, curveList *list);
	curveElem * curve(curveElem *current, curveList *header, float threshold, float inner1 = cos(PI / 4.0f), float inner2 = 0);
};

