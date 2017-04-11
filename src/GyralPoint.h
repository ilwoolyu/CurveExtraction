/*************************************************
*	GyralPoint.h
*
*	Release: Oct 2014
*	Update: Apr 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include "Mesh.h"
#include "Slicer.h"

using namespace std;

class GyralPoint
{
private:
	struct contourElem
	{
		float v[3];
		contourElem *next;
		contourElem *prev;
		bool deleted;
		bool selected;
		bool isVertex;
		bool isEndpoint;
		int id;
	};

	struct contourList
	{
		contourElem *data;
		contourList *next;
		int nElem;
		int sid;
		int nRidge;
		bool hasSVertex;
	};

	const Mesh *m_mesh;
	Slicer *m_slicer;
	float *m_curvature;
	float *m_direction;
	float *m_shapeIndex;
	float *m_tangent;
	float m_away;
	int m_nRidge;
	bool *m_seed;
	bool *m_ridge;

public:
	GyralPoint(void);
	GyralPoint(const Mesh *mesh, int iterTensor);
	~GyralPoint(void);
	void setSeed(const float threshold = -0.055f);
	void saveRidgePoint(const char *filename);
	float localMaxima(int idx, vector<float *> *list, vector<bool *> *ridge, vector<int> *size, const float normal[3] = NULL);
	float away(void);
	int nRidge(void);
	bool isGPoint(int idx);
	const bool *getSeed(void);
	const float *getDirection(void);

private:
	void setCurvature(int iterTensor);
	void convex(contourElem *list, const Vector n);
	void findCandidate(contourList *list, float w, float h, const float normal[3], contourElem *begin = NULL, contourElem *end = NULL);
	void deleteContour(contourList *list);
	int maxDeviation(contourElem *begin, contourElem *end, float w, float h);
	float longestAxis(contourElem *list, contourElem **begin, contourElem **end);
	void tanAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3], const float normal[3]);
	void projAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3]);
	void arbitraryAxis(contourElem *list, contourElem **begin, contourElem **end, const float direction[3], const float center[3] = NULL);
	contourList * contour(int idx, const float normal[3]);
};

