/*************************************************
*	SulcalPoint.h
*
*	Release: Nov 2013
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

class SulcalPoint
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
		int nValley;
		bool hasSVertex;
	};

	const Mesh *m_mesh;
	Slicer *m_slicer;
	float *m_curvature;
	float *m_direction;
	float *m_shapeIndex;
	float *m_tangent;
	float m_away;
	float m_threshold;
	int m_nValley;
	bool *m_seed;
	bool *m_valley;

public:
	SulcalPoint(void);
	SulcalPoint(const Mesh *mesh, int iterTensor, float threshold);
	~SulcalPoint(void);
	void setSeed(const float threshold = -0.055f);
	void saveValleyPoint(const char *filename);
	float localMaxima(int idx, vector<float *> *list, vector<bool *> *valley, vector<int> *size, const float normal[3] = NULL);
	float away(void);
	int nValley(void);
	bool isSPoint(int idx);
	const bool *getSeed(void);
	const float *getDirection(void);

private:
	void setCurvature(int iterTensor);
	void concave(contourElem *list, const Vector n);
	void findCandidate(contourList *list, float w, float h, const float normal[3], contourElem *begin = NULL, contourElem *end = NULL);
	void deleteContour(contourList *list);
	int maxDeviation(contourElem *begin, contourElem *end, float w, float h);
	float longestAxis(contourElem *list, contourElem **begin, contourElem **end);
	void tanAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3], const float normal[3]);
	void projAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3]);
	void arbitraryAxis(contourElem *list, contourElem **begin, contourElem **end, const float direction[3], const float center[3] = NULL);
	contourList * contour(int idx, const float normal[3]);
};

