/*************************************************
*	SulcalPoint.cpp
*
*	Release: Nov 2013
*	Update: Apr 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <cstring>
#include <stdio.h>
#include <float.h>
#include "SurfaceUtil.h"
#include "SulcalPoint.h"

using namespace std;

SulcalPoint::SulcalPoint(void)
{
	m_mesh = NULL;
	m_curvature = NULL;
	m_direction = NULL;
	m_tangent = NULL;
	m_seed = NULL;
	m_valley = NULL;
}

SulcalPoint::SulcalPoint(const Mesh *mesh, int iterTensor, float threshold)
{
	m_mesh = mesh;
	const int nVertex = m_mesh->nVertex();
	m_curvature = new float[nVertex];
	m_shapeIndex = new float[nVertex];
	m_direction = new float[nVertex * 3];
	m_tangent = new float[nVertex * 3];
	m_seed = new bool[nVertex];
	m_valley = new bool[nVertex]; memset(m_valley, 0, sizeof(bool) * nVertex);
	m_slicer = new Slicer(m_mesh);
	setCurvature(iterTensor);
	m_threshold = threshold;
}

SulcalPoint::~SulcalPoint(void)
{
	delete m_slicer;
	delete [] m_curvature;
	delete [] m_seed;
	delete [] m_direction;
	delete [] m_tangent;
	delete [] m_valley;
	delete [] m_shapeIndex;
}

void SulcalPoint::setSeed(const float threshold)
{
	int nVertex = m_mesh->nVertex();
	int nSeed = 0;

	for (int i = 0; i < nVertex; i++)
	{
		m_seed[i] = (m_curvature[i] < threshold && m_shapeIndex[i] < -0.1);
		//m_seed[i] = (m_curvature[i] < threshold);
		if (m_seed[i])
		{
			nSeed++;
			/*Vertex *v = (Vertex *)m_mesh->vertex(i);
			v->setColor(1, 0, 0);*/
		}
	}

	cout << "# of Seed Points: " << nSeed << endl;
}

void SulcalPoint::setCurvature(int iterTensor)
{
	int nVertex = m_mesh->nVertex();

	float *cmin = new float[nVertex];
	float *cmax = new float[nVertex];
	float **umin = new float*[nVertex];
	float **umax = new float*[nVertex];
	float *u = new float[nVertex * 3 * 2];
	for (int i = 0; i < nVertex; i++)
	{
		umin[i] = &u[i * 6];
		umax[i] = &u[i * 6 + 3];
	}
	SurfaceUtil::curvature(m_mesh, cmin, cmax, umin, umax, iterTensor);

	for (int i = 0; i < nVertex; i++)
	{
		m_tangent[i * 3] = umin[i][0];
		m_tangent[i * 3 + 1] = umin[i][1];
		m_tangent[i * 3 + 2] = umin[i][2];

		m_direction[i * 3] = umax[i][0];
		m_direction[i * 3 + 1] = umax[i][1];
		m_direction[i * 3 + 2] = umax[i][2];

		m_curvature[i] = cmin[i];
		if (cmax[i] == cmin[i]) m_shapeIndex[i] = (cmax[i] + cmin[i] > 0)? 1: -1;
		else m_shapeIndex[i] = 2.0f / PI * atan((cmax[i] + cmin[i]) / (cmax[i] - cmin[i]));
	}

	for (int i = 0; i < nVertex; i++)
	{
		Vector N = Vector(&m_direction[i * 3]).cross(Vector(&m_tangent[i * 3])); N.unit();
		if (N * Vector(m_mesh->normal(i)->fv()) < 0)
		{
			m_tangent[i * 3 + 0] *= -1;
			m_tangent[i * 3 + 1] *= -1;
			m_tangent[i * 3 + 2] *= -1;
		}
	}

	delete [] cmin;
	delete [] cmax;
	delete [] umin;
	delete [] umax;
	delete [] u;
}

const bool * SulcalPoint::getSeed(void)
{
	return m_seed;
}

const float * SulcalPoint::getDirection(void)
{
	return m_direction;
}

int SulcalPoint::nValley(void)
{
	return m_nValley;
}

float SulcalPoint::away(void)
{
	return m_away;
}

float SulcalPoint::localMaxima(int idx, vector<float *> *list, vector<bool *> *valley, vector<int> *size, const float normal[3])
{
	m_away = FLT_MAX;
	m_nValley = 0;
	float planeNormal[3];

	if (normal != NULL)
		memcpy(planeNormal, normal, sizeof(float) * 3);
	else
		memcpy(planeNormal, &m_direction[idx * 3], sizeof(float) * 3);

	//cout << "Contour Delineation.. ";
	Vector dn(planeNormal); dn.unit();
	contourList *cList = contour(idx, dn.fv());

	//cout << "Done" << endl;

	int g = 0;
	for (contourList *liter = cList; liter != NULL; liter = liter->next, g++)
	{
		if (!liter->hasSVertex) continue;

		int n = liter->nElem;
		//cout << "Contour Size: " << n << endl;
		float *lElem = (list != NULL) ? new float[n * 3]: NULL;
		bool *vElem = (valley != NULL) ? new bool[n]: NULL;
		if (valley != NULL) memset(vElem, 0, sizeof(bool) * n);

		int index = 0;
		contourElem *iter = liter->data;
		contourElem *query;
		int i = 0;
		do
		{
			const float *id = m_slicer->getEdgeID(iter->id, g);
			if (id[0] == liter->sid && id[1] == liter->sid)
			{
				index = i;
				query = iter;
			}
			if (list != NULL) memcpy(&lElem[i * 3], iter->v, sizeof(float) * 3);
			iter = iter->next;
			i++;
		}
		while (iter != liter->data);

		for (int k = 0; k < 2; k++)
		{
			contourElem *begin, *end;
			if (k == 0) longestAxis(liter->data, &begin, &end);
			else if (k == 1) tanAxis(liter->data, &begin, &end, &m_tangent[liter->sid * 3], &m_direction[liter->sid * 3]);
			//else projAxis(liter->data, &begin, &end, &m_tangent[liter->sid * 3]);
			//cout << "Valley Point Detection.. ";
			//findCandidate(liter, 0.0f, 2.0f, dn.fv(), begin, end);		// freesurfer mid smooth
			//findCandidate(liter, 0.0f, 1.0f, dn.fv(), begin, end);		// freesurfer pial smooth
			//findCandidate(liter, 0.0f, 1.5f, dn.fv(), begin, end);		// freesurfer pial smooth s20
			//findCandidate(liter, 0.0f, 2.0f, dn.fv(), begin, end);		// civet mid
			findCandidate(liter, 0.0f, m_threshold, dn.fv(), begin, end);
			//cout << "Done" << endl;
			if (query->selected) break;
		}

		if (query->selected)
		{
			m_away = 0;
		}
		else
		{
			iter = liter->data;
			do
			{
				if (!iter->deleted && iter->selected)
				{
					liter->nValley++;
					float dist = Vector(iter->v, m_mesh->vertex(liter->sid)->fv()).norm();
					if (dist < m_away) m_away = dist;
					//if (valley != NULL) vElem[i] = true;
				}
				iter = iter->next;
			}
			while (iter != liter->data);
		}

		//if (m_away < 0.15f)	// freesurfer mid/pial smooth
		if (m_away < 0.5f)	// freesurfer pial smooth s20
		{
			m_valley[idx] = true;
			if (valley != NULL) vElem[index] = true;
		}
		m_nValley += liter->nValley;
		if (size != NULL) size->push_back(i);
		if (list != NULL) list->push_back(lElem);
		if (valley != NULL) valley->push_back(vElem);
	}
	
	const Vertex &v = *m_mesh->vertex(idx);
	Vector p(v.fv());
	float dist = p * Vector(planeNormal);

	deleteContour(cList);

	return dist;
}

void SulcalPoint::deleteContour(contourList *list)
{
	int g = 0;
	for (contourList *liter = list; liter != NULL; liter = liter->next) g++;
	
	contourList *liter = list;
	for (int i = 0; i < g; i++)
	{
		int n = liter->nElem;
		contourList *next = liter->next;

		contourElem *iter = liter->data;
		for (int j = 0; j < n; j++)
		{
			contourElem *nextItem = iter->next;
			delete iter;
			iter = nextItem;
		}
		delete liter;
		liter = next;
	}
}

SulcalPoint::contourList * SulcalPoint::contour(int idx, const float normal[3])
{
	const Vertex &v = *m_mesh->vertex(idx);
	Vector p(v.fv()), dp(normal);
	float d = -p * normal;

	int nGroup = m_slicer->slicing(normal[0], normal[1], normal[2], d);

	contourList *cList = NULL, *cPrev = NULL;

	//cout << "# of Contours: " << nGroup << endl;

	for (int g = 0; g < nGroup; g++)
	{
		contourList *list = new contourList();
		list->next = NULL;
		list->hasSVertex = false;
		list->sid = -1;
		list->nValley = 0;
		if (cList == NULL) cList = list;
		if (cPrev != NULL) cPrev->next = list;

		int n = m_slicer->size(g);
		float *slice = new float[n * 3];
		m_slicer->getSlice(slice, g);
		contourElem *prev = NULL, *next = NULL, *begin = NULL;

		for (int i = 0; i < n; i++)
		{
			contourElem *elem = new contourElem();
			memcpy(elem->v, &slice[i * 3], sizeof(float) * 3);
			elem->deleted = false;
			elem->selected = false;
			if (begin == NULL) begin = elem;
			else prev->next = elem;
			elem->prev = prev;
			elem->next = NULL;
			elem->id = i;
			elem->isVertex = false;
			elem->isEndpoint = false;
			const float *id = m_slicer->getEdgeID(i, g);
			if (id[0] == id[1] && id[0] == idx)
			{
				list->hasSVertex = true;
				list->sid = idx;
				list->data = elem;
				elem->isVertex = true;
			}
			prev = elem;
		}
		prev->next = begin;
		begin->prev = prev;

		list->nElem = n;
		if (!list->hasSVertex) list->data = begin;
		cPrev = list;

		delete [] slice;
	}

	return cList;
}

void SulcalPoint::findCandidate(contourList *list, float w, float h, const float normal[3], contourElem *begin, contourElem *end)
{
	//if (list->hasSVertex)
	//	arbitraryAxis(list->data, &begin, &end, &m_tangent[list->sid * 3]);
	//	//tanAxis(list->data, &begin, &end, &m_tangent[list->sid * 3], &m_direction[list->sid * 3]);
	//else	// the two farthest points if fails
	//	longestAxis(list->data, &begin, &end);
	//	//cout << "Axis Length: " << axisLen << endl;
	
	if (begin == NULL || end == NULL) longestAxis(list->data, &begin, &end);
	//tanAxis(list->data, &begin, &end, &m_tangent[list->sid * 3], &m_direction[list->sid * 3]);

	// max deviations
	int count = 0;

	for (int i = 0; i < 2; i++)
	{
		if (i == 1) swap(begin, end); // need the opposite direction
		for (contourElem *iter = begin->next; iter != end; iter = iter->next)
		{
			if (iter->isVertex)
			{
				count = maxDeviation(begin, end, w, h);
				/*for (contourElem *iter = begin->next; iter != end; iter = iter->next)
					if (!iter->selected) iter->deleted = true;*/
				break;
					/*else iter->selected = false;
				if (count > 0)
				{
					count = maxDeviation(begin, end, 10, 2);
					for (contourElem *iter = begin->next; iter != end; iter = iter->next)
						if (!iter->selected) iter->deleted = true;
				}*/
			}
		}
	}

	// concave
	//if (count > 0) concave(list->data, -Vector(normal));

	// a new starting point
	list->data = end;
	begin->isEndpoint = true;
	end->isEndpoint = true;
}

int SulcalPoint::maxDeviation(contourElem *begin, contourElem *end, float w, float h)
{
	Vector v0 = begin->v;
	Vector v1 = end->v;
	float maxLen = (v1 - v0).norm();
	Vector axis = (v1 - v0) / maxLen;

	if (maxLen < w) return 0;

	float height, maxHeight = 0;

	contourElem *maxDev;
	for (contourElem *iter = begin->next; iter != end; iter = iter->next)
	{
		Vector v = iter->v;
		if (iter->deleted) continue;
		height = (v0 - v).cross(axis).norm();
		if (maxHeight < height)
		{
			maxHeight = height;
			maxDev = iter;
		}
	}
	if (maxHeight < h) return 0;

	maxDev->selected = true;

	if (maxDev->isVertex) return 1;

	return 1 + maxDeviation(begin, maxDev, w, h) + maxDeviation(maxDev, end, w, h);
}

void SulcalPoint::tanAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3], const float normal[3])
{
	*begin = *end = NULL;

	// ortho tangent vector
	Vector T(tangent), N(normal); T.unit(); N.unit();
	Vector B = T.cross(N); B.unit();
	arbitraryAxis(list, begin, end, B.fv());
	
	Vector P1((*end)->v), P2((*begin)->v);
	Vector B1 = P1 - P2;
	float len = B1.norm() * 0.5f;
	B1.unit();

	Vector P = Vector(list->v) + B1 * len;

	// make sure P is inside of contour
	if ((P2 - P) * (P1 - P) > 0)
		P = Vector(list->v) - B1 * len;

	arbitraryAxis(list, begin, end, tangent, P.fv());
}

void SulcalPoint::arbitraryAxis(contourElem *list, contourElem **begin, contourElem **end, const float direction[3], const float center[3])
{
	*begin = *end = NULL;
	contourElem *iter;

	// tangent vector
	Vector T(direction);
	T.unit();

	// sign changes
	Vector P;
	if (center == NULL)
	{
		P = list->v;

		// end points
		float maxDist = 0;
		for (iter = list; iter->next != list; iter = iter->next)
		{
			if (!iter->deleted)
			{
				Vector Q(P.fv(), iter->v);
				float dist = fabs(Q * T);
				if (dist >= maxDist)
				{
					maxDist = dist;
					*end = iter;
				}
			}
		}
		*begin = list;
		if (*end == NULL) *end = list;
	}
	else
	{
		P = center;

		float minDist = FLT_MAX, maxDist = -FLT_MAX;
		vector<contourElem *> cand;
		iter = list;
		do 
		{
			if (!iter->deleted)
			{
				Vector D1(P.fv(), iter->v);
				Vector D2(P.fv(), iter->next->v);
				float sign = T.cross(D1) * T.cross(D2);

				if (sign <= 0) cand.push_back(iter);
			}
			iter = iter->next;
		}
		while (iter != list);

		for (int i = 0; i < cand.size(); i++)
		{
			float dist = Vector(P.fv(), cand[i]->v).norm();
			if (minDist > dist)
			{
				minDist = dist;
				*begin = cand[i];
			}
			else if (maxDist < dist)
			{
				maxDist = dist;
				*end = cand[i];
			}
		}

		if (*begin == NULL) *begin = list;
		if (*end == NULL) *end = list;
	}
	//MathVector P;
	//if (center == NULL)
	//{
	//	P = list->v;

	//	iter = list;
	//	do 
	//	{
	//		if (!iter->deleted)
	//		{
	//			MathVector D1 = MathVector(iter->v) - P;
	//			MathVector D2 = MathVector(iter->next->v) - P;
	//			float sign = T.cross(D1) * T.cross(D2);

	//			if (sign < 0) cand.push_back(iter);
	//		}
	//		iter = iter->next;
	//	}
	//	while (iter != list);

	//	// end points
	//	float maxDist = 0;
	//	for (int i = 0; i < cand.size(); i++)
	//	{
	//		MathVector Q = MathVector(cand[i]->v) - P;
	//		float dist = Q * T;
	//		if (dist >= maxDist)
	//		{
	//			maxDist = dist;
	//			*end = cand[i];
	//		}
	//	}
	//	*begin = list;
	//}
	//else
	//{
	//	P = center;

	//	iter = list;
	//	do 
	//	{
	//		if (!iter->deleted)
	//		{
	//			MathVector D1 = MathVector(iter->v) - P;
	//			MathVector D2 = MathVector(iter->next->v) - P;
	//			float sign = T.cross(D1) * T.cross(D2);

	//			if (sign <= 0) cand.push_back(iter);
	//		}
	//		iter = iter->next;
	//	}
	//	while (iter != list);

	//	// end points
	//	float minDist = FLT_MAX, maxDist = -FLT_MAX;
	//	for (int i = 0; i < cand.size(); i++)
	//	{
	//		MathVector Q = MathVector(cand[i]->v) - P;
	//		float dist = Q * T;
	//		if (dist <= minDist)
	//		{
	//			minDist = dist;
	//			*begin = cand[i];
	//		}
	//		if (dist >= maxDist)
	//		{
	//			maxDist = dist;
	//			*end = cand[i];
	//		}
	//	}
	//	if (*begin == NULL) *begin = list;
	//	if (*end == NULL) *end = list;
	//}
}

void SulcalPoint::projAxis(contourElem *list, contourElem **begin, contourElem **end, const float tangent[3])
{
	*begin = *end = NULL;
	contourElem *iter;

	// tangent vector
	Vector T(tangent);
	T.unit();

	// end points
	float minDist = FLT_MAX, maxDist = -FLT_MAX;
	iter = list;
	do 
	{
		if (iter->deleted) continue;
		Vector Q(list->v, iter->v);
		float dist = Q * T;
		if (dist <= minDist)
		{
			minDist = dist;
			*begin = iter;
		}
		if (dist >= maxDist)
		{
			maxDist = dist;
			*end = iter;
		}
	}
	while (iter != list);
	if (*begin == NULL) *begin = list;
	if (*end == NULL) *end = list;
}

float SulcalPoint::longestAxis(contourElem *list, contourElem **begin, contourElem **end)
{
	contourElem *iter;
	Vector v0 = list->v, v1;
	float len, maxLen = -1;
	*begin = list;

	iter = list;
	do 
	{
		if (iter->deleted)
		{
			iter = iter->next;
			continue;
		}
		len = (v0 - iter->v).norm();
		if (maxLen < len)
		{
			maxLen = len;
			v1 = iter->v;
			*end = iter;
		}
		iter = iter->next;
	}
	while (iter != list);

	iter = list;
	do
	{
		if (iter->deleted)
		{
			iter = iter->next;
			continue;
		}
		len = (v1 - iter->v).norm();
		if (maxLen < len)
		{
			maxLen = len;
			*begin = iter;
		}
		iter = iter->next;
	}
	while (iter != list);

	return maxLen;
}

void SulcalPoint::concave(contourElem *list, const Vector n)
{
	Vector v[5];
	contourElem *iter, *prev, *next;
	contourElem *start = NULL;
	iter = list;
	while (iter->deleted) iter = iter->next;

	while (iter != start)
	{
		prev = iter->prev;
		next = iter->next;
		while (prev->deleted) prev = prev->prev;
		while (next->deleted) next = next->next;
		if (prev == next)
		{
			iter->deleted = true;
			iter->selected = false;
			prev->deleted = true;
			prev->selected = false;
			break;
		}

		v[0] = prev->v;
		v[1] = iter->v;
		v[2] = next->v;
		v[3] = v[0] - v[1]; v[3].unit();
		v[4] = v[2] - v[1]; v[4].unit();

		if (v[3] * v[4] > -1) // an angle between two points
		{
			if (v[3].cross(v[4]) * n < 0)
				iter->selected = true; // concave
			else
				iter->selected = false; // convex
			if (start == NULL) start = iter;
			iter = next;
		}
		else
		{
			iter->deleted = true;
			iter->selected = false;
			// backtracking
			iter = prev;
			if (iter == start) start = NULL;
			//iter = next;
		}
	}
}

bool SulcalPoint::isSPoint(int idx)
{
	return m_valley[idx];
}

void SulcalPoint::saveValleyPoint(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		fprintf(fp, "%d\n", (int)m_valley[i]);
	}
	fclose(fp);
}
