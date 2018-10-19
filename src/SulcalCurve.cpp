/*************************************************
*	SulcalCurve.cpp
*
*	Release: Nov 2013
*	Update: Aug 2017
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <algorithm>
#include <cstring>
#include <float.h>
#include "Geom.h"
#include "SulcalCurve.h"

SulcalCurve::SulcalCurve(void)
{
}

SulcalCurve::SulcalCurve(const char *mesh, const bool *valley)
{
}

SulcalCurve::SulcalCurve(const Mesh *mesh, const char *valley, const float *curvature, const float *likelihood)
{
	m_mesh = mesh;
	m_geodesic = new Geodesic(m_mesh);
	int n = m_mesh->nVertex();
	m_list = NULL;

	bool *isValley = new bool[n];
	FILE *fp = fopen(valley, "r");
	for (int i = 0; i < n; i++)
	{
		int flag;
		fscanf(fp, "%d", &flag);
		isValley[i] = (bool)flag;
	}
	fclose(fp);
	
	vector<const float *> point, normal;
	vector<int> id;
	for (int i = 0; i < n; i++)
	{
		if (isValley[i])
		{
			point.push_back(m_mesh->vertex(i)->fv());
			normal.push_back(m_mesh->normal(i)->fv());
			id.push_back(i);
		}
	}

	m_sulcalPoint = new bool[n]; memset(m_sulcalPoint, 0, sizeof(bool) * n);

	m_nPoints = point.size();

	m_candEndPoint = new bool[m_nPoints]; memset(m_candEndPoint, 0, sizeof(bool) * m_nPoints);
	m_curveElem = new curveElem*[m_nPoints];

	for (int i = 0; i < m_nPoints; i++)
	{
		m_curveElem[i] = new curveElem;
		m_curveElem[i]->id = i;
		m_curveElem[i]->vid = id[i];
		m_curveElem[i]->header = NULL;
		m_curveElem[i]->v = point[i];
		m_curveElem[i]->normal = normal[i];
		m_curveElem[i]->isEndPoint = false;
		m_curveElem[i]->isJunction = false;
		m_curveElem[i]->nJunctions = 0;
		m_curveElem[i]->traced = false;
		m_curveElem[i]->deleted = false;
		m_curveElem[i]->direction = (curvature != NULL) ? &curvature[id[i] * 3]: NULL;
		m_curveElem[i]->likelihood = (likelihood != NULL) ? likelihood[id[i]]: 0;
		m_curveElem[i]->orientation = 0.0f;
	}

	delete [] isValley;
}

SulcalCurve::SulcalCurve(const Mesh *mesh, const bool *valley, const float *curvature, const float *likelihood)
{
	m_mesh = mesh;
	m_geodesic = new Geodesic(m_mesh);
	int n = m_mesh->nVertex();
	m_list = NULL;
	
	vector<const float *> point, normal;
	vector<int> id;
	for (int i = 0; i < n; i++)
	{
		if (valley[i])
		{
			point.push_back(m_mesh->vertex(i)->fv());
			normal.push_back(m_mesh->normal(i)->fv());
			id.push_back(i);
		}
	}

	m_sulcalPoint = new bool[n]; memset(m_sulcalPoint, 0, sizeof(bool) * n);

	m_nPoints = point.size();

	m_candEndPoint = new bool[m_nPoints]; memset(m_candEndPoint, 0, sizeof(bool) * m_nPoints);
	m_curveElem = new curveElem*[m_nPoints];

	for (int i = 0; i < m_nPoints; i++)
	{
		m_curveElem[i] = new curveElem;
		m_curveElem[i]->id = i;
		m_curveElem[i]->vid = id[i];
		m_curveElem[i]->header = NULL;
		m_curveElem[i]->v = point[i];
		m_curveElem[i]->normal = normal[i];
		m_curveElem[i]->isEndPoint = false;
		m_curveElem[i]->isJunction = false;
		m_curveElem[i]->nJunctions = 0;
		m_curveElem[i]->traced = false;
		m_curveElem[i]->deleted = false;
		m_curveElem[i]->direction = (curvature != NULL) ? &curvature[id[i] * 3]: NULL;
		m_curveElem[i]->likelihood = (likelihood != NULL) ? likelihood[id[i]]: 0;
		m_curveElem[i]->orientation = 0.0f;
	}
}

SulcalCurve::~SulcalCurve(void)
{
	delete [] m_candEndPoint;
	delete [] m_sulcalPoint;
	delete [] m_curveElem;
	for (int i = 0; i < m_nPoints; i++)
		delete [] m_dist[i];
	delete [] m_dist;
	delete m_geodesic;

	curveList *iter = m_list;
	while (iter != NULL)
	{
		curveList *next = iter->next;
		deleteCurves(iter);
		iter = next;
	}
}

void SulcalCurve::setThreshold(float threshold1, float threshold2, float threshold3, float threshold4)
{
	m_threshold1 = threshold1;
	m_threshold2 = threshold2;
	m_threshold3 = threshold3;
	m_threshold4 = threshold4;
}


void SulcalCurve::run(void)
{
	m_gamma = 1;
	cout << "Geodesic distance.. ";
	fflush(stdout);
	detectNearestPoints(m_threshold1);
	cout << "Done" << endl;
	fflush(stdout);

	/*FILE *fp=fopen("left.ep","w");
	detectEndPoints(m_threshold3, -0);
	for (int i = 0; i < m_nPoints; i++)
	if(m_candEndPoint[i])fprintf(fp, "%d\n",m_curveElem[i]->vid);*/
	
	cout << "Grouping.. ";
	fflush(stdout);
	grouping(m_threshold1, m_threshold2, m_threshold3);
	cout << "Done" << endl;
	fflush(stdout);
	
	cout << "Refining.. ";
	fflush(stdout);
	refineCurves(m_threshold1, m_threshold4);
	cout << "Done" << endl;
	fflush(stdout);
	
	/*for (int i = 0; i < m_nPoints; i++)
		m_curveElem[i]->deleted = false;
	grouping(m_threshold1, m_threshold2, m_threshold3);*/
	//refineCurves(m_threshold1);
	
	cout << "Integrity test.. ";
	fflush(stdout);
	if (testIntegrity()) cout << "Passed" << endl;
	else cout << "Failed" << endl;
	fflush(stdout);

	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		for (int i = 0; i < iter->item.size(); i++)
			m_sulcalPoint[iter->item[i]->id] = true;
}

void SulcalCurve::grouping(float threshold1, float threshold2, float threshold3)
{
	/*do
	{
		group = current;
		deleteNearestPoints(1.5f);
		current = delineation(treshold);
	}
	while (group != current);*/

	int nCurvePrev = 0, nCurve = 0;
	do
	{
		nCurvePrev = nCurve;
		delineation(threshold1, threshold2, threshold3);
		deleteNearestPoints(threshold2);
		nCurve = nCurves();
	} while (nCurvePrev != nCurve);
}

void SulcalCurve::deleteNearestPoints(float threshold)
{
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		deleteNearestPoints(iter, threshold);
}

void SulcalCurve::deleteNearestPoints(curveList *list, float threshold)
{
	for (int i = 0; i < list->item.size(); i++)
	{
		for (int j = 0; j < m_nPoints; j++)
		{
			if (m_curveElem[j]->deleted || m_curveElem[j]->header != NULL) continue;
			if (m_dist[list->item[i]->id][j] < threshold)
			{
				deleteCurveElem(m_curveElem[j]);
				m_candEndPoint[j] = false;
			}
		}
	}
}

int SulcalCurve::delineation(float threshold1, float threshold2, float threshold3)
{
	int n = m_nPoints;

	for (int i = 0; i < n; i++)
		if (!m_curveElem[i]->deleted && m_curveElem[i]->header == NULL)
			m_curveElem[i]->traced = false;

	int nEndPoints;
	float inner = 0;
	do
	{
		nEndPoints = detectEndPoints(threshold3, inner);
		inner -= 0.1;
	}
	while (nEndPoints == 0 && inner >= -1.0f);
	
	int ii=0;
	vector<curveElem> endPointList;
	for (int i = 0; i < n; i++)
	{
		if (!m_candEndPoint[i]) continue;
		//if (m_curveElem[i]->vid == 20667 || m_curveElem[i]->vid == 19008) {ii=i;m_candEndPoint[ii]=false;}
		//if (m_curveElem[i]->vid == 94236 || m_curveElem[i]->vid == 93185) {ii=i;m_candEndPoint[ii]=false;}
		//if (m_curveElem[i]->vid == 88829) continue;
		endPointList.push_back(*m_curveElem[i]);
	}
	sort(endPointList.begin(), endPointList.end());

	curveList *prev;
	for (prev = m_list; prev != NULL && prev->next != NULL; prev = prev->next);

	//for (int i = 0; i < n; i++)
	//for (int i = n - 1; i >= 0; i--)
	for (int c = 0; c < endPointList.size(); c++)
	//for (int c = endPointList.size() - 1; c >= 0; c--)
	{
		int i = endPointList[c].id;
		if (m_curveElem[i]->traced || !m_candEndPoint[i] || m_curveElem[i]->deleted || m_curveElem[i]->header != NULL) continue;
		curveList *list = new curveList;
		curveElem *endPoint = curve_dijkstra(m_curveElem[i], list, threshold1);
		//curveElem *endPoint = curve(m_curveElem[i], list, threshold1);
		if (list->item.size() == 1)
		{
			deleteCurveElem(m_curveElem[i]);
			list->item.clear();
			delete list;
			continue;
		}
		list->next = NULL;
		list->prev = prev;
		list->length = 0;
		if (m_list == NULL) m_list = list;
		else prev->next = list;
		prev = list;
		m_curveElem[i]->isEndPoint = true;

		if (endPoint->header != list) joinCurves(list->item[list->item.size() -1], endPoint);

		//deleteNearestPoints(list, threshold2);
	}

	return nCurves();
}

void SulcalCurve::deleteCurves(curveList *list)
{
	for (int i = 0; i < list->item.size(); i++)
	{
		if (list->item[i]->isJunction) deleteJunction(list->item[i], list);
		else deleteCurveElem(list->item[i]);
	}
			
	if (list->next != NULL) list->next->prev = list->prev;
	if (list->prev != NULL) list->prev->next = list->next;
	else m_list = list->next;

	list->item.clear();
	delete list;
}

bool SulcalCurve::pruneCurves(float threshold)
{
	updateCurveLength();
	float minLen = FLT_MAX;
	curveList *candidate = NULL; 
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		if (iter->item[0]->isEndPoint || iter->item[iter->item.size() - 1]->isEndPoint)
		{
			/*for (int i = 0; i < iter->item.size(); i++)
				m_candEndPoint[iter->item[i]->id] = true;*/
			curveElem *elem = NULL;
			if (iter->item[0]->isJunction) elem = iter->item[0];
			if (iter->item[iter->item.size() - 1]->isJunction) elem = iter->item[iter->item.size() - 1];
			float w = 0;
			if (elem != NULL)
			{
				for (curveList *iter1 = m_list; iter1 != NULL; iter1 = iter1->next)
				{
					if (elem == iter1->item[0])
						w = max(w, exp(fabs(Vector(elem->direction) * Vector(iter1->item[0]->direction)) / 2));
					if (elem == iter1->item[iter1->item.size() - 1])
						w = max(w, exp(fabs(Vector(elem->direction) * Vector(iter1->item[iter1->item.size() - 1]->direction)) / 2));
				}
			}
			else w = 1;
			if (minLen > iter->length * w)
			{
				minLen = iter->length * w;
				candidate = iter;
			}
		}
	}
	if (minLen < threshold)
	{
		deleteCurves(candidate);
		updateCurveLength();
	}
	
	return (minLen < threshold);
}

void SulcalCurve::updateCurveLength(void)
{
	// update length of the curves
	//vector<float> hist;
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		iter->length = 0;
		for (int i = 1; i < iter->item.size(); i++)
			iter->length += m_dist[iter->item[i]->id][iter->item[i - 1]->id];
		/*if (iter->item[0]->isEndPoint || iter->item[iter->item.size() - 1]->isEndPoint)
			hist.push_back(iter->length);*/
	}
	/*sort(hist.begin(), hist.end());
	float lThreshold = hist[(int)(hist.size() * 0.1)];*/
}

void SulcalCurve::deleteJunction(curveElem *elem, curveList *parent)
{
	elem->nJunctions--;
	if (elem->nJunctions == 0)
	{
		for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		{
			if (iter == parent) continue;
			if (elem == iter->item[0] || elem == iter->item[iter->item.size() - 1])
			{
				elem->isJunction = false;
				elem->isEndPoint = true;
				elem->header = iter;
				break;
			}
		}
	}
	else if (elem->nJunctions == 1)
	{
		curveList *list1 = NULL, *list2 = NULL;
		for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		{
			if (iter == parent) continue;
			if (elem == iter->item[0] || elem == iter->item[iter->item.size() - 1])
			{
				if (list1 == NULL) list1 = iter;
				else
				{
					list2 = iter;
					break;
				}
			}
		}
		elem->isJunction = false;
		elem->isEndPoint = false;
		elem->nJunctions = 0;

		if (list1 != NULL && list2 != NULL)
		{
			if (list1->item[list1->item.size() - 1] != elem)
				reverseCurveOrder(list1);
			if (list2->item[0] != elem)
				reverseCurveOrder(list2);
			elem->header = list1;
			
			for (int i = 1; i < list2->item.size(); i++)
			{
				list1->item.push_back(list2->item[i]);
				if (list2->item[i]->header == list2)
					list2->item[i]->header = list1;
			}
			if (list2->next != NULL) list2->next->prev = list2->prev;
			if (list2->prev != NULL) list2->prev->next = list2->next;
			else m_list = list2->next;

			list2->item.clear();
			delete list2;
		}
	}
	else
	{
		for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		{
			if (iter == parent) continue;
			if (elem == iter->item[0] || elem == iter->item[iter->item.size() - 1])
			{
				elem->header = iter;
				break;
			}
		}
	}
}

void SulcalCurve::deleteCurveElem(curveElem *elem)
{
	elem->header = NULL;
	elem->deleted = true;
	elem->isEndPoint = false;
	elem->isJunction = false;
}

void SulcalCurve::detectNearestPoints(float threshold)
{
	m_dist = new float*[m_nPoints];
	for (int i = 0; i < m_nPoints; i++)
		m_dist[i] = new float[m_nPoints];
		
	for (int i = 0; i < m_nPoints; i++)
	{
		m_geodesic->perform_front_propagation(m_curveElem[i]->vid, (double)threshold);
		for (int j = i; j < m_nPoints; j++)
		{
			// Geodesic distance
			float dist = (float)m_geodesic->dist()[m_curveElem[j]->vid];
			
			// Euclidean distance
			//float dist = Vector(m_curveElem[i]->v, m_curveElem[j]->v).norm();

			m_dist[i][j] = (dist <= threshold)? dist: FLT_MAX;
			m_dist[j][i] = m_dist[i][j];
		}
	}
}

int SulcalCurve::detectEndPoints(float threshold, float inner)
{
	int c = 0;
	int n = m_nPoints;
	for (int i = 0; i < n; i++)
	{
		m_candEndPoint[i] = false;
		if (m_curveElem[i]->deleted || m_curveElem[i]->header != NULL) continue;
		Vector V0(m_curveElem[i]->v);

		m_candEndPoint[i] = true;
		for (int j = 0; j < n && m_candEndPoint[i]; j++)
		{
			if (m_dist[i][j] < threshold && !m_curveElem[j]->deleted && m_curveElem[j]->header == NULL)
			{
				Vector V = (Vector(m_curveElem[j]->v) - V0).unit();
				for (int k = j + 1; k < n && m_candEndPoint[i]; k++)
				{
					if (m_dist[i][k] < threshold && !m_curveElem[k]->deleted && m_curveElem[k]->header == NULL)
					{
						float innerProd = V * (Vector(m_curveElem[k]->v) - V0).unit();
						if (innerProd < inner) m_candEndPoint[i] = false;
					}
				}
			}
		}
		if (m_candEndPoint[i]) c++;
	}
	return c;
}

void SulcalCurve::extendCurves(curveElem *elem, float threshold, float inner1, float inner2)
{
	if (!elem->isEndPoint) return;
	int n = m_nPoints;
	float minDist = FLT_MAX;
	int closestID = -1;
	Vector V(elem->v);

	// reverse
	if (elem->header->item[0] == elem) reverseCurveOrder(elem->header);

	bool adjGroup = false;
	for (int i = 0; i < n; i++)
	{
		if (m_curveElem[i]->header == elem->header) continue;

		Vector V1(m_curveElem[i]->v);
		float dist = m_dist[elem->id][i];

		if (dist < threshold)
		{
			if (m_curveElem[i]->header != NULL && m_curveElem[i]->header != elem->header)
			{
				if (!adjGroup)
				{
					minDist = dist;
					adjGroup = true;
				}
			}
			else if (adjGroup) continue;

			float dev = elem->orientation * (V1 - V).unit();
			if (((!adjGroup && dev > inner1) || (adjGroup && dev > inner2)) && minDist >= dist)
			{
				minDist = dist;
				closestID = i;
			}
		}
	}
	
	if (closestID != -1)
	{
		//m_candEndPoint[closestID] = true;
		if (m_curveElem[closestID]->header != NULL && m_curveElem[closestID]->header != elem->header)
		{
			joinCurves(elem, m_curveElem[closestID], 0);
			/*if (m_curveElem[closestID]->isEndPoint)
			{
				curveElem *end = elem->header->item[elem->header->item.size() - 1];
				if (end != elem && !end->isJunction)
					extendCurves(end, threshold);
			}*/
		}
		else
		{
			elem->header->item.push_back(m_curveElem[closestID]);
			elem->isEndPoint = false;
			m_curveElem[closestID]->isEndPoint = true;
			m_curveElem[closestID]->header = elem->header;
			updateOrientation(elem->header, m_curveElem[closestID]);
			m_curveElem[closestID]->deleted = false;
			m_curveElem[closestID]->traced = true;

			extendCurves(m_curveElem[closestID], threshold);
		}
	}
}

void SulcalCurve::extendCurves_dijkstra(curveElem *elem, float threshold, float inner1, float inner2)
{
	if (!elem->isEndPoint) return;
	int n = m_nPoints;
	
	// reverse
	if (elem->header->item[0] == elem) reverseCurveOrder(elem->header);
	
	// init Dijkstra path
	vector<int> Q;
	for (int i = 0; i < n; i++)
	{
		m_curveElem[i]->dijNode = -1;
		m_curveElem[i]->dijDist = FLT_MAX;
		m_curveElem[i]->dijEnd = false;
	}
	m_curveElem[elem->id]->dijNode = elem->id;
	m_curveElem[elem->id]->dijDist = 0;
	Q.push_back(elem->id);
		
	while (!Q.empty())
	{
		float minDist = FLT_MAX;
		int activeID = -1;
		for (int i = 0; i < Q.size(); i++)
		{
			if (m_curveElem[Q[i]]->dijDist < minDist)
			{
				minDist = m_curveElem[Q[i]]->dijDist;
				activeID = Q[i];
			}
		}
		Q.erase(std::remove(Q.begin(), Q.end(), activeID), Q.end()); // remove active point

		Vector V(m_curveElem[activeID]->v);
		Vector B1 = m_curveElem[activeID]->direction;
		if (B1 * m_curveElem[activeID]->orientation < 0) B1 *= -1;

		for (int i = 0; i < n; i++)
		{
			if (activeID == i) continue;

			Vector V1(m_curveElem[i]->v);
			float wdist = m_curveElem[activeID]->dijDist;
			float dist = m_dist[activeID][i];

			if (dist < threshold)
			{
				Vector B2 = m_curveElem[i]->direction;
				
				if (B1 * B2 < 0) B2 *= -1;
				float v[3];
				Coordinate::sphmean(B1.fv(), B2.fv(), v);
				
				float w = exp(m_gamma * ((V1 - V).unit()).cross(v).norm());
				wdist += dist * w;
				bool adjGroup = (m_curveElem[i]->header != NULL && m_curveElem[i]->header != elem->header);
				
				float dev = Vector(v) * (V1 - V).unit();

				if ((!adjGroup && dev > inner1) || (adjGroup && dev > inner2))
				{
					if (m_curveElem[i]->dijDist > wdist)
					{
						m_curveElem[i]->dijDist = wdist;
						m_curveElem[i]->dijNode = activeID;
						if (!m_curveElem[i]->traced)
						{
							Q.erase(std::remove(Q.begin(), Q.end(), i), Q.end());
							Q.push_back(i);
							m_curveElem[i]->orientation = v;
							if ((V1 - V) * Vector(v) < 0) m_curveElem[i]->orientation *= -1;
						}
					}
				}
			}
		}
	}
	
	// check if there exist valid points
	int nValidEndPoints = 0;
	for (int i = 0; i < n; i++)
	{
		if (m_curveElem[i]->header == elem->header || i == elem->id || !m_curveElem[i]->traced || m_curveElem[i]->dijNode == -1) continue;

		m_curveElem[i]->dijEnd = true;
		nValidEndPoints++;
	}
	
	// find the shortest path
	if (nValidEndPoints > 0)
	{
		int target = elem->id;
		float minDist = FLT_MAX;
		for (int i = 0; i < n; i++)
		{
			if (m_curveElem[i]->dijEnd && m_curveElem[i]->dijDist < minDist)
			{
				target = i;
				minDist = m_curveElem[i]->dijDist;
			}
		}
		
		for (curveList *iter = m_list; iter != NULL; iter = iter->next) iter->visit = false;
		if (testLoop(elem->header, m_curveElem[target]->header)) return;

		int id = target;
		vector<int> trace;
		while (id != elem->id)
		{
			trace.push_back(id);
			id = m_curveElem[id]->dijNode;
		}
		std::reverse(trace.begin(), trace.end());
		for (int i = 0; i < trace.size(); i++)
		{
			id = trace[i];
			if (m_curveElem[id]->header != NULL)
			{
				elem->header->item[elem->header->item.size() - 1]->isEndPoint = true;
				joinCurves(elem->header->item[elem->header->item.size() - 1], m_curveElem[id], 0);
				/*if (m_curveElem[id]->isEndPoint)
				{
					curveElem *end = m_curveElem[id]->header->item[m_curveElem[id]->header->item.size() - 1];
					if (end != m_curveElem[id] && !end->isJunction)
						extendCurves_dijkstra(end, threshold);
				}*/
			}
			else
			{
				elem->isEndPoint = false;
				elem->header->item.push_back(m_curveElem[id]);
				m_curveElem[id]->isEndPoint = (i == trace.size() - 1);
				m_curveElem[id]->header = elem->header;
				m_curveElem[id]->deleted = false;
				m_curveElem[id]->traced = true;
			}
		}
		updateOrientation(elem->header);
	}
}

bool SulcalCurve::testLoop(curveList *from, curveList *to)
{
	from->visit = true;
	if (to == NULL) return false;
	if (from == to) return true;
	
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		if (iter->visit) continue;
		if (iter->item[0] == from->item[0] || iter->item[0] == from->item[from->item.size() - 1] || iter->item[iter->item.size() - 1] == from->item[0] || iter->item[iter->item.size() - 1] == from->item[from->item.size() - 1])
			if (testLoop(iter, to)) return true;
	}
	
	return false;
}

void SulcalCurve::updateOrientation(curveList *list, curveElem *_elem)
{
	for (int i = 0; i < list->item.size(); i++)
	{
		if (_elem != NULL && _elem != list->item[i]) continue;

		curveElem *elem = list->item[i];

		if (list->item[i]->isJunction && list->item[i]->header != list) continue;

		Vector V;
		if (i == 0)
		{
			V = Vector(elem->v, list->item[1]->v).unit();
			if (elem->direction != NULL)
			{
				elem->orientation = elem->direction;
				if (elem->orientation * V < 0)
					elem->orientation *= -1;
			}
		}
		else
		{
			V = Vector(list->item[i - 1]->v, elem->v).unit();
			if (elem->direction != NULL)
			{
				elem->orientation = elem->direction;
				if (elem->orientation * V < 0)
					elem->orientation *= -1;
			}
		}

		if (elem->direction == NULL) elem->orientation = V;
		else
		{
			// spherical mean
			float v[3];
			Coordinate::sphmean(elem->orientation.fv(), V.unit().fv(), v);
			elem->orientation = v;
		}
	}
}

void SulcalCurve::reverseCurveOrder(curveList *list)
{
	reverse(list->item.begin(), list->item.end());
	updateOrientation(list);
}

SulcalCurve::curveElem * SulcalCurve::curve(curveElem *current, curveList *header, float threshold, float inner1, float inner2)
{
	int n = m_nPoints;

	int closestID = -1;
	float minDist = FLT_MAX;
	current->traced = true;
	current->header = header;
	
	Vector V(current->v);
	header->item.push_back(current);

	if (header->item.size() > 1) updateOrientation(header, current);
	if (header->item.size() == 2) updateOrientation(header, header->item[0]);
	if (header->item.size() == 1) current->orientation = current->direction;
	
	bool adjGroup = false;
	for (int i = 0; i < n; i++)
	//for (int i = n - 1; i >= 0; i--)
	{
		// prevent a "large" u-turn curve
		/*if (m_curveElem[i]->header == header)
		{
			if (distCurveElem(m_curveElem[i], current) > threshold * 5.0)
			{
				closestID = -1;
				break;
			}
		}*/
		if (m_curveElem[i]->deleted) continue;

		Vector V1(m_curveElem[i]->v);
		float dist = m_dist[current->id][i];

		if (dist < threshold)
		{
			float wdist;
			if (m_curveElem[i]->header != NULL && m_curveElem[i]->header != header)
			{
				wdist = dist;
				if (!adjGroup) // first detection
				{
					minDist = dist;
					adjGroup = true;
				}
			}
			else
			{
				if (adjGroup) continue;
				else wdist = dist * ((V1 - V).unit()).cross(current->orientation).norm();
				//wdist = dist * (V1 - V).unit().cross(current->orientation).norm() * (m_curveElem[i]->likelihood * 2 + 1);
			}

			float dev = current->orientation * (V1 - V).unit();
			if (header->item.size() == 1) dev = fabs(dev);

			if ((!adjGroup && minDist >= wdist && dev > inner1) || (adjGroup && minDist >= wdist && dev > inner2))
			{
				minDist = wdist;
				closestID = i;
			}
		}
	}
	
	if (closestID != -1)
	{
		if (!m_curveElem[closestID]->traced)
		{
			m_candEndPoint[closestID] = false;
			return curve(m_curveElem[closestID], header, threshold);
		}
		else
		{
			//if (m_curveElem[closestID]->header == header)
			current->isEndPoint = true;
			return m_curveElem[closestID];
		}
	}
	
	current->isEndPoint = true;
	return current;
}

SulcalCurve::curveElem * SulcalCurve::curve_dijkstra(curveElem *current, curveList *header, float threshold, float inner1, float inner2)
{
	int n = m_nPoints;

	current->traced = true;
	current->header = header;
	current->orientation = current->direction;
	
	// init Dijkstra path
	vector<int> Q;
	for (int i = 0; i < n; i++)
	{
		m_curveElem[i]->dijNode = -1;
		m_curveElem[i]->dijDist = FLT_MAX;
		m_curveElem[i]->dijEnd = false;
	}
	m_curveElem[current->id]->dijNode = current->id;
	m_curveElem[current->id]->dijDist = 0;
	Q.push_back(current->id);
		
	while (!Q.empty())
	{
		float minDist = FLT_MAX;
		int activeID = -1;
		for (int i = 0; i < Q.size(); i++)
		{
			if (m_curveElem[Q[i]]->dijDist < minDist)
			{
				minDist = m_curveElem[Q[i]]->dijDist;
				activeID = Q[i];
			}
		}
		Q.erase(std::remove(Q.begin(), Q.end(), activeID), Q.end()); // remove active point

		Vector V(m_curveElem[activeID]->v);
		Vector B1 = m_curveElem[activeID]->direction;
		if (B1 * m_curveElem[activeID]->orientation < 0) B1 *= -1;

		for (int i = 0; i < n; i++)
		{
			if (activeID == i) continue;
			if (m_curveElem[i]->deleted) continue;

			Vector V1(m_curveElem[i]->v);
			float wdist = m_curveElem[activeID]->dijDist;
			float dist = m_dist[activeID][i];
			bool adjGroup;

			if (dist < threshold)
			{
				Vector B2 = m_curveElem[i]->direction;
				
				if (B1 * B2 < 0) B2 *= -1;
				float v[3];
				Coordinate::sphmean(B1.fv(), B2.fv(), v);
				
				float w = exp(m_gamma * ((V1 - V).unit()).cross(v).norm());
				//float w = (V1 - V).unit() * Vector(v); w = (w == 0)? FLT_MAX: exp(gamma / w);
				wdist += dist * w;
				adjGroup = (m_curveElem[i]->header != NULL && m_curveElem[i]->header != header);
				
				//float dev = m_curveElem[activeID]->orientation * (V1 - V).unit();
				float dev = Vector(v) * (V1 - V).unit();
				if (activeID == current->id) dev = fabs(dev);

				if ((!adjGroup && dev > inner1) || (adjGroup && dev > inner2))
				{
					if (m_curveElem[i]->dijDist > wdist)
					{
						m_curveElem[i]->dijDist = wdist;
						m_curveElem[i]->dijNode = activeID;
						if (!m_curveElem[i]->traced)
						//if (!m_curveElem[i]->isEndPoint)
						{
							Q.erase(std::remove(Q.begin(), Q.end(), i), Q.end());
							Q.push_back(i);
							m_curveElem[i]->orientation = v;
							if ((V1 - V) * Vector(v) < 0) m_curveElem[i]->orientation *= -1;
						}
					}
				}
			}
		}
	}

	// check if there exist traced points
	int nValidEndPoints = 0;
	int target = current->id;
	for (int i = 0; i < n; i++)
	{
		if (m_curveElem[i]->dijNode == -1 || current->id == i || !m_curveElem[i]->traced) continue;

		m_curveElem[i]->dijEnd = true;
		nValidEndPoints++;
	}
	if (nValidEndPoints > 0)
	{
		// find the maximum distance of the shortest path
		float minDist = FLT_MAX;
		for (int i = 0; i < n; i++)
		{
			if (m_curveElem[i]->dijEnd && m_curveElem[i]->dijDist < minDist)
			{
				target = i;
				minDist = m_curveElem[i]->dijDist;
			}
		}
	}
	// if not, find valid endpoints
	else
	{
		// global endpoints
		for (int i = 0; i < n; i++)
		{
			if (m_curveElem[i]->dijNode == -1) continue;

			m_curveElem[i]->dijEnd = true;
			for (int j = 0; j < n && m_curveElem[i]->dijEnd; j++)
			{
				if (i == j) continue;
				if (m_curveElem[j]->dijNode != -1 && m_dist[i][j] < threshold)
					m_curveElem[i]->dijEnd = (m_curveElem[i]->dijDist > m_curveElem[j]->dijDist);
			}
			if (m_curveElem[i]->dijEnd && (m_candEndPoint[i] || m_curveElem[i]->isEndPoint)) nValidEndPoints++;
		}
		// local endpoints
		if (nValidEndPoints == 0)
		{
			for (int i = 0; i < n; i++)
			{
				if (m_curveElem[i]->dijNode == -1 || !m_candEndPoint[i] || !m_curveElem[i]->isEndPoint) continue;

				m_curveElem[i]->dijEnd = true;
				nValidEndPoints++;
			}
		}
		// if not, find local maximum
		if (nValidEndPoints == 0)
		{
			for (int i = 0; i < n; i++)
			{
				if (m_curveElem[i]->dijNode == -1) continue;

				m_curveElem[i]->dijEnd = true;
				for (int j = 0; j < n && m_curveElem[i]->dijEnd; j++)
				{
					if (i == j) continue;
					if (m_curveElem[j]->dijNode != -1 && m_dist[i][j] < threshold)
						m_curveElem[i]->dijEnd = (m_curveElem[i]->dijDist > m_curveElem[j]->dijDist);
				}
			}
		}
		// find the maximum distance of the shortest path
		float maxDist = 0;
		for (int i = 0; i < n; i++)
		{
			if (m_curveElem[i]->dijEnd && m_curveElem[i]->dijDist > maxDist)
			{
				target = i;
				maxDist = m_curveElem[i]->dijDist;
			}
		}
	}
	
	bool isEndPointMarked = false;
	int id = target;
	while (id != current->id)
	{
		if (!m_curveElem[id]->traced)
		{
			if (!isEndPointMarked)
			{
				isEndPointMarked = true;
				m_curveElem[id]->isEndPoint = true;
			}
			header->item.push_back(m_curveElem[id]);
			m_curveElem[id]->header = header;
			m_curveElem[id]->traced = true;
			m_candEndPoint[id] = false;
		}
		id = m_curveElem[id]->dijNode;
	}
	header->item.push_back(current);
	std::reverse(header->item.begin(), header->item.end());
	m_candEndPoint[current->id] = false;
	
	if (header->item.size() > 1) updateOrientation(header);

	return m_curveElem[target];
}

void SulcalCurve::refineCurves(float threshold1, float threshold2)
{
	memset(m_candEndPoint, 0, sizeof(bool) * m_nPoints);
	int current, previous;
	
	// extend endpoints
	current = nCurves();
	do
	{
		previous = current;
		// separate branches
		for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		{
			//if (iter->item[0]->vid == 82053)
			//if (iter->item[0]->vid == 116928)
			/*{
			curveElem *elem1 = iter->item[0];
			if (elem1->isEndPoint) extendCurves_dijkstra(elem1, threshold1;
			separateBranch();
			curveElem *elem2 = iter->item[iter->item.size() - 1];
			if (elem2->isEndPoint) extendCurves_dijkstra(elem2, threshold1);
			separateBranch();
			}
			else{*/
			curveElem *elem1 = iter->item[iter->item.size() - 1];
			if (elem1->isEndPoint) extendCurves_dijkstra(elem1, threshold1);
			separateBranch();
			curveElem *elem2 = iter->item[0];
			if (elem2->isEndPoint) extendCurves_dijkstra(elem2, threshold1);
			separateBranch();
			//}
		}
		current = nCurves();
	}
	while (previous != current);

	// connect closest endpoints
	//joinCurves(threshold);

	// separate branches
	separateBranch();

	// delete a group if its arclength is below a threshold
	while (pruneCurves(threshold2));
}

void SulcalCurve::joinCurves(float threshold)
{
	int previous;
	int current = nCurves();
	do
	{
		previous = current;
		for (curveList *iter = m_list; iter != NULL;)
		{
			curveList *next = iter->next;
			curveElem *elem1 = iter->item[iter->item.size() - 1];
			curveElem *elem2 = iter->item[0];
			if (elem1->isEndPoint || elem2->isEndPoint)
			{
				curveElem *cand1 = NULL, *cand2 = NULL;
				curveElem *elem = NULL, *cand = NULL;
				float minDist = FLT_MAX;
				for (curveList *subIter = m_list; subIter != NULL; subIter = subIter->next)
				{
					if (iter != subIter)
					{
						if (elem1->isEndPoint) cand1 = closestPoint(elem1, subIter);
						if (elem2->isEndPoint) cand2 = closestPoint(elem2, subIter);
						float dist1 = FLT_MAX;
						if (cand1 != NULL) dist1 = m_dist[cand1->id][elem1->id];
						float dist2 = FLT_MAX;
						if (cand2 != NULL) dist2 = m_dist[cand2->id][elem2->id];

						if (minDist < dist1 && dist1 < dist2)
						{
							minDist = dist1;
							elem = elem1;
							cand = cand1;
						}
						else if (minDist < dist2 && dist1 < dist2)
						{
							minDist = dist2;
							elem = elem2;
							cand = cand2;
						}
					}
				}
				if (elem != NULL && elem->isEndPoint && minDist < threshold)
					joinCurves(elem, cand, 0);
			}
			iter = next;
		}
		current = nCurves();
	}
	while (previous != current);
}

void SulcalCurve::joinCurves(curveElem *elem1, curveElem *elem2, float inner)
{
	if (!elem1->isEndPoint) return;
	curveList *list1 = elem1->header;
	curveList *list2 = elem2->header;
	
	// reverse
	if (elem1 == list1->item[0]) reverseCurveOrder(list1);

	if (elem2->isEndPoint)
	{
		// reverse
		if (elem2 == list2->item[list2->item.size() - 1])
			reverseCurveOrder(list2);

		// merge
		if (elem2->orientation * Vector(elem1->v, elem2->v).unit() > inner &&
			Vector(elem2->v, list2->item[1]->v).unit() * Vector(elem1->v, elem2->v).unit() > inner)
		{
			elem1->isEndPoint = false;
			elem2->isEndPoint = false;

			//list->item.insert(list->item.end(), closeHead->item.begin(), closeHead->item.end());
			for (int j = 0; j < list2->item.size(); j++)
			{
				list1->item.push_back(list2->item[j]);
				if (!list2->item[j]->isJunction || (list2->item[j]->isJunction && list2->item[j]->header == list2))
					list2->item[j]->header = list1;
			}
			updateOrientation(list1, elem2);
			if (list2->next != NULL) list2->next->prev = list2->prev;
			if (list2->prev != NULL) list2->prev->next = list2->next;
			else m_list = list2->next;

			list2->item.clear();
			delete list2;
		}
	}
	else	// branch - junction point
	{
		elem2->isJunction = true;
		elem2->nJunctions++;
		list1->item.push_back(elem2);
		elem1->isEndPoint = false;
	}
}

SulcalCurve::curveElem * SulcalCurve::closestPoint(curveElem *elem, curveList *list)
{
	if (!elem->isEndPoint) return NULL;
	curveElem *elem_opposite;
	if (elem->header->item[0] == elem)
		elem_opposite = elem->header->item[elem->header->item.size() - 1];
	else
		elem_opposite = elem->header->item[0];
	
	if (list->item[0]->header == elem_opposite->header || list->item[list->item.size() - 1]->header == elem_opposite->header)
		return NULL;

	float minDist = FLT_MAX;
	curveElem *elem2 = NULL;

	for (int i = 0; i < list->item.size(); i++)
	{
		float dist = m_dist[elem->id][list->item[i]->id];
		if (minDist > dist)
		{
			minDist = dist;
			elem2 = list->item[i];
		}
	}
	return elem2;
}

void SulcalCurve::separateBranch(void)
{
	curveList *prev;
	for (prev = m_list; prev != NULL && prev->next != NULL; prev = prev->next);
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 1; i < iter->item.size() - 1; i++)
		{
			if (iter->item[i]->isJunction)
			{
				iter->item[i]->nJunctions++;
				curveList *list = new curveList;
				list->next = NULL;
				list->prev = prev;
				list->length = 0;
				prev->next = list;

				int size = iter->item.size();
				for (int j = i; j < size; j++)
				{
					list->item.push_back(iter->item[j]);
					if (iter->item[j]->header == iter)
						iter->item[j]->header = list;
				}

				for (int j = i + 1; j < size; j++)
					iter->item.pop_back();
				
				prev = list;
				break;
			}
		}
	}
}

int SulcalCurve::nCurves(void)
{
	int n = 0;
	for (curveList *iter = m_list; iter != NULL; iter = iter->next) n++;
	
	return n;
}

void SulcalCurve::getSeedPoint(bool *isValley)
{
	for (int i = 0; i < m_nPoints; i++)
		isValley[m_curveElem[i]->vid] = true;
}

void SulcalCurve::saveSulcalPoint(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		fprintf(fp, "%d\n", (int)m_sulcalPoint[i]);
	}
	fclose(fp);
}

void SulcalCurve::SimplifyCurves(float threshold, float inner)
{
	curveList *iter = m_list;
	while (iter != NULL)
	{
		curveList *next = iter->next;
		float arc = 0;
		int i = 0;
		while (i < iter->item.size())
		{
			if (i > 0 && i < iter->item.size() - 1 && !iter->item[i]->isJunction)
			{
				Vector V1(iter->item[i - 1]->v, iter->item[i]->v);
				Vector V2(iter->item[i]->v, iter->item[i + 1]->v);
				arc += m_dist[i - 1][i];
				if (V1.unit() * V2.unit() > inner && arc < threshold)
				{
					deleteCurveElem(iter->item[i]);
					iter->item.erase(iter->item.begin() + i);
					m_sulcalPoint[iter->item[i]->id] = false;
					continue;
				}
				else arc = 0;
			}
			i++;
		}
		if (iter->item.size() < 1) deleteCurves(iter);
		iter = next;
	}
}

void SulcalCurve::saveSulcalCurves(const char *filename, bool incJunc)
{
	FILE *fp = fopen(filename, "w");
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		int nPoint = 0;
		for (int i = 0; i < iter->item.size(); i++)
		{
			if (iter->item[i]->isJunction && iter->item[i]->header != iter)
				if (!incJunc) continue;	// junction
			fprintf(fp, "%d ", iter->item[i]->vid);
			nPoint++;
		}
		if (nPoint > 0) fprintf(fp, "\n");
	}
	fclose(fp);
}

void SulcalCurve::saveGeodesicPath(const char *filename, bool barycentric)
{
	GeodesicPath gp(m_geodesic->dist(), m_mesh);
	vector<float *> path;
	FILE *fp = fopen(filename, "w");
	fprintf(fp, "%d\n", nCurves());
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 0; i < path.size(); i++) delete [] path[i];
		path.clear();
		for (int i = 0; i < iter->item.size(); i++)
		{
			int target = min(i + 1, (int)iter->item.size() - 1);
			m_geodesic->perform_front_propagation(iter->item[i]->vid, iter->item[target]->vid);
			gp.computeGeodesicPath(iter->item[target]->vid);
			for (int j = gp.size() - 1; j > 0; j--)
			{
				const float *p = (!barycentric) ? gp.getPoint(j): gp.getBarycentricPoint(j);
				float *v = new float[3];
				v[0] = p[0]; v[1] = p[1]; v[2] = p[2];
				path.push_back(v);
			}
		}
		fprintf(fp, "%d\n", path.size());
		for (int i = 0; i < path.size(); i++)
		{
			const float *v = path[i];
			if (barycentric)
				fprintf(fp, "%d %d %f\n", (int)v[0], (int)v[1], v[2]);
			else
				fprintf(fp, "%f %f %f\n", v[0], v[1], v[2]);
		}
	}
	fclose(fp);
}

void SulcalCurve::saveVTK(const char *filename)
{
	((Mesh *)m_mesh)->saveFile(filename, "vtk", false);
	FILE *fp = fopen(filename, "a");
	int nElem = 0;
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		nElem += iter->item.size();
	}
	int nCurve = nCurves();
	fprintf(fp, "LINES %d %d\n", nCurve, nCurve + nElem);
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		fprintf(fp, "%d ", iter->item.size());
		for (int i = 0; i < iter->item.size(); i++)
		{
			fprintf(fp, "%d ", iter->item[i]->vid);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	int *label = new int[m_mesh->nVertex()];
	memset(label, 0, sizeof(int) * m_mesh->nVertex());
	fprintf(fp, "POINT_DATA %d\nSCALARS Labels int\nLOOKUP_TABLE default\n", m_mesh->nVertex());
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 0; i < iter->item.size(); i++)
		{
			if (!iter->item[i]->isJunction)
				label[iter->item[i]->vid] = 1;
			else
				label[iter->item[i]->vid] = iter->item[i]->nJunctions + 1;
		}
	}
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		fprintf(fp, "%d ", label[i]);
		if ((i + 1) % 256 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fclose(fp);
	
	delete [] label;
}

void SulcalCurve::showInfo(void)
{
	int nPoints = 0, nCurves = 0, nJunctions = 0;
	memset(m_sulcalPoint, 0, sizeof(bool) * m_nPoints);
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 0; i < iter->item.size(); i++)
		{
			nPoints++;

			if (m_curveElem[iter->item[i]->id]->isJunction)
				nJunctions++;
		}
		nCurves++;
	}
	cout << "# of the sulcal points: " << m_nPoints << endl;
	cout << "# of the connected points: " << nPoints << endl;
	cout << "# of the curves: " << nCurves << endl;
	cout << "# of the junction points: " << nJunctions << endl;
}

float SulcalCurve::distCurveElem(curveElem *elem1, curveElem *elem2)
{
	if (elem1->header != elem2->header) return FLT_MAX;

	curveList *header = elem1->header;
	float dist = FLT_MAX;
	for (int i = 0; i < header->item.size(); i++)
	{
		if (dist < FLT_MAX)
		{
			dist += Vector(header->item[i - 1]->v, header->item[i]->v).norm();
			if (header->item[i] == elem1 || header->item[i] == elem2)
				break;
		}
		else
			if (header->item[i] == elem1 || header->item[i] == elem2)
				dist = 0;
	}
	return dist;
}

bool SulcalCurve::testIntegrity(void)
{
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		if ((iter->item[0]->isEndPoint && iter != iter->item[0]->header) // starting point
			|| (iter->item[0]->isJunction && iter->item[0]->header == NULL)
			|| (!iter->item[0]->isEndPoint && !iter->item[0]->isJunction))
		{
			cout << "Error: Code1 at " << iter->item[0]->id << endl;
			cout << "# of points: " << iter->item.size() << endl;
			cout << "EndPoint: " << iter->item[0]->isEndPoint << endl;
			cout << "Junction: " << iter->item[0]->isJunction << endl;
			cout << "NULL Header: " << (iter->item[0]->header == NULL) << endl;
			return false;
		}
		if ((iter->item[iter->item.size() - 1]->isEndPoint && iter != iter->item[iter->item.size() - 1]->header)
			|| (iter->item[iter->item.size() - 1]->isJunction && iter->item[iter->item.size() - 1]->header == NULL)
			|| (!iter->item[iter->item.size() - 1]->isEndPoint && !iter->item[iter->item.size() - 1]->isJunction))
		{
			cout << "Error: Code2 at " << iter->item[iter->item.size() - 1]->id << endl;
			cout << "# of points: " << iter->item.size() << endl;
			cout << "EndPoint: " << iter->item[iter->item.size() - 1]->isEndPoint << endl;
			cout << "Junction: " << iter->item[iter->item.size() - 1]->isJunction << endl;
			cout << "NULL Header: " << (iter->item[iter->item.size() - 1]->header == NULL) << endl;
			return false;
		}
		for (int i = 1; i < iter->item.size() - 1; i++)	// intermediate points
			if (iter->item[i]->isEndPoint || iter != iter->item[i]->header)
			{
				cout << "Error: Code3 at " << iter->item[i]->id << endl;
				cout << "Index: " << i + 1 << "/" << iter->item.size() << endl;
				cout << "EndPoint: " << iter->item[i]->isEndPoint << endl;
				cout << "Invalid Header: " << (iter != iter->item[i]->header) << endl;
				return false;
			}
		for (int i = 0; i < iter->item.size(); i++)	// junction points
			if ((iter->item[i]->isEndPoint && iter->item[i]->isJunction)
				|| (iter->item[i]->nJunctions == 0 && iter->item[i]->isJunction)
				|| (iter->item[i]->nJunctions > 0 && !iter->item[i]->isJunction)
				|| (iter->item[i]->header == NULL && iter->item[i]->isJunction))
			{
				cout << "Error: Code";
				if (iter->item[i]->isEndPoint && iter->item[i]->isJunction)
					cout << "4-1";
				if (iter->item[i]->nJunctions == 0 && iter->item[i]->isJunction)
					cout << "4-2";
				if (iter->item[i]->nJunctions > 0 && !iter->item[i]->isJunction)
					cout << "4-3";
				if (iter->item[i]->header == NULL && iter->item[i]->isJunction)
					cout << "4-4";
				cout << " at " << iter->item[i]->id << endl;
				return false;
			}
		for (int i = 0; i < iter->item.size(); i++)	// junction points
			if (iter->item[i]->deleted)
			{
				cout << "Error: Code5 at " << iter->item[i]->id << endl;
				return false;
			}
		if (iter->item.size() == 1)	// not a curve
		{
			cout << "Error: Code6 - Not a curve" << endl;
			return false;
		}
	}
	return true;
}
