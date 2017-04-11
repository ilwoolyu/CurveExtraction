/*************************************************
*	GyralCurve.cpp
*
*	Release: Oct 2014
*	Update: Oct 2016
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
#include "GyralCurve.h"


GyralCurve::GyralCurve(void)
{
}

GyralCurve::GyralCurve(const char *mesh, const bool *ridge)
{
}

GyralCurve::GyralCurve(const Mesh *mesh, const char *ridge, const float *curvature, const float *likelihood)
{
	m_mesh = mesh;
	int n = m_mesh->nVertex();
	m_list = NULL;

	bool *isRidge = new bool[n];
	FILE *fp = fopen(ridge, "r");
	for (int i = 0; i < n; i++)
	{
		int flag;
		fscanf(fp, "%d", &flag);
		isRidge[i] = (bool)flag;
	}
	fclose(fp);
	
	vector<const float *> point, normal;
	vector<int> id;
	for (int i = 0; i < n; i++)
	{
		if (isRidge[i])
		{
			point.push_back(m_mesh->vertex(i)->fv());
			normal.push_back(m_mesh->normal(i)->fv());
			id.push_back(i);
		}
	}

	m_gyralPoint = new bool[n]; memset(m_gyralPoint, 0, sizeof(bool) * n);

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

	delete [] isRidge;
}

GyralCurve::GyralCurve(const Mesh *mesh, const bool *ridge, const float *curvature, const float *likelihood)
{
	m_mesh = mesh;
	int n = m_mesh->nVertex();
	m_list = NULL;

	vector<const float *> point, normal;
	vector<int> id;
	for (int i = 0; i < n; i++)
	{
		if (ridge[i])
		{
			point.push_back(m_mesh->vertex(i)->fv());
			normal.push_back(m_mesh->normal(i)->fv());
			id.push_back(i);
		}
	}

	m_gyralPoint = new bool[n]; memset(m_gyralPoint, 0, sizeof(bool) * n);

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

GyralCurve::~GyralCurve(void)
{
	delete [] m_candEndPoint;
	delete [] m_gyralPoint;
	delete [] m_curveElem;

	curveList *iter = m_list;
	while (iter != NULL)
	{
		curveList *next = iter->next;
		deleteCurves(iter);
		iter = next;
	}
}

void GyralCurve::run(void)
{
	cout << "Grouping.. ";
	fflush(stdout);
	grouping();
	cout << "Done" << endl;
	fflush(stdout);

	cout << "Refining.. ";
	fflush(stdout);
	refineCurves();
	cout << "Done" << endl;
	fflush(stdout);

	cout << "Integrity test.. ";
	fflush(stdout);
	if (testIntegrity()) cout << "Passed" << endl;
	else cout << "Failed" << endl;
	fflush(stdout);

	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		for (int i = 0; i < iter->item.size(); i++)
			m_gyralPoint[iter->item[i]->id] = true;
}

void GyralCurve::grouping(float threshold1, float threshold2, float threshold3)
{
	/*do
	{
		group = current;
		deleteNearestPoints(1.5f);
		current = delineation(treshold);
	}
	while (group != current);*/

	bool cand;
	do
	{
		cand = false;
		deleteNearestPoints(threshold2);
		delineation(threshold1, threshold2, threshold3);
		for (int i = 0; i < m_nPoints && !cand; i++)
		{
			cand = m_candEndPoint[i];
		}
	}
	while (cand);
}

void GyralCurve::deleteNearestPoints(float threshold)
{
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 0; i < iter->item.size(); i++)
		{
			for (int j = 0; j < m_nPoints; j++)
			{
				if (m_curveElem[j]->deleted || m_curveElem[j]->header != NULL) continue;
				if (Vector(iter->item[i]->v, m_curveElem[j]->v).norm() < threshold)
				{
					deleteCurveElem(m_curveElem[j]);
					m_candEndPoint[j] = false;
				}
			}
		}
	}
}

void GyralCurve::deleteNearestPoints(curveList *list, float threshold)
{
	for (int i = 0; i < list->item.size(); i++)
	{
		for (int j = 0; j < m_nPoints; j++)
		{
			if (m_curveElem[j]->deleted || m_curveElem[j]->header != NULL) continue;
			if (Vector(list->item[i]->v, m_curveElem[j]->v).norm() < threshold)
			{
				deleteCurveElem(m_curveElem[j]);
				m_candEndPoint[j] = false;
			}
		}
	}
}

int GyralCurve::delineation(float threshold1, float threshold2, float threshold3)
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
	
	vector<curveElem> endPointList;
	for (int i = 0; i < n; i++)
	{
		if (!m_candEndPoint[i]) continue;
		endPointList.push_back(*m_curveElem[i]);
	}
	sort(endPointList.begin(), endPointList.end());

	curveList *prev;
	for (prev = m_list; prev != NULL && prev->next != NULL; prev = prev->next);

	//for (int i = 0; i < n; i++)
	//for (int i = n - 1; i >= 0; i--)
	for (int c = 0; c < endPointList.size(); c++)
	{
		int i = endPointList[c].id;
		if (m_curveElem[i]->traced || !m_candEndPoint[i] || m_curveElem[i]->deleted || m_curveElem[i]->header != NULL) continue;
		curveList *list = new curveList;
		curveElem *endPoint = curve(m_curveElem[i], list, threshold1);

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

		deleteNearestPoints(list, threshold2);
	}

	return nCurves();
}

void GyralCurve::deleteCurves(curveList *list)
{
	if (list->next != NULL) list->next->prev = list->prev;
	if (list->prev != NULL) list->prev->next = list->next;
	else m_list = list->next;

	list->item.clear();
	delete list;
}

void GyralCurve::deleteCurves(float threshold)
{
	for (curveList *iter = m_list; iter != NULL;)
	{
		curveList *next = iter->next;
		if (threshold > iter->length && (iter->item[0]->isEndPoint || iter->item[iter->item.size() - 1]->isEndPoint))
		{
			/*for (int i = 0; i < iter->item.size(); i++)
				m_candEndPoint[iter->item[i]->id] = true;*/
			if (iter->item[0]->isJunction) iter->item[0]->nJunctions--;
			if (iter->item[iter->item.size() - 1]->isJunction) iter->item[iter->item.size() - 1]->nJunctions--;
			for (int i = 0; i < iter->item.size(); i++)
			{
				if (iter->item[i]->header == iter)
				{
					if (iter->item[i]->isJunction)
						iter->item[i]->header = NULL;
					else deleteCurveElem(iter->item[i]);
				}
			}
			deleteCurves(iter);
		}
		iter = next;
	}

	for (int i = 0; i < m_nPoints; i++)
	{
		if (m_curveElem[i]->isJunction && m_curveElem[i]->header == NULL && m_curveElem[i]->nJunctions > 0)
		{
			findHeader(m_curveElem[i]);
		}
	}

	for (int i = 0; i < m_nPoints; i++)
	{
		if (m_curveElem[i]->nJunctions == 0)
		{
			for (curveList *iter = m_list; iter != NULL; iter = iter->next)
			{
				if (m_curveElem[i] == iter->item[0] || m_curveElem[i] == iter->item[iter->item.size() - 1])
				{
					m_curveElem[i]->isJunction = false;
					m_curveElem[i]->isEndPoint = true;
					m_curveElem[i]->header = iter;
					break;
				}
			}
		}
		else if (m_curveElem[i]->nJunctions == 1)
		{
			curveList *list1 = NULL, *list2 = NULL;
			for (curveList *iter = m_list; iter != NULL; iter = iter->next)
			{
				if (m_curveElem[i] == iter->item[0] || m_curveElem[i] == iter->item[iter->item.size() - 1])
				{
					if (list1 == NULL) list1 = iter;
					else
					{
						list2 = iter;
						break;
					}
				}
			}
			m_curveElem[i]->isJunction = false;
			m_curveElem[i]->isEndPoint = true;
			m_curveElem[i]->nJunctions = 0;
			
			if (list1 == NULL || list2 == NULL) continue;

			if (list1->item.size() > list2->item.size())
			{
				if (list1->item[list1->item.size() - 1] != m_curveElem[i])
					reverseCurveOrder(list1);
				if (list2->item[0] != m_curveElem[i])
					reverseCurveOrder(list2);
				m_curveElem[i]->header = list2;
			}
			else
			{
				if (list1->item[0] != m_curveElem[i])
					reverseCurveOrder(list1);
				if (list2->item[list2->item.size() - 1] != m_curveElem[i])
					reverseCurveOrder(list2);
				m_curveElem[i]->header = list1;
				swap(list1, list2);
			}
			if (list1->item.size() == 2)
			{
				curveElem *elem = NULL;
				if (list1->item[0]->isJunction)
				{
					list1->item[0]->nJunctions--;
					if (list1->item[0]->header == list1)
					{
						list1->item[0]->header = NULL;
						elem = list1->item[0];
					}
					if (list1->item[0]->id < i) i = list1->item[0]->id;
				}
				else deleteCurveElem(list1->item[0]);
				deleteCurves(list1);
				if (elem != NULL) findHeader(elem);
			}
			else
			{
				list1->item.pop_back();
				list1->item[list1->item.size() - 1]->isEndPoint = true;
				joinCurves(list1->item[list1->item.size() - 1], list2->item[0], 0);
			}
		}
	}
}

void GyralCurve::findHeader(curveElem *elem)
{
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		if (elem == iter->item[0] || elem == iter->item[iter->item.size() - 1])
		{
			elem->header = iter;
			break;
		}
	}
}

void GyralCurve::deleteCurveElem(curveElem *elem)
{
	elem->header = NULL;
	elem->deleted = true;
	elem->isEndPoint = false;
	elem->isJunction = false;
}

int GyralCurve::detectEndPoints(float threshold, float inner)
{
	int c = 0;
	int n = m_nPoints;
	vector<int> cand;
	for (int i = 0; i < n; i++)
	{
		m_candEndPoint[i] = false;
		if (m_curveElem[i]->deleted || m_curveElem[i]->header != NULL) continue;
		cand.clear();
		Vector V0(m_curveElem[i]->v);
		for (int j = 0; j < n; j++)
		{
			if (i == j || m_curveElem[j]->deleted || m_curveElem[j]->header != NULL) continue;

			Vector V1(m_curveElem[j]->v);
			float dist = (V1 - V0).norm();
			if (dist < threshold)
			{
				cand.push_back(j);
			}
		}

		m_candEndPoint[i] = true;
		for (int j = 0; j < cand.size() && m_candEndPoint[i]; j++)
		{
			Vector V1(m_curveElem[cand[j]]->v);
			Vector V = (V1 - V0).unit();
			for (int k = j + 1; k < cand.size(); k++)
			{
				Vector V2(m_curveElem[cand[k]]->v);
				float innerProd = V * (V2 - V0).unit();
				if (innerProd < inner) m_candEndPoint[i] = false;
			}
		}
		if (m_candEndPoint[i]) c++;
	}
	return c;
}

void GyralCurve::extendCurves(curveElem *elem, float threshold, float inner1, float inner2)
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
		float dist = (V1 - V).norm();

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
			else if(adjGroup) continue;

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
			if (m_curveElem[closestID]->isEndPoint)
			{
				curveElem *end = elem->header->item[elem->header->item.size() - 1];
				if (end != elem && !end->isJunction)
					extendCurves(end, threshold);
			}
		}
		else
		{
			elem->header->item.push_back(m_curveElem[closestID]);
			elem->isEndPoint = false;
			m_curveElem[closestID]->isEndPoint = true;
			m_curveElem[closestID]->header = elem->header;
			updateOrientation(elem->header, m_curveElem[closestID]);
			m_curveElem[closestID]->deleted = false;

			extendCurves(m_curveElem[closestID], threshold);
		}
	}
}

void GyralCurve::updateOrientation(curveList *list, curveElem *elem)
{
	for (int i = 0; i < list->item.size(); i++)
	{
		if (elem != NULL && elem != list->item[i]) continue;

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
			float v[3];
			Coordinate::sphmean(elem->orientation.fv(), V.unit().fv(), v);
			elem->orientation = v;
		}
	}
}

void GyralCurve::reverseCurveOrder(curveList *list)
{
	reverse(list->item.begin(), list->item.end());
	updateOrientation(list);
}

GyralCurve::curveElem * GyralCurve::curve(curveElem *current, curveList *header, float threshold, float inner1, float inner2)
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
		float dist = (V1 - V).norm();

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
				else wdist = (V1 - V).cross(current->orientation).norm();
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
			current->isEndPoint = true;
			return m_curveElem[closestID];
		}
	}
	
	current->isEndPoint = true;
	return current;
}

void GyralCurve::refineCurves(float threshold)
{
	memset(m_candEndPoint,0,sizeof(bool) * m_nPoints);
	int current, previous;

	// extend endpoints
	do
	{
		for (curveList *iter = m_list; iter != NULL; iter = iter->next)
		{
			curveElem *elem1 = iter->item[iter->item.size() - 1];
			curveElem *elem2 = iter->item[0];
			if (elem1->isEndPoint) extendCurves(elem1, threshold);
			if (elem2->isEndPoint) extendCurves(elem2, threshold);
		}
	}
	while (previous != current);
	
	// connect closest endpoints
	joinCurves(threshold);

	// separate branches
	separateBranch();

	// update length of the curves
	vector<float> hist;
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		iter->length = 0;
		for (int i = 1; i < iter->item.size(); i++)
			iter->length += (Vector(iter->item[i]->v) - Vector(iter->item[i - 1]->v)).norm();
		if (iter->item[0]->isEndPoint || iter->item[iter->item.size() - 1]->isEndPoint)
			hist.push_back(iter->length);
	}
	/*sort(hist.begin(), hist.end());
	float lThreshold = hist[(int)(hist.size() * 0.1)];*/

	// delete a group if its arclength is below a threshold
	current = nCurves();
	do
	{
		previous = current;
		deleteCurves(1.0f);
	}
	while (previous != current);
}

void GyralCurve::joinCurves(float threshold)
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
						if (cand1 != NULL) dist1 = Vector(cand1->v, elem1->v).norm();
						float dist2 = FLT_MAX;
						if (cand2 != NULL) dist2 = Vector(cand2->v, elem2->v).norm();

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

void GyralCurve::joinCurves(curveElem *elem1, curveElem *elem2, float inner)
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
			deleteCurves(list2);
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

GyralCurve::curveElem * GyralCurve::closestPoint(curveElem *elem, curveList *list)
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
	curveElem *elem2;

	for (int i = 0; i < list->item.size(); i++)
	{
		Vector V = Vector(elem->v, list->item[i]->v);
		float dist = V.norm();
		if (minDist > dist)
		{
			minDist = dist;
			elem2 = list->item[i];
		}
	}
	return elem2;
}

void GyralCurve::separateBranch(void)
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

int GyralCurve::nCurves(void)
{
	int n = 0;
	for (curveList *iter = m_list; iter != NULL; iter = iter->next) n++;
	
	return n;
}

void GyralCurve::getSeedPoint(bool *isRidge)
{
	for (int i = 0; i < m_nPoints; i++)
		isRidge[m_curveElem[i]->vid] = true;
}

void GyralCurve::saveGyralPoint(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		fprintf(fp, "%d\n", (int)m_gyralPoint[i]);
	}
	fclose(fp);
}

void GyralCurve::SimplifyCurves(float threshold, float inner)
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
				float len = V1.norm();
				arc += len;
				if (V1.unit() * V2.unit() > inner && arc < threshold)
				{
					deleteCurveElem(iter->item[i]);
					iter->item.erase(iter->item.begin() + i);
					m_gyralPoint[iter->item[i]->id] = false;
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

void GyralCurve::saveGyralCurves(const char *filename, bool incJunc)
{
	FILE *fp = fopen(filename, "w");
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		for (int i = 0; i < iter->item.size(); i++)
		{
			if (iter->item[i]->isJunction && iter->item[i]->header != iter)
				if (!incJunc) continue;	// junction
			fprintf(fp, "%d ", iter->item[i]->vid);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void GyralCurve::saveVTK(const char *filename)
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
		fprintf(fp, "%d ", label[i]);
	fprintf(fp, "\n");
	fclose(fp);
	
	delete [] label;
}

void GyralCurve::showInfo(void)
{
	int nPoints = 0, nCurves = 0, nJunctions = 0;
	memset(m_gyralPoint, 0, sizeof(bool) * m_nPoints);
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
	cout << "# of the gyral points: " << m_nPoints << endl;
	cout << "# of the connected points: " << nPoints << endl;
	cout << "# of the curves: " << nCurves << endl;
	cout << "# of the junction points: " << nJunctions << endl;
}

float GyralCurve::distCurveElem(curveElem *elem1, curveElem *elem2)
{
	if (elem1->header != elem2->header) return FLT_MAX;

	curveList *header = elem1->header;
	float dist = FLT_MAX;
	for (int i = 0; i < header->item.size(); i++)
	{
		if (dist < FLT_MAX)
			dist += Vector(header->item[i - 1]->v, header->item[i]->v).norm();
			if (header->item[i] == elem1 || header->item[i] == elem2)
				break;
		else
			if (header->item[i] == elem1 || header->item[i] == elem2)
				dist = 0;
	}
	return dist;
}

bool GyralCurve::testIntegrity(void)
{
	for (curveList *iter = m_list; iter != NULL; iter = iter->next)
	{
		if ((iter->item[0]->isEndPoint && iter != iter->item[0]->header) // starting point
			|| (iter->item[0]->isJunction && iter->item[0]->header == NULL)
			|| (!iter->item[0]->isEndPoint && !iter->item[0]->isJunction))
		{
			cout <<"Error: Code1 at " << iter->item[0]->id << endl;
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
			cout <<"Error: Code2 at " << iter->item[iter->item.size() - 1]->id << endl;
			cout << "# of points: " << iter->item.size() << endl;
			cout << "EndPoint: " << iter->item[iter->item.size() - 1]->isEndPoint << endl;
			cout << "Junction: " << iter->item[iter->item.size() - 1]->isJunction << endl;
			cout << "NULL Header: " << (iter->item[iter->item.size() - 1]->header == NULL) << endl;
			return false;
		}
		for (int i = 1; i < iter->item.size() - 1; i++)	// intermediate points
			if (iter->item[i]->isEndPoint || iter != iter->item[i]->header)
			{
				cout <<"Error: Code3 at " << iter->item[i]->id << endl;
				return false;
			}
		for (int i = 0; i < iter->item.size(); i++)	// junction points
			if ((iter->item[i]->isEndPoint && iter->item[i]->isJunction)
				|| (iter->item[i]->nJunctions == 0 && iter->item[i]->isJunction)
				|| (iter->item[i]->nJunctions > 0 && !iter->item[i]->isJunction)
				|| (iter->item[i]->header == NULL && iter->item[i]->isJunction))
			{
				cout <<"Error: Code " << iter->item[i]->id;
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
