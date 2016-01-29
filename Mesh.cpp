// --------------------------------------------------------------------------
// Copyright(C) 2009-2015
// Tamy Boubekeur
//
// All rights reserved.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
// for more details.
// --------------------------------------------------------------------------

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

void Mesh::loadOFF (const std::string & filename) {
	ifstream in (filename.c_str ());
    if (!in)
        exit (1);
	string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    V.resize (sizeV);
    T.resize (sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> V[i].p;
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> T[i].v[j];
    }
    in.close ();
    centerAndScaleToUnit ();
    recomputeNormals ();
}

void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < V.size (); i++)
        V[i].n = Vec3f (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < T.size (); i++) {
        Vec3f e01 = V[T[i].v[1]].p -  V[T[i].v[0]].p;
        Vec3f e02 = V[T[i].v[2]].p -  V[T[i].v[0]].p;
        Vec3f n = cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            V[T[i].v[j]].n += n;
    }
    for (unsigned int i = 0; i < V.size (); i++)
        V[i].n.normalize ();
}

void Mesh::smooth () {
		std::vector< vector<Vec3f> > voisinage;
		for (size_t i = 0; i < V.size(); i++) {
			vector<Vec3f> vosin;
			for (size_t k = 0; k < T.size(); k++) {
				Vec3f point  = V[i].p;
				Vec3f e00 = V[T[k].v[0]].p;
				Vec3f e01 = V[T[k].v[1]].p;
				Vec3f e02 = V[T[k].v[2]].p;
				if (point == e00) {
					vosin.push_back(e01);
					vosin.push_back(e02);
					continue;
				}
				if (point == e01) {
					vosin.push_back(e00);
					vosin.push_back(e02);
					continue;
				}
				if (point == e02) {
					vosin.push_back(e00);
					vosin.push_back(e01);
					continue;
				}
				}
				voisinage.push_back(vosin);
			}
			for (size_t i = 0; i < V.size(); i++) {
				Vec3f sum(0.0, 0.0, 0.0);
				for (size_t j = 0;j < voisinage[i].size();j++) {
					sum += voisinage[i][j];
				}
				sum = sum / ( voisinage[i].size());
				V[i].p = 	 sum;
			}
			recomputeNormals();
		}

void Mesh::centerAndScaleToUnit () {
    Vec3f c;
    for  (unsigned int i = 0; i < V.size (); i++)
        c += V[i].p;
    c /= V.size ();
    float maxD = dist (V[0].p, c);
    for (unsigned int i = 0; i < V.size (); i++){
        float m = dist (V[i].p, c);
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < V.size (); i++)
        V[i].p = (V[i].p - c) / maxD;
}

void Mesh::splitEdges (float l ) {
	float coff = l * 4.0 / 3.0;
	std::cout << "The coff is " << coff << std::endl;
	for (size_t i = 0; i < T.size(); i++) {
		/* code */
		int currentTriangle = i;
	 Vertex & v0 = V[T[i].v[0]];
	 Vertex & v1 = V[T[i].v[1]];
	 Vertex & v2 = V[T[i].v[2]];

	Vec3f point0 = Vec3f(v0.p);
	Vec3f point1 = Vec3f(v1.p);
	Vec3f point2 = Vec3f(v2.p);

	Vec3f length0 = point1 - point0;
	Vec3f length1 = point2 - point1;
	Vec3f length2 = point2 - point0;
 if (length0.length() >= coff ) {

	splitEdgesHander (currentTriangle, T[currentTriangle].v[1], T[currentTriangle].v[0],T[currentTriangle].v[2] );
	i += 2;
}

if (length1.length() >= coff) {
	/* code */

		splitEdgesHander ( currentTriangle, T[currentTriangle].v[2], T[currentTriangle].v[1],T[currentTriangle].v[0]);
		i += 2;
}

if (length2.length() >= coff) {
	/* code */
	splitEdgesHander (currentTriangle, T[currentTriangle].v[0], T[currentTriangle].v[2],T[currentTriangle].v[1]);
	T.erase(T.begin() + currentTriangle);
	i += 1;

}


}
return;
}

void Mesh::splitEdgesHander (int numberOfTriangle, int nbPointA, int nbPointB, int nbPointC) {
	Vec3f pointA = Vec3f (V[nbPointA].p);
	Vec3f pointB = Vec3f (V[nbPointB].p);

	Vec3f middlePoint = (pointA + pointB) / 2.0;
	Vec3f middleNormale = Vec3f ((V[nbPointA].n + V[nbPointB].n) / 2.0 );
	Vertex middlePointAdd =  Vertex(middlePoint,middleNormale);
	V.push_back(middlePointAdd);
	//T.erase(T.begin() + numberOfTriangle);

	Triangle triangle0 = Triangle(nbPointA, nbPointC,V.size() -1);
	T.insert(T.begin() + numberOfTriangle + 1 ,triangle0);
	triangle0 = Triangle(V.size() -1,nbPointC, nbPointB); ;
	T.insert(T.begin() + numberOfTriangle + 2, triangle0);
	return;
}
