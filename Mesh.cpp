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

//void Mesh::splitEdges (float l ) {
//	float coff = l * 4.0 / 3.0;
//	std::cout << "The coff is " << coff << std::endl;
//	for (size_t i = 0; i < T.size(); i++) {
//		/* code */
//		int currentTriangle = i;
//	 Vertex & v0 = V[T[i].v[0]];
//	 Vertex & v1 = V[T[i].v[1]];
//	 Vertex & v2 = V[T[i].v[2]];
//	Edge e0 = Edge(T[i].v[0], T[i].v[1]);
//	Edge e1 = Edge(T[i].v[1], T[i].v[2]);
//	Edge e2 = Edge(T[i].v[2], T[i].v[0]);
//
//	std::vector<Edge> edgesWaiting;
//
//
//	Vec3f point0 = Vec3f(v0.p);
//	Vec3f point1 = Vec3f(v1.p);
//	Vec3f point2 = Vec3f(v2.p);
//
//	Vec3f length0 = point1 - point0;
//	Vec3f length1 = point2 - point1;
//	Vec3f length2 = point2 - point0;
// if (length0.length() >= coff ) {
//	edgesWaiting.push_back(e0);
//
//	//splitEdgesHander (currentTriangle, T[currentTriangle].v[1], T[currentTriangle].v[0],T[currentTriangle].v[2] );
//	//i += 2;
//}
//
//if (length1.length() >= coff) {
//	/* code */
//	edgesWaiting.push_back(e1);
//	//edgesToDeal++;
//	//	splitEdgesHander ( currentTriangle, T[currentTriangle].v[2], T[currentTriangle].v[1],T[currentTriangle].v[0]);
//		//i += 2;
//}
//
//if (length2.length() >= coff) {
//	/* code */
//	edgesWaiting.push_back(e2);
//	//edgesToDeal++
//	//splitEdgesHander (currentTriangle, T[currentTriangle].v[0], T[currentTriangle].v[2],T[currentTriangle].v[1]);
//	//T.erase(T.begin() + currentTriangle);
//	//i += 1;
//
//}
//
//switch (edgesWaiting.size()) {
//	case 0:
//	break;
//	case 1:
//		splitEdgesHanderOne (edgesWaiting,currentTriangle);
//		break;
//	case 2:
////		splitEdgesHanderTwo (edgesWaiting,currentTriangle);
//		break;
//	case 3:
//		//splitEdgesHanderThree (edgesWaiting,currentTriangle);
//	break;
//
//}
//
//
//}
//return;
//}

void Mesh::createEdgeList() {
    E.clear();
    for (int i = 0; i < V.size(); i++) {
        V[i].edge.clear(); //clear all the edges
    }

		std::cout << "Creating list" << std::endl;
   // unsigned int id = 0;
    for (int i  = 0; i < T.size(); i++) {
        Edge edge0 = Edge(T[i].v[0], T[i].v[1]);
        Edge edge1 = Edge(T[i].v[1], T[i].v[2]);
        Edge edge2 = Edge(T[i].v[2], T[i].v[0]);

        int edgeIndex0 = -1;
        int edgeIndex1 = -1;
        int edgeIndex2 = -1;

        for (int j  = 0; j < E.size(); j++) {
					//check whether this edge is already in the list
            if ( edge0 == E[j]) {
                edgeIndex0 = j;
            }

            if (edge1 == E[j]) {
                edgeIndex1 = j;
            }

            if (edge2 == E[j]) {
                edgeIndex2 = j;
            }
        }

        // if edgeIndex == -1 add new edge
        if (edgeIndex0 == -1) {
            edge0.t.push_back(i);
            V[T[i].v[0]].edge.push_back(E.size());
            V[T[i].v[1]].edge.push_back(E.size());
            T[i].e[0] = E.size();
            E.push_back(edge0);
        }
        else {
            E[edgeIndex0].t.push_back(i);
            V[T[i].v[0]].edge.push_back(edgeIndex0);
            V[T[i].v[1]].edge.push_back(edgeIndex0);
            T[i].e[0] = edgeIndex0;

        }

        if (edgeIndex1 == -1) {
            edge1.t.push_back(i);
            V[T[i].v[1]].edge.push_back(E.size());
            V[T[i].v[2]].edge.push_back(E.size());
            T[i].e[1] = E.size();
            E.push_back(edge1);

        }
        else {

            E[edgeIndex1].t.push_back(i);
            V[T[i].v[1]].edge.push_back(edgeIndex1);
            V[T[i].v[2]].edge.push_back(edgeIndex1);
            T[i].e[1] = edgeIndex1;


        }

        if (edgeIndex2 == -1) {
            edge2.t.push_back(i);
            V[T[i].v[2]].edge.push_back(E.size());
            V[T[i].v[0]].edge.push_back(E.size());
            T[i].e[2] = E.size();
            E.push_back(edge2);
            }

        else {
            E[edgeIndex2].t.push_back(i);
            V[T[i].v[2]].edge.push_back(edgeIndex2);
            V[T[i].v[0]].edge.push_back(edgeIndex2);
            T[i].e[2] = edgeIndex2;


        }

    }

		std::cout << "E.SIZE() is " << E.size() << std::endl;
		return;

}
//void Mesh::splitEdgesHanderOne (std::vector<Edge> edgesWaiting, int numberOfTriangle) {
//	std::vector<Edge>::iterator iter = edgesWaiting.begin();
//	Vec3f pointA = Vec3f (V[iter->v[0]].p);
//	Vec3f pointB = Vec3f (V[iter->v[1]].p);
//
//	Vec3f middlePoint = (pointA + pointB) / 2.0;
//	Vec3f middleNormale = Vec3f ((V[iter->v[0]].n + V[iter->v[0]].n) / 2.0 );
//	Vertex middlePointAdd =  Vertex(middlePoint,middleNormale);
//	V.push_back(middlePointAdd);
//
//	int nbPointC = 0;
//	while (iter->contains(T[numberOfTriangle].v[nbPointC]) == false) {
//		/* code */
//		nbPointC++;
//	}
//
//	T.erase(T.begin() + numberOfTriangle);
//	Triangle triangle0 = Triangle(iter->v[0], nbPointC,V.size() -1);
//	T.insert(T.begin() + numberOfTriangle ,triangle0);
//	triangle0 = Triangle(V.size() -1,nbPointC, iter->v[1]); ;
//	T.insert(T.begin() + numberOfTriangle + 1, triangle0);
//	return;
//}
//
//Mesh::splitEdgesHanderTwo (std::vector<Edge> edgesWaiting, int numberOfTriangle) {
//
//}


void Mesh::splitEdges(float l) {

    float coff = l * 4.0 / 3.0;

    int originalSize = T.size();

		std::cout << "The coff is " << coff << std::endl;
		std::cout << "function splitEdges called " << std::endl;
    for (int i = 0; i < E.size(); i++) {
        if (E[i].traiter == true) {
            continue;
        }
        Vertex pointA = V[E[i].v[0]];
        Vertex pointB = V[E[i].v[1]]; // two points
        Vertex pointC;
        Vertex pointD;
        int nbPointA = E[i].v[0];
        int nbPointB = E[i].v[1];
        int nbPointC = 0;
        int nbPointD = 0;

        Vertex middlePoint;
        Vec3f length = pointA.p - pointB.p;
				//std::cout << length.length() << "THe number of tiangle "<<  E[i].t.size() <<  std::endl;
        if (length.length() >= coff && E[i].t.size() == 2 ) {
            //split do nothing with frontier

            int triangleAdj0 = E[i].t[0];
            int triangleAdj1 = E[i].t[1];
            if (T[triangleAdj0].willBeDelete == true || T[triangleAdj1].willBeDelete == true) {
                    continue;
            }

            //search for another point
            for (int k = 0; k < 3; k++) {
                if ( E[i].contains(T[triangleAdj0].v[k]) == false) {
                    nbPointC = T[triangleAdj0].v[k];
                    pointC = V[T[triangleAdj0].v[k]];
                }
                if ( E[i].contains(T[triangleAdj1].v[k]) == false) {
                    nbPointD = T[triangleAdj1].v[k];
                    pointD = V[T[triangleAdj1].v[k]];
                }

            }

            Triangle newTriangle0 = Triangle(nbPointA,nbPointC,V.size());
            Triangle newTriangle1 = Triangle(nbPointC,nbPointB,V.size());
            Triangle newTriangle2 = Triangle(nbPointB,nbPointD,V.size());
            Triangle newTriangle3 = Triangle(nbPointD,nbPointA,V.size());
            std::cout << "add triangles" << std::endl;

            T[triangleAdj0].willBeDelete = true;
            T[triangleAdj1].willBeDelete = true;
            middlePoint.p = (pointB.p + pointA.p) / 2;
            middlePoint.n = (pointB.n + pointA.n) / 2;

            //add two edges and four triangles



            V.push_back(middlePoint);
            T.push_back(newTriangle0);
            T.push_back(newTriangle1);
            T.push_back(newTriangle2);
            T.push_back(newTriangle3);

        }

    }
    // for (std::vector<Triangle>::iterator iter = T.begin(); iter < T.begin() + originalSize; iter++) {
    //         if (iter->willBeDelete == true) //delete them
    //         {
    //             iter->willBeDelete = false;
    //             T.erase(iter);
		//
    //         }
    // }

    createEdgeList();
    return;

}

void Mesh::flipEdges(){

	std::cout << "function flipEdges called " << std::endl;
	int originalSize = T.size();

	for (int i = 0; i < E.size(); i++) {
			//if (E[i].traiter == true) {
			//		continue;
			//}
			Vertex pointA = V[E[i].v[0]];
			Vertex pointB = V[E[i].v[1]];
			Vertex pointC;
			Vertex pointD;
			int nbPointA = E[i].v[0];
			int nbPointB = E[i].v[1];
			int nbPointC = 0;
			int nbPointD = 0;

			if (E[i].t.size() == 2 ) {
					//On exclue les frontières

					int triangleAdj0 = E[i].t[0];
					int triangleAdj1 = E[i].t[1];

					//Si l'un des deux va être effacé, on passe à l'itération suivante
					if (T[triangleAdj0].willBeDelete == true || T[triangleAdj1].willBeDelete == true) {
									continue;
					}

					//On cherche les points qui appartiennet aux triangles adjacents
					//mais pas à l'edge étudié
					for (int k = 0; k < 3; k++) {
							if ( E[i].contains(T[triangleAdj0].v[k]) == false) {
									nbPointC = T[triangleAdj0].v[k];
									pointC = V[T[triangleAdj0].v[k]];
							}
							if ( E[i].contains(T[triangleAdj1].v[k]) == false) {
									nbPointD = T[triangleAdj1].v[k];
									pointD = V[T[triangleAdj1].v[k]];
							}
					}

					Triangle newTriangle0 = Triangle(nbPointA,nbPointC,nbPointD);
					Triangle newTriangle1 = Triangle(nbPointB,nbPointC,nbPointD);

					std::cout << "add triangles" << std::endl;

					T[triangleAdj0].willBeDelete = true;
					T[triangleAdj1].willBeDelete = true;

					T.push_back(newTriangle0);
					T.push_back(newTriangle1);
					//T.erase(triangleAdj0);
					//T.erase(triangleAdj1);

				}
	}

	for (std::vector<Triangle>::iterator iter = T.begin(); iter < T.begin() + originalSize; iter++) {
	       if (iter->willBeDelete == true) //delete them
	       {
	           iter->willBeDelete = false;
	           T.erase(iter);
         }
	}

	createEdgeList();
	return;

}

float Mesh::computeL() {

	float sum = 0;

	//Faisable avec les edges directement plutôt que les triangles
	for (size_t i = 0; i < T.size(); i++) {
		Vertex & v0 = V[T[i].v[0]];
		Vertex & v1 = V[T[i].v[1]];
		Vertex & v2 = V[T[i].v[2]];

		Vec3f point0 = Vec3f(v0.p);
		Vec3f point1 = Vec3f(v1.p);
		Vec3f point2 = Vec3f(v2.p);

		Vec3f length0 = point1 - point0;
		Vec3f length1 = point2 - point1;
		Vec3f length2 = point2 - point0;

		sum += length0.length() + length1.length() + length2.length();
	}

	float average = sum / (3*T.size());
	std::cout << "Average edges length : " << average << std::endl;

	return average * (0.78); //We take a target l equals to 78% of the average edges length
}

//void Mesh::colapseEdges(float l){
//	float coff = l * 4.0 / 5.0;
//	int originalSize = E.size();
//	std::cout << "The coff is " << coff << std::endl;
//	std::cout << "function colapseEdges called " << std::endl;
//
//	for (int i = 0; i < originalSize; i++) {
//				//if (E[i].traiter == true) {
//					//	continue;
//				//}
//
//				Vertex pointA = V[E[i].v[0]];
//				Vertex pointB = V[E[i].v[1]];
//				int nbPointA = E[i].v[0];
//				int nbPointB = E[i].v[1];
//
//				Vertex middlePoint;
//				Vec3f length = pointA.p - pointB.p;
//
//				if (length.length() < coff ) {
//
//						 middlePoint.p = (pointB.p + pointA.p) / 2;
//						 middlePoint.n = (pointB.n + pointA.n) / 2;
//						// il reste à ajouter les edges du middlePoint qui sont en fait les edges de A et B (pas d'intersection entre les deux hormis E[i])
//						 V.push_back(middlePoint);
//						 int	valenceA = pointA.edge.size();
//						 int  valenceB = pointB.edge.size();
//						 int trianglesAdj;
//						 for (int j = 0; j < valenceA; j++) {
//								trianglesAdj = E[pointA.edge[j]].t.size();
//								for(int k = 0; k < trianglesAdj; k++){
//										// We verify if the triangle has both pointA and pointB as vertices in that case we'll delete without creating new triangle
//										if(T[E[pointA.edge[j]].t[k]].contains(nbPointB)){
//											std::cout << "no triangle added" << std::endl;
//											T[E[pointA.edge[j]].t[k]].willBeDelete = true;
//										}
//										// The case where the triangle doesn't have pointB in his vertices, we create a new triangle with the midpoint and delete the old one
//										// don't forget to make a clean delete and put the right data in the redunduncy
//										else{
//												unsigned int otherPointThanA[2];
//												if(T[E[pointA.edge[j]].t[k]].v[0] = nbPointB){
//													otherPointThanA[0] = T[E[pointA.edge[j]].t[k]].v[1];
//													otherPointThanA[1] = T[E[pointA.edge[j]].t[k]].v[2];
//												}
//
//												if(T[E[pointA.edge[j]].t[k]].v[1] = nbPointB){
//													otherPointThanA[0] = T[E[pointA.edge[j]].t[k]].v[0];
//													otherPointThanA[1] = T[E[pointA.edge[j]].t[k]].v[2];
//												}
//
//												if(T[E[pointA.edge[j]].t[k]].v[2] = nbPointB){
//													otherPointThanA[0] = T[E[pointA.edge[j]].t[k]].v[0];
//													otherPointThanA[1] = T[E[pointA.edge[j]].t[k]].v[1];
//												}
//												int nbPointC = otherPointThanA[0];
//												int nbPointD = otherPointThanA[1];
//												T[E[pointA.edge[j]].t[k]].willBeDelete = true;
//											  Triangle newTriangle = Triangle(nbPointC, nbPointD, V.size());
//												std::cout << "add triangles A" << std::endl;
//												T.push_back(newTriangle);
//										}
//								}
//						 }
//
//						for (int j = 0; j < valenceB; j++) {
//							trianglesAdj = E[pointB.edge[j]].t.size();
//							for(int k = 0; k < trianglesAdj; k++){
//									// We verify if the triangle has both pointA and pointB as vertices in that case we'll delete without creating new triangle
//									if(T[E[pointB.edge[j]].t[k]].contains(nbPointA)){
//										std::cout << "no triangle added" << std::endl;
//										//T[E[pointB.edge[j]].t[k]].willBeDelete = true;
//									}
//									// The case where the triangle doesn't have pointA in his vertices, we create a new triangle with the midpoint and delete the old one
//									// don't forget to make a clean delete and put the right data in the redunduncy
//									else{
//											unsigned int otherPointThanB[2];
//											if(T[E[pointB.edge[j]].t[k]].v[0] = nbPointA){
//												otherPointThanB[0] = T[E[pointB.edge[j]].t[k]].v[1];
//												otherPointThanB[1] = T[E[pointB.edge[j]].t[k]].v[2];
//											}
//
//											if(T[E[pointB.edge[j]].t[k]].v[1] = nbPointA){
//												otherPointThanB[0] = T[E[pointB.edge[j]].t[k]].v[0];
//												otherPointThanB[1] = T[E[pointB.edge[j]].t[k]].v[2];
//											}
//
//											if(T[E[pointB.edge[j]].t[k]].v[2] = nbPointA){
//												otherPointThanB[0] = T[E[pointB.edge[j]].t[k]].v[0];
//												otherPointThanB[1] = T[E[pointB.edge[j]].t[k]].v[1];
//											}
//											int nbPointC = otherPointThanB[0];
//											int nbPointD = otherPointThanB[1];
//											T[E[pointB.edge[j]].t[k]].willBeDelete = true;
//											Triangle newTriangle = Triangle(nbPointC, nbPointD, V.size());
//											std::cout << "add triangles B" << std::endl;
//											T.push_back(newTriangle);
//									}
//							}
//						}
//
//
//					}
//			}
//
//			for (std::vector<Triangle>::iterator iter = T.begin(); iter < T.end(); iter++) {
//	            if (iter->willBeDelete == true) //delete them
//	            {
//	                T.erase(iter);
//	            }
//	    }
//				std::cout << "E.SIZE() IS" << E.size() << std::endl;
//			//createEdgeList();
//	    return;
//}
