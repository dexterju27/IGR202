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

#pragma once
#include <cmath>
#include <vector>
#include "Vec3.h"
#include "Edge.h"

/// A simple vertex class storing position() and normal
class Vertex {
public:
    inline Vertex () {}
    inline Vertex (const Vec3f & p, const Vec3f & n) : p (p), n (n) {}
    inline virtual ~Vertex () {}
    Vec3f p;
    Vec3f n;
    std::vector<unsigned long> edge;
};

/// A Triangle class expressed as a triplet of indices (over an external vertex list)
class Triangle {
public:
    inline Triangle () {
        v[0] = v[1] = v[2] = 0;
        willBeDelete = false;
    }
    inline Triangle (const Triangle & t) {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
        willBeDelete = t.willBeDelete;
    }
    inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        willBeDelete = false;
    }
    inline virtual ~Triangle () {}
    inline Triangle & operator= (const Triangle & t) {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
        return (*this);
    }
    unsigned int v[3];
    unsigned int e[3];
    bool willBeDelete;
};

/// A Mesh class, storing a list of vertices and a list of triangles indexed over it.
class Mesh {
public:
	std::vector<Vertex> V;
	std::vector<Triangle> T;
    std::vector<Edge>  E; // vector for edges


    /// Loads the mesh from a <file>.off
	void loadOFF (const std::string & filename);

    /// Compute smooth per-vertex normals
    void recomputeNormals ();

    /// scale to the unit cube and center at original
    void centerAndScaleToUnit ();

    void smooth ();

//    void simplifyMesh (unsigned int r);

    float max (float x, float y);

  void splitEdges (float l );
    
    void createEdgeList ();

//  void splitEdgesHanderOne (std::vector<Edge> edgesWaiting, int numberOfTriangle);
};
