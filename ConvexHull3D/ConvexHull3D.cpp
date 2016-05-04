#include <stdio.h>
#include <stdlib.h>
#include <unordered_set>
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Timer.h"
#include "Util/Ply.h"
#include <queue>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <iostream>
#include <array>
#include <unordered_set>
#include <limits>
////////////////////////////
// Command line parsing info

enum {
    ALGORITHM_GIFT_WRAP,
    ALGORITHM_INCREMENTAL,
    ALGORITHM_COUNT
};
const char* AlgorithmNames[] = {"Gift-Wrap", "Incremental"};

cmdLineParameter< char* > Out("out");
cmdLineParameter< int > Count("count"), AlgorithmType("aType", ALGORITHM_GIFT_WRAP), RandomSeed("sRand", 0), Resolution("res", 1024);
cmdLineReadable ForVisualization("viewable");
cmdLineReadable* params[] = {&Count, &Out, &Resolution, &AlgorithmType, &RandomSeed, &ForVisualization, NULL};

void ShowUsage(const char* ex) {
    printf("Usage %s:\n", ex);
    printf("\t --%s <input vertex count>\n", Count.name);
    printf("\t[--%s <output 3D triangulation>]\n", Out.name);
    printf("\t[--%s <algorithm type>=%d]\n", AlgorithmType.name, AlgorithmType.value);
    for (int i = 0; i < ALGORITHM_COUNT; i++) printf("\t\t%d] %s\n", i, AlgorithmNames[i]);
    printf("\t[--%s <random seed>=%d]\n", RandomSeed.name, RandomSeed.value);
    printf("\t[--%s <grid resolution>=%d]\n", Resolution.name, Resolution.value);
    printf("\t[--%s]\n", ForVisualization.name);
}
// Command line parsing info
////////////////////////////

using namespace Geometry;

void RandomPoints(std::vector< Point3i >& points, int count, int seed, int res) {
    const int MAX_EXP = (sizeof (long long) * 8) / 3;
    if (res > (1 << MAX_EXP)) fprintf(stderr, "[WARNING] Maximum resolution is: %d\n", 1 << MAX_EXP), res = 1 << MAX_EXP;

    srand(seed);

    // Add distinct random points in a disk
    points.reserve(count);
    Point3i center;
    center[0] = center[1] = center[2] = res / 2;
    long long r = res / 2;
    std::unordered_set< long long > usedPoints;
    while (points.size() < count) {
        Point3i p;
        p[0] = rand() % res, p[1] = rand() % res, p[2] = rand() % res;
        {
            long long d[] = {center[0] - p[0], center[1] - p[1], center[2] - p[2]};
            if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > r * r) continue;
        }
        long long key = (((long long) p[0]) << (2 * MAX_EXP)) | (((long long) p[1]) << MAX_EXP) | ((long long) p[2]);
        if (usedPoints.find(key) == usedPoints.end()) points.push_back(p), usedPoints.insert(key);
    }
}

long long SquaredArea(Point<3, int> &p0, Point<3, int> &p1, Point<3, int> &p2) {
    long long a = 0;
    a += ((long long) (p1[0] + p0[0])) * (p1[1] - p0[1]);
    a += ((long long) (p2[0] + p1[0])) * (p2[1] - p1[1]);
    a += ((long long) (p0[0] + p2[0])) * (p0[1] - p2[1]);
    return a*a;
}

long SignedVolume(Point<3, int> &p0, Point<3, int> &p1, Point<3, int> &p2, Point<3, int> &p3) {
    long long vol;
    long ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;

    bx = p1[0] - p0[0];
    cx = p2[0] - p0[0];
    dx = p3[0] - p0[0];

    by = p1[1] - p0[1];
    cy = p2[1] - p0[1];
    dy = p3[1] - p0[1];

    bz = p1[2] - p0[2];
    cz = p2[2] - p0[2];
    dz = p3[2] - p0[2];

    vol = dx * (by * cz - bz * cy) + dy * (bz * cx - bx * cz) + dz * (bx * cy - by * cx);

    return vol;
}

// this is the incremental algorithm implementation using custom design half-edge data structure from the textbook
// this implementation is closely following the textbook in the cleanup procedure
// this implementation is for grading.
namespace IncrementalAlg{
    typedef struct vertexStruct *vertex_ptr;
    typedef struct halfEdgeStruct *halfEdge_ptr;
    typedef struct faceStruct *face_ptr;
    
    struct halfEdgeStruct{
        vertex_ptr vert;
        halfEdge_ptr pair;
        face_ptr face;
        halfEdge_ptr next,prev;
        
        bool deleteFlag;
        bool mark; /* T iff point already processed. */        
        face_ptr newFace;
    };
    
    struct vertexStruct{
        Point<3, int> v;
        halfEdge_ptr edge;
        vertex_ptr next, prev;
        
        int idx;
        bool onhull;
        bool mark;
        halfEdge_ptr newHalfEdge;
    };
    
    struct faceStruct{
        halfEdge_ptr edge[3];
        face_ptr next;
        face_ptr prev;
        vertex_ptr vertex[3];
        bool visible;
    };
    template<typename T>
    void SWAP(T& t, T& x, T& y) {
        t = x;x = y;y = t;
    }

    template<typename T>
    void fREE(T& p) {
        if (p) { free(p); p = NULL; }
    }

    template<typename T>
    void ADD(T& head, T& p) {
        if (head){
        p->next = head;
        p->prev = head->prev;
        head->prev = p;
        p->prev->next = p;} else{
            head = p;
            head->next = head->prev = p;
        }
    }

    template<typename T>
    void DELETE(T& head, T& p) {
        if (head) {
            if (head == head->next)
                head = NULL;
            else if (p == head)
                head = head->next;
            p->next->prev = p->prev;
            p->prev->next = p->next;
            FREE(p);
        }
    }
    bool Collinear(vertex_ptr a, vertex_ptr b, vertex_ptr c) {
        return (SquaredArea(a->v, b->v, c->v) == 0);
    }
    vertex_ptr vertices = NULL;
    halfEdge_ptr edges = NULL;
    face_ptr faces = NULL;
    
    
    vertex_ptr MakeNewVertex() {
        vertex_ptr v = new vertexStruct;
        v->newHalfEdge = NULL;
        v->onhull = false;
        v->mark = false;
        ADD(vertices, v);
        return v;
    }

    halfEdge_ptr MakeNewEdge() {
        halfEdge_ptr e = new halfEdgeStruct;
        e->face = e->newFace = NULL;
        e->vert = NULL;
        e->deleteFlag = false;
        ADD(edges, e);
        return e;
    }

    face_ptr MakeNewFace() {
        face_ptr f = new faceStruct;
        for (int i = 0; i < 3; ++i) {
            f->edge[i] = NULL;
            f->vertex[i] = NULL;
        }
        f->visible = false;
        ADD(faces, f);
        return f;
    }    
    
    
    face_ptr MakeFace(vertex_ptr v0, vertex_ptr v1, vertex_ptr v2) {
        face_ptr f;
        halfEdge_ptr e0, e1, e2;

        
        e0 = MakeNewEdge();e1 = MakeNewEdge();e2 = MakeNewEdge();
        
        e0->vert = v0;e1->vert = v1;e2->vert = v2;

        /* Create face for triangle. */
        f = MakeNewFace();
        f->edge[0] = e0;
        f->edge[1] = e1;
        f->edge[2] = e2;
        f->vertex[0] = v0;
        f->vertex[1] = v1; 
        f->vertex[2] = v2;
        /* Link edges to face. */
        e0->face = e1->face = e2->face = f;
        return f;
    }

    
    void BuildInitialHull() {
        vertex_ptr v0, v1, v2, v3, t;
        face_ptr f0, f1 = NULL;
        halfEdge_ptr e0, e1, e2, s;
        
        /* Find 3 noncollinear points. */
        v0 = vertices;
        while (Collinear(v0, v0->next, v0->next->next)){
            if ((v0 = v0->next) == vertices)
                printf("All points are Collinear!\n"), exit(0);
        }
        v1 = v0->next;
        v2 = v1->next;

        /* Mark the vertices as processed. */
        v0->mark = true;
        v1->mark = true;
        v2->mark = true;

        /* Create the two "twin" faces. */
        f0 = MakeFace(v0, v1, v2);
        f1 = MakeFace(v2, v1, v0);

        
        // link half edge pairs
        
        f0->edge[0]->pair = f1->edge[1];
        f0->edge[1]->pair = f1->edge[0];
        f0->edge[2]->pair = f1->edge[2];
        
        f1->edge[0]->pair = f0->edge[1];
        f1->edge[1]->pair = f0->edge[0];
        f1->edge[2]->pair = f0->edge[2];
        
        
        
        /* find a fourth, noncoplanar point to form tetrahedron. */
        v3 = v2->next;
        auto vol = SignedVolume(f0->vertex[0]->v, f0->vertex[1]->v, f0->vertex[2]->v, v3->v);
        while (!vol) {
            if ((v3 = v3->next) == v0)
                printf("DoubleTriangle:  All points are coplanar!\n"), exit(0);
            vol = SignedVolume(f0->vertex[0]->v, f0->vertex[1]->v, f0->vertex[2]->v, v3->v);
        }
        vertices = v3;
    }
    
    face_ptr MakeConeFace(halfEdge_ptr e, vertex_ptr p) {
        halfEdge_ptr new_edge[2];
        face_ptr new_face;
        /* Make two new edges (if don't already exist). */

        new_edge[0] = MakeNewEdge();
        new_edge[0]->vert = e->pair->vert;

        new_edge[1] = MakeNewEdge();
        new_edge[1]->vert = p;

        // here need to add pair information, and new half-edge info to vertex
        if (e->vert->newHalfEdge == NULL) {
            e->vert->newHalfEdge = new_edge[1];
        } else {
            new_edge[1]->pair = e->vert->newHalfEdge;
            e->vert->newHalfEdge->pair = new_edge[1];
        }

        if (e->pair->vert->newHalfEdge == NULL) {
            e->pair->vert->newHalfEdge = new_edge[0];
        } else {
            new_edge[0]->pair = e->pair->vert->newHalfEdge;
            e->pair->vert->newHalfEdge->pair = new_edge[0];
        }
        
        /* Make the new face. */
        new_face = MakeNewFace();
        new_face->edge[0] = e;new_face->edge[1] = new_edge[0];new_face->edge[2] = new_edge[1];
        new_face->vertex[0] = e->vert;new_face->vertex[1] = new_edge[0]->vert;new_face->vertex[2] = p;

        /* Set the newfaces for new edges*/       
        new_edge[0]->face = new_face;
        new_edge[1]->face = new_face;
        
        return new_face;
    }
    
    
        bool AddPointToHull(vertex_ptr p) {
        face_ptr f;
        halfEdge_ptr e, temp;
        bool vis = false;

        /* Mark faces visible from p. */
        f = faces;
        do {
            auto vol = SignedVolume(f->vertex[0]->v, f->vertex[1]->v, f->vertex[2]->v, p->v);
            if (vol < 0) {
                f->visible = true;
                vis = true;
            }
            f = f->next;
        } while (f != faces);

        /* If no faces are visible from p, then p is inside the hull. */
        if (!vis) {
            p->onhull = false;
            return false;
        }

        // iterate through all half edges
        e = edges;
        do {
            temp = e->next;
            if (e->face->visible && e->pair->face->visible) {// if both visible mark for deletion
                                
                e->deleteFlag = true;
            }else if (e->face->visible && !(e->pair->face->visible)){
                e->newFace = MakeConeFace(e, p);           
            }
            e = temp;
        } while (e != edges);
        return true;
    }
        
        void CleanEdges(void) {
        halfEdge_ptr e; /* Primary index into edge list. */
        halfEdge_ptr t; /* Temporary edge pointer. */

        e = edges;
        do {
            if (e->newFace) {
                if (e->face->visible)
                    e->face = e->newFace;
                e->newFace = NULL;
            }
            e = e->next;
        } while (e != edges);

        /* Delete any edges where its face and its pair's face are visible marked for deletion. */
        while (edges && edges->deleteFlag) {
            e = edges;
            DELETE(edges, e);
        }
        
        e = edges->next;
        do {
            if (e->deleteFlag) {
                t = e;
                e = e->next;
                DELETE(edges, t);
            } else e = e->next;
        } while (e != edges);
    }

    void CleanFaces(void) {
        face_ptr f; /* Primary pointer into face list. */
        face_ptr t; 
        while (faces && faces->visible) {
            f = faces;
            DELETE(faces, f);
        }
        f = faces->next;
        do {
            if (f->visible) {
                t = f;
                f = f->next;
                DELETE(faces, t);
            } else f = f->next;
        } while (f != faces);
    }

    void CleanVertices(vertex_ptr *pvnext) {
        halfEdge_ptr e;
        vertex_ptr v, t;

        /* Mark all vertices incident to some undeleted edge as on the hull. */
        e = edges;
        do {
//            e->vert->onhull = e->pair->vert->onhull = true;
            e->vert->onhull = true;
            e = e->next;
        } while (e != edges);

        /* Delete all vertices that have been processed but
           are not on the hull. */
        while (vertices && vertices->mark && !vertices->onhull) {
            /* If about to delete vnext, advance it first. */
            v = vertices;
            if (v == *pvnext)
                *pvnext = v->next;
            DELETE(vertices, v);
        }
        v = vertices->next;
        do {
            if (v->mark && !v->onhull) {
                t = v;
                v = v->next;
                if (t == *pvnext)
                    *pvnext = t->next;
                DELETE(vertices, t);
            } else v = v->next;
        } while (v != vertices);

        /* Reset flags. */
        v = vertices;
        do {
            v->newHalfEdge = NULL;
            v->onhull = false;
            v = v->next;
        } while (v != vertices);
    }

    void CleanUp(vertex_ptr *pvnext) {
        CleanEdges();
        CleanFaces();
        CleanVertices(pvnext);
    }

    void ConstructHull() {
        vertex_ptr v, vnext;
        bool changed; /* T if addition changes hull; not used. */

        v = vertices;
        do {
            vnext = v->next;
            if (!v->mark) {
                v->mark = true;
                changed = AddPointToHull(v);
 //               PrintOut(v);
                CleanUp(&vnext); /* Pass down vnext in case it gets deleted. */                
 //               PrintOut(v);
            }
            v = vnext;
        }while (v != vertices);
    }        
}


// this is the incremental algorithm implementation using full edge data structure from the textbook
// this implementation is closely following the textbook
// this implementation is simply for comparison, not for grading.
namespace IncrementalAlg2 {
    typedef struct vertexStruct *vertex_ptr;
    typedef struct edgeStruct *edge_ptr;
    typedef struct faceStruct *face_ptr;
    void PrintOut( vertex_ptr v );
    struct vertexStruct {
        Point<3, int> v;
        int idx;
        edge_ptr duplicate; /* pointer to incident cone edge (or NULL) */
        bool onhull; /* T iff point on hull. */
        bool mark; /* T iff point already processed. */
        vertex_ptr next, prev;
    };

    struct edgeStruct {
        face_ptr adjface[2];
        vertex_ptr endpts[2];
        face_ptr newface; /* pointer to incident cone face. */
        bool deleteFlag; /* T iff edge should be delete. */
        edge_ptr next, prev;
    };

    struct faceStruct {
        edge_ptr edge[3];
        vertex_ptr vertex[3];
        bool visible; /* T iff face visible from new point. */
        face_ptr next, prev;
    };

    template<typename T>
    void SWAP(T& t, T& x, T& y) {
        t = x;x = y;y = t;
    }

    template<typename T>
    void fREE(T& p) {
        if (p) { free(p); p = NULL; }
    }

    template<typename T>
    void ADD(T& head, T& p) {
        if (head){
        p->next = head;
        p->prev = head->prev;
        head->prev = p;
        p->prev->next = p;} else{
            head = p;
            head->next = head->prev = p;
        }
    }

    template<typename T>
    void DELETE(T& head, T& p) {
        if (head) {
            if (head == head->next)
                head = NULL;
            else if (p == head)
                head = head->next;
            p->next->prev = p->prev;
            p->prev->next = p->next;
            FREE(p);
        }
    }
    bool Collinear(vertex_ptr a, vertex_ptr b, vertex_ptr c) {
        return (SquaredArea(a->v, b->v, c->v) == 0);
    }
    vertex_ptr vertices = NULL;
    edge_ptr edges = NULL;
    face_ptr faces = NULL;

    vertex_ptr MakeNewVertex() {
        vertex_ptr v = new vertexStruct;
        v->duplicate = NULL;
        v->onhull = false;
        v->mark = false;
        ADD(vertices, v);
        return v;
    }

    edge_ptr MakeNewEdge() {
        edge_ptr e = new edgeStruct;
        e->adjface[0] = e->adjface[1] = e->newface = NULL;
        e->endpts[0] = e->endpts[1] = NULL;
        e->deleteFlag = false;
        ADD(edges, e);
        return e;
    }

    face_ptr MakeNewFace() {
        face_ptr f = new faceStruct;
        for (int i = 0; i < 3; ++i) {
            f->edge[i] = NULL;
            f->vertex[i] = NULL;
        }
        f->visible = false;
        ADD(faces, f);
        return f;
    }

    face_ptr MakeFace(vertex_ptr v0, vertex_ptr v1, vertex_ptr v2, face_ptr fold) {
        face_ptr f;
        edge_ptr e0, e1, e2;

        /* Create edges of the initial triangle. */
        if (!fold) {
            e0 = MakeNewEdge();
            e1 = MakeNewEdge();
            e2 = MakeNewEdge();
        } else { /* Copy from fold, in reverse order. */
            e0 = fold->edge[2];
            e1 = fold->edge[1];
            e2 = fold->edge[0];
        }
        e0->endpts[0] = v0;
        e0->endpts[1] = v1;
        e1->endpts[0] = v1;
        e1->endpts[1] = v2;
        e2->endpts[0] = v2;
        e2->endpts[1] = v0;

        /* Create face for triangle. */
        f = MakeNewFace();
        f->edge[0] = e0;
        f->edge[1] = e1;
        f->edge[2] = e2;
        f->vertex[0] = v0;
        f->vertex[1] = v1; 
        f->vertex[2] = v2;

        /* Link edges to face. */
        e0->adjface[0] = e1->adjface[0] = e2->adjface[0] = f;

        return f;

    }

    void MakeCcw(face_ptr f, edge_ptr e, vertex_ptr p) {
        face_ptr fv; /* The visible face adjacent to e */
        int i; /* Index of e->endpoint[0] in fv. */
        edge_ptr s; /* Temporary, for swapping */

        if (e->adjface[0]->visible)
            fv = e->adjface[0];
        else fv = e->adjface[1];

        /* Set vertex[0] & [1] of f to have the same orientation
           as do the corresponding vertices of fv. */
        for (i = 0; fv->vertex[i] != e->endpts[0]; ++i)
            ;
        /* Orient f the same as fv. */
        if (fv->vertex[ (i + 1) % 3 ] != e->endpts[1]) {
            f->vertex[0] = e->endpts[1];
            f->vertex[1] = e->endpts[0];
        } else {
            f->vertex[0] = e->endpts[0];
            f->vertex[1] = e->endpts[1];
            SWAP(s, f->edge[1], f->edge[2]);
        }
        /* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
           edge[2] on endpt[1].  So if e is oriented "forwards," we
           need to move edge[1] to follow [0], because it precedes. */

        f->vertex[2] = p;
    }

    void BuildInitialHull() {
        vertex_ptr v0, v1, v2, v3, t;
        face_ptr f0, f1 = NULL;
        edge_ptr e0, e1, e2, s;
        
        /* Find 3 noncollinear points. */
        v0 = vertices;
        while (Collinear(v0, v0->next, v0->next->next)){
            if ((v0 = v0->next) == vertices)
                printf("DoubleTriangle:  All points are Collinear!\n"), exit(0);
        }
        v1 = v0->next;
        v2 = v1->next;

        /* Mark the vertices as processed. */
        v0->mark = true;
        v1->mark = true;
        v2->mark = true;

        /* Create the two "twin" faces. */
        f0 = MakeFace(v0, v1, v2, f1);
        f1 = MakeFace(v2, v1, v0, f0);

        /* Link adjacent face fields. */
        f0->edge[0]->adjface[1] = f1;
        f0->edge[1]->adjface[1] = f1;
        f0->edge[2]->adjface[1] = f1;
        f1->edge[0]->adjface[1] = f0;
        f1->edge[1]->adjface[1] = f0;
        f1->edge[2]->adjface[1] = f0;

        /* Find a fourth, noncoplanar point to form tetrahedron. */
        v3 = v2->next;
        auto vol = SignedVolume(f0->vertex[0]->v, f0->vertex[1]->v, f0->vertex[2]->v, v3->v);
        while (!vol) {
            if ((v3 = v3->next) == v0)
                printf("DoubleTriangle:  All points are coplanar!\n"), exit(0);
            vol = SignedVolume(f0->vertex[0]->v, f0->vertex[1]->v, f0->vertex[2]->v, v3->v);
        }
        /* Insure that v3 will be the first added. */
        vertices = v3;
//        PrintOut(v3);
    }

    face_ptr MakeConeFace(edge_ptr e, vertex_ptr p) {
        edge_ptr new_edge[2];
        face_ptr new_face;
        /* Make two new edges (if don't already exist). */
        for (int i = 0; i < 2; ++i)
            /* If the edge exists, copy it into new_edge. */
            if (!(new_edge[i] = e->endpts[i]->duplicate)) {
                /* Otherwise (duplicate is NULL), MakeNullEdge. */
                new_edge[i] = MakeNewEdge();
                new_edge[i]->endpts[0] = e->endpts[i];
                new_edge[i]->endpts[1] = p;
                e->endpts[i]->duplicate = new_edge[i];
            }

        /* Make the new face. */
        new_face = MakeNewFace();
        new_face->edge[0] = e;new_face->edge[1] = new_edge[0];new_face->edge[2] = new_edge[1];
 //       PrintOut(vertices);
        MakeCcw(new_face, e, p);
 //       PrintOut(vertices);
        /* Set the adjacent face pointers. */
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                /* Only one NULL link should be set to new_face. */
                if (!(new_edge[i]->adjface[j])) {
                    new_edge[i]->adjface[j] = new_face;
                    break;
                }
            }
        }
        return new_face;
    }

    bool AddPointToHull(vertex_ptr p) {
        face_ptr f;
        edge_ptr e, temp;
       bool vis = false;

        /* Mark faces visible from p. */
        f = faces;
        do {
            auto vol = SignedVolume(f->vertex[0]->v, f->vertex[1]->v, f->vertex[2]->v, p->v);
            if (vol > 0) {
                f->visible = true;
                vis = true;
            }
            f = f->next;
        } while (f != faces);

        /* If no faces are visible from p, then p is inside the hull. */
        if (!vis) {
            p->onhull = false;
            return false;
        }

        e = edges;
        do {
            temp = e->next;
            if (e->adjface[0]->visible && e->adjface[1]->visible) {// if both visible
                /* mark for deletion. */
                e->deleteFlag = true;
            }
            else if (e->adjface[0]->visible || e->adjface[1]->visible){ // else if only one is visible
                /* e border: make a new face. */
                e->newface = MakeConeFace(e, p);
 //               PrintOut(vertices);
            }
            e = temp;
        } while (e != edges);
        return true;

    }

    void CleanEdges(void) {
        edge_ptr e; /* Primary index into edge list. */
        edge_ptr t; /* Temporary edge pointer. */

        /* Integrate the newface's into the data structure. */
        /* Check every edge. */
        e = edges;
        do {
            if (e->newface) {
                if (e->adjface[0]->visible)
                    e->adjface[0] = e->newface;
                else e->adjface[1] = e->newface;
                e->newface = NULL;
            }
            e = e->next;
        } while (e != edges);

        /* Delete any edges marked for deletion. */
        while (edges && edges->deleteFlag) {
            e = edges;
            DELETE(edges, e);
        }
        e = edges->next;
        do {
            if (e->deleteFlag) {
                t = e;
                e = e->next;
                DELETE(edges, t);
            } else e = e->next;
        } while (e != edges);
    }

    void CleanFaces(void) {
        face_ptr f; /* Primary pointer into face list. */
        face_ptr t; /* Temporary pointer, for deleting. */

        while (faces && faces->visible) {
            f = faces;
            DELETE(faces, f);
        }
        f = faces->next;
        do {
            if (f->visible) {
                t = f;
                f = f->next;
                DELETE(faces, t);
            } else f = f->next;
        } while (f != faces);
    }

    void CleanVertices(vertex_ptr *pvnext) {
        edge_ptr e;
        vertex_ptr v, t;

        /* Mark all vertices incident to some undeleted edge as on the hull. */
        e = edges;
        do {
            e->endpts[0]->onhull = e->endpts[1]->onhull = true;
            e = e->next;
        } while (e != edges);

        /* Delete all vertices that have been processed but
           are not on the hull. */
        while (vertices && vertices->mark && !vertices->onhull) {
            /* If about to delete vnext, advance it first. */
            v = vertices;
            if (v == *pvnext)
                *pvnext = v->next;
            DELETE(vertices, v);
        }
        v = vertices->next;
        do {
            if (v->mark && !v->onhull) {
                t = v;
                v = v->next;
                if (t == *pvnext)
                    *pvnext = t->next;
                DELETE(vertices, t);
            } else v = v->next;
        } while (v != vertices);

        /* Reset flags. */
        v = vertices;
        do {
            v->duplicate = NULL;
            v->onhull = false;
            v = v->next;
        } while (v != vertices);
    }

    void CleanUp(vertex_ptr *pvnext) {
        CleanEdges();
        CleanFaces();
        CleanVertices(pvnext);
    }

    void ConstructHull() {
        vertex_ptr v, vnext;
        bool changed; /* T if addition changes hull; not used. */

        v = vertices;
        do {
            vnext = v->next;
            if (!v->mark) {
                v->mark = true;
                changed = AddPointToHull(v);
 //               PrintOut(v);
                CleanUp(&vnext); /* Pass down vnext in case it gets deleted. */                
 //               PrintOut(v);
            }
            v = vnext;
        }while (v != vertices);

    }

}
namespace GiftWrap {

    typedef std::array<Point<3, int>, 2> EdgePoint;

    int PivotOnEdge(EdgePoint &edge, std::vector< Point < 3, int>>&points) {
        int p = 0;
        int np = points.size();
//        std::vector<double> volSet;
        long long area2 = SquaredArea(edge[0], edge[1], points[p]);
        for (int i = 1; i < np; i++) {
            long long volume6 = SignedVolume(edge[0], edge[1], points[p], points[i]);
//            volSet.push_back(volume6);

            if (volume6 < 0) {
                p = i; // this step is comparing p and i to  see whether i is on the outer side of the face e0,e1,p, if negative, then it is
                //this step can always ensure this is the first encountered point during the pivoting process
            } else if (volume6 == 0) {
                long long _area2 = SquaredArea(edge[0], edge[1], points[i]);
                if (_area2 > area2) {
                    area2 = _area2; // this step can always eliminate points inside a triangle
                    p = i;
                }

            }

        }

        return p;
    }

    Edge FindEdgeOnHull(std::vector< Point < 3, int >>&points) {
        // first find the bottom most left most back most point
        int p = 0;
        for (int i = 1; i < points.size(); i++) {
            if (points[i][0] < points[p][0]
                    || (points[i][0] == points[p][0] && points[i][1] < points[p][1])
                    || (points[i][0] == points[p][0] && points[i][1] == points[p][1] && points[i][2] < points[p][2])) {
                p = i;
            }
        }
        int q = p;
        for (int i = 0; i < points.size(); i++) {
            if ((points[i][0] == points[q][0] && points[i][1] == points[q][1] && points[i][2] > points[q][2])) {
                q = i;
            }
        }

        EdgePoint e;
        e[0] = points[p];
        e[1] = points[q];
        if (q == p) {
            e[1][2] = e[1][2] + 1;
        }

        q = PivotOnEdge(e, points);
        Edge res;
        res[0] = p;
        res[1] = q;

        return res;
    }

    Triangle FindTriangleOnHull(std::vector< Point < 3, int >>&points) {

        Edge edge = FindEdgeOnHull(points);
        EdgePoint edgeP = {points[edge[0]], points[edge[1]]};
        int r = PivotOnEdge(edgeP, points);
        Triangle tri;
        tri[0] = edge[0];
        tri[1] = edge[1];
        tri[2] = r;
        return tri;
    }

    struct Edgehash {

        std::size_t operator()(const Edge &x) const {
            return std::hash<int>()(x[0])+std::hash<int>()(x[1]);
        }
    };

    struct EdgeEqual {

        bool operator()(const Edge &x1, const Edge &x2) const {
            return ((x1[0]==x2[0] && x1[1] ==x2[1])||(x1[0]==x2[1] && x1[1] ==x2[0]));
        }
    };

    struct Trihash {

        std::size_t operator()(const Triangle &x) const {
            return std::hash<int>()(x[0])+std::hash<int>()(x[1])+std::hash<int>()(x[2]);
        }
    };

    struct TriEqual {

        bool operator()(const Triangle &x1, const Triangle &x2) const {
            bool case1 = (x1[0]==x2[0] && x1[1] ==x2[1] && x1[2]==x2[2]);
            bool case2 = (x1[0]==x2[1] && x1[1] ==x2[2] && x1[2]==x2[0]);
            bool case3 = (x1[0]==x2[2] && x1[1] ==x2[0] && x1[2]==x2[1]);
            return (case1 || case2 || case3);
       }
    };

    
    struct EdgeInfo {
        std::unordered_map<Edge, int, Edgehash, EdgeEqual> mark;

        bool NotProcessed(Edge e) {
            if (mark.find(e) == mark.end()) {
                return true; // if not found, then not processed
            }
            return false;
        }

        void MarkProcessedEdges(Edge e) {
            mark[e] = 1;
        }
    };


}

template< class CType >
void GiftWrapAlgorithm(std::vector< Point< 3, CType > >& points, std::vector< Point< 3, CType > >& hullV, std::vector< Triangle >& hullF) {
   
    Triangle t = GiftWrap::FindTriangleOnHull(points);
    std::queue<Edge> Q;
    Edge e0, e1, e2;
    e0[0] = t[1];
    e0[1] = t[0];
    e1[0] = t[2];
    e1[1] = t[1];
    e2[0] = t[0];
    e2[1] = t[2];
    Q.push(e0);
    Q.push(e1);
    Q.push(e2);
    hullF.push_back(t);
    GiftWrap::EdgeInfo edgeInfo;
    std::unordered_set<Triangle,GiftWrap::Trihash,GiftWrap::TriEqual> tri_set;
    tri_set.insert(t);
    while (!Q.empty()) {
        Edge e = Q.front();
        Q.pop();
        if (edgeInfo.NotProcessed(e)) {
            GiftWrap::EdgePoint edgeP = {points[e[0]], points[e[1]]};
            int q = GiftWrap::PivotOnEdge(edgeP, points);
            t[0] = e[0];
            t[1] = e[1];
            t[2] = q;
            // if it is new triangle
            if (tri_set.find(t) == tri_set.end()){
                tri_set.insert(t);
                hullF.push_back(t);
            // hullV.push_back(points[t[0]]);hullV.push_back(points[t[1]]);hullV.push_back(points[t[2]]);
                e0[0] = t[1];
                e0[1] = t[0];
                e1[0] = t[2];
                e1[1] = t[1];
                e2[0] = t[0];
                e2[1] = t[2];
                Q.push(e0);
                Q.push(e1);
                Q.push(e2);
            }
            edgeInfo.MarkProcessedEdges(e);
        }
    }
    std::unordered_set<int> V;
    for (auto &t : hullF) {
        V.insert(t[0]);
        V.insert(t[1]);
        V.insert(t[2]);
    }
    int count = 0;
    std::unordered_map<int,int> v_map;
    for (auto &v : V) {
        hullV.push_back(points[v]);
        v_map[v] = count++;
//        std::cout << v << std::endl;
    }
    for (auto &t : hullF) {
        t[0] = v_map[t[0]];
        t[1] = v_map[t[1]];
        t[2] = v_map[t[2]];
    }
    
    std::ofstream os;
    os.open("faces.dat");
    
    for (auto &t: hullF){
        os << t[0] << "\t" << t[1] << "\t" << t[2] << std::endl;    
    }
    os.close();
}

template< class CType >
void IncrementalAlgorithm(std::vector< Point< 3, CType > >& points, std::vector< Point< 3, CType > >& hullV, std::vector< Triangle >& hullF) {

    int count = 0;
    for(auto &p: points){
        auto v = IncrementalAlg::MakeNewVertex();
        v->v = p;
        v->idx = count++;
    }
    
    IncrementalAlg::BuildInitialHull();
    IncrementalAlg::ConstructHull();
    
    IncrementalAlg::vertex_ptr v = IncrementalAlg::vertices;
    count = 0;
    std::unordered_map<int,int> v_map;
    do{
        hullV.push_back(v->v);
        v_map[v->idx] = count++;
        v = v->next;
    }while(v!=IncrementalAlg::vertices);
    
    IncrementalAlg::face_ptr f = IncrementalAlg::faces;
    do{
        Triangle tri;
        tri[0] = v_map[f->vertex[0]->idx];tri[1] = v_map[f->vertex[1]->idx];tri[2] = v_map[f->vertex[2]->idx];
        hullF.push_back(tri);
        f = f->next;
    }while(f!=IncrementalAlg::faces);
    
          std::ofstream os;
    os.open("faces.dat");
    
    for (auto &t: hullF){
        os << t[0] << "\t" << t[1] << "\t" << t[2] << std::endl;    
    }
    os.close(); 
}

int main(int argc, char* argv[]) {
    PlyProperty Point3iProperties[] ={
        { "x", PLY_INT, PLY_INT, int( offsetof(Point3i, coordinates[0])), 0, 0, 0, 0},
        { "y", PLY_INT, PLY_INT, int( offsetof(Point3i, coordinates[1])), 0, 0, 0, 0},
        { "z", PLY_INT, PLY_INT, int( offsetof(Point3i, coordinates[2])), 0, 0, 0, 0}
    };
    PlyProperty Point3fProperties[] ={
        { "x", PLY_FLOAT, PLY_FLOAT, int( offsetof(Point3f, coordinates[0])), 0, 0, 0, 0},
        { "y", PLY_FLOAT, PLY_FLOAT, int( offsetof(Point3f, coordinates[1])), 0, 0, 0, 0},
        { "z", PLY_FLOAT, PLY_FLOAT, int( offsetof(Point3f, coordinates[2])), 0, 0, 0, 0}
    };
//    std::cout << std::numeric_limits<long long>::min() << std::endl;

    cmdLineParse(argc - 1, argv + 1, params);
    if (!Count.set) {
        ShowUsage(argv[0]);
        return EXIT_FAILURE;
    }

    std::vector< Point< 3, int > > points, hullVertices;
    std::vector< Triangle > hullFaces;
    {
        Timer t;
        RandomPoints(points, Count.value, RandomSeed.value, Resolution.value);
        printf("Got random points: %.2f (s)\n", t.elapsed());
    }

 
    FILE * fp;

   fp = fopen ("inputpoints4.txt", "r");
   points.clear();
   int x,y,z;
   while ( fscanf (fp,"%d %d %d", &x, &y, &z ) != EOF )  {
       Point<3,int> p;
       p[0] = x; p[1] = y; p[2]=z;
       points.push_back(p);
      
   }

    std::ofstream os;
    os.open("points.txt");

    for (auto &p : points) {
        os << p[0] << "\t" << p[1] << "\t" << p[2] << std::endl;
    }
    os.close();
    {
        Timer t;
        switch (AlgorithmType.value) {
            case ALGORITHM_GIFT_WRAP: GiftWrapAlgorithm(points, hullVertices, hullFaces);
                break;
            case ALGORITHM_INCREMENTAL: IncrementalAlgorithm(points, hullVertices, hullFaces);
                break;
            default: fprintf(stderr, "[ERROR] Unrecognized algorithm type: %d\n", AlgorithmType.value), exit(0);
        }
        printf("Computed hull %d -> %d in %.2f (s)\n", Count.value, (int) hullVertices.size(), t.elapsed());
//         printf("%d %d  %.8f \n", Count.value, (int) hullVertices.size(), t.elapsed());
    }
    os.open("outpoints.txt");

    for (auto &p : hullVertices) {
        os << p[0] << "\t" << p[1] << "\t" << p[2] << std::endl;
    }
    os.close();
    if (Out.set) {
        if (ForVisualization.set) {
            std::vector< std::vector< unsigned int > > _hullFaces(hullFaces.size());
            std::vector< Point3f > _hullVertices(hullVertices.size());
            for (int i = 0; i < hullFaces.size(); i++) {
                _hullFaces[i].resize(3);
                for (int j = 0; j < 3; j++) _hullFaces[i][j] = hullFaces[i][j];
            }
            for (int i = 0; i < hullVertices.size(); i++) {
                _hullVertices[i][0] = (float) hullVertices[i][0];
                _hullVertices[i][1] = (float) hullVertices[i][1];
                _hullVertices[i][2] = (float) hullVertices[i][2];
            }
            PLY::Write(Out.value, _hullVertices, _hullFaces, Point3fProperties, 3, PLY_ASCII);
        } else PLY::Write(Out.value, hullVertices, NULL, &hullFaces, NULL, Point3iProperties, 3, PLY_ASCII);
    }

    return EXIT_SUCCESS;
}

