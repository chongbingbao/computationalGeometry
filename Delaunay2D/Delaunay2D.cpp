#include <stdio.h>
#include <stdlib.h>
#include <unordered_set>
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Timer.h"
#include "Util/Ply.h"
#include <fstream>
#include <iostream>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <queue>
////////////////////////////

cmdLineParameter< char* > Out( "out" );
cmdLineParameter< int > Count( "count" ) , RandomSeed( "sRand" , 0 ) , Resolution( "res" , 1024 );
cmdLineReadable ForVisualization( "viewable" );
cmdLineReadable* params[] = { &Count , &Out , &Resolution , &RandomSeed , &ForVisualization , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input vertex count>\n" , Count.name );
	printf( "\t[--%s <output 2D triangulation>]\n" , Out.name );
	printf( "\t[--%s <random seed>=%d]\n" , RandomSeed.name , RandomSeed.value );
	printf( "\t[--%s <grid resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s]\n" , ForVisualization.name );
}
// Command line parsing info
////////////////////////////

using namespace Geometry;

void RandomPoints( std::vector< Point2i >& points , int count , int seed , int res )
{
	srand( seed );

	// Add distinct random points in a disk
	points.reserve( count );
	Point2i center;
	center[0] = center[1] = res/2;
	long long r = res / 2;
	std::unordered_set< long long > usedPoints;
	while( points.size()<count )
	{
		Point2i p;
		p[0] = rand() % res , p[1] = rand() % res;
		{
			long long d[] = { center[0] - p[0] , center[1] - p[1] };
			if( d[0]*d[0] + d[1]*d[1]>r*r ) continue;
		}
		long long key = ( ( long long )p[0] ) << 32 | ( ( long long )p[1] );
		if( usedPoints.find( key )==usedPoints.end() ) points.push_back( p ) , usedPoints.insert( key );
	}
}

long long SquaredArea(Point<3, int> &p0, Point<3, int> &p1, Point<3, int> &p2) {
    long long a = 0;
    a += ((long long) (p1[0] + p0[0])) * (p1[1] - p0[1]);
    a += ((long long) (p2[0] + p1[0])) * (p2[1] - p1[1]);
    a += ((long long) (p0[0] + p2[0])) * (p0[1] - p2[1]);
    return a*a;
}

long long SignedVolume(Point<3, int> &p0, Point<3, int> &p1, Point<3, int> &p2, Point<3, int> &p3) {
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
//    count = 0;
//    std::unordered_map<int,int> v_map;
    do{
        hullV.push_back(v->v);
//        v_map[v->idx] = count++;
        v = v->next;
    }while(v!=IncrementalAlg::vertices);
    
    IncrementalAlg::face_ptr f = IncrementalAlg::faces;
    do{
        Triangle tri;
        tri[0] = f->vertex[0]->idx;tri[1] = f->vertex[1]->idx;tri[2] = f->vertex[2]->idx;
        hullF.push_back(tri);
        f = f->next;
    }while(f!=IncrementalAlg::faces);
    
//          std::ofstream os;
//    os.open("faces.dat");
    
//    for (auto &t: hullF){
//        os << t[0] << "\t" << t[1] << "\t" << t[2] << std::endl;    
//    }
//    os.close(); 
}

template< class CType >
void Delaunay( std::vector< Point< 2 , CType > >& points , std::vector< Point< 2 , CType > >& dVertices , std::vector< Triangle >& dTriangles )
{
	// first lift to 3d points
	std::vector<Point<3,CType>> points_3d;
	for(auto &p: points){
		Point<3,CType> p_3d;
		p_3d[0] = p[0];p_3d[1] = p[1]; p_3d[2] = p[0]*p[0] + p[1]*p[1];
		points_3d.push_back(p_3d);
	}

	// construct the convex hull for the lifted points
	std::vector< Point< 3, CType > > hullV;
	std::vector< Triangle > hullF;
       IncrementalAlgorithm(points_3d, hullV, hullF);

        dVertices = points;
        // now project lower triangles, we use visibility test
        for (auto &t: hullF){
            Point<3,CType> testP;
            testP[0] = points_3d[t[0]][0];testP[1] = points_3d[t[0]][1];testP[2] = -1;
            auto vol = SignedVolume(points_3d[t[0]],points_3d[t[1]],points_3d[t[2]],testP);
            if (vol < 0) {
            // if visible, then the triangle is the lower hull
                dTriangles.push_back(t);
            }
        }
//    std::ofstream os;
//    os.open("faces.dat");
    
//    for (auto &t: dTriangles){
//        os << t[0] << "\t" << t[1] << "\t" << t[2] << std::endl;    
//    }
//    os.close(); 
        
        
}

int main( int argc , char* argv[] )
{
	PlyProperty Point2iProperties[] =
	{
		{ "x" , PLY_INT , PLY_INT , int( offsetof( Point2i , coordinates[0] ) ) , 0 , 0 , 0 , 0 } ,
		{ "y" , PLY_INT , PLY_INT , int( offsetof( Point2i , coordinates[1] ) ) , 0 , 0 , 0 , 0 }
	};
	PlyProperty Point3fProperties[] =
	{
		{ "x" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Point3f , coordinates[0] ) ) , 0 , 0 , 0 , 0 } ,
		{ "y" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Point3f , coordinates[1] ) ) , 0 , 0 , 0 , 0 } ,
		{ "z" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Point3f , coordinates[2] ) ) , 0 , 0 , 0 , 0 }
	};


	cmdLineParse( argc-1 , argv+1 , params );
	if( !Count.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	std::vector< Point< 2 , int > > points , dVertices;
	std::vector< Triangle > dTriangles;
	{
		Timer t;
		RandomPoints( points , Count.value , RandomSeed.value , Resolution.value );
		printf( "Got random points: %.2f(s)\n" , t.elapsed() );
	}

//    std::ofstream os;
//    os.open("points.txt");

//    for (auto &p : points) {
//        os << p[0] << "\t" << p[1]  << std::endl;
//    }
//    os.close();

	{
		Timer t;
		Delaunay( points , dVertices , dTriangles );
		printf( "Computed Deluanay Triangulation %d -> %d in %.2f (s)\n" , Count.value , (int)dVertices.size() , t.elapsed() );
	}

	if( Out.set )
	{
		if( ForVisualization.set )
		{
			std::vector< Point3f > _dVertices( dVertices.size() );
			std::vector< std::vector< unsigned int > > _dTriangles( dTriangles.size() );
			for( int i=0 ; i<dTriangles.size() ; i++ )
			{
				_dTriangles[i].resize( 3 );
				for( int j=0 ; j<3 ; j++ ) _dTriangles[i][j] = dTriangles[i][j];
			}
			for( int i=0 ; i<dVertices.size() ; i++ )
			{
				_dVertices[i][0] = (float) dVertices[i][0];
				_dVertices[i][1] = (float) dVertices[i][1];
				_dVertices[i][2] = 0.f;
			}
			PLY::Write( Out.value , _dVertices , _dTriangles , Point3fProperties , 3 , PLY_ASCII );
		}
		else PLY::Write( Out.value , dVertices , NULL , &dTriangles , NULL , Point2iProperties , 2 , PLY_ASCII );
	}

	return EXIT_SUCCESS;
}

