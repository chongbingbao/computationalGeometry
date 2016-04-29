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

namespace GiftWrap {

    typedef std::array<Point<3, int>, 2> EdgePoint;

    int PivotOnEdge(EdgePoint &edge, std::vector< Point < 3, int>>&points) {
        int p = 0;
        int np = points.size();
        std::vector<double> volSet;
        long long area2 = SquaredArea(edge[0], edge[1], points[p]);
        for (int i = 1; i < np; i++) {
            long long volume6 = SignedVolume(edge[0], edge[1], points[p], points[i]);
            volSet.push_back(volume6);

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
    for (auto &v : V) {
        hullV.push_back(points[v]);
    }

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
	GiftWrapAlgorithm(points_3d, hullV, hullF);
        
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

    std::ofstream os;
    os.open("points.txt");

    for (auto &p : points) {
        os << p[0] << "\t" << p[1] << "\t" << p[2] << std::endl;
    }
    os.close();

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

