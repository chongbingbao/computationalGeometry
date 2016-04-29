#include <stdio.h>
#include <stdlib.h>
#include <unordered_set>
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Timer.h"
#include "Util/Ply.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
////////////////////////////
// Command line parsing info
enum
{
	ALGORITHM_NAIVE ,
	ALGORITHM_QUICK_HULL ,
	ALGORITHM_GRAHAM ,
	ALGORITHM_COUNT
};
const char* AlgorithmNames[] = { "Naive" , "QuickHull" , "Graham" };

cmdLineParameter< char* > Out( "out" );
cmdLineParameter< int > Count( "count" ) , AlgorithmType( "aType" , ALGORITHM_NAIVE ) , RandomSeed( "sRand" , 0 ) , Resolution( "res" , 1024 );
cmdLineReadable ForVisualization( "viewable" );
cmdLineReadable* params[] = { &Count , &Out , &Resolution , &AlgorithmType , &RandomSeed , &ForVisualization , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input vertex count>\n" , Count.name );
	printf( "\t[--%s <output 2D triangulation>]\n" , Out.name );
	printf( "\t[--%s <algorithm type>=%d]\n" , AlgorithmType.name , AlgorithmType.value );
	for( int i=0 ; i<ALGORITHM_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , AlgorithmNames[i] );
	printf( "\t[--%s <random seed>=%d]\n" , RandomSeed.name , RandomSeed.value );
	printf( "\t[--%s <grid resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s]\n" , ForVisualization.name );
}
// Command line parsing info
////////////////////////////

using namespace Geometry;

long long Area2( const Point2i p[3] )
{
	long long a = 0;
	a += ( (long long)( p[1][0] + p[0][0] ) ) * ( p[1][1] - p[0][1] );
	a += ( (long long)( p[2][0] + p[1][0] ) ) * ( p[2][1] - p[1][1] );
	a += ( (long long)( p[0][0] + p[2][0] ) ) * ( p[0][1] - p[2][1] );
	return a;
}

long long Area2( Point2i p1 , Point2i p2 , Point2i p3 )
{
	Point2i p[] = { p1 , p2 , p3 };
	return Area2( p );
}

long long Dist2( Point2i p1 , Point2i p2 )
{
	long long dx = p1[0] - p2[0] , dy = p1[1] - p2[1];
	return dx*dx + dy*dy;
}

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

template< class CType >
void NaiveAlgorithm( std::vector< Point< 2 , CType > >& points , std::vector< Point< 2 , CType > >& hull )
{
	// Find the bottom-most (left-most) point
	int idx = 0;
	for( int i=1 ; i<points.size() ; i++ ) if( points[i][1]<points[idx][1] || ( points[i][1]==points[idx][1] && points[i][0]<points[idx][0] ) ) idx = i;

	// Remove it from the list of points
	hull.push_back( points[idx] ) , points[idx] = points.back() , points.pop_back();

	// Add the vertex making the right-most angle until either all the points are used, or you circle back around
	while( points.size() )
	{
		idx = 0;
		for( int i=1 ; i<points.size() ; i++ )
		{
			long long a = Area2( hull.back() , points[idx] , points[i] );
			if( a<0 || ( a==0 && Dist2( hull.back() , points[i] )>Dist2( hull.back() , points[idx] ) ) ) idx = i;
		}
		long long a = Area2( hull.back() , points[idx] , hull[0] );
		if( a<0 || ( a==0 && Dist2( hull.back() , hull[0] )>Dist2( hull.back() , points[idx] ) ) ) break;
		hull.push_back( points[idx] ) , points[idx] = points.back() , points.pop_back();
	}
}

template< class CType >
std::vector< Point<2, CType> > RightOf(const std::vector< Point< 2, CType > >& S, const Point< 2, CType > &a, const Point< 2, CType > &b) {
   std::vector< Point<2, int> > res;

    for (int i = 0; i < S.size(); i++) {
        if (Area2(S[i],a, b) < 0) {
            res.push_back(S[i]);
        }
    }
    return res;
}

template< class CType >
std::vector< Point<2, CType> > QuickHalfHull(const std::vector< Point< 2, CType > >& S, const Point< 2, CType > &a, const Point< 2, CType > &b) {
    if (S.size() == 0) {
        return S;
    } else {
        // the furthest will have the largest absolute area and (smallest distance to a in order to distinguish colinear case)
        int areaMax = -1;
		long long distance = -1;
        Point<2, CType> c;
        for (int i = 0; i < S.size(); i++) {
            if (abs(Area2(S[i], a, b)) > areaMax || (abs(Area2(S[i], a, b)) == areaMax && Dist2(a,S[i]) > distance)) {
				distance = Dist2(a,S[i]);              
				c = S[i];
                areaMax = abs(Area2(S[i], a, b));
//				if (areaMax == 0) std::cout << "area zero!" << std::endl;
            }
        }


        std::vector< Point<2, CType> > A = RightOf(S, a, c);
        std::vector< Point<2, CType> > B = RightOf(S, c, b);
        std::vector< Point<2, CType> > QA = QuickHalfHull(A, a, c);
        std::vector< Point<2, CType> > QB = QuickHalfHull(B, c, b);
        
		std::vector< Point<2, CType> > res;
        for (int i = 0; i < QA.size(); i++) {
            res.push_back(QA[i]);
        }
        res.push_back(c);
        for (int i = 0; i < QB.size(); i++) {
            res.push_back(QB[i]);
        }
        
        return res;
    }
}

template< class CType >
void QuickHullAlgorithm(std::vector< Point< 2, CType > >& points, std::vector< Point< 2, CType > >& hull) {
    // first find the horizonal(vertical) extremes
    Point<2, CType> a, b;
    a = points[0];
    b = points[0];
    for (int i = 0; i < points.size(); i++) {
        if (points[i][0] < a[0]) a = points[i];
        if (points[i][0] > b[0]) b = points[i];
    }

    for (int i = 0; i < points.size(); i++) {
        if (points[i][0] == a[0] && points[i][1] < a[1]) a = points[i];
        if (points[i][0] == b[0] && points[i][1] > b[1]) b = points[i];
    }


    std::vector< Point<2, CType> > A = RightOf(points, a, b);
    std::vector< Point<2, CType> > B = RightOf(points, b, a);
    std::vector< Point<2, CType> > QA = QuickHalfHull(A, a, b);
    std::vector< Point<2, CType> > QB = QuickHalfHull(B, b, a);

    hull.push_back(a);
    for (int i = 0; i < QA.size(); i++) {
        hull.push_back(QA[i]);
    }
    
    hull.push_back(b);
    for (int i = 0; i < QB.size(); i++) {
        hull.push_back(QB[i]);
    }
    
    
}

struct AngleLengthInfo{
	double angle;
	long long dist2;
	int idx;
	AngleLengthInfo(double angle0,long long dist0, int idx0){
		angle = angle0;
		dist2 = dist0;
		idx = idx0;
	}
};

bool AngleLengthComparator ( const AngleLengthInfo& l, const AngleLengthInfo& r)
{ return (l.angle < r.angle)||(l.angle==r.angle && l.dist2 < r.dist2); }


template< class CType >
void GrahamsAlgorithm( std::vector< Point< 2 , CType > >& points , std::vector< Point< 2 , CType > >& hull )
{
	// Find the bottom-most (left-most) point
	int idx = 0;
	for( int i=1 ; i<points.size() ; i++ ) if( points[i][1]<points[idx][1] || ( points[i][1]==points[idx][1] && points[i][0]<points[idx][0] ) ) idx = i;
	// now sort by the angle and lengthpoints
	std::vector<AngleLengthInfo> sortAngleLength;
	for (int i = 0; i < points.size(); i++){
		if (i!=idx){
			double angle = atan2(points[i][1]-points[idx][1],points[i][0]-points[idx][0]);
			long long dist2 = Dist2(points[i],points[idx]);
			sortAngleLength.push_back(AngleLengthInfo(angle,dist2,i));
		
		}
	}
	std::sort(sortAngleLength.begin(), sortAngleLength.end(), AngleLengthComparator);
	// now remove points having the same angle but smaller distance
	std::list<Point<2, int> > uniqueSortAngleLength;
	
	uniqueSortAngleLength.push_front(points[sortAngleLength.back().idx]);
	
	for (int i = sortAngleLength.size()-2; i >= 0;i-- ){
            // checking collinear
		if(Area2(points[idx],points[sortAngleLength[i].idx],points[sortAngleLength[i+1].idx])!=0){
			uniqueSortAngleLength.push_front(points[sortAngleLength[i].idx]);
                } 
//                else {
//                    std::cout << "colinear happen!" << std::endl;
//                }
	}

	hull.push_back(points[idx]);
	hull.push_back(uniqueSortAngleLength.front());
	
        uniqueSortAngleLength.pop_front();


	while(!uniqueSortAngleLength.empty()){
		if(Area2(uniqueSortAngleLength.front(),hull[hull.size()-2],hull[hull.size()-1])>0){
			hull.push_back(uniqueSortAngleLength.front());
			
                    uniqueSortAngleLength.pop_front();    
		
		}else{
			hull.pop_back();
		}
		
	}

}

	/////////////////////////////
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

	std::vector< Point< 2 , int > > points , hullVertices;
	{
		Timer t;
		RandomPoints( points , Count.value , RandomSeed.value , Resolution.value );
		printf( "Got random points: %.2f(s)\n" , t.elapsed() );
	}
//        std::ofstream is;
//        is.open("input.dat");
//        for(int i = 0; i < points.size(); i++){
//            is << points[i][0] << "\t" << points[i][1] << std::endl;
//        }

	{
		Timer t;
		switch( AlgorithmType.value )
		{
			case ALGORITHM_NAIVE:          NaiveAlgorithm( points , hullVertices ) ; break;
			case ALGORITHM_QUICK_HULL: QuickHullAlgorithm( points , hullVertices ) ; break;
			case ALGORITHM_GRAHAM:       GrahamsAlgorithm( points , hullVertices ) ; break;
			default: fprintf( stderr , "[ERROR] Unrecognized algorithm type: %d\n" , AlgorithmType.value ) , exit( 0 );
		}
		printf( "Computed hull %d -> %d in %.12f (s)\n" , Count.value , (int)hullVertices.size() , t.elapsed() );
//		printf( "%d  %d  %.12f \n" , Count.value , (int)hullVertices.size() , t.elapsed() );

	}

	if( Out.set )
	{
		if( ForVisualization.set )
		{
			std::vector< std::vector< unsigned int > > polygons(1);
			std::vector< Point3f > _hullVertices( hullVertices.size() );
			polygons[0].resize( hullVertices.size() );
			for( int i=0 ; i<hullVertices.size() ; i++ ) polygons[0][i] = i;
			for( int i=0 ; i<hullVertices.size() ; i++ )
			{
				_hullVertices[i][0] = (float) hullVertices[i][0];
				_hullVertices[i][1] = (float) hullVertices[i][1];
				_hullVertices[i][2] = 0.f;
			}
			PLY::Write( Out.value , _hullVertices , polygons , Point3fProperties , 2 , PLY_ASCII );
		}
		else PLY::Write( Out.value , hullVertices , NULL , NULL , NULL , Point2iProperties , 2 , PLY_ASCII );
	}

	return EXIT_SUCCESS;
}

