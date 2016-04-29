#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "Util/CmdLineParser.h"
#include "Util/Geometry.h"
#include "Util/Ply.h"


////////////////////////////
// Command line parsing info
cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineReadable ForVisualization( "viewable" );
cmdLineReadable* params[] = { &In , &Out , &ForVisualization , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input 2D polygon vertices>\n" , In.name );
	printf( "\t[--%s <output 2D triangulation>]\n" , Out.name );
	printf( "\t[--%s]\n" , ForVisualization.name );
}
// Command line parsing info
////////////////////////////

using namespace Geometry;


struct PVertex{
    Point< 2, int> p;
    bool isEar;
    unsigned int idx;
    PVertex *prev, *next;
    PVertex( Point<2, int> _p, unsigned int idx0);
    PVertex& addBefore( Point< 2, int> p, unsigned int idx0);
    unsigned int size() const;
    long long area2() const;
    
    static PVertex* Remove(PVertex* v);
};

PVertex::PVertex(Point<2, int> _p, unsigned int idx0){
    idx = idx0;
    p=_p;
    prev=next=this;
}

PVertex& PVertex::addBefore(Point<2,int> p, unsigned int idx0){
    PVertex* v = new PVertex(p,idx0);
    v->prev = this->prev;
    v->next = this;
    this->prev->next = v;
    this->prev = v;
    
    return *v;
}

PVertex* PVertex::Remove(PVertex* v){
    PVertex* temp = v->prev;
    
    temp->next = v->next;
    v->next->prev = temp;
    
    delete v;
    return temp==v? NULL:temp;

}

unsigned int PVertex::size() const{
    unsigned int s = 0;
    for (const PVertex* i=this;;i=i->next){
        s++;
        if (i->next == this) break;
    }
    return s;
}


long long Integral2(const Point<2,int> p[3]){
    long long a = 0;
    a += ((long long)(p[1][0]+p[0][0])) * (p[1][1] - p[0][1]);
    a += ((long long)(p[2][0]+p[1][0])) * (p[2][1] - p[1][1]);
    a += ((long long)(p[0][0]+p[2][0])) * (p[0][1] - p[2][1]);
    return a;
}

long long Area2(const Point<2,int> p,const Point<2,int> q,const Point<2,int> r){
    Point<2,int> triangle[3];
    triangle[0] = p;
    triangle[1] = q;
    triangle[2] = r;
    return Integral2(triangle);
}


long long PVertex::area2() const{
    Point<2, int> p[3];
    p[0][0] = 0;
    p[0][1] = 0;
    long long a = 0;
    for (const PVertex* i = this;;i = i->next){
        p[1] = i->p;
        p[2] = i->next->p;
        a += Integral2(p);
        if (i->next == this) break;
    }
    return a;
}

PVertex* BuildList(const std::vector< Point< 2 , int> >& polygon){
    PVertex* polygonList = new PVertex(polygon[0],0);
    for (int i = 1; i < polygon.size(); i++){
        polygonList->addBefore(polygon[i],i);
    }
    
    for (const PVertex* p = polygonList;;p = p->next){
        std::cout << p->idx << "\t" << p->p[0] << "\t" << p->p[1] << std::endl;
        
        if(p->next == polygonList) break;
    }
    
    return polygonList;
}

bool Left(Point<2,int> p, Point<2,int> q, Point<2,int> r){ 
    return Area2(p,q,r)>0;
}
bool LeftOn(Point<2,int> p, Point<2,int> q, Point<2,int> r){ 
    return Area2(p,q,r)>=0;
}
bool Collinear(Point<2,int> p, Point<2,int> q, Point<2,int> r){ 
    return Area2(p,q,r)==0;
}

bool InCone(Point<2,int> p, Point<2,int> q, Point<2,int> r, Point<2,int> s){
    if (LeftOn(p,q,r)){
        return (Left(q,s,p)&&Left(s,q,r));
    }
    
    return !(LeftOn(q,s,r)&&LeftOn(s,q,p));
}

bool InCones(PVertex *q, PVertex *s){

    return  InCone(q->prev->p,q->p, q->next->p, s->p)&&InCone(s->prev->p,s->p,s->next->p,q->p);
            
}



bool IsectProper(Point<2,int> p, Point<2,int> q, Point<2,int> r, Point<2,int> s){

    if(Collinear(p,q,r) || Collinear(p,q,s)) return false;
    if(Collinear(r,s,p) || Collinear(r,s,q)) return false;
    if(Left(p,q,r) == Left(p,q,s)) return false;
    if(Left(r,s,p) == Left(r,s,q)) return false;
    
    return true;
}

bool Between(Point<2,int> p, Point<2,int> q, Point<2,int> r){
    if (!Collinear(p,q,r)) return false;
    
    unsigned int dir = (p[0]!=q[0]?0:1);
    return 
        (p[dir] <= r[dir] && r[dir] <= q[dir]) ||
        (q[dir] <= r[dir] && r[dir] <= p[dir]);
}

bool Isect(Point<2,int> p, Point<2,int> q, Point<2,int> r, Point<2,int> s){
    return
            IsectProper(p,q,r,s) ||
            Between(p,q,r) || Between(p,q,s)||Between(r,s,p) ||Between(r,s,q);
}

bool DiagonalIsect(PVertex *r, PVertex *s){
    for( const PVertex* i=r;;i=i->next){
        if(i->prev!=r && i->prev !=s && i!=r && i!=s){
            if(Isect(r->p,s->p, i->prev->p, i->p)) return true;
        
        }
    
        if (i->next == r) break;
    }

    return false;
}

bool IsDiagonal(PVertex *r, PVertex *s){
    return  InCones(r,s)&&(!DiagonalIsect(r,s));
}

bool InitEars(PVertex *poly){
    for (PVertex *i=poly;;i=i->next){
        i->isEar = IsDiagonal(i->prev, i->next);
//        std::cout << i->idx << "\t" <<i->isEar << std::endl;
        if (i->next == poly) break;
        
    }
}

PVertex* ProcessEar(PVertex *e, std::vector< Triangle >& triangles){
    
    Triangle tri;
    tri[0] = e->prev->idx;
    tri[1] = e->idx;
    tri[2] = e->next->idx;
    
    e->prev->isEar = IsDiagonal(e->prev->prev,e->next);
    e->next->isEar = IsDiagonal(e->prev,e->next->next);
    
    triangles.push_back(tri);

    return PVertex::Remove(e);
}


template< class CType >
void Triangulate( const std::vector< Point< 2 , CType > >& polygon , std::vector< Triangle >& triangles )
{
    // first build a circularly linked list for the polygon
    PVertex* poly = BuildList(polygon);
    // initiallize the ears
    InitEars(poly);
    // do the ear removal algorithm 
    unsigned int sz = poly->size();
    int count = 0;
    while (sz > 3){
        count++;
        for (PVertex *v = poly;;v=v->next){
            if (v->isEar){
                poly = ProcessEar(v,triangles);
                sz--;
                break;
            }
            if (v->next == poly){break;}
            
            
        }
        
    }
    // process the last triangle
    ProcessEar(poly,triangles);
    
}



int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	// Information for how the vertices are described in a PLY file
	PlyProperty PointProperties[] =
	{
		{ "x" , PLY_INT , PLY_INT , int( offsetof( Point2i , coordinates[0] ) ) , 0 , 0 , 0 , 0 } ,
		{ "y" , PLY_INT , PLY_INT , int( offsetof( Point2i , coordinates[1] ) ) , 0 , 0 , 0 , 0 }
	};

	std::vector< Point< 2 , int > > polygonVertices;
	std::vector< Triangle > triangles;
	int file_type;
	PLY::Read( In.value , polygonVertices , NULL , NULL , NULL , PointProperties , NULL , 2 , file_type );

        if (polygonVertices.size() == 0){
            printf( "polygon is empty");
            exit(1);
        } else {
        
            std::cout << "polygon size: " << polygonVertices.size() << std::endl; 
            for (int i = 0;i < polygonVertices.size(); i++){
                std::cout << i << "\t" << polygonVertices[i][0] << "\t" << polygonVertices[i][1] << std::endl;
            }
        }
        
	Triangulate( polygonVertices , triangles );

	if( Out.set )
	{
		if( ForVisualization.set )
		{
			std::vector< std::vector< unsigned int > > polygons( triangles.size() );
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				polygons[i].resize(3);
				for( int j=0 ; j<3 ; j++ ) polygons[i][j] = triangles[i][j];
			}
			PLY::Write( Out.value , polygonVertices , polygons , PointProperties , 2 , file_type );
		}
		else PLY::Write( Out.value , polygonVertices , NULL , &triangles , NULL , PointProperties , 2 , file_type );
	}

	return EXIT_SUCCESS;
}

