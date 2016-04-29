#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED

namespace Geometry
{
	template< unsigned int Dim , class CType=int >
	struct Point
	{
		CType coordinates[Dim];

		CType& operator[] ( unsigned int i )       { return coordinates[i]; }
		CType  operator[] ( unsigned int i ) const { return coordinates[i]; }
	};

	template< unsigned int Dim >
	struct Simplex
	{
		unsigned int indices[Dim+1];

		unsigned int& operator[] ( unsigned int i )       { return indices[i]; }
		unsigned int  operator[] ( unsigned int i ) const { return indices[i]; }
	};

	typedef Point< 2 , int > Point2i;
	typedef Point< 3 , int > Point3i;
	typedef Point< 2 , float > Point2f;
	typedef Point< 3 , float > Point3f;
	typedef Point< 2 , double > Point2d;
	typedef Point< 3 , double > Point3d;
	typedef Simplex< 1 > Edge;
	typedef Simplex< 2 > Triangle;
	typedef Simplex< 3 > Tet;

};

#endif // GEOMETRY_INCLUDED