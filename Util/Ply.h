/*

   Header for PLY polygon files.

   - Greg Turk, March 1994

   A PLY file contains a single polygonal _object_.

   An object is composed of lists of _elements_.  Typical elements are
   vertices, faces, edges and materials.

   Each type of element for a given object has one or more _properties_
   associated with the element type.  For instance, a vertex element may
   have as properties three floating-point values x,y,z and three unsigned
   chars for red, green and blue.

   ---------------------------------------------------------------

   Copyright (c) 1994 The Board of Trustees of The Leland Stanford
   Junior University.  All rights reserved.   

   Permission to use, copy, modify and distribute this software and its   
   documentation for any purpose is hereby granted without fee, provided   
   that the above copyright notice and this permission notice appear in   
   all copies of this software and that you do not sell the software.   

   THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
   EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
   WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

#ifndef PLY_INCLUDED
#define PLY_INCLUDED
#include "Util/Geometry.h"
#include "Util/PlyFile.h"


namespace PLY
{
	int DefaultFileType( void ){ return PLY_ASCII; }

	// Note, this code resets the vertex, edge, triangle, and tet lists.
	template< class Vertex >
	bool Read
	(
		const char* fileName ,
		std::vector< Vertex >& vertices ,
		std::vector< Geometry::Edge >* edges ,
		std::vector< Geometry::Triangle >* triangles ,
		std::vector< Geometry::Tet >* tets ,
		PlyProperty* properties , bool* propertyFlags , int propertyNum ,
		int& file_type ,
		char*** comments=NULL , int* commentNum=NULL
	);
	template<class Vertex>
	bool Write
	(
		const char* fileName ,
		const std::vector< Vertex >& vertices ,
		const std::vector< Geometry::Edge >* edges ,
		const std::vector< Geometry::Triangle >* triangles ,
		const std::vector< Geometry::Tet >* tets ,
		PlyProperty* properties , int propertyNum ,
		int file_type ,
		char** comments=NULL , int commentNum=0
	);
	template<class Vertex>
	bool Write
	(
		const char* fileName ,
		const std::vector< Vertex >& vertices ,
		const std::vector< std::vector< unsigned int > >& faces ,
		PlyProperty* properties , int propertyNum ,
		int file_type ,
		char** comments=NULL , int commentNum=0
	);
};
#include "Ply.inl"

#endif // PLY_INCLUDED
