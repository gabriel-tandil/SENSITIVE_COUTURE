/*
DelFEM (Finite Element Analysis)
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*! @file
@brief the basical definition refered from whole system including mesh and fem
@author Nobuyuki Umetani
*/

#if !defined(CAD_COM_H)
#define CAD_COM_H

#ifdef __VISUALC__
	#pragma warning( disable : 4786 )
#endif

namespace Cad{


/*!
@ingroup CAD
@brief the type of component 
*/
enum CAD_ELEM_TYPE{
  NOT_SET,  //!< not setted(for the error handling)
  VERTEX,   //!< vertex
  EDGE,     //!< edge
  LOOP,     //!< loop
	SOLID,		//!< solid
};

//! iterator go around loop
class IItrLoop
{
public:
  virtual void Begin() = 0;	//!< back to initial point of current use-loop
  virtual void operator++() = 0; //!< move to next edge
  virtual void operator++(int n)= 0;	//!< dummy operator (for ++)
  //! return current edge id and whether if this edge is same dirrection as loop
  virtual bool GetIdEdge(unsigned int& id_e, bool& is_same_dir) const = 0;	
  virtual bool ShiftChildLoop() = 0;	//!< move to next use-loop in this loop
  virtual bool IsEndChild() const = 0;	//!< return true if iterator go around
  virtual unsigned int GetIdVertex() const = 0;	//!< return current vertex id
  virtual unsigned int GetIdVertex_Ahead()  const = 0;	//!< return next vertex
  virtual unsigned int GetIdVertex_Behind() const = 0;	//!< return previous vertex
  virtual bool IsEnd() const = 0;	//!< return true if iterator go around
};
  
class IItrVertex
{
public:		
  virtual void operator++() = 0; //!< go around (cc-wise) loop around vertex 
  virtual void operator++(int n) = 0; //!< dummy operator (for ++)		
  //! cc-wise ahead  edge-id and its direction(true root of edge is this vertex)
  virtual bool GetIdEdge_Ahead( unsigned int& id_e, bool& is_same_dir) const = 0;
  //! cc-wise behind edge-id and its direction(true root of edge is this vertex)
  virtual bool GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const = 0;
  virtual unsigned int GetIdLoop() const = 0; //!< get loop-id		
  virtual bool IsEnd() const = 0; //!< return true if iterator go around
};     

}

#endif
