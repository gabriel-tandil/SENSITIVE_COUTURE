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


#if !defined(TRI_ARY_TOPOLOGY_H)
#define TRI_ARY_TOPOLOGY_H

class CTriAryTopology
{
public:
  class CItr
  {
    friend class CTriAryTopology;
  public:
    CItr(const CItr& ie):icur(ie.icur), iend(ie.iend), edge(ie.edge){}
    void operator++(int){ icur++; }
    bool IsEnd(){    
      if( icur==iend ) return true;
      return false;
    }
    unsigned int operator*(){ return edge[icur]; }    
  private:
    CItr(){}
    CItr(unsigned int icur, unsigned int iend, unsigned int* pe):edge(pe),icur(icur),iend(iend){}
    unsigned int icur;
    unsigned int iend;
    const unsigned int* edge;
  };
public:
  CTriAryTopology();
  ~CTriAryTopology();
  CItr GetItrEdge(unsigned int ino) const{ 
    assert(ino<nnode_);
    return CItr(edge_ind_[ino],edge_ind_[ino+1],edge_);
  }
  CItr GetItrElSuP(unsigned int ino) const{ 
    return CItr(elsup_ind_[ino],elsup_ind_[ino+1],elsup_);
  }  
  int GetElSuEl(unsigned int itri,unsigned int iedtri) const{ 
    assert( itri < ntri_ );
    assert( iedtri < 3 );
    return elsuel[itri*3+iedtri]; 
  }
  void SetTriAry(unsigned int ntri, unsigned int* atri, unsigned int nnode);  
private:
  void Clear();  
  void MakeEdge();
  void MakeElSurPo();
  void MakeElSuEl();    
private:
  unsigned int nnode_;
  
  unsigned int ntri_;
  unsigned int* aTri_;
  
  int* elsuel;
  
  unsigned int nedge_;
  unsigned int* edge_ind_;
  unsigned int* edge_;
  ////
  unsigned int  nelsup_;
  unsigned int* elsup_ind_;
  unsigned int* elsup_;
};


#endif