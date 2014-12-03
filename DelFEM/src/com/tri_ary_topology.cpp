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

#include <iostream>
#include <assert.h>

#include "delfem/tri_ary_topology.h"

CTriAryTopology::CTriAryTopology()
{
  this->ntri_ = 0;
  this->nnode_ = 0;
  this->aTri_ = 0;
  
  this->nedge_ = 0;
  this->edge_ = 0;
  this->edge_ind_ = 0;
  
  this->nelsup_ = 0;
  this->elsup_ind_ = 0;
  this->elsup_ = 0;
  
  this->elsuel = 0;
}

CTriAryTopology::~CTriAryTopology()
{
  this->Clear();
}


void CTriAryTopology::SetTriAry(unsigned int ntri, unsigned int* atri, unsigned int nnode)
{
  this->Clear();
  this->ntri_ = ntri;
  this->nnode_ = nnode;
  if( aTri_ != 0 ){ delete[] aTri_; }
  aTri_ = new unsigned int [ntri_*3];
  for(unsigned int i=0;i<ntri*3;i++){ aTri_[i] = atri[i]; }
  this->MakeElSurPo();
  this->MakeElSuEl();
  this->MakeEdge();
}


void CTriAryTopology::MakeElSuEl()
{
  if( elsuel != 0 ){ delete[] elsuel; }
  elsuel = new int [ntri_*3];
  for(unsigned int i=0;i<ntri_*3;i++){ elsuel[i] = -1; }
  unsigned int edge_tbl[3][2] = { {1,2},{2,0},{0,1} };
  for(unsigned int itri=0;itri<ntri_;itri++){
    for(unsigned int inotri=0;inotri<3;inotri++){
      if( elsuel[itri*3+inotri] != -1 ) continue;
      //      unsigned int ino0 = aTri_[itri*3+inotri];
      unsigned int ino1 = aTri_[itri*3+edge_tbl[inotri][0]];
      unsigned int ino2 = aTri_[itri*3+edge_tbl[inotri][1]];
      for(unsigned int ielsup=elsup_ind_[ino1];ielsup<elsup_ind_[ino1+1];ielsup++){
        unsigned int jtri1 = elsup_[ielsup];
        if( itri == jtri1 ) continue;
        unsigned int* jno = aTri_ + jtri1*3;
        for(unsigned int jnotri=0;jnotri<3;jnotri++){
          unsigned int jno1 = jno[edge_tbl[jnotri][0]];
          unsigned int jno2 = jno[edge_tbl[jnotri][1]];
          if( ino1 == jno2 && ino2 == jno1 ){
            elsuel[itri *3+inotri] = jtri1;
            elsuel[jtri1*3+jnotri] = itri;
            break;
          }
        }
        if( elsuel[itri*3+inotri] != -1 ) break;
      }
      if( elsuel[itri*3+inotri] == -1 ){
//        std::cout << "Open Edge : " << itri << " " << inotri << "   " << ino1 << " " << ino2 << std::endl;
      }
    }
  }
}

void CTriAryTopology::MakeElSurPo()
{
  if( elsup_ind_ != 0 ){ delete[] elsup_ind_; }
  elsup_ind_ = new unsigned int [nnode_+1];
  for(unsigned int ino=0;ino<nnode_+1;ino++){ elsup_ind_[ino] = 0; }
  for(unsigned int itri=0;itri<ntri_;itri++){
    for(unsigned int inotri=0;inotri<3;inotri++){
      unsigned int ino1 = aTri_[itri*3+inotri];
      elsup_ind_[ino1+1] += 1;
    }
  }
  for(unsigned int ino=0;ino<nnode_;ino++){
    elsup_ind_[ino+1] += elsup_ind_[ino];
  }
  nelsup_ = elsup_ind_[nnode_];
  if( elsup_ != 0 ){ delete[] elsup_; }
  elsup_ = new unsigned int [nelsup_];
  for(unsigned int itri=0;itri<ntri_;itri++){
    for(unsigned int inotri=0;inotri<3;inotri++){
      unsigned int ino1 = aTri_[itri*3+inotri];
      unsigned int ind1 = elsup_ind_[ino1];
      elsup_[ind1] = itri;
      elsup_ind_[ino1] += 1;
    }
  }  
  for(int ino=nnode_;ino>=1;ino--){
    elsup_ind_[ino] = elsup_ind_[ino-1];
  }  
  elsup_ind_[0] = 0;
  /*
   for(unsigned int ino=0;ino<nnode_;ino++){
   std::cout << ino << " ---> ";
   for(unsigned int ielsup=elsup_ind[ino];ielsup<elsup_ind[ino+1];ielsup++){
   unsigned int itri1 = elsup[ielsup];      
   std::cout << itri1 << " ";
   }
   std::cout << std::endl;
   }
   */
}

void CTriAryTopology::MakeEdge()
{  
  assert( elsup_ind_ != 0 );
  assert( elsup_ != 0 );
  unsigned int* aflg = new unsigned int [nnode_];
  for(unsigned int ino=0;ino<nnode_;ino++){ aflg[ino] = 0; }
  edge_ind_ = new unsigned int [nnode_+1];
  edge_ind_[0] = 0;
  for(unsigned int ino=0;ino<nnode_;ino++){
    edge_ind_[ino+1] = edge_ind_[ino];
    aflg[ino] = ino;
    for(unsigned int ielsup=elsup_ind_[ino];ielsup<elsup_ind_[ino+1];ielsup++){
      unsigned int itri1 = elsup_[ielsup];      
      for(unsigned int inotri=0;inotri<3;inotri++){
        unsigned int ino1 = aTri_[itri1*3+inotri];
        if( aflg[ino1] == ino ) continue;
        edge_ind_[ino+1]++;
        aflg[ino1] = ino;
      }
    }
  }
  nedge_ = edge_ind_[nnode_];  
//  std::cout << "nedge : " << nedge_ << std::endl;∫
  edge_ = new unsigned int [nedge_];
  for(unsigned int ino=0;ino<nnode_;ino++){ aflg[ino] = 0; }
  unsigned int iedge = 0;
  for(unsigned int ino=0;ino<nnode_;ino++){
    assert( edge_ind_[ino] == iedge );
    aflg[ino] = ino;
    for(unsigned int ielsup=elsup_ind_[ino];ielsup<elsup_ind_[ino+1];ielsup++){
      unsigned int itri1 = elsup_[ielsup];      
      for(unsigned int inotri=0;inotri<3;inotri++){
        unsigned int ino1 = aTri_[itri1*3+inotri];
        if( aflg[ino1] == ino ) continue;
        edge_[iedge] = ino1;
        iedge++;
        aflg[ino1] = ino;
      }
    }
  }
  assert( iedge == nedge_ );
  
  delete[] aflg;
  /*
   for(unsigned int ino=0;ino<nnode_;ino++){
   std::cout << ino << " ---> ";
   for(unsigned int ielsup=elsup_ind[ino];ielsup<elsup_ind[ino+1];ielsup++){
   unsigned int itri1 = elsup[ielsup];      
   std::cout << itri1 << " ";
   }
   std::cout << std::endl;
   }
   */
  /*
   for(unsigned int ino=0;ino<nnode_;ino++){
   std::cout << ino << "  -->   ";
   for(unsigned int iedge=edge_ind_[ino];iedge<edge_ind_[ino+1];iedge++){
   unsigned int ino1 = edge_[iedge];
   std::cout << ino1 << " ";
   }
   std::cout << std::endl;
   }
   */
}

void CTriAryTopology::Clear()
{
  nedge_ = 0;    
  nelsup_ = 0;
  if( aTri_      != 0 ){ delete[] aTri_;      aTri_      = 0; }
  if( elsuel     != 0 ){ delete[] elsuel;     elsuel     = 0; }
  ////
  if( elsup_ind_ != 0 ){ delete[] elsup_ind_; elsup_ind_ = 0; }
  if( elsup_     != 0 ){ delete[] elsup_;     elsup_     = 0; }
  ////
  if( edge_ind_  != 0 ){ delete[] edge_ind_;  edge_ind_  = 0; }
  if( edge_      != 0 ){ delete[] edge_;      edge_      = 0; }
}


