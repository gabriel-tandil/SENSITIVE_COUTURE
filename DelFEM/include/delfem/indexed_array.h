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
@brief the interface and implementation of jagged-array class (Com::CIndexedArray)
@author Nobuyuki Umetani
*/


#if !defined(INDEXED_ARRAY_H)
#define INDEXED_ARRAY_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <assert.h>
#include <iostream>

namespace Com{
	
//! Jagged array
class CIndexedArray
{
public:
	CIndexedArray(){}

  void InitializeSize(unsigned int size){
    index.clear();
    index.resize(size+1,0);
    array.clear();
  }
  unsigned int Size() const{
    if( index.size() == 0 ){ return 0; }
    return index.size()-1;
  }

	//! initialize with transpose pattern
	void SetTranspose(unsigned int size, const CIndexedArray& crs){	
    index.clear();
		index.resize( size+1, 0 );
		for(unsigned int icrs=0;icrs<crs.array.size();icrs++){
			const unsigned int jno = crs.array[icrs];
			if( jno < size ){ index[jno+1]++; }
		}
		for(unsigned int j=0;j<size;j++){ index[j+1] += index[j]; }
		array.resize(index[size]);
		for(unsigned int i=0;i<crs.Size();i++){
		for(unsigned int icrs=crs.index[i];icrs<crs.index[i+1];icrs++){
			const unsigned int jno = crs.array[icrs];
			const unsigned int jcrs = index[jno];
			array[jcrs] = i;
			index[jno]++;
		}
		}
		for(int k=size-1;k>=0;k--){ index[k+1] = index[k]; }
		index[0] = 0;
/*
		////////////////
		std::cout << array.size() << std::endl;
		for(unsigned int l=0;l<size;l++){
			std::cout << l << " --> ";
			for(unsigned int icrs=index[l];icrs<index[l+1];icrs++){
				const unsigned int jno = array[icrs];
				std::cout << jno << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		////////////////
		std::cout << crs.array.size() << std::endl;
		for(unsigned int m=0;m<crs.size;m++){
			std::cout << m << " --> ";
			for(unsigned int icrs=crs.index[m];icrs<crs.index[m+1];icrs++){
				const unsigned int jno = crs.array[icrs];
				std::cout << jno << " ";
			}
			std::cout << std::endl;
		}
*/
	}
public:
  //! initialize as a dense matrix
	void Fill(unsigned int ncol, unsigned int nrow){
		index.resize(ncol+1);
		index[0] = 0;
		array.resize(ncol*nrow);
		for(unsigned int i=0;i<ncol;i++){
			index[i+1] = (i+1)*nrow;
			for(unsigned int j=0;j<nrow;j++){
				array[i*nrow+j] = j;
			}
		}
	}
	//! sort that each sub array have incremental odering
	void Sort(){
    if( index.size() == 0 ) return;
    const unsigned int size = index.size()-1;
		for(unsigned int ipoin=0;ipoin<size;ipoin++){
			const unsigned int is = index[ipoin  ];
			const unsigned int ie = index[ipoin+1];
			if( is == ie ) continue;
			assert( is < ie );
			int itmp;
			for(unsigned int i=is;i<ie-1;i++){
				for(unsigned int j=ie-1;j>i;j--){
					if( array[j] < array[j-1] ){
						itmp = array[j];
						array[j] = array[j-1];
						array[j-1] = itmp;
					}
				}
			}
		}
	}
	bool CheckValid() const {
    unsigned int size = index.size()-1;    
    if( index[size] > array.size() ){ return false; }
		return true;
	}
public:
	std::vector<unsigned int> index;
	std::vector<unsigned int> array;
};

}	// end namespace Fem

#endif
