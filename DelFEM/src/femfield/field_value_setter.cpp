/*
 *  field_value_setter.cpp
 *  dfm_core
 *
 *  Created by Nobuyuki Umetani on 1/20/11.
 *  Copyright 2011 The University of Tokyo. All rights reserved.
 *
 */

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/femeqn/ker_emat_hex.h"
#include "delfem/femeqn/ker_emat_tet.h"

#include "delfem/field_value_setter.h"


//! set constant value to the field
bool Fem::Field::SetFieldValue_Constant
(unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
 Fem::Field::CFieldWorld& world,
 double val)
{
  if( !world.IsIdField(id_field_to) ){ return false; }
  Fem::Field::CField& field = world.GetField(id_field_to);
  
	if( idofns >= field.GetNLenValue() ) return false;
	if( !(field.GetFieldDerivativeType() & fdt) ) return false;
  
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
  
  if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
    const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(CORNER);      
    assert( world.IsIdNA(nsna.id_na_va) );
    CNodeAry& na = world.GetNA(nsna.id_na_va);			
    unsigned int id_ns = 0;	// target node segment;
    {
      if(      fdt & VALUE        ) id_ns = nsna.id_ns_va;
      else if( fdt & VELOCITY     ) id_ns = nsna.id_ns_ve;
      else if( fdt & ACCELERATION ) id_ns = nsna.id_ns_ac;
      else{ assert(0); }
      assert( na.IsSegID(id_ns) );
      const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns);
      assert( idofns < ns.Length() );
    }
    /*
     if( m_aElemIntp.size() == 0 || !this->IsPartial() ){	// 剛体の場合
     const unsigned int nnode = na.Size();
     CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns);
     for(unsigned int inode=0;inode<nnode;inode++){
     ns.SetValue(inode,idofns,const_val);
     }
     }*/
    {
      for(unsigned int iei=0;iei<aIdEA.size();iei++){
        const unsigned int id_ea = aIdEA[iei];
        assert( world.IsIdEA(id_ea) );
        const CElemAry& ea = world.GetEA(id_ea);
        unsigned int id_es_c_va = field.GetIdElemSeg(id_ea,CORNER,true,world);
        assert( ea.IsSegID(id_es_c_va) );
        na.SetValueToNodeSegment(ea,id_es_c_va,id_ns,idofns,val);
      }
    }
  }
  /*		if( m_na_e.id_na_va != 0 ){
   assert( world.IsIdNA(m_na_e.id_na_va) );
   CNodeAry& na = world.GetNA(m_na_e.id_na_va);
   const CNodeAry::CNodeSeg& ns_val = na.GetSeg(m_na_e.id_ns_va);
   for(unsigned int iei=0;iei<m_aElemIntp.size();iei++){
   unsigned int id_ea = m_aElemIntp[iei].id_ea;
   assert( world.IsIdEA(id_ea) );
   const CElemAry& ea = world.GetEA(id_ea);
   assert( ea.IsSegID( m_aElemIntp[iei].id_es_e_va ) );
   na.SetValueToNodeSegment(ea,m_aElemIntp[iei].id_es_e_va,
   m_na_e.id_ns_va,idofns,const_val);
   }
   }*/
  if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 ){
    const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(BUBBLE);
    assert( world.IsIdNA(nsna.id_na_va) );
    CNodeAry& na = world.GetNA(nsna.id_na_va);
    unsigned int id_ns = 0;	// target node segment;
    {
      if(      fdt & VALUE        ) id_ns = nsna.id_ns_va;
      else if( fdt & VELOCITY     ) id_ns = nsna.id_ns_ve;
      else if( fdt & ACCELERATION ) id_ns = nsna.id_ns_ac;
      else{ assert(0); }
      assert( na.IsSegID(id_ns) );
      const CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns);
      assert( idofns < ns.Length() );
    }
    //			const CNodeAry::CNodeSeg& ns_val = na.GetSeg(id_ns);
    for(unsigned int iea=0;iea<aIdEA.size();iea++){
      const unsigned int id_ea = aIdEA[iea];
      assert( world.IsIdEA(id_ea) );
      const CElemAry& ea = world.GetEA(id_ea);
      unsigned int id_es_b_va = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
      assert( ea.IsSegID(id_es_b_va) );
      na.SetValueToNodeSegment(ea,id_es_b_va,id_ns,idofns,val);
    }
  }  
  return true;
}

//! set mathematical expression to the field
bool Fem::Field::SetFieldValue_MathExp
(unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
 Fem::Field::CFieldWorld& world,
 std::string math_exp, double t)
{
  Fem::Field::CField& field = world.GetField(id_field_to);
  
	if( idofns >= field.GetNLenValue() ) return false;
	if( !(field.GetFieldDerivativeType() & fdt) ) return false;
  
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();  
  
  if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
    CEval eval;
    {
      eval.SetKey("x",0.0 );
      eval.SetKey("y",0.0 );
      eval.SetKey("z",0.0 );
      eval.SetKey("t",t   );
      if( !eval.SetExp(math_exp) ) return false;
    }
    const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(CORNER);
    assert( world.IsIdNA(nsna.id_na_va) );
    CNodeAry& na_va = world.GetNA(nsna.id_na_va);	
    unsigned int id_ns_va;
    {
      if(      fdt & VALUE        ) id_ns_va = nsna.id_ns_va;
      else if( fdt & VELOCITY     ) id_ns_va = nsna.id_ns_ve;
      else if( fdt & ACCELERATION ) id_ns_va = nsna.id_ns_ac;
      else{ assert(0); }
      assert( na_va.IsSegID(id_ns_va) );
      const CNodeAry::CNodeSeg& ns = na_va.GetSeg(id_ns_va);
      assert( idofns < ns.Length() );
    }
    assert( na_va.IsSegID(id_ns_va) );
    CNodeAry::CNodeSeg& ns_va = na_va.GetSeg(id_ns_va);
    assert( world.IsIdNA(nsna.id_na_co) );
    CNodeAry& na_co = world.GetNA(nsna.id_na_co);
    unsigned int id_ns_co = nsna.id_ns_co;
    assert( na_co.IsSegID(id_ns_co) );
    const CNodeAry::CNodeSeg& ns_co = na_co.GetSeg(id_ns_co);
    const unsigned int ndim = field.GetNDimCoord();
    assert( ns_co.Length() == ndim );
    assert( ndim <= 3 );
    double coord[3];
    if( !field.IsPartial() ){	// 親フィールドなら節点を全部参照している。
      for(unsigned int inode=0;inode<na_va.Size();inode++){
        unsigned int inode_co = field.GetMapVal2Co(inode);
        ns_co.GetValue(inode_co,coord);
        switch(ndim){
          case 3: eval.SetKey("z",coord[2]);
          case 2: eval.SetKey("y",coord[1]);
          case 1: eval.SetKey("x",coord[0]);
            break;
          default:
            assert(0);
            break;
        }
        double val = eval.Calc();
        ns_va.SetValue(inode,idofns,val);
      }
    }
    else{	// 要素に参照される節点だけを指している。
      for(unsigned int iei=0;iei<aIdEA.size();iei++){
        unsigned int id_ea = aIdEA[iei];          
        CElemAry& ea = world.GetEA(id_ea);
        unsigned int id_es_c_va = field.GetIdElemSeg(id_ea,CORNER,true,world);
        assert( ea.IsSegID(id_es_c_va) );
        const CElemAry::CElemSeg& es_c_va = ea.GetSeg(id_es_c_va);
        const unsigned int nnoes = es_c_va.Length();
        unsigned int noes[16];
        for(unsigned int ielem=0;ielem<ea.Size();ielem++){
          es_c_va.GetNodes(ielem,noes);
          for(unsigned int inoes=0;inoes<nnoes;inoes++){
            const unsigned int inode0 = noes[inoes];
            unsigned int inode_co0 = field.GetMapVal2Co(inode0);
            ns_co.GetValue(inode_co0,coord);
            switch(ndim){
              case 3: eval.SetKey("z",coord[2]);
              case 2: eval.SetKey("y",coord[1]);
              case 1: eval.SetKey("x",coord[0]);
                break;
              default:
                assert(0);
                break;
            }
            double val = eval.Calc();
            ns_va.SetValue(inode0,idofns,val);
          }
        }
      }
    }
  }
  ////////////////////////////////
  if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 ){		
    CEval eval;
    {
      eval.SetKey("x",0.0 );
      eval.SetKey("y",0.0 );
      eval.SetKey("z",0.0 );
      eval.SetKey("t",t );
      if( !eval.SetExp(math_exp) ) return false;
    }
    const CField::CNodeSegInNodeAry& nsna_b = field.GetNodeSegInNodeAry(BUBBLE);
    assert( world.IsIdNA(nsna_b.id_na_va) );
    CNodeAry& na_va = world.GetNA(nsna_b.id_na_va);	
    unsigned int id_ns;
    {
      if(      fdt & VALUE        ) id_ns = nsna_b.id_ns_va;
      else if( fdt & VELOCITY     ) id_ns = nsna_b.id_ns_ve;
      else if( fdt & ACCELERATION ) id_ns = nsna_b.id_ns_ac;
      else{ assert(0); }
    }
    assert( na_va.IsSegID(id_ns) );
    CNodeAry::CNodeSeg& ns_va = na_va.GetSeg(id_ns);
    if( !field.IsPartial() && nsna_b.id_na_co ){ // the bubble ns have coordinate and this is not partial
      std::cout << "Error!-->Not Implimented" << std::endl;
      assert(0);
      for(unsigned int inode=0;inode<na_va.Size();inode++){
        double val = eval.Calc();
        ns_va.SetValue(inode,idofns,val);
      }
    }
    else{
      assert( field.IsNodeSeg(CORNER,false,world) );
      const CNodeAry::CNodeSeg& ns_c_co = field.GetNodeSeg(CORNER,false,world);
      const unsigned int ndim = ns_c_co.Length();      
      for(unsigned int iei=0;iei<aIdEA.size();iei++){
        unsigned int id_ea = aIdEA[iei];
        assert( world.IsIdEA(id_ea) );
        CElemAry& ea = world.GetEA(id_ea);
        ////////////////
        unsigned int id_es_b_va = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
        assert( ea.IsSegID(id_es_b_va) );
        const CElemAry::CElemSeg& es_b_va = ea.GetSeg(id_es_b_va);
        assert( es_b_va.Length() == 1);
        ////////////////
        const unsigned int id_es_c_co = field.GetIdElemSeg(id_ea,CORNER, false, world);
        assert( ea.IsSegID(id_es_c_co) );
        const CElemAry::CElemSeg& es_c_co = ea.GetSeg(id_es_c_co);
        const unsigned int nnoes = es_c_co.Length();
        ////////////////
        unsigned int noes_c[16];
        double coord[3], coord_cnt[3];
        unsigned int inoes_b;
        for(unsigned int ielem=0;ielem<ea.Size();ielem++){
          // calc the position of element center
          es_c_co.GetNodes(ielem,noes_c);
          es_b_va.GetNodes(ielem,&inoes_b);
          for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] = 0.0; }
          for(unsigned int inoes=0;inoes<nnoes;inoes++){
            const unsigned int inode0 = noes_c[inoes];
            ns_c_co.GetValue(inode0,coord);
            for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] += coord[idim]; }
          }
          for(unsigned int idim=0;idim<ndim;idim++){ coord_cnt[idim] /= nnoes; }
          ////////////////
          switch(ndim){
            case 3: eval.SetKey("z",coord_cnt[2]);
            case 2: eval.SetKey("y",coord_cnt[1]);
            case 1: eval.SetKey("x",coord_cnt[0]);
              break;
            default:
              assert(0);
              break;
          }
          double val = eval.Calc();
          ns_va.SetValue(inoes_b,idofns,val);
//          std::cout << inoes_b << " " << idofns << " " << val << std::endl;
        }
      }
    }
  }
  if( field.GetNodeSegInNodeAry(EDGE).id_na_va != 0 ){
    std::cout << "Error!-->Not Implimented" << std::endl;
    assert(0);
    getchar();
  }  
  return true;
}

//! set random field to the field
void Fem::Field::SetFieldValue_Random
(unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
 Fem::Field::CFieldWorld& world,
 double ave, double range)
{
  if( !world.IsIdField(id_field_to) ) return;
  Fem::Field::CField field_to = world.GetField(id_field_to);
	srand(0);
	if( field_to.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
    const Fem::Field::CField::CNodeSegInNodeAry& m_na_c = field_to.GetNodeSegInNodeAry(CORNER);
		assert( world.IsIdNA(m_na_c.id_na_va) );
		CNodeAry& na = world.GetNA(m_na_c.id_na_va);
		if( !na.IsSegID(m_na_c.id_ns_va) ){
			std::cout << "Valueセグメントがない（速度場に設定しようとしている)" << std::endl;
			std::cout << "そのうちValueセグメントを追加で作る関数を加える" << std::endl;
			assert(0);
		}
		assert( na.IsSegID(m_na_c.id_ns_va) );
		CNodeAry::CNodeSeg& ns_val = na.GetSeg(m_na_c.id_ns_va);
		const unsigned int nnode = na.Size();
		for(unsigned int inode=0;inode<nnode;inode++){
			ns_val.SetValue(inode,0,rand()*2.0/(1.0+RAND_MAX)-1.0);
		}
	}
  if( field_to.GetNodeSegInNodeAry(EDGE).id_na_va != 0 ){
    const Fem::Field::CField::CNodeSegInNodeAry& m_na_e = field_to.GetNodeSegInNodeAry(EDGE);
		assert( world.IsIdNA(m_na_e.id_na_va) );
		CNodeAry& na = world.GetNA(m_na_e.id_na_va);
		assert( na.IsSegID(m_na_e.id_ns_va) );
		CNodeAry::CNodeSeg& ns_val = na.GetSeg(m_na_e.id_ns_va);
		const unsigned int nnode = na.Size();
		for(unsigned int inode=0;inode<nnode;inode++){
			ns_val.SetValue(inode,0,rand()*2.0/(1.0+RAND_MAX)-1.0);
		}
	}
  if( field_to.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 ){
    const Fem::Field::CField::CNodeSegInNodeAry& m_na_b = field_to.GetNodeSegInNodeAry(BUBBLE);  
		assert( world.IsIdNA(m_na_b.id_na_va) );
		CNodeAry& na = world.GetNA(m_na_b.id_na_va);
		assert( na.IsSegID(m_na_b.id_ns_va) );
		CNodeAry::CNodeSeg& ns_val = na.GetSeg(m_na_b.id_ns_va);
		const unsigned int nnode = na.Size();
		for(unsigned int inode=0;inode<nnode;inode++){
			ns_val.SetValue(inode,0,rand()*2.0/(1.0+RAND_MAX)-1.0);
		}
	}  
}

//! copy value to the field
void Fem::Field::SetFieldValue_Copy
(unsigned int id_field_to, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
 Fem::Field::CFieldWorld& world,
 unsigned int id_field_from)
{
} 


static double TriArea(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

static void TriDlDx(double dldx[][2], double const_term[],
                    const double p0[], const double p1[], const double p2[]){
  
	const double area = TriArea(p0,p1,p2);
	const double tmp1 = 0.5 / area;
  
	const_term[0] = tmp1*(p1[0]*p2[1]-p2[0]*p1[1]);
	const_term[1] = tmp1*(p2[0]*p0[1]-p0[0]*p2[1]);
	const_term[2] = tmp1*(p0[0]*p1[1]-p1[0]*p0[1]);
  
	dldx[0][0] = tmp1*(p1[1]-p2[1]);
	dldx[1][0] = tmp1*(p2[1]-p0[1]);
	dldx[2][0] = tmp1*(p0[1]-p1[1]);
  
	dldx[0][1] = tmp1*(p2[0]-p1[0]);
	dldx[1][1] = tmp1*(p0[0]-p2[0]);
	dldx[2][1] = tmp1*(p1[0]-p0[0]);
}

//! set gradient value to the field
bool Fem::Field::SetFieldValue_Gradient
(unsigned int id_field_to, Fem::Field::CFieldWorld& world,
 unsigned int id_field_from)
{  
	if( !world.IsIdField(id_field_from) ) return false;
	Fem::Field::CField& field_from = world.GetField(id_field_from);
  
	if( !world.IsIdField(id_field_to) ) return false;
	Fem::Field::CField& field_to = world.GetField(id_field_to);  
  
	if( field_to.GetAryIdEA().size() != 1 ){
		std::cout << "Error!-->Not Implimented" << std::endl;
		getchar();
		assert(0);
	}
  
	Fem::Field::INTERPOLATION_TYPE type_from, type_to;
	{
		const std::vector<unsigned int>& aIdEA_from = field_from.GetAryIdEA();
		const std::vector<unsigned int>& aIdEA_to   = field_to.GetAryIdEA();
		if( aIdEA_from.size() != aIdEA_to.size() ) return false;
		const unsigned int niea = aIdEA_from.size();
		assert( niea == 1 );
		if( aIdEA_from[0] != aIdEA_to[0] ) return false;
		const unsigned int id_ea = aIdEA_from[0];
		type_from = field_from.GetInterpolationType(id_ea,world);
		type_to   = field_to.GetInterpolationType(id_ea,world);
	}
  
	unsigned int nnoes, ndim;
	if( type_from==HEX11 && type_to==HEX1001 ){
		nnoes = 8; ndim = 3;
	}
	else if( type_from==TET11 && type_to==TET1001 ){
		nnoes = 4; ndim = 3;
	}
	else if( type_from==TRI11 && type_to==TRI1001 ){
		nnoes = 3; ndim = 2;
	}
	else{
		std::cout << "NotImplimented!" << std::endl;
		assert(0);
		getchar();
	}
  
	unsigned int id_ea = field_to.GetAryIdEA()[0];
	const CElemAry& ea = world.GetEA(id_ea);
	const CElemAry::CElemSeg& es_c_co = field_from.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_c_va = field_from.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_b_va = field_to.GetElemSeg(id_ea,BUBBLE,true, world);
  
	const Fem::Field::CField::CNodeSegInNodeAry& nsna_c = field_from.GetNodeSegInNodeAry(CORNER);
	assert( world.IsIdNA(nsna_c.id_na_co) );
	assert( world.IsIdNA(nsna_c.id_na_va) );
	const CNodeAry& na_c_co = world.GetNA(nsna_c.id_na_co);
	const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(nsna_c.id_ns_co);
	const CNodeAry& na_c_va = world.GetNA(nsna_c.id_na_va);
	const CNodeAry::CNodeSeg& ns_c_va = na_c_va.GetSeg(nsna_c.id_ns_va);
	unsigned int id_na_b_va = field_to.GetNodeSegInNodeAry(BUBBLE).id_na_va;
	unsigned int id_ns_b_va = field_to.GetNodeSegInNodeAry(BUBBLE).id_ns_va;
  
	assert( world.IsIdNA(id_na_b_va) );
	CNodeAry& na_b_va = world.GetNA(id_na_b_va);
	CNodeAry::CNodeSeg& ns_b_va = na_b_va.GetSeg(id_ns_b_va);
  
	double coord[16][3];
	double value[16];
	double grad[3];
  
  //	const unsigned int nnoes_c = es_c_co.GetSizeNoes();
	unsigned int noes[64];
  
	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		{	// 座標(coord)と値(value)を作る
			es_c_co.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				unsigned int ipoi0 = noes[inoes];
				assert( ipoi0 < na_c_co.Size() );
				ns_c_co.GetValue(ipoi0,coord[inoes]);
			}
      /*			if( id_es_c_va == id_es_c_co ){	
       for(unsigned int inoes=0;inoes<nnoes;inoes++){
       unsigned int ipoi0 = noes[inoes];
       assert( ipoi0 < na_c_va.Size() );
       na_c_va.GetValueFromNode(ipoi0,id_ns_c_va,0,val);
       value[inoes] = val;
       }
       }
       else{*/
      es_c_va.GetNodes(ielem,noes);
      for(unsigned int inoes=0;inoes<nnoes;inoes++){
        unsigned int ipoi0 = noes[inoes];
        assert( ipoi0 < na_c_va.Size() );
        ns_c_va.GetValue(ipoi0,&value[inoes]);
      }
      //			}
		}
		if( type_from == HEX11 ){
			double dndx[8][3];
			double an[8];
			double detjac;
			ShapeFunc_Hex8(0.0,0.0,0.0, coord, detjac,dndx,an);
			for(unsigned int idim=0;idim<ndim;idim++){ grad[idim] = 0.0; }
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
        for(unsigned int idim=0;idim<ndim;idim++){
          grad[idim] += value[inoes]*dndx[inoes][idim];
        }
			}
		}
		else if( type_from == TET11 ){
			double dldx[4][3];
			double const_term[4];
			TetDlDx(dldx,const_term, coord[0],coord[1],coord[2],coord[3]);
			for(unsigned int idim=0;idim<ndim;idim++){ grad[idim] = 0.0; }
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
        for(unsigned int idim=0;idim<ndim;idim++){
          grad[idim] += value[inoes]*dldx[inoes][idim];
        }
			}
		}
		else if( type_from == TRI11 ){
			double dldx[3][2];
			double const_term[3];
			TriDlDx(dldx,const_term, coord[0],coord[1],coord[2]);
			for(unsigned int idim=0;idim<ndim;idim++){ grad[idim] = 0.0; }
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
        for(unsigned int idim=0;idim<ndim;idim++){
          grad[idim] += value[inoes]*dldx[inoes][idim];
        }
			}
			grad[2] = 0.0;
		}
    //		std::cout << grad[0] << " " << grad[1] << " " << grad[2] << std::endl;
		{
			double norm = sqrt(grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2]);
			norm *=0.2;
      //			norm = 15.0;
			grad[0] /= norm;
			grad[1] /= norm;
			grad[2] /= norm;
		}
		{
			unsigned int noes[16];
			es_b_va.GetNodes(ielem,noes);
			unsigned int ipoi0 = noes[0];
			assert( ipoi0 < na_b_va.Size() );
			for(unsigned int idim=0;idim<ndim;idim++){
				ns_b_va.SetValue(ipoi0,idim,grad[idim]);
			}
		}
	}  
}

Fem::Field::CFieldValueSetter::CFieldValueSetter
(unsigned int id_field, Fem::Field::CFieldWorld& world)
{
  if( !world.IsIdField(id_field) ){ return; }
  Fem::Field::CField& field = world.GetField(id_field);
  const unsigned int nlen = field.GetNLenValue();
  id_field_ = id_field;
  aValueFieldDof_.resize(nlen*3);
}

void Fem::Field::CFieldValueSetter::SetMathExp
(const std::string& math_exp, unsigned int idof, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
 Fem::Field::CFieldWorld& world)
{
  if( !world.IsIdField(id_field_) ){ return; }
  Fem::Field::CField& field = world.GetField(id_field_);
  const unsigned int nlen = field.GetNLenValue();
  if( idof >= nlen ) return;
  if( fdt & Fem::Field::VALUE ){
    aValueFieldDof_[idof+nlen*0].SetValue(math_exp);
    return;
  }
  if( fdt & Fem::Field::VELOCITY ){
    aValueFieldDof_[idof+nlen*1].SetValue(math_exp);
    return;
  }
  if( fdt & Fem::Field::ACCELERATION ){
    aValueFieldDof_[idof+nlen*2].SetValue(math_exp);
    return;
  }    
}

void Fem::Field::CFieldValueSetter::SetConstant
(double val, unsigned int idof, Fem::Field::FIELD_DERIVATION_TYPE fdt,
 Fem::Field::CFieldWorld& world)
{
  if( !world.IsIdField(id_field_) ){ return; }
  Fem::Field::CField& field = world.GetField(id_field_);
  const unsigned int nlen = field.GetNLenValue();
  if( idof >= nlen ) return;
  if( fdt & Fem::Field::VALUE ){
    aValueFieldDof_[idof+nlen*0].SetValue(val);
    return;
  }
  if( fdt & Fem::Field::VELOCITY ){
    aValueFieldDof_[idof+nlen*1].SetValue(val);
    return;
  }
  if( fdt & Fem::Field::ACCELERATION ){
    aValueFieldDof_[idof+nlen*2].SetValue(val);
    return;
  }        
}


void Fem::Field::CFieldValueSetter::SetGradient
(unsigned int id_field_from, Fem::Field::CFieldWorld& world)
{
  if( !world.IsIdField(id_field_from) ){ return; }
  this->id_field_gradient_ = id_field_from;
}


bool Fem::Field::CFieldValueSetter::ExecuteValue
(double cur_time, Fem::Field::CFieldWorld& world)
{
  if( id_field_gradient_ != 0 ){
    SetFieldValue_Gradient(id_field_, world, id_field_gradient_);
  }
  if( !world.IsIdField(id_field_) ){ return false; }  
  Fem::Field::CField& field = world.GetField(id_field_);
  const unsigned int nlen = field.GetNLenValue();
  for(unsigned int ilen=0;ilen<nlen;ilen++){
    if(      aValueFieldDof_[ilen+nlen*0].itype == 1 ){
      SetFieldValue_Constant(id_field_,ilen,VALUE,world,aValueFieldDof_[ilen+nlen*0].val);      
    }
    else if( aValueFieldDof_[ilen+nlen*0].itype == 2 ){
      SetFieldValue_MathExp(id_field_,ilen,VALUE,world,aValueFieldDof_[ilen+nlen*0].math_exp,cur_time);
    }
  }
  for(unsigned int ilen=0;ilen<nlen;ilen++){
    if(      aValueFieldDof_[ilen+nlen*1].itype == 1 ){
      SetFieldValue_Constant(id_field_,ilen,VELOCITY,world,aValueFieldDof_[ilen+nlen*1].val);  
    }
    else if( aValueFieldDof_[ilen+nlen*1].itype == 2 ){
      SetFieldValue_MathExp(id_field_,ilen,VELOCITY,world,aValueFieldDof_[ilen+nlen*1].math_exp,cur_time);
    }
  }
  for(unsigned int ilen=0;ilen<nlen;ilen++){
    if(      aValueFieldDof_[ilen+nlen*2].itype == 1 ){
      SetFieldValue_Constant(id_field_,ilen,ACCELERATION,world,aValueFieldDof_[ilen+nlen*2].val);      
    }
    else if( aValueFieldDof_[ilen+nlen*2].itype == 2 ){
      SetFieldValue_MathExp(id_field_,ilen,ACCELERATION,world,aValueFieldDof_[ilen+nlen*2].math_exp,cur_time);
    }
  }    
}

