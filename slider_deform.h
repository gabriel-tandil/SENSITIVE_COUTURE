/*
 *  slider_deform.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 12/20/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if !defined(SLIDER_DEFORM_H)
#define SLIDER_DEFORM_H


#include <iostream>
#include <vector>

#include "delfem/cad_obj2d_move.h"

class CSliderDeform
{
public:
  void Clear(){
    aSlider.clear();
    aLoopDeform.clear();    
  }
  unsigned int AddSlider(const std::string& str, double val, double min, double max)  // register slider
  {
    unsigned int isl = aSlider.size();
    aSlider.push_back( CSlider(str,val,min,max) );
    return isl;    
  }
  unsigned int count() const { return aSlider.size(); }
  void AddSliderParamToLoop(unsigned int id_l, unsigned int islider,
                            unsigned int ixy, unsigned int idir)  // set loop which slider it belongs and slider's direction
  {
    aLoopDeform.push_back( CLoopDeform(id_l,islider,ixy,idir) );    
  }
  void ChangeSlider(double val_pos, double val_pre);  // deform the loop with parameter change
  void GetLamVtx(unsigned int islider, Cad::CCadObj2D_Move& cad_2d, std::vector<double>& aLamXY_Vtx)  // set loop the parameter
  {    
    {
      unsigned int max_id_v = 0;    
      const std::vector<unsigned int>& aIdV = cad_2d.GetAryElemID(Cad::VERTEX);
      for(unsigned int iiv=0;iiv<aIdV.size();iiv++){
        max_id_v = ( aIdV[iiv] > max_id_v ) ? aIdV[iiv] : max_id_v;
      }
      aLamXY_Vtx.clear();
      aLamXY_Vtx.resize((max_id_v+1)*4,0);    
    }  
    if( islider >= aSlider.size() ) return;
    double val = aSlider[islider].GetValue();
    for(unsigned int ild=0;ild<aLoopDeform.size();ild++){
      aLoopDeform[ild].SetLambda(islider,val,cad_2d,aLamXY_Vtx);
    }        
  }
  void Draw(); // draw FFD structure
  double GetValueSlider(unsigned int islider, double& min, double& max) const {
    if( islider < aSlider.size() ){
      aSlider[islider].GetMinMax(min,max);
      return aSlider[islider].GetValue();
    }
    return 0;
  }
  void GetSliderProperty(unsigned int islider, std::string& name){
    if( islider < aSlider.size() ){
      name = aSlider[islider].name;
      return;
    }
    return;     
  }
  void SetLoopCenter(unsigned int id_l, double cnt_x, double cnt_y){
    for(unsigned int ild=0;ild<aLoopDeform.size();ild++){
      if( id_l == aLoopDeform[ild].id_l ){
        aLoopDeform[ild].cnt_x = cnt_x;
        aLoopDeform[ild].cnt_y = cnt_y;
      }
    }
  }
  void SetValueSlider(unsigned int islider, double val0){    
    if( islider >= aSlider.size() ) return;
    aSlider[islider].SetValue(val0);
  }
  bool MoveCad(unsigned int islider, double val_pos, double val_pre, 
               Cad::CCadObj2D_Move& cad_2d)
  {    
    bool res = false;
    for(unsigned int ild=0;ild<aLoopDeform.size();ild++){
      bool res0 = aLoopDeform[ild].MoveCad2(islider,val_pos,val_pre,cad_2d);
      if( res0 ){ res = true; }
    }  
    return res;    
  }
private:
  class CSlider{
  public:
    CSlider(const std::string& str, double val, double min, double max){
      name = str;
      this->val = val;
      this->min = min;
      this->max = max;      
    }
    double GetValue() const { return val; }
    void SetValue(double val){ 
      if(      val < min ){ val = min; }
      else if( val > max ){ val = max; }
      this->val = val; 
    }
    void GetMinMax(double& min, double& max) const {
      min = this->min;
      max = this->max;
    }
  public:
    std::string name;
  private:
    double min;
    double max;
    double val;
  };
  class CLoopDeform{
  public:
    CLoopDeform(unsigned int id_l, unsigned int islider, unsigned int ixy, unsigned int idir){
      this->id_l = id_l;
      this->islider = islider;
      this->ixy = ixy;
      this->idir = idir;
      this->cnt_x = 0;
      this->cnt_y = 0;
    }
    bool MoveCad2(unsigned int islider, double val_pos, double val_pre, 
                 Cad::CCadObj2D_Move& cad_2d)
    {      
      if( this->islider != islider ) return false;
      if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return false;
      std::vector< std::pair<unsigned int, Com::CVector2D> > aIdDist;  
      for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl.ShiftChildLoop()){
        for(itrl.Begin();!itrl.IsEnd();itrl++){        
          unsigned int id_v = itrl.GetIdVertex();  
          Com::CVector2D vec0 = cad_2d.GetVertexCoord(id_v);
          vec0.x -= cnt_x;
          vec0.y -= cnt_y;      
          //      std::cout << id_v << " " << ixy << " " << idir << " " << vec0.x << " " << vec0.y << std::endl;
          if(      ixy == 0 ){        
            if(   (idir == 0 && vec0.x > 0)
               || (idir == 1 && vec0.x < 0)
               || (idir == 2 && vec0.y > 0)
               || (idir == 3 && vec0.y < 0)
               ||  idir == 4 ){
              double x1 = vec0.x*(val_pos+1.0)/(val_pre+1.0);
              aIdDist.push_back( std::make_pair(id_v,Com::CVector2D(x1+cnt_x,vec0.y+cnt_y)) );
            }        
          }
          else if( ixy == 1 ){        
            if(   (idir == 0 && vec0.x > 0)
               || (idir == 1 && vec0.x < 0)
               || (idir == 2 && vec0.y > 0)
               || (idir == 3 && vec0.y < 0)
               ||  idir == 4 ){
              double y1 = vec0.y*(val_pos+1.0)/(val_pre+1.0);
              aIdDist.push_back( std::make_pair(id_v,Com::CVector2D(vec0.x+cnt_x,y1+cnt_y)) );
            }
          }
        }  
      }
      return cad_2d.MoveVertex(aIdDist);      
    }
    void SetLambda(unsigned int islider, double val,
                   const Cad::CCadObj2D_Move& cad_2d,
                   std::vector<double>& aLamXY_Vtx)
    {      
      if( this->islider != islider ) return;
      if( !cad_2d.IsElemID(Cad::LOOP,id_l) ) return;
      std::vector< std::pair<unsigned int, Com::CVector2D> > aIdDist;  
      for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl.ShiftChildLoop()){
        for(itrl.Begin();!itrl.IsEnd();itrl++){        
          unsigned int id_v = itrl.GetIdVertex();  
          Com::CVector2D vec0 = cad_2d.GetVertexCoord(id_v);
          vec0.x -= cnt_x;
          vec0.y -= cnt_y;      
          if(      ixy == 0 ){        
            if(   (idir == 0 && vec0.x > 0)
               || (idir == 1 && vec0.x < 0)
               || (idir == 2 && vec0.y > 0)
               || (idir == 3 && vec0.y < 0)
               ||  idir == 4 ){
              double s0 = vec0.x/(val+1.0);
              aLamXY_Vtx[id_v*4+0] = s0;
            }
          }
          else if( ixy == 1 ){
            if(   (idir == 0 && vec0.x > 0)
               || (idir == 1 && vec0.x < 0)   
               || (idir == 2 && vec0.y > 0) 
               || (idir == 3 && vec0.y < 0) 
               ||  idir == 4 ){
              double s0 = vec0.y/(val+1.0);
              aLamXY_Vtx[id_v*4+2] = s0;
            }
          }
        }  
      }          
    }
  public:
    double cnt_x, cnt_y;
    unsigned int id_l;
    unsigned int islider;
    unsigned int ixy; // 0:x, 0:y
    unsigned int idir; // 0:yposi 1:ynega 2:xposi 3:xnega 4:all
  };
private:
  std::vector<CSlider>  aSlider;
  std::vector<CLoopDeform> aLoopDeform;
};

#endif