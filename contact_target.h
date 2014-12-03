#if !defined(CONTACT_TARGET_H)
#define CONTACT_TARGET_H


#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/spatial_hash_grid3d.h"

class CContactTarget3D{
public:
	virtual void Draw() const = 0;
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const = 0;
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const = 0;
  virtual void GetMesh(std::vector<unsigned int>& aTri,
                       std::vector<double>& aXYZ,
                       double elen) const = 0;
};

class CContactTarget3D_Mesh : public CContactTarget3D
{
public:
	CContactTarget3D_Mesh();
	~CContactTarget3D_Mesh();
	virtual void Draw() const;	
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const;
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const;
//	void Load_Off(const std::string& fname);
//	void Load_Gmv(const std::string& fname);
//  void Load_Ply(const std::string& fname);
  void GetCenterWidth(double& cx, double& cy, double& cz,  double& wx, double& wy, double& wz);
	void Translate(double x, double y, double z);	
	void BuildBoxel();
	void SetHole(bool is_hole){	this->is_hole = is_hole; }		
  virtual void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const;
  void SetMesh(const std::vector<unsigned int>& aTri, const std::vector<double>& aXYZ);  
private:private:
	double Distance_Mesh
	(double px, double py, double pz,
	 double n[3]) const;
	double Distance_Mesh_Boxel
	(double px, double py, double pz,
	 double n[3]) const;
	
	unsigned int FindInOut_IntersectionRay
	(double px, double py, double pz,
	 const double dir[3]) const;
	unsigned int FindInOut_IntersectionRay_Boxel
	(double px, double py, double pz,
	 const double dir[3]) const;
	virtual unsigned int FindInOut(double px, double py, double pz) const;
	virtual unsigned int FindInOut_Boxel
	(double px, double py, double pz) const;
private:
	bool is_hole;
	unsigned int nnode_;
	double* pXYZs_;
	unsigned int ntri_;	
	unsigned int* aTri_;
	CSpatialHash_Grid3D* pBoxel_;		
	mutable std::vector<unsigned int> aFlgTriUsed;
	mutable std::vector<unsigned int> aIndTriCand;
};

class CContactTarget3D_AdaptiveDistanceField3D : public CContactTarget3D
{
public:
	CContactTarget3D_AdaptiveDistanceField3D();
    ~CContactTarget3D_AdaptiveDistanceField3D();
	void SetUp(const CContactTarget3D& ct, double bb[6]);
	void Draw() const;
  void SetFaceColor(double r, double g, double b){ color_[0] = r; color_[1] = g; color_[2] = b; }
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const;
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double dir[3]) const { return false; }   
  virtual void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const{}
  /////
	void BuildIsoSurface_MarchingCube();
	void BuildMarchingCubeEdge();
	void SetShowCage(bool is_show){ this->is_show_cage = is_show; }
private:
	class CNode
	{
	public:
		CNode();
		CNode(const CNode& no);
		void SetCornerDist(const CContactTarget3D& ct);
		void Draw_Wire() const;
		void DrawThisAndChild_Wire(const std::vector<CNode>& aNo) const ;
		void MakeChildTree(const CContactTarget3D& ct, std::vector<CNode>& aNo, double min_hw, double max_hw);
		double FindDistNormal
		(double px, double py, double pz,
		 double n[3],
		 const std::vector<CNode>& aNo) const;
		void GenerateIsoSurface
		(std::vector<double>& aTri,
		 const std::vector<CNode>& aNo) const;		
	public:
		double cent_[3];
		double hw_;
		int ichilds_[8];
		double dists_[8];
	};
private:
	std::vector<CNode> aNode;
	double dist_min, dist_max;
	unsigned int nIsoTri_;
	double* aIsoTri_;
	double* aIsoEdge_;
	bool is_show_cage;
  double color_[3];
};

class CContactTarget3D_Sphere : public CContactTarget3D
{
public:
	CContactTarget3D_Sphere(double rad, double cent[3], bool is_out){        
    cent_[0] = cent[0];
    cent_[1] = cent[1];
    cent_[2] = cent[2];
    radius_ = rad;
    this->is_out_ = is_out;    
  }
	virtual void Draw() const{    
    const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
    const bool is_texture = ::glIsEnabled(GL_TEXTURE_2D);  
    ::glDisable(GL_LIGHTING);
    ::glDisable(GL_TEXTURE_2D);  
    ////
    ::glLineWidth(1);
    ::glColor3d(1,0,0);
    ::glMatrixMode(GL_MODELVIEW);
    ::glPushMatrix();
    ::glTranslated(cent_[0],cent_[1],cent_[2]);
    
    const unsigned int nlg = 32;
    const unsigned int nlt = 18;	
    const double rlg = 6.28/nlg;	// longtitude
    const unsigned int ndiv = 32;
    const double rdiv = 6.28/ndiv;
    for(unsigned int ilg=0;ilg<nlg;ilg++){
      ::glBegin(GL_LINE_LOOP);
      for(unsigned int idiv=0;idiv<ndiv;idiv++){
        ::glVertex3d(radius_*cos(idiv*rdiv)*cos(ilg*rlg), 
                     radius_*cos(idiv*rdiv)*sin(ilg*rlg),
                     radius_*sin(idiv*rdiv) );
      }		
      ::glEnd();
    }
    for(unsigned int ilt=0;ilt<nlt;ilt++){
      const double d = ((double)ilt/nlt-0.5)*radius_*2.0;
      const double r0 = sqrt(radius_*radius_-d*d);
      ::glBegin(GL_LINE_LOOP);
      for(unsigned int idiv=0;idiv<ndiv;idiv++){
        ::glVertex3d(r0*cos(idiv*rdiv), 
                     r0*sin(idiv*rdiv),
                     d);
      }		
      ::glEnd();
    }
    //   ::glutWireSphere(radius_,32,32);
    ::glPopMatrix();
    ////
    if(is_lighting){ glEnable(GL_LIGHTING); }
    if(is_texture ){ glEnable(GL_TEXTURE_2D); }    
  }
  
	virtual double Projection
	(double px, double py, double pz,
	 double n[3]) const {    
    double dir[3] = { px-cent_[0], py-cent_[1], pz-cent_[2] };
    const double len = sqrt( dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2] );		
    const double invlen = 1.0/len;
    if( !is_out_ ){
      n[0] = -dir[0]*invlen;
      n[1] = -dir[1]*invlen;		
      n[2] = -dir[2]*invlen;
      return +len-radius_;
    }
    n[0] = dir[0]*invlen;
    n[1] = dir[1]*invlen;		
    n[2] = dir[2]*invlen;	
    return radius_-len;    
  }
	virtual unsigned int FindInOut(double px, double py, double pz) const{    
    double n[3];
    double pd = this->Projection(px, py, pz, n);	
    if( !is_out_ ) pd *= -1.0;
    if( pd > 0 ){ return 0; }
    return 1;
  }
  
  virtual bool IntersectionPoint
  (double p[3], 
   const double org[3], const double d[3]) const{    
    const double q[3] = { org[0]-cent_[0], org[1]-cent_[1], org[2]-cent_[2] };
    const double a = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    const double b = q[0]*d[0] + q[1]*d[1] + q[2]*d[2];
    const double c = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] - radius_*radius_;
    const double det = b*b-a*c;
    if( det < 0 )  return false;
    const double t = (-b+sqrt(det))/a;
    p[0] = org[0] + t*d[0];
    p[1] = org[1] + t*d[1];
    p[2] = org[2] + t*d[2];  
    return true;
  }
    
  
  virtual void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const
  {    
    double pi = 3.1415;
    const unsigned int nlg = 32;
    const unsigned int nlt = 18;	
    const double rlg = 2*pi/nlg;
    const double rlt = 1*pi/nlt;
    aXYZ.resize(nlg*(nlt-1)*3+2*3);
    for(unsigned int ilg=0;ilg<nlg;ilg++){    
      for(unsigned int ilt=0;ilt<(nlt-1);ilt++){
        double alt = ilt*rlt+rlt-pi*0.5;
        aXYZ[(ilg*(nlt-1)+ilt)*3+0] = cent_[0]+radius_*cos(alt)*sin(ilg*rlg);
        aXYZ[(ilg*(nlt-1)+ilt)*3+1] = cent_[1]+radius_*cos(alt)*cos(ilg*rlg);
        aXYZ[(ilg*(nlt-1)+ilt)*3+2] = cent_[2]+radius_*sin(alt);
      }
    }
    const unsigned int iu1 = nlg*(nlt-1);
    {
      aXYZ[iu1*3+0] = cent_[0];
      aXYZ[iu1*3+1] = cent_[1];
      aXYZ[iu1*3+2] = cent_[2]+radius_;
    }
    const unsigned int iu2 = nlg*(nlt-1)+1;  
    {    
      aXYZ[iu2*3+0] = cent_[0];
      aXYZ[iu2*3+1] = cent_[1];
      aXYZ[iu2*3+2] = cent_[2]-radius_;   
    }
    ////
    aTri.resize(nlg*(nlt-2)*2*3+nlg*2*3);
    for(unsigned int ilg=0;ilg<nlg;ilg++){        
      for(unsigned int ilt=0;ilt<nlt-2;ilt++){
        unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;
        aTri[(ilg*(nlt-2)+ilt)*6+0] = ilg*(nlt-1)+ilt;
        aTri[(ilg*(nlt-2)+ilt)*6+1] = iug*(nlt-1)+ilt+1;
        aTri[(ilg*(nlt-2)+ilt)*6+2] = iug*(nlt-1)+ilt;
        ////
        aTri[(ilg*(nlt-2)+ilt)*6+3] = ilg*(nlt-1)+ilt;
        aTri[(ilg*(nlt-2)+ilt)*6+4] = ilg*(nlt-1)+ilt+1;
        aTri[(ilg*(nlt-2)+ilt)*6+5] = iug*(nlt-1)+ilt+1;
      }
    }
    const unsigned int itri1 = nlg*(nlt-2)*2;
    for(unsigned int ilg=0;ilg<nlg;ilg++){        
      unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
      aTri[(itri1+ilg)*3+0] = nlt-2+(nlt-1)*ilg;
      aTri[(itri1+ilg)*3+1] = iu1;
      aTri[(itri1+ilg)*3+2] = nlt-2+(nlt-1)*iug;
    }
    const unsigned int itri2 = nlg*(nlt-2)*2+nlg;
    for(unsigned int ilg=0;ilg<nlg;ilg++){        
      unsigned int iug = ( ilg == nlg-1 ) ? 0 : ilg+1;    
      aTri[(itri2+ilg)*3+0] = (nlt-1)*iug;
      aTri[(itri2+ilg)*3+1] = iu2;
      aTri[(itri2+ilg)*3+2] = (nlt-1)*ilg;    
    }        
  }
private:
	double cent_[3];
	double radius_;
	bool is_out_;	
};


#endif
