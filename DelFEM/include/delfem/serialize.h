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
@brief シリアライゼーションのためのクラス(Com::Serializer)のインターフェイス
@author Nobuyuki Umetani
*/

#if !defined(SERIALIZE_H)
#define SERIALIZE_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

//#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <cstring> //(strlen)

namespace Com{

//! class for serialization
class CSerializer
{
public:
	CSerializer(std::string fname, bool is_loading, bool isnt_binary = true)
		: m_buff_size(512)
	{
		m_file_name = fname;
		m_is_loading = is_loading;
		m_is_open = false;
		m_isnt_binary = isnt_binary;
		m_idepth = 1;
		////////////////
		m_pos = 0;
		if( !is_loading ){
			fp = fopen(m_file_name.c_str(),"w");
			fclose(fp);
		}
	}
	~CSerializer(){
		this->Close();
	}
	bool IsLoading(){ return m_is_loading; }
//	bool IsOpen(){ return m_is_open; }
	void Get(const char* format,...){
		assert( m_is_loading );
		if( !m_is_open ){
			if( m_isnt_binary ){ fp = fopen(m_file_name.c_str(),"r");  }
			else{                fp = fopen(m_file_name.c_str(),"rb"); }
			::fseek(fp,m_pos,SEEK_SET);
			m_is_open = true;
		}
		fgets(m_buffer,m_buff_size,fp);
		va_list ap;
		va_start(ap, format);
		{
			const unsigned int nform = strlen(format);
			const unsigned int nbuff = strlen(m_buffer);
			unsigned int jstat=0, jend=0;
			for(unsigned int i=0;i<nform;i++){
				if( format[i] != '%' ) continue;
				for(;jstat<nbuff;jstat++){
					if( m_buffer[jstat] != ' ' ) break;
				}
				for(jend=jstat;jend<nbuff;jend++){
					if( m_buffer[jend] == ' ' ) break;
				}
				if( format[i+1] == 's' ){
					sscanf(m_buffer+jstat,"%s", va_arg(ap,char*));
					i+=1;
				}
				else if( format[i+1] == 'd' ){
					sscanf(m_buffer+jstat,"%d", va_arg(ap,int*));
					i+=1;
				}
				else if( format[i+1] == 'l' && format[i+2] == 'f' ){
					sscanf(m_buffer+jstat,"%lf",va_arg(ap,double*));
					i+=2;
				}
				else{ assert(0); }
				jstat = jend+1;
				if( jstat >= nbuff ) break;
			}
		}
		va_end(ap);
	}
	void GetLine(char* buffer, unsigned int buff_size){
		assert( m_is_loading );
		if( !m_is_open ){
			if( m_isnt_binary ){ fp = fopen(m_file_name.c_str(),"r");  }
			else{                fp = fopen(m_file_name.c_str(),"rb"); }
			::fseek(fp,m_pos,SEEK_SET);
			m_is_open = true;
		}
		fgets(buffer,buff_size,fp);
    }
	void Out(const char* format,...){
		assert( !m_is_loading );
		if( !m_is_open ){
			if( m_isnt_binary ){ fp = fopen(m_file_name.c_str(),"a");  }
			else{                fp = fopen(m_file_name.c_str(),"ab"); }
			m_is_open = true;
		}
		va_list ap;
		va_start(ap, format);
		vfprintf(fp,format,ap);
		va_end(ap);
	}
	void ReadDepthClassName(char* class_name, unsigned int buff_size){
		assert( m_is_loading );
		if( !m_is_open ){
			if( m_isnt_binary ){ fp = fopen(m_file_name.c_str(),"r");  }
			else{                fp = fopen(m_file_name.c_str(),"rb"); }
			::fseek(fp,m_pos,SEEK_SET);
			m_is_open = true;
		}
		fgets(m_buffer,m_buff_size,fp);
		assert( m_buffer[0] == '#' );
		const unsigned int idepth0 = atoi(m_buffer+1);
		fgets(class_name,buff_size,fp);
		assert( m_idepth == idepth0 );
	}
	void WriteDepthClassName(const char* class_name){
		assert( !m_is_loading );
		if( !m_is_open ){
			if( m_isnt_binary ){ fp = fopen(m_file_name.c_str(),"a");  }
			else{                fp = fopen(m_file_name.c_str(),"ab"); }
			m_is_open = true;
		}
		fprintf(fp,"#%d\n",m_idepth);
		fprintf(fp,"%s\n",class_name);
	}
	void ShiftDepth( bool is_add ){
		if( is_add ){ m_idepth++; return; }
		assert( m_idepth > 1 );
		m_idepth -= 1;
	}
	unsigned int GetDepth(){ return m_idepth; }
	void Close(){
		if( m_is_open ){ fclose(fp); }
		m_is_open = false;
	}
private:
	std::string m_file_name;
	bool m_is_loading;
	bool m_is_open;
	bool m_isnt_binary;
	////////////////
	FILE* fp;
	long m_pos;
	////////////////
	const int m_buff_size;	// init as 512 in the constructor
	char m_buffer[512];
	////////////////
	unsigned int m_idepth;
};

}	// end namespace Com;

#endif
