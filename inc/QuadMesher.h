//**************************************************************************
// File:     QuadMesher.h
// Author:   Zhang Jun
// Date:     2016-05-05
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: ��߽ṹ����������Ƭ��Ź���(��ԵĶ���ͱߵı��һ��  ��ʱ�뷽��)
//**************************************************************************

#ifndef _QMORPHMESHER_H
#define _QMORPHMESHER_H

#include "USE_CGAL.h"
#include "toolsdefine.h"
#include "sldelaunay2d.h"

#define EPS_QMorph 1.0e-4

/*class Vector3D;*/
class QuadMeshGenerator;
struct DelaunayAFT2DInterface;
class QuadMesher
{
public:
	QuadMesher(HDS_mesh &backgroundMesh);
	~QuadMesher() {};

	/** \brief ��inp�ļ��ж�ȡ����������
	*    ���ݽڵ������ж�����ά��
	   \param file_name  �ļ���(.inp��ʽ)
	*/
	bool read_file(std::string file_name);

	bool read_tri_mesh(DelaunayAFT2DInterface& tri_mesh);

	bool mesh(const char* method);

	/** \brief �������
	*
		\param file_name �ļ���
	*/
	void write_my_mesh(std::string file_name);

	/** \brief nas��ʽ�������ļ�
	*
		\param file_name �ļ���
	*/
	void write_my_mesh_nas(std::string file_name);

private:

	/** \brief ͨ�������������������α߽��н�������߲�������������������������ں����ı������������
	*/
	bool localReconnect();

	//----------------------------------------------------------------------------//

	void  print_node_long_format_nas(int id,double x,double y,double z,char* strBuffer);

private:
	std::string _file_name; ///< inp�ļ���

	HDS_mesh &_backgroundMesh; ///< ��������

	///< �ı������������㷨ʵ��
	QuadMeshGenerator* _meshGenerator;
};

class QMorphInput : public DelaunayAFT2DInterface::Input
{
public:
	QMorphInput(unsigned nnode,const double* nodes,
		unsigned  nedge,const int* edges);

	void      node_begin(){_ni=0;}
	bool      node_next(int& id,double xy[2]);
	void      edge_begin(){_ei=0;}
	bool      edge_next(int nodes[2]);
private:
	unsigned  _nnode,_nedge,_ni,_ei;
	const double* _nodes;
	const int*    _edges;
};
#endif