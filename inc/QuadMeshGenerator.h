//**************************************************************************
// File:     QuadMeshGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-07
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: �ı������������㷨ʵ��(������)
//**************************************************************************
#ifndef QUAD_MESH_IMP_H
#define QUAD_MESH_IMP_H
#include "USE_CGAL.h"

class QuadMeshGenerator
{
public:
	virtual void mesh() = 0;

	virtual void set_SourceMesh(HDS_mesh *_source_mesh) = 0;

	/** \brief ִ�������Ż�����
	* ���ʷֽ�����г��ȼ�Ȩ���Ϲ�˳�������Ż�
	*/
	//virtual bool quad_mesh_OPT() = 0;

	virtual ~QuadMeshGenerator() {}
};
#endif