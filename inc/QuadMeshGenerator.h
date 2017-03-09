//**************************************************************************
// File:     QuadMeshGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-07
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 四边形网格生成算法实现(抽象类)
//**************************************************************************
#ifndef QUAD_MESH_IMP_H
#define QUAD_MESH_IMP_H
#include "USE_CGAL.h"

class QuadMeshGenerator
{
public:
	virtual void mesh() = 0;

	virtual void set_SourceMesh(HDS_mesh *_source_mesh) = 0;

	/** \brief 执行网格优化操作
	* 对剖分结果进行长度加权拉氏光顺及拓扑优化
	*/
	//virtual bool quad_mesh_OPT() = 0;

	virtual ~QuadMeshGenerator() {}
};
#endif