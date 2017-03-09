//**************************************************************************
// File:     TriMergeGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-07
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 四边形网格生成算法实现(基于三角形合并的算法)
//**************************************************************************
#ifndef TRI_MERGE_IMP_H
#define TRI_MERGE_IMP_H
#include "QuadMeshGenerator.h"

class MeshOptimizer;

class TriMergeGenerator: public QuadMeshGenerator
{
public:
	TriMergeGenerator();
	virtual void mesh();
	virtual void set_SourceMesh(HDS_mesh *_source_mesh);
	virtual ~TriMergeGenerator();
private:
	void   deleteFromUnactiveFront(HDS_mesh::Halfedge_handle eh);

	bool   get_Front_Merge(HDS_mesh::Halfedge_handle& eh);///< 从活动前沿中获取前沿边

	/** \brief 将前沿边从相应前沿链表中删除
	    \param eh 待删除半边的手柄
	*/
	void   deleteFromFront(HDS_mesh::Halfedge_handle eh);
	HDS_mesh::Face_handle merge_Tri(HDS_mesh::Face_handle& tri1, HDS_mesh::Face_handle& tri2);
private:
	HDS_mesh* _source_mesh;
	MeshOptimizer* _optimizer;
	std::list<HDS_mesh::Halfedge_handle> _active_front_list; ///< 活动前沿边链表
	std::list<HDS_mesh::Halfedge_handle> _unactive_front_list; ///< 非活动前沿边链表
};
#endif