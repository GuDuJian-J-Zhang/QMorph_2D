//**************************************************************************
// File:     TriMergeGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-07
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: �ı������������㷨ʵ��(���������κϲ����㷨)
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

	bool   get_Front_Merge(HDS_mesh::Halfedge_handle& eh);///< �ӻǰ���л�ȡǰ�ر�

	/** \brief ��ǰ�رߴ���Ӧǰ��������ɾ��
	    \param eh ��ɾ����ߵ��ֱ�
	*/
	void   deleteFromFront(HDS_mesh::Halfedge_handle eh);
	HDS_mesh::Face_handle merge_Tri(HDS_mesh::Face_handle& tri1, HDS_mesh::Face_handle& tri2);
private:
	HDS_mesh* _source_mesh;
	MeshOptimizer* _optimizer;
	std::list<HDS_mesh::Halfedge_handle> _active_front_list; ///< �ǰ�ر�����
	std::list<HDS_mesh::Halfedge_handle> _unactive_front_list; ///< �ǻǰ�ر�����
};
#endif