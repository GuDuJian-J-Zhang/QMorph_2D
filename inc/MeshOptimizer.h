//**************************************************************************
// File:     MeshOptimizer.h
// Author:   Zhang Jun
// Date:     2016-11-10
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: ���������Ż�
//**************************************************************************
#ifndef MESH_OPTIMIZER_H
#define MESH_OPTIMIZER_H

#include "USE_CGAL.h"
const double EPS_QuadMesh = 1.0e-4;
class MeshOptimizer
{
public:
	MeshOptimizer() { _source_mesh = nullptr; }

	void set_SourceMesh(HDS_mesh *source_mesh);

	bool mesh_OPT(); ///< ʵ�ָýӿ�, ����������Ż�

	/** \brief ���������ȫ�ֹ�˳
	    \return n ��˳����
	*/
	bool global_smooth(int n);

	/** \brief �ֲ���˳
	* ÿ���²�����Ԫ����Ҫ���оֲ���˳\n
	* ��˳����Ϊ�µ�Ԫ���ĸ����㼰���һ���ڽӶ���\n
	* ��˳�㷨Ϊ���ڳ��ȵ�Laplacian��˳
		\param[int] smooth_node_list ��ִ�й�˳�Ľڵ㼯
	*/
	void localSmoothing(std::list<HDS_mesh::Vertex_handle>& smooth_node_list);

	/** \brief ���������ε�Ԫ������
	    \param v0~v2 ��������������
	*/
	double tri_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2);

	/** \brief ���������ε�Ԫ������
	    \param v0~v2 ��������������
	*/
	double tri_quality(Point_3& p0, Point_3&  p1, Point_3&  p2);

	/** \breif �������ĵ�����ɵ��ıߵ�����ϵ��(0~1)
	    \param �ı��ε��ĸ�����, ����ʱ��˳�����
	*/
	double quad_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
		                        HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3);

	/** \breif �Ľ����ı�������ϵ�����㹫ʽ
	*   ���ڻ����Ż��Ĺ�˳����
	    \param �ı��ε��ĸ�����, ����ʱ��˳�����
	*/
	double quad_quality2(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
		                         HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3);

	/** \breif �Ľ����ı�������ϵ�����㹫ʽ
	*   ���ڻ����Ż��Ĺ�˳����
	    \param �ı��ε��ĸ�����, ����ʱ��˳�����
	*/
	double quad_quality2(Point_3& p0, Point_3& p1, 
		                         Point_3& p2, Point_3& p3);

	/** \brief ��ǰ�ص����Laplacian��˳
	* localSmoothing()�ĸ�������
	    \param vh ����˳�ڵ�
	*/
	void smoothFrontNode(HDS_mesh::Vertex_handle& vh);

	int smoothInteriorNode(HDS_mesh::Vertex_handle& vh);

	int smoothAdjustment(HDS_mesh::Vertex_handle& vh, Point_3& new_point);

	bool localReconnect();

private:

	/** \brief �ж���������ɵ���������Ƭ�Ƿ�Ϸ�
	    \return true �Ϸ� ���������ʱ������
		\return false ���Ϸ� �������˳ʱ������
	*/
	bool isValid(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2);

	/** brief ʹ�ǶȽ���120~240֮��ı߽��Ķ�Ϊ2
	*/
	bool boundary_unform();

	/** brief ���һ�������ڽӵ�Ԫ����С��4��ֻ�������ڽ�������,�����ִ�з��Ѳ���
	*/
	bool split_vertex();

	/** brief ���һ���ı��ε����Ա߾������������ڣ��������
	* ���ı��η��Ѳ��������������ηֱ���ԭ�������������κϲ����������µ��ı���
	*/
	bool clean_tri_zj();

	/** \brief ִ���ı��������������
	* ���ڽ��������ڽӵ�Ԫ�Ľڵ�, ɾ���õ㲢���������ڽӵ�Ԫ�ϲ���һ���ı��ε�Ԫ
	*/
	bool quad_mesh_clean();

	/** \brief ִ���ı���������Ƭ�������
	* ���һ���ı��ε�Ԫ���������Զ�����ڽӵ�Ԫ����Ϊ��, ��������ϲ�
	*/
	bool quad_face_clean();

	/** \brief quad_face_clean()��������
	* �ϲ�e_base->vertex()��e_base->next()->next()->vertex()����
	*/
	bool face_clean_aux(HDS_mesh::Halfedge_handle& e_base);

	/** \brief ִ���ı���������Ƭ�������
	* �����������ڽڵ�ȶ�Ϊ�������
	*/
	bool quad_face_clean_aux1();

	/** \brief ִ���ı�������ڵ��������
	* ���ڽ��������ڽӵ�Ԫ�Ľڵ�, ɾ���õ㲢���������ڽӵ�Ԫ�ϲ���һ���ı��ε�Ԫ
	*/
	bool quad_bivalent_vertex_clean();

	/** \brief ִ���ı���������Ƭ�������
	* �����ڽǽӽ�180�ȵ��ı��ε�Ԫ���������ڽӵ�Ԫ�ϲ�����ù̶�ģ��������
	*/
	bool quad_face_clean_aux2();

private:
	HDS_mesh* _source_mesh;
};

class Vertex_Comp
{
public:  
	bool operator()(const HDS_mesh::Vertex_handle& v1, const HDS_mesh::Vertex_handle& v2)  
	{          
		return v1->getIndex() < v2->getIndex();     
	}  
};

/** \brief ���������������qp��qr�ļн�
	\return ����ֵ��λΪ��
*/
double angle_evaluate(const Point_3 &p, const Point_3 &q, const Point_3 &r);

/** \brief ��������A������B�нǵ�����ֵ
    \param p A���������
    \param q A�������յ�
    \param r B���������
    \param s B�������յ�
*/
double angle_evaluate2(const Point_3 &p, const Point_3 &q, 
							  const Point_3 &r, const Point_3 &s);

Point_3 operator+ (Point_3&p, Point_3& q);
Point_3 operator- (Point_3&p);
Point_3 operator* (double t, Point_3&p);
Point_3 operator* (Point_3&p, double t);
Point_3 operator/ (Point_3&p, double t);
#endif