//**************************************************************************
// File:     QMorphGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-09
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: �ı������������㷨ʵ��(����QMorph���㷨)
//**************************************************************************
#ifndef QMORPH_GENERATOR_H
#define QMORPH_GENERATOR_H

#include "QuadMeshGenerator.h"
class MeshOptimizer;

const double state_tolerence = 135; ///< �Ƕ�ֵ(��λ����)
const double side_tolerence = 30; ///< �Ƕ�ֵ(��λ����)

const double PAVING_SIDE = 240; ///< Paving�㷨�б߽ڵ㴦ǰ�ر߼нǵ�����
const double PAVINE_CORNER = 300; ///< Paving�㷨�нǽڵ㴦ǰ�ر߼нǵ�����
const double SEAM_ANGLE = 45; ///< ��Ͻ�(������ǰ�ؼн�С�ڸ�ֵʱִ��ǰ�ط�ϲ���)

class QMorphGenerator: public QuadMeshGenerator
{
public:
	QMorphGenerator();
	virtual void mesh();
	virtual void set_SourceMesh(HDS_mesh *_source_mesh);
	virtual ~QMorphGenerator();

private:
	void classify_Front_Edge(); ///< ����ǰ�ر�״̬��ǰ�ر߷�Ϊ����

	void set_Front_Node_State(HDS_mesh::Vertex_handle vh); ///< ����ǰ�ص��״ֵ̬

	/** \brief ����ǰ��
		\param influenced_front_edge_list ���ڹ�˳����Ӱ��ı�����
	*/
	void updateFront(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list);

	/** \brief ��ȡ��һ��ǰ�ر�
	* �ڱ�Ҫʱ����ǰ��
	*/
	bool get_Next_Front(HDS_mesh::Halfedge_handle& eh);

	/** \brief ֱ�Ӵ��ĸ�ǰ�ر������л�ȡ��һ��ǰ�ر�
	    \return True ��ȡ�ɹ�
		\return False ��ȡʧ��
	*/
	bool get_Front_Aux(HDS_mesh::Halfedge_handle& eh);

	/** \breif ��ȡ���߽��
	    \param front ��ǰǰ��
	    \param eh ����ȡ�İ��
	    \return true ��ȡ�ɹ�
		\return false ��ȡʧ��
	*/
	bool get_Side_Left(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh);

	/** \breif ��ȡ�Ҳ�߽��
	    \param front ��ǰǰ��
	    \param eh ����ȡ�İ��
	    \return true ��ȡ�ɹ�
		\return false ��ȡʧ��
	*/
	bool get_Side_Right(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh);

	/** \brief ��ȡ����
	    \param vl ������˵���ֱ�
		\param vr �����Ҷ˵���ֱ�
		\return ���ض��ߵ��ֱ�
	*/
	HDS_mesh::Halfedge_handle  get_Top(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr);
	
	/** \breif ʹ�������ߴ����ı��ε�Ԫ
	*   ���ı��ε�Ԫ�Ĺ������Ϊ��ǰǰ�ر�
	    \param e_down �ı����²��(����ǰǰ�ر�)
		\param e_right �ı����Ҳ��
		\param e_up �ı����ϲ��
		\param e_left �ı�������
	*/
	HDS_mesh::Face_handle constructNewFace(HDS_mesh::Halfedge_handle e_down, HDS_mesh::Halfedge_handle e_right, 
		                                   HDS_mesh::Halfedge_handle e_up,   HDS_mesh::Halfedge_handle e_left);

	/** \breif ��face_base��ʼ���ù�������������ռ����д�ɾ��������Ƭ
	*   ��������
	    \param[in] face_base �Ӹ���Ƭ��ʼִ�й����������
		\param[out] face_list ��ɾ����Ƭ����
		\param[out] old_edge_list ��ɾ��������
		\param[out] old_vertex_list ��ɾ���ڵ�����
	*/
	void getFaceToDelete(HDS_mesh::Face_handle face_base, 
		                 std::list<HDS_mesh::Face_handle>& old_face_list,
		                 std::list<HDS_mesh::Halfedge_handle>& old_edge_list,
		                 std::list<HDS_mesh::Vertex_handle>& old_vertex_list);

	/** \brief ��ȡ����˳�Ľڵ�
	* localSmoothing()�ĸ�������
	    \param quad_face �²������ı�����Ƭ
	*/
	void getVertexToSmooth(HDS_mesh::Face_handle quad_face, 
		                   std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list);

	/** \brief ��ȡ���ڹ�˳���ܵ�Ӱ���ǰ�ر�
	    \param[in] vertex_smooth_list ����˳�ĵ�
		\param[out] influenced_front_edge_list ��Ӱ���ǰ�ر�
	*/
	void getInfluencedFrontEdge(std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list,
		                        std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list);

	/** \brief ��ǰ�رߴ���Ӧǰ��������ɾ��
	    \param ed ��ɾ����ߵ��ֱ�
	*/
	void deleteFromFront(HDS_mesh::Halfedge_handle eh);

	/** \brief ��ǰ�رߴӷǼ���ǰ��������ɾ��
	    \param ed ��ɾ����ߵ��ֱ�
	*/
	void deleteFromUnactiveFront(HDS_mesh::Halfedge_handle eh);

	/** \brief �жϱ߽���ǰ�������Ƿ��С
	*/
	bool MinMaxSwap2D(Point_3& pt0,Point_3& pt1,Point_3& pt2,Point_3& pt3);

	/** \brief ִ�б߷�ת����
	    \param eh �����ѱ�
		\return �����±ߵ��ֱ�
	*/
    bool flip_edge(HDS_mesh::Halfedge_handle& eh, bool check_angle = false);

	/** \brief �ָ���vl�͵�vr֮��ı�
	    \param[in] vl ���ָ��ߵ����˵�
		\param[in] vr ���ָ��ߵ��Ҷ˵�
		\param[out] edge ���ָ��ı�
		\return false �ָ�ʧ��
		\return true �ָ��ɹ�
	*/
	bool recover_edge(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr, HDS_mesh::Halfedge_handle& edge);

	/** \brief �۵���
	* ���ѱߵ������˵�ϲ�
	    \param t ����tָ���¶�����ԭ���ϵ�λ��(0:��ʾ����˵㣬1:��ʾ���Ҷ˵�)
	*/
	bool merge_edge(HDS_mesh::Halfedge_handle& edge, double t);

	/** \brief �ֲ���˳
	* ÿ���²�����Ԫ����Ҫ���оֲ���˳\n
	* ��˳����Ϊ�µ�Ԫ���ĸ����㼰���һ���ڽӶ���\n
	* ��˳�㷨Ϊ���ڳ��ȵ�Laplacian��˳
		\param[int] smooth_node_list ��ִ�й�˳�Ľڵ㼯
	*/
	void localSmoothing(std::list<HDS_mesh::Vertex_handle>& smooth_node_list);

	/** \brief ��ǰ�ص����Laplacian��˳
	* localSmoothing()�ĸ�������
	    \param vh ����˳�ڵ�
	*/
	void smoothFrontNode(HDS_mesh::Vertex_handle& vh);

	/** \brief �����Ż��Ĺ�˳
	    \param vh ����˳�ڵ�
	*/
	bool OptBasedSmooth(HDS_mesh::Vertex_handle& vh);

	/** \brief ��ǰ�����²�ص���Ӧǰ�ض��е�β��
	*/
	void rePushFront(HDS_mesh::Halfedge_handle eh);

	bool has_node(HDS_mesh::Face_const_handle tri, HDS_mesh::Vertex_handle vh);

	/** \brief ����ǰ��������6ʱ����true�����򷵻�false
	* 
	    \param e_left �����ǰ�ر�
		\param e_right �����ǰ�ر�
	*/
	bool canSeam(HDS_mesh::Halfedge_handle& e_left);

	/** \brief ִ��ǰ�ط�ϲ���
	* ��������ǰ�ؼн�С�ڸ���ֵʱִ��ǰ�غϲ�����
	    \param e_left �����ǰ�ر�
		\param e_right �����ǰ�ر�
		\return ���ط�ϵ������ǰ�ر��ֱ�
	*/
	HDS_mesh::Halfedge_handle seamFront(HDS_mesh::Halfedge_handle& e_left, 
		                                HDS_mesh::Halfedge_handle& e_right, 
		                                bool is_top = false);

	/** \brief ִ��ǰ�ط�ϲ���
	* seamFront()�������������ڽϳ����ϲ�����϶̴��ϲ��߳�֮��С�ڵ���2.5�����
	*/
	HDS_mesh::Halfedge_handle seamFrontNormal(HDS_mesh::Halfedge_handle& e_left,
		                                      HDS_mesh::Halfedge_handle& e_right,
		                                      bool is_top = false);

	/** \brief ִ��ǰ�ط�ϲ���
	* seamFront()�������������ڽϳ����ϲ�����϶̴��ϲ��߳�֮�ȴ���2.5�����
	*/
	HDS_mesh::Halfedge_handle seamFrontTrans(HDS_mesh::Halfedge_handle& e_biger, 
		                                     HDS_mesh::Halfedge_handle& e_shorter, 
		                                     bool is_top = false);

	/** \brief ����Ƿ���ǰ����Ҫ���
	* ��ÿ�������ǰ�أ��ж�����ǰ��ǰ�ؼ����ǰ��֮��ļнǣ����ݽǶ�ֵ��ǰ�ص�Ļ��ϼ۾����Ƿ���\n
	* һ��ִ�з�ϲ����������ڷ�ϲ������ܵ�Ӱ����ռ�����
	    \param[in] influenced_front_edge_list ������ǰ��
		\param[out] node_need_smooth ���ڷ�ϲ�������Ӱ��Ľڵ�
		\return true ִ�й���ϲ���
		\return false δִ�й���ϲ��� (��ʱ����Ҫ������˳����)
	*/
	bool checkFrontToSeam(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list,
		                  std::list<HDS_mesh::Vertex_handle>& node_need_smooth,
		                  int n);

	/** \brief ִ�����������������
	* ���ڽ��������ڽӵ�Ԫ�Ľڵ�, ɾ���õ㲢���������ڽӵ�Ԫ�ϲ���һ�������ε�Ԫ
	*/
	bool tri_mesh_clean();
private:
	HDS_mesh* _source_mesh;
	MeshOptimizer* _optimizer;

	std::list<HDS_mesh::Halfedge_handle> _active_front_list; ///< �ǰ�ر�����
	std::list<HDS_mesh::Halfedge_handle> _unactive_front_list; ///< �ǻǰ�ر�����

	std::list<HDS_mesh::Halfedge_handle> _front_edge_11; ///< ״ֵ̬Ϊ11��ǰ�ر�����
	std::list<HDS_mesh::Halfedge_handle> _front_edge_10; ///< ״ֵ̬Ϊ10��ǰ�ر�����
	std::list<HDS_mesh::Halfedge_handle> _front_edge_01; ///< ״ֵ̬Ϊ01��ǰ�ر�����
	std::list<HDS_mesh::Halfedge_handle> _front_edge_00; ///< ״ֵ̬Ϊ00��ǰ�ر�����

	int _max_index;
	double _reatio; ///< ԭʼ�߽����ߴ�����С�ߴ�ı�ֵ(����ǰ�ص��˳)
	double _deta; ///< ���ڻ����Ż��Ĺ�˳����
};
#endif 