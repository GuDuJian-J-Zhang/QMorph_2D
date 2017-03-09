//**************************************************************************
// File:     QMorphGenerator.h
// Author:   Zhang Jun
// Date:     2016-11-09
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 四边形网格生成算法实现(基于QMorph的算法)
//**************************************************************************
#ifndef QMORPH_GENERATOR_H
#define QMORPH_GENERATOR_H

#include "QuadMeshGenerator.h"
class MeshOptimizer;

const double state_tolerence = 135; ///< 角度值(单位：度)
const double side_tolerence = 30; ///< 角度值(单位：度)

const double PAVING_SIDE = 240; ///< Paving算法中边节点处前沿边夹角的上限
const double PAVINE_CORNER = 300; ///< Paving算法中角节点处前沿边夹角的上限
const double SEAM_ANGLE = 45; ///< 缝合角(当相邻前沿夹角小于该值时执行前沿缝合操作)

class QMorphGenerator: public QuadMeshGenerator
{
public:
	QMorphGenerator();
	virtual void mesh();
	virtual void set_SourceMesh(HDS_mesh *_source_mesh);
	virtual ~QMorphGenerator();

private:
	void classify_Front_Edge(); ///< 根据前沿边状态将前沿边分为四类

	void set_Front_Node_State(HDS_mesh::Vertex_handle vh); ///< 设置前沿点的状态值

	/** \brief 更新前沿
		\param influenced_front_edge_list 由于光顺而受影响的边链表
	*/
	void updateFront(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list);

	/** \brief 获取下一条前沿边
	* 在必要时更新前沿
	*/
	bool get_Next_Front(HDS_mesh::Halfedge_handle& eh);

	/** \brief 直接从四个前沿边链表中获取下一条前沿边
	    \return True 获取成功
		\return False 获取失败
	*/
	bool get_Front_Aux(HDS_mesh::Halfedge_handle& eh);

	/** \breif 获取左侧边界边
	    \param front 当前前沿
	    \param eh 待获取的半边
	    \return true 获取成功
		\return false 获取失败
	*/
	bool get_Side_Left(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh);

	/** \breif 获取右侧边界边
	    \param front 当前前沿
	    \param eh 待获取的半边
	    \return true 获取成功
		\return false 获取失败
	*/
	bool get_Side_Right(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh);

	/** \brief 获取顶边
	    \param vl 顶边左端点的手柄
		\param vr 顶边右端点的手柄
		\return 返回顶边的手柄
	*/
	HDS_mesh::Halfedge_handle  get_Top(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr);
	
	/** \breif 使用四条边创建四边形单元
	*   新四边形单元的关联半边为当前前沿边
	    \param e_down 四边形下侧边(即当前前沿边)
		\param e_right 四边形右侧边
		\param e_up 四边形上侧边
		\param e_left 四边形左侧边
	*/
	HDS_mesh::Face_handle constructNewFace(HDS_mesh::Halfedge_handle e_down, HDS_mesh::Halfedge_handle e_right, 
		                                   HDS_mesh::Halfedge_handle e_up,   HDS_mesh::Halfedge_handle e_left);

	/** \breif 从face_base开始采用广度优先搜索来收集所有待删除三角面片
	*   辅助函数
	    \param[in] face_base 从该面片开始执行广度优先搜索
		\param[out] face_list 待删除面片链表
		\param[out] old_edge_list 待删除边链表
		\param[out] old_vertex_list 待删除节点链表
	*/
	void getFaceToDelete(HDS_mesh::Face_handle face_base, 
		                 std::list<HDS_mesh::Face_handle>& old_face_list,
		                 std::list<HDS_mesh::Halfedge_handle>& old_edge_list,
		                 std::list<HDS_mesh::Vertex_handle>& old_vertex_list);

	/** \brief 获取待光顺的节点
	* localSmoothing()的辅助函数
	    \param quad_face 新产生的四边形面片
	*/
	void getVertexToSmooth(HDS_mesh::Face_handle quad_face, 
		                   std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list);

	/** \brief 获取由于光顺而受到影响的前沿边
	    \param[in] vertex_smooth_list 被光顺的点
		\param[out] influenced_front_edge_list 受影响的前沿边
	*/
	void getInfluencedFrontEdge(std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list,
		                        std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list);

	/** \brief 将前沿边从相应前沿链表中删除
	    \param ed 待删除半边的手柄
	*/
	void deleteFromFront(HDS_mesh::Halfedge_handle eh);

	/** \brief 将前沿边从非激活前沿链表中删除
	    \param ed 待删除半边的手柄
	*/
	void deleteFromUnactiveFront(HDS_mesh::Halfedge_handle eh);

	/** \brief 判断边交换前后，最大角是否减小
	*/
	bool MinMaxSwap2D(Point_3& pt0,Point_3& pt1,Point_3& pt2,Point_3& pt3);

	/** \brief 执行边翻转操作
	    \param eh 待分裂边
		\return 返回新边的手柄
	*/
    bool flip_edge(HDS_mesh::Halfedge_handle& eh, bool check_angle = false);

	/** \brief 恢复点vl和点vr之间的边
	    \param[in] vl 待恢复边的做端点
		\param[in] vr 待恢复边的右端点
		\param[out] edge 待恢复的边
		\return false 恢复失败
		\return true 恢复成功
	*/
	bool recover_edge(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr, HDS_mesh::Halfedge_handle& edge);

	/** \brief 折叠边
	* 即把边的两个端点合并
	    \param t 参数t指定新顶点在原边上的位置(0:表示在左端点，1:表示在右端点)
	*/
	bool merge_edge(HDS_mesh::Halfedge_handle& edge, double t);

	/** \brief 局部光顺
	* 每次新产生单元，都要进行局部光顺\n
	* 光顺对象为新单元的四个顶点及其第一级邻接顶点\n
	* 光顺算法为基于长度的Laplacian光顺
		\param[int] smooth_node_list 被执行光顺的节点集
	*/
	void localSmoothing(std::list<HDS_mesh::Vertex_handle>& smooth_node_list);

	/** \brief 对前沿点进行Laplacian光顺
	* localSmoothing()的辅助函数
	    \param vh 待光顺节点
	*/
	void smoothFrontNode(HDS_mesh::Vertex_handle& vh);

	/** \brief 基于优化的光顺
	    \param vh 待光顺节点
	*/
	bool OptBasedSmooth(HDS_mesh::Vertex_handle& vh);

	/** \brief 把前沿重新插回到相应前沿队列的尾部
	*/
	void rePushFront(HDS_mesh::Halfedge_handle eh);

	bool has_node(HDS_mesh::Face_const_handle tri, HDS_mesh::Vertex_handle vh);

	/** \brief 环上前沿数大于6时返回true，否则返回false
	* 
	    \param e_left 待缝合前沿边
		\param e_right 待缝合前沿边
	*/
	bool canSeam(HDS_mesh::Halfedge_handle& e_left);

	/** \brief 执行前沿缝合操作
	* 当相邻两前沿夹角小于给定值时执行前沿合并操作
	    \param e_left 待缝合前沿边
		\param e_right 待缝合前沿边
		\return 返回缝合点的入射前沿边手柄
	*/
	HDS_mesh::Halfedge_handle seamFront(HDS_mesh::Halfedge_handle& e_left, 
		                                HDS_mesh::Halfedge_handle& e_right, 
		                                bool is_top = false);

	/** \brief 执行前沿缝合操作
	* seamFront()辅助函数。用于较长待合并边与较短待合并边长之比小于等于2.5的情况
	*/
	HDS_mesh::Halfedge_handle seamFrontNormal(HDS_mesh::Halfedge_handle& e_left,
		                                      HDS_mesh::Halfedge_handle& e_right,
		                                      bool is_top = false);

	/** \brief 执行前沿缝合操作
	* seamFront()辅助函数。用于较长待合并边与较短待合并边长之比大于2.5的情况
	*/
	HDS_mesh::Halfedge_handle seamFrontTrans(HDS_mesh::Halfedge_handle& e_biger, 
		                                     HDS_mesh::Halfedge_handle& e_shorter, 
		                                     bool is_top = false);

	/** \brief 检查是否有前沿需要缝合
	* 对每条待检查前沿，判断它与前驱前沿及后继前沿之间的夹角，根据角度值及前沿点的化合价决定是否缝合\n
	* 一旦执行缝合操作，则将由于缝合操作而受到影响边收集起来
	    \param[in] influenced_front_edge_list 待检查的前沿
		\param[out] node_need_smooth 由于缝合操作而受影响的节点
		\return true 执行过缝合操作
		\return false 未执行过缝合操作 (此时不需要后续光顺操作)
	*/
	bool checkFrontToSeam(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list,
		                  std::list<HDS_mesh::Vertex_handle>& node_need_smooth,
		                  int n);

	/** \brief 执行三角网格清理操作
	* 对于仅有三个邻接单元的节点, 删除该点并将其三个邻接单元合并成一个三角形单元
	*/
	bool tri_mesh_clean();
private:
	HDS_mesh* _source_mesh;
	MeshOptimizer* _optimizer;

	std::list<HDS_mesh::Halfedge_handle> _active_front_list; ///< 活动前沿边链表
	std::list<HDS_mesh::Halfedge_handle> _unactive_front_list; ///< 非活动前沿边链表

	std::list<HDS_mesh::Halfedge_handle> _front_edge_11; ///< 状态值为11的前沿边链表
	std::list<HDS_mesh::Halfedge_handle> _front_edge_10; ///< 状态值为10的前沿边链表
	std::list<HDS_mesh::Halfedge_handle> _front_edge_01; ///< 状态值为01的前沿边链表
	std::list<HDS_mesh::Halfedge_handle> _front_edge_00; ///< 状态值为00的前沿边链表

	int _max_index;
	double _reatio; ///< 原始边界最大尺寸与最小尺寸的比值(用于前沿点光顺)
	double _deta; ///< 用于基于优化的光顺操作
};
#endif 