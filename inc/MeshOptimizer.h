//**************************************************************************
// File:     MeshOptimizer.h
// Author:   Zhang Jun
// Date:     2016-11-10
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 负责网格优化
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

	bool mesh_OPT(); ///< 实现该接口, 对网格进行优化

	/** \brief 对网格进行全局光顺
	    \return n 光顺次数
	*/
	bool global_smooth(int n);

	/** \brief 局部光顺
	* 每次新产生单元，都要进行局部光顺\n
	* 光顺对象为新单元的四个顶点及其第一级邻接顶点\n
	* 光顺算法为基于长度的Laplacian光顺
		\param[int] smooth_node_list 被执行光顺的节点集
	*/
	void localSmoothing(std::list<HDS_mesh::Vertex_handle>& smooth_node_list);

	/** \brief 计算三角形单元的质量
	    \param v0~v2 三角形三个顶点
	*/
	double tri_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2);

	/** \brief 计算三角形单元的质量
	    \param v0~v2 三角形三个顶点
	*/
	double tri_quality(Point_3& p0, Point_3&  p1, Point_3&  p2);

	/** \breif 计算有四点所组成的四边的质量系数(0~1)
	    \param 四边形的四个顶点, 按逆时针顺序给出
	*/
	double quad_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
		                        HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3);

	/** \breif 改进的四边形质量系数计算公式
	*   用于基于优化的光顺操作
	    \param 四边形的四个顶点, 按逆时针顺序给出
	*/
	double quad_quality2(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
		                         HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3);

	/** \breif 改进的四边形质量系数计算公式
	*   用于基于优化的光顺操作
	    \param 四边形的四个顶点, 按逆时针顺序给出
	*/
	double quad_quality2(Point_3& p0, Point_3& p1, 
		                         Point_3& p2, Point_3& p3);

	/** \brief 对前沿点进行Laplacian光顺
	* localSmoothing()的辅助函数
	    \param vh 待光顺节点
	*/
	void smoothFrontNode(HDS_mesh::Vertex_handle& vh);

	int smoothInteriorNode(HDS_mesh::Vertex_handle& vh);

	int smoothAdjustment(HDS_mesh::Vertex_handle& vh, Point_3& new_point);

	bool localReconnect();

private:

	/** \brief 判断由三点组成的三角形面片是否合法
	    \return true 合法 即三点成逆时针走向
		\return false 不合法 即三点成顺时针走向
	*/
	bool isValid(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2);

	/** brief 使角度介于120~240之间的边界点的度为2
	*/
	bool boundary_unform();

	/** brief 如果一个顶点邻接单元数不小于4且只有两个邻接三角形,则对其执行分裂操作
	*/
	bool split_vertex();

	/** brief 如果一个四边形的两对边均与三角形相邻，则将其分裂
	* 用四边形分裂产生的两个三角形分别与原来的两个三角形合并生成两个新的四边形
	*/
	bool clean_tri_zj();

	/** \brief 执行四边形网格清理操作
	* 对于仅有两个邻接单元的节点, 删除该点并将其两个邻接单元合并成一个四边形单元
	*/
	bool quad_mesh_clean();

	/** \brief 执行四边形网格面片清理操作
	* 如果一个四边形单元中有两个对顶点的邻接单元数均为三, 则将这两点合并
	*/
	bool quad_face_clean();

	/** \brief quad_face_clean()辅助函数
	* 合并e_base->vertex()和e_base->next()->next()->vertex()两点
	*/
	bool face_clean_aux(HDS_mesh::Halfedge_handle& e_base);

	/** \brief 执行四边形网格面片清理操作
	* 负责处理两相邻节点度都为三的情况
	*/
	bool quad_face_clean_aux1();

	/** \brief 执行四边形网格节点清理操作
	* 对于仅有两个邻接单元的节点, 删除该点并将其两个邻接单元合并成一个四边形单元
	*/
	bool quad_bivalent_vertex_clean();

	/** \brief 执行四边形网格面片清理操作
	* 对于内角接近180度的四边形单元，将其与邻接单元合并后采用固定模板进行填充
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

/** \brief 计算两共点的向量qp和qr的夹角
	\return 返回值单位为度
*/
double angle_evaluate(const Point_3 &p, const Point_3 &q, const Point_3 &r);

/** \brief 计算向量A与向量B夹角的余弦值
    \param p A向量的起点
    \param q A向量的终点
    \param r B向量的起点
    \param s B向量的终点
*/
double angle_evaluate2(const Point_3 &p, const Point_3 &q, 
							  const Point_3 &r, const Point_3 &s);

Point_3 operator+ (Point_3&p, Point_3& q);
Point_3 operator- (Point_3&p);
Point_3 operator* (double t, Point_3&p);
Point_3 operator* (Point_3&p, double t);
Point_3 operator/ (Point_3&p, double t);
#endif