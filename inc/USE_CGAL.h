//**************************************************************************
// File:     USE_CGAL.h
// Author:   Zhang Jun
// Date:     2016-06-19
// Email:    zhangjun_dg@mail.dlut.edu.cn
//**************************************************************************

#ifndef _USE_CGAL_H
#define _USE_CGAL_H

#define PI 3.141592653589793238462643
#define TOLERENCE 1.0e-6

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include "foundation.h"

///< 半边数据结构相关
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

///< 仿射变换相关
#include <CGAL/Aff_transformation_3.h>

///< 相交计算相关
#include <CGAL/intersections.h>

///< 通用计算相关
#include <CGAL/Cartesian/function_objects.h>

template < class Refs, class P>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>, public Marker
{
public:
	enum State{BOUNDARY = 1, ///< 边界节点
			   VISITED,      ///< 访问标记
			   STATE,        ///< 前沿点状态(只用于前沿点的标记)
			   OBTUSE,       ///< 标记共享该点的两前沿形成钝角
			   FrontNode,    ///< 标记为前沿点
			   QuadNode,     ///< 标记节点为四边形节点
			   HASCHOICED,   ///< 防止了四边形两侧边共点
			   LOOPCUT,      ///< 环切割点(用于环分裂操作)
			   SEAMCHECK,    ///< 用于前沿缝合
			   ERASED,
			   SIDENODE,    ///< Paving方法边节点
			   LOCKED       ///< 被锁定的节点不能被光顺
		};
public:
	My_vertex(): CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>() {}
	My_vertex(const P& pp): CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>(pp){}
	public:
		using Marker::setDead;
		using Marker::isDead;
		/** \brief 检查点是否被做过某种标记
		*
			\return true  是
			\return false 否
		*/
		using Marker::testMark;

		/** \brief 给点做某种标记
		*
			即将相应标记位置1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< 返回点的索引号
		using Marker::setIndex;
		void         set_angle(double t)      { _angle = t; } ///< 设置该点处两个前沿边之间的夹角
		double       get_angle()              { return _angle; }
	private:
		//int     some_additional_index;
		double  _angle; ///< 该点处两个前沿边之间的夹角
};

template <class Refs>
class My_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>, public Marker
{
public:
    typedef typename Refs::Halfedge_handle                Halfedge_handle;
	typedef typename CGAL::HalfedgeDS_halfedge_base<Refs>::Base_base  Base_base;
	public:
	enum State{BOUNDARY = 1, ///< 几何边界
			   VISITED,      ///< 访问标记
			   FRONT,        ///< 前沿标记
			   QuadEdge,     ///< 四边形边
			   UNACTIVE,     ///< 未激活前沿
			   LOOPBEGIN,    ///< 前沿环起始边
			   ISSPLIT,     ///< 用于分裂缝合
			   ISDOWN,      ///< 用于三角形合并
	          };
	public:
		using Marker::setDead;
		using Marker::isDead;

		/** \brief 检查点是否被做过某种标记
		*
			\return true  是
			\return false 否
		*/
		using Marker::testMark;

		/** \brief 给点做某种标记
		*
			即将相应标记位置1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< 返回点的索引号
		using Marker::setIndex;
	public:
		///< 前沿推进相关
		Halfedge_handle           loopNext()                               { return _loop_nxt;}
		void                      set_loopNext(Halfedge_handle h) { _loop_nxt = h;} ///< 设置所在环后继半边

		Halfedge_handle           loopPrev()                      { return _loop_prv; }
		void                      set_loopPrev(Halfedge_handle h) { _loop_prv = h;    } ///< 设置所在环前驱半边

		void                      set_levelNum(int n)             { _level_num = n; } ///< 设置前沿所在层数
		int                       get_levelNum()                  { return _level_num; }

		void                      set_FrontState(int t)           { _state = t; } ///< 设置前沿边的状态值
		int                       get_FrontState()                { return _state; }

		void                      set_Key(int t)                  { _key = t; }
		int                       get_Key()                       { return _key; }

		void                      set_SquardLength(double t)      { _squard_length = t; }
		double                    get_SquardLength()              { return _squard_length; }

		void                      set_opposite( Halfedge_handle h)  { Base_base::set_opposite(h);}
	private:
		Halfedge_handle  _loop_nxt; ///< 所在环后继半边
		Halfedge_handle  _loop_prv; ///< 所在环前驱半边

		double _squard_length;///< 边长,用于前沿排序
		int _level_num; ///< 前沿所处的层数
		int _state; ///< 前沿的状态值
		int _key; ///< 边的键值(用于快速查找)
};

template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs>, public Marker
{
	enum State{
		       VISITED = 1,   ///< 访问标记
			   TRI,           ///< 三角形
			   QUAD,          ///< 四边形
	          };
		using Marker::setDead;
		using Marker::isDead;
		/** \brief 检查点是否被做过某种标记
		*
			\return true  是
			\return false 否
		*/
		using Marker::testMark;

		/** \brief 给点做某种标记
		*
			即将相应标记位置1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< 返回点的索引号
		using Marker::setIndex;
};

struct My_items : public CGAL::Polyhedron_items_3
{
	template < class Refs, class Traits>
	struct Vertex_wrapper
	{
		typedef typename Traits::Point_3 Point;
		typedef My_vertex<Refs, Point> Vertex;
	};

	template < class Refs, class Traits>
	struct Halfedge_wrapper
	{
		typedef My_halfedge<Refs> Halfedge;
	};

	template <class Refs, class Traits>
	struct Face_wrapper
	{
		typedef My_face<Refs> Face;
	};
};
typedef CGAL::Simple_cartesian<double>                     MeshKernel;
typedef CGAL::Polyhedron_3<MeshKernel, My_items>           HDS_mesh;
typedef HDS_mesh::HalfedgeDS                               HalfedgeDS;
//typedef HDS_mesh::Halfedge_handle                          Halfedge_handle;
typedef MeshKernel::Point_2                                Point_2;
typedef MeshKernel::Point_3                                Point_3;
typedef MeshKernel::Vector_3                               Vector_3;
typedef MeshKernel::Line_3                                 Line_3;
typedef CGAL::Aff_transformation_3<MeshKernel>             Transformation_3;
typedef CGAL::Direction_3<MeshKernel>                      Direction_3;
typedef MeshKernel::Compute_scalar_product_3               Compute_scalar_product_3;
///< 用于二维相交判断
// #include <CGAL/Arr_segment_traits_2.h>
// #include <CGAL/Sweep_line_2_algorithms.h>
// #include <CGAL/intersections.h>
//
// #include <CGAL/intersections.h>

///< 用于二维相交判断
// typedef CGAL::Exact_predicates_exact_constructions_kernel                                   Kernel;
// typedef CGAL::Arr_segment_traits_2<Kernel>                                                  Intersect_Traits_2;
// typedef Intersect_Traits_2::Curve_2                                                         Intersect_Segement_2;
// typedef Kernel::Segment_2                                                                   Segment_2;
// typedef Kernel::Point_2                                                                     Point_2;
// typedef Kernel::Vector_2                                                                    Vector_2;
// typedef Kernel::Line_2                                                                      Line_2;
// typedef Kernel::Intersect_2                                                                 Intersect_2;

#endif
