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

///< ������ݽṹ���
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

///< ����任���
#include <CGAL/Aff_transformation_3.h>

///< �ཻ�������
#include <CGAL/intersections.h>

///< ͨ�ü������
#include <CGAL/Cartesian/function_objects.h>

template < class Refs, class P>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>, public Marker
{
public:
	enum State{BOUNDARY = 1, ///< �߽�ڵ�
			   VISITED,      ///< ���ʱ��
			   STATE,        ///< ǰ�ص�״̬(ֻ����ǰ�ص�ı��)
			   OBTUSE,       ///< ��ǹ���õ����ǰ���γɶ۽�
			   FrontNode,    ///< ���Ϊǰ�ص�
			   QuadNode,     ///< ��ǽڵ�Ϊ�ı��νڵ�
			   HASCHOICED,   ///< ��ֹ���ı�������߹���
			   LOOPCUT,      ///< ���и��(���ڻ����Ѳ���)
			   SEAMCHECK,    ///< ����ǰ�ط��
			   ERASED,
			   SIDENODE,    ///< Paving�����߽ڵ�
			   LOCKED       ///< �������Ľڵ㲻�ܱ���˳
		};
public:
	My_vertex(): CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>() {}
	My_vertex(const P& pp): CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, P>(pp){}
	public:
		using Marker::setDead;
		using Marker::isDead;
		/** \brief �����Ƿ�����ĳ�ֱ��
		*
			\return true  ��
			\return false ��
		*/
		using Marker::testMark;

		/** \brief ������ĳ�ֱ��
		*
			������Ӧ���λ��1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< ���ص��������
		using Marker::setIndex;
		void         set_angle(double t)      { _angle = t; } ///< ���øõ㴦����ǰ�ر�֮��ļн�
		double       get_angle()              { return _angle; }
	private:
		//int     some_additional_index;
		double  _angle; ///< �õ㴦����ǰ�ر�֮��ļн�
};

template <class Refs>
class My_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>, public Marker
{
public:
    typedef typename Refs::Halfedge_handle                Halfedge_handle;
	typedef typename CGAL::HalfedgeDS_halfedge_base<Refs>::Base_base  Base_base;
	public:
	enum State{BOUNDARY = 1, ///< ���α߽�
			   VISITED,      ///< ���ʱ��
			   FRONT,        ///< ǰ�ر��
			   QuadEdge,     ///< �ı��α�
			   UNACTIVE,     ///< δ����ǰ��
			   LOOPBEGIN,    ///< ǰ�ػ���ʼ��
			   ISSPLIT,     ///< ���ڷ��ѷ��
			   ISDOWN,      ///< ���������κϲ�
	          };
	public:
		using Marker::setDead;
		using Marker::isDead;

		/** \brief �����Ƿ�����ĳ�ֱ��
		*
			\return true  ��
			\return false ��
		*/
		using Marker::testMark;

		/** \brief ������ĳ�ֱ��
		*
			������Ӧ���λ��1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< ���ص��������
		using Marker::setIndex;
	public:
		///< ǰ���ƽ����
		Halfedge_handle           loopNext()                               { return _loop_nxt;}
		void                      set_loopNext(Halfedge_handle h) { _loop_nxt = h;} ///< �������ڻ���̰��

		Halfedge_handle           loopPrev()                      { return _loop_prv; }
		void                      set_loopPrev(Halfedge_handle h) { _loop_prv = h;    } ///< �������ڻ�ǰ�����

		void                      set_levelNum(int n)             { _level_num = n; } ///< ����ǰ�����ڲ���
		int                       get_levelNum()                  { return _level_num; }

		void                      set_FrontState(int t)           { _state = t; } ///< ����ǰ�رߵ�״ֵ̬
		int                       get_FrontState()                { return _state; }

		void                      set_Key(int t)                  { _key = t; }
		int                       get_Key()                       { return _key; }

		void                      set_SquardLength(double t)      { _squard_length = t; }
		double                    get_SquardLength()              { return _squard_length; }

		void                      set_opposite( Halfedge_handle h)  { Base_base::set_opposite(h);}
	private:
		Halfedge_handle  _loop_nxt; ///< ���ڻ���̰��
		Halfedge_handle  _loop_prv; ///< ���ڻ�ǰ�����

		double _squard_length;///< �߳�,����ǰ������
		int _level_num; ///< ǰ�������Ĳ���
		int _state; ///< ǰ�ص�״ֵ̬
		int _key; ///< �ߵļ�ֵ(���ڿ��ٲ���)
};

template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs>, public Marker
{
	enum State{
		       VISITED = 1,   ///< ���ʱ��
			   TRI,           ///< ������
			   QUAD,          ///< �ı���
	          };
		using Marker::setDead;
		using Marker::isDead;
		/** \brief �����Ƿ�����ĳ�ֱ��
		*
			\return true  ��
			\return false ��
		*/
		using Marker::testMark;

		/** \brief ������ĳ�ֱ��
		*
			������Ӧ���λ��1
		*/
		using Marker::setMark;

		using Marker::resetMark;
		using Marker::clearMark;
		using Marker::getIndex; ///< ���ص��������
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
///< ���ڶ�ά�ཻ�ж�
// #include <CGAL/Arr_segment_traits_2.h>
// #include <CGAL/Sweep_line_2_algorithms.h>
// #include <CGAL/intersections.h>
//
// #include <CGAL/intersections.h>

///< ���ڶ�ά�ཻ�ж�
// typedef CGAL::Exact_predicates_exact_constructions_kernel                                   Kernel;
// typedef CGAL::Arr_segment_traits_2<Kernel>                                                  Intersect_Traits_2;
// typedef Intersect_Traits_2::Curve_2                                                         Intersect_Segement_2;
// typedef Kernel::Segment_2                                                                   Segment_2;
// typedef Kernel::Point_2                                                                     Point_2;
// typedef Kernel::Vector_2                                                                    Vector_2;
// typedef Kernel::Line_2                                                                      Line_2;
// typedef Kernel::Intersect_2                                                                 Intersect_2;

#endif
