//**************************************************************************
// File:     MeshOptimizer.cpp
// Author:   Zhang Jun
// Date:     2016-11-10
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 负责网格优化
//**************************************************************************
#include "MeshOptimizer.h"

double angle_evaluate2(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
{
	const double dx1 = q.x() - p.x();
	const double dx2 = s.x() - r.x();

	const double dy1 = q.y() - p.y();
	const double dy2 = s.y() - r.y();

	const double dz1 = q.z() - p.z();
	const double dz2 = s.z() - r.z();

	double ans = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1) * std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2));

	if (ans >= 0)
		return ans < 1 ? ans : 1;
	else
		return ans < -1 ? -1 : ans;
}

double angle_evaluate(const Point_3 &p, const Point_3 &q, const Point_3 &r)
{
	double cosAngle = angle_evaluate2(q, p, q, r);
	double angle_temp = std::acos(cosAngle) * 180.0 / PI;

	CGAL::Orientation ori = CGAL::orientation(Point_2(p.x(), p.y()), Point_2(q.x(), q.y()), Point_2(r.x(), r.y()));
	if (ori == CGAL::RIGHT_TURN)
		angle_temp = 360 - angle_temp;
	return angle_temp;
}

Point_3 operator+ (Point_3&p, Point_3& q)
{
	return Point_3(p.x()+q.x(), p.y()+q.y(), p.z()+q.z());
}
Point_3 operator- (Point_3&p)
{
	return Point_3(-p.x(), -p.y(), -p.z());
}
Point_3 operator* (double t, Point_3&p)
{
	return Point_3(p.x() * t, p.y() * t, p.z() * t);
}
Point_3 operator* (Point_3&p, double t)
{
	return Point_3(p.x() * t, p.y() * t, p.z() * t);
}
Point_3 operator/ (Point_3&p, double t)
{
	assert(t != 0);
	return Point_3(p.x() / t, p.y() / t, p.z() / t);
}



void MeshOptimizer::set_SourceMesh(HDS_mesh *source_mesh)
{
	_source_mesh = source_mesh;
}

bool MeshOptimizer::mesh_OPT()
{
	assert(_source_mesh);

	global_smooth(3); ///< 全局加权光顺
	//return true;
	boundary_unform();

	global_smooth(3); ///< 全局加权光顺
	//return true;
	split_vertex(); ///< 通过点分裂操作来消除三角形

	global_smooth(3); ///< 全局加权光顺

	clean_tri_zj();

	global_smooth(3); ///< 全局加权光顺
	//return true;
	quad_mesh_clean(); ///< 拓扑优化

	global_smooth(3); ///< 全局加权光顺
	//return true;
	quad_face_clean_aux2();

	global_smooth(3); ///< 全局加权光顺
	//boundary_unform();

	//global_smooth(3); ///< 全局加权光顺
	return true;
}

bool MeshOptimizer::global_smooth(int n)
{
	HDS_mesh::Vertex_iterator it_v;
	for (int i = 0; i < n; ++i)
	{
		for (it_v = _source_mesh->vertices_begin(); it_v != _source_mesh->vertices_end(); ++it_v)
		{
			if (it_v->testMark(HDS_mesh::Vertex::BOUNDARY))
				continue;
			if (it_v->testMark(HDS_mesh::Vertex::FrontNode))
				continue;
			smoothInteriorNode(it_v);
		}
	}
	return true;
}

int MeshOptimizer::smoothInteriorNode(HDS_mesh::Vertex_handle& vh)
{
	if (vh->testMark(HDS_mesh::Vertex::LOCKED))
		return 1;
	///< 首先获取点vh的邻接节点链表
	Point_3 old_point = vh->point();
	Point_3 new_point(0,0,0);
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh->vertex_begin();
	double len_sum = 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	do 
	{
		double len = CGAL::sqrt(CGAL::squared_distance(old_point, edge_begin->prev()->vertex()->point()));
		x += len * (edge_begin->prev()->vertex()->point().x() - old_point.x());
		y += len * (edge_begin->prev()->vertex()->point().y() - old_point.y());
		z += len * (edge_begin->prev()->vertex()->point().z() - old_point.z());
		len_sum += len;
		++edge_begin;
	} while (edge_begin != edge_end);
	x = x / len_sum + old_point.x();
	y = y / len_sum + old_point.y();
	z = z / len_sum + old_point.z();
	new_point = Point_3(x, y, z);
	return smoothAdjustment(vh, new_point);
	return 1;
}

int MeshOptimizer::smoothAdjustment(HDS_mesh::Vertex_handle& vh, Point_3& new_point)
{
	Point_3 old_point = vh->point();
	Vector_3 vec_temp(new_point, old_point);
	vh->point() = new_point;
	Point_3 p_temp = new_point;
	int num = 0;///< 统计调整次数
	///< 遍历节点vh的入射边表
	bool flag = false;
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh->vertex_begin();
	do 
	{
		if (!edge_begin->is_border())
		{
			HDS_mesh::Face_handle face_temp = edge_begin->face();
			if (face_temp->testMark(HDS_mesh::Face::TRI))
			{
				Point_3 pa = face_temp->halfedge()->prev()->vertex()->point();
				Point_3 pb = face_temp->halfedge()->vertex()->point();
				Point_3 pc = face_temp->halfedge()->next()->vertex()->point();
				CGAL::Orientation ori = CGAL::orientation(Point_2(pa.x(), pa.y()), Point_2(pb.x(), pb.y()), Point_2(pc.x(), pc.y()));
				if (ori == CGAL::RIGHT_TURN)
				{
					++num;
					p_temp = p_temp + 0.1 * vec_temp;
					vh->point() = p_temp;
					flag = true;
					if (num == 10)
						break;
					continue;
				}
				else
					flag = false;
			}
		}
		++edge_begin;
	} while (edge_begin != edge_end || flag);
	return 0;
}

double MeshOptimizer::tri_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2)
{
	return tri_quality(v0->point(), v1->point(), v2->point());
}

double MeshOptimizer::tri_quality(Point_3& p0, Point_3&  p1, Point_3&  p2)
{
	double d0 = CGAL::sqrt(CGAL::squared_distance(p0, p1));
	double d1 = CGAL::sqrt(CGAL::squared_distance(p0, p2));
	double d2 = CGAL::sqrt(CGAL::squared_distance(p1, p2));

	double min = d0 < d1 ? d0 : d1;
	min = min < d2 ? min : d2;
	min *= 0.1;
	min = 1.0 / min;
	d0 *= min;
	d1 *= min;
	d2 *= min;
	double ditor = 2 * sqrt(3.0);///< modified by ZhangJun 2016.06.07
	double bottom = d0 * d0 + d1 * d1 + d2 * d2;
	double det = (d0 + d1 + d2) * ( d0 + d1 - d2) * ( d0 + d2 - d1) * ( d1 + d2 - d0);
	if(det < EPS_QuadMesh) return 0.0;
	double up = sqrt(det) * 0.5;
	int I = 1;
	if (CGAL::orientation(Point_2(p0.x(), p0.y()), 
		Point_2(p1.x(), p1.y()), 
		Point_2(p2.x(), p2.y())) == CGAL::RIGHT_TURN)
		I = -1;
	return I * (ditor * up / bottom);
}

double MeshOptimizer::quad_quality(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
									   HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3)
{
	std::vector<double> q_tri;
	q_tri.push_back(tri_quality(v0, v1, v2));
	q_tri.push_back(tri_quality(v0, v2, v3));
	q_tri.push_back(tri_quality(v0, v1, v3));
	q_tri.push_back(tri_quality(v1, v2, v3));
	std::sort(q_tri.begin(), q_tri.end());

	return (q_tri[0] * q_tri[1] / (q_tri[2] * q_tri[3]));
}

double MeshOptimizer::quad_quality2(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, 
										HDS_mesh::Vertex_handle v2, HDS_mesh::Vertex_handle v3)
{
	return quad_quality2(v0->point(), v1->point(), v2->point(), v3->point());
}

double MeshOptimizer::quad_quality2(Point_3& p0, Point_3& p1,Point_3& p2, Point_3& p3)
{
	std::vector<double> q_tri;
	q_tri.push_back(tri_quality(p0, p1, p2));
	q_tri.push_back(tri_quality(p0, p2, p3));
	q_tri.push_back(tri_quality(p0, p1, p3));
	q_tri.push_back(tri_quality(p1, p2, p3));
	std::sort(q_tri.begin(), q_tri.end());
	int negval = 0;
	int n_negative = 0;
	for (int i = 0; i < 4; ++i)
	{
		if (q_tri[i] < 0)
			++n_negative;
	}
	if (2 == n_negative)
		negval = 1;
	if (3 == n_negative)
		negval = 2;
	else if (4 == n_negative)
		negval = 3;
	return (q_tri[0] - negval);
}

bool MeshOptimizer::boundary_unform()
{
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	HDS_mesh::Vertex_iterator it_v;
	std::list<HDS_mesh::Vertex_handle> v_list;
	for (it_v = _source_mesh->vertices_begin(); it_v != _source_mesh->vertices_end(); ++ it_v)
	{
		if (it_v->testMark(HDS_mesh::Vertex::BOUNDARY)
			&& it_v->testMark(HDS_mesh::Vertex::SIDENODE))
		{
			v_list.push_back(it_v);
			if (it_v->degree() != 3)
			{
				HDS_mesh::Halfedge_around_vertex_circulator edge_begin = it_v->vertex_begin();
				HDS_mesh::Halfedge_around_vertex_circulator edge_end = it_v->vertex_begin();
				do 
				{
					if (!edge_begin->is_border())
					{
						if (!edge_begin->prev()->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY)
							&& !edge_begin->next()->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY))
						{
							break;
						}
					}
					++edge_begin;
				} while (edge_begin != edge_end);
				if (edge_begin->opposite()->is_border())
					continue;
				Point_3 p_temp = (edge_begin->prev()->vertex()->point() 
					              + edge_begin->next()->vertex()->point()) / 2.0;
				HDS_mesh::Halfedge_handle edgeprev = edge_begin->prev();
				HDS_mesh::Halfedge_handle edgenext = edge_begin->next();
				HDS_mesh::Halfedge_handle edge_temp = decorator.split_face(edgeprev, edgenext);
				decorator.join_face(edge_begin);
				decorator.join_face(edgeprev);
				edge_temp = decorator.join_vertex(edge_temp);
				edge_temp->vertex()->point() = p_temp;
				//edge_temp->vertex()->setMark(HDS_mesh::Vertex::LOCKED);
				//return true;
			}
		}
	}
	std::list<HDS_mesh::Vertex_handle>::iterator it_VL;
	for (it_VL = v_list.begin(); it_VL != v_list.end(); ++it_VL)
	{
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = (*it_VL)->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = (*it_VL)->vertex_begin();
		do 
		{
			edge_begin->prev()->vertex()->setMark(HDS_mesh::Vertex::LOCKED);
			++edge_begin;
		} while (edge_begin != edge_end);
	}
	return true;
}

bool MeshOptimizer::split_vertex()
{
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	std::list<HDS_mesh::Vertex_iterator> vertex_list;
	HDS_mesh::Vertex_iterator it_v;
	for (it_v = _source_mesh->vertices_begin(); it_v != _source_mesh->vertices_end(); ++ it_v)
	{
		if (!it_v->testMark(HDS_mesh::Vertex::BOUNDARY))
		{
			if (it_v->is_bivalent() || it_v->is_trivalent())
				continue;

			int n_tri = 0;
			int n_quad = 0;
			HDS_mesh::Halfedge_around_vertex_circulator edge_begin = it_v->vertex_begin();
			HDS_mesh::Halfedge_around_vertex_circulator edge_end = it_v->vertex_begin();
			do 
			{
				if (edge_begin->is_triangle())
					++n_tri;
				else
					++n_quad;
				++edge_begin;
			} while (edge_begin != edge_end);
			assert(n_tri + n_quad >= 4);
			if (n_tri < 2)
				continue;
			vertex_list.push_back(it_v);
		}
	}
	std::list<HDS_mesh::Vertex_iterator>::iterator i_v;
	for (i_v = vertex_list.begin(); i_v != vertex_list.end(); ++i_v)
	{
		std::vector<HDS_mesh::Halfedge_handle> edge_vector;
		std::vector<Point_3> point_vector;
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = (*i_v)->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = (*i_v)->vertex_begin();
		do 
		{
			if (edge_begin->is_triangle())
				edge_vector.push_back(edge_begin);
			++edge_begin;
		} while (edge_begin != edge_end);
		if (edge_vector.size() != 2)
			continue;
		assert(edge_vector.size() == 2);
		_source_mesh->split_vertex(edge_vector[0], edge_vector[1]);
	}

	return true;
}

bool MeshOptimizer::clean_tri_zj()
{
	std::vector<std::pair<HDS_mesh::Face_handle, int>> face_vector;
	HDS_mesh::Face_iterator it_f;
	int n_merge = 0;
	for (it_f = _source_mesh->facets_begin(); it_f != _source_mesh->facets_end();)
	{
		if (!it_f->testMark(HDS_mesh::Face::TRI))
		{
			++it_f;
			continue;
		}
		HDS_mesh::Halfedge_handle edge_begin = it_f->halfedge();
		HDS_mesh::Halfedge_handle edge_end = it_f->halfedge();
		int eid = 0;
		do 
		{
			HDS_mesh::Halfedge_handle edge_temp = edge_begin->opposite()->next()->next()->opposite();
			if (!edge_temp->is_border())
			{
				HDS_mesh::Face_handle face_temp = edge_temp->face();
				if (face_temp->is_triangle())
				{
					face_vector.push_back(std::make_pair(it_f, eid));
					break;
				}
			}
			++eid;
			edge_begin = edge_begin->next();
		} while (edge_begin != edge_end);
		++it_f;
	}
	std::vector<std::pair<HDS_mesh::Face_handle, int>>::iterator it_tri;
	for (it_tri = face_vector.begin(); it_tri != face_vector.end(); ++it_tri)
	{
		if (!(*it_tri).first->is_triangle())
			continue;
		HDS_mesh::Face_handle tri1_temp = (*it_tri).first;
		HDS_mesh::Halfedge_handle edge_temp = tri1_temp->halfedge();
		for (int i = 0; i < (*it_tri).second; ++i)
		{
			edge_temp = edge_temp->next();
		}
		HDS_mesh::Face_handle quad_temp = edge_temp->opposite()->face();
		HDS_mesh::Face_handle tri2_temp =  edge_temp->opposite()->next()->next()->opposite()->face();
		assert(tri1_temp->is_triangle());
		if (!tri2_temp->is_triangle())
			continue;
		assert(tri2_temp->is_triangle());
		if (!quad_temp->is_quad())
			continue;
		assert(quad_temp->is_quad());

		HDS_mesh::Halfedge_handle v_split_edge1 = edge_temp->opposite()->prev()->opposite()->prev();
		HDS_mesh::Halfedge_handle v_split_edge2 = edge_temp->opposite()->next()->next()->opposite()->prev();

		_source_mesh->join_vertex(edge_temp->opposite()->prev());
		_source_mesh->join_facet(edge_temp);
		_source_mesh->split_vertex(v_split_edge1, v_split_edge2);
	}
	return true;
}

bool MeshOptimizer::quad_mesh_clean()
{
	//write_my_mesh_nas("clean.nas");
	// 	quad_face_clean();
	// 	quad_bivalent_vertex_clean();
	while (1)
	{
		bool flag1 = quad_face_clean();
		bool flag3 = quad_face_clean_aux1();
		bool flag2 = quad_bivalent_vertex_clean();

		if (!flag1 && !flag2 && !flag3)
			break;
	}
	// 	quad_face_clean();
	// 	quad_bivalent_vertex_clean();
	//global_smooth(3);
	//quad_vertex_clean();
	return true;
}

bool MeshOptimizer::quad_face_clean()
{
	HDS_mesh::Face_iterator it_f;
	int n = 0;
	int n_merge = 0;
	bool flag = false;
	do 
	{
		n_merge = 0;
		for (it_f = _source_mesh->facets_begin(); it_f != _source_mesh->facets_end();)
		{
			if (!it_f->testMark(HDS_mesh::Face::QUAD))
			{
				++it_f;
				continue;
			}
			HDS_mesh::Vertex_handle v1 = it_f->halfedge()->vertex();
			HDS_mesh::Vertex_handle v2 = it_f->halfedge()->next()->vertex();
			HDS_mesh::Vertex_handle v3 = it_f->halfedge()->next()->next()->vertex();
			HDS_mesh::Vertex_handle v4 = it_f->halfedge()->prev()->vertex();
			if (!v1->testMark(HDS_mesh::Vertex::BOUNDARY) && v1->degree() == 3
				&& !v3->testMark(HDS_mesh::Vertex::BOUNDARY) && v3->degree() == 3)
			{
				++n;
				HDS_mesh::Halfedge_handle e_base =  it_f->halfedge();

				it_f = ++it_f;
				face_clean_aux(e_base);
				++n_merge;
				flag = true;
				continue;
			}
			if (!v2->testMark(HDS_mesh::Vertex::BOUNDARY) && v2->degree() == 3
				&& !v4->testMark(HDS_mesh::Vertex::BOUNDARY) && v4->degree() == 3)
			{
				++n;
				HDS_mesh::Halfedge_handle e_base =  it_f->halfedge()->prev();
				it_f = ++it_f;
				face_clean_aux(e_base);
				++n_merge;
				flag = true;
				continue;
			}
			++it_f;
		}
	} while (n_merge > 0);

	return flag;
}

bool MeshOptimizer::face_clean_aux(HDS_mesh::Halfedge_handle& e_base)
{
	assert(e_base->face()->testMark(HDS_mesh::Face::QUAD));
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());

	HDS_mesh::Vertex_handle v = e_base->vertex();
	// 
	// 	HDS_mesh::Vertex_handle vo = e_base->next()->next()->vertex();
	// 
	// 	HDS_mesh::Halfedge_handle e_base_naxe_next = e_base->next()->next();
	// 
	// 	HDS_mesh::Halfedge_handle eo__naxe_next = e_base->opposite()->next()->next();
	// 
	// 	double x = (v->point().x() + vo->point().x()) / 2.0;
	// 	double y = (v->point().y() + vo->point().y()) / 2.0;
	// 	double z = (v->point().z() + vo->point().z()) / 2.0;
	// 	vo->point() = Point_3(x, y, z);

	// 	_source_mesh->erase_center_vertex(e_base);
	// 
	// 	_source_mesh->split_facet(e_base_naxe_next, eo__naxe_next);
	HDS_mesh::Vertex_handle vo = e_base->next()->next()->vertex();
	HDS_mesh::Halfedge_handle e_base_next = e_base->next();
	HDS_mesh::Halfedge_handle e_base_prev = e_base->prev();
	HDS_mesh::Halfedge_handle e_base_naxe_next = e_base_next->next();

	HDS_mesh::Halfedge_handle e_base_o = e_base->opposite();
	HDS_mesh::Halfedge_handle e_base_next_o = e_base_next->opposite();
	HDS_mesh::Halfedge_handle e_base_prev_o = e_base_prev->opposite();
	HDS_mesh::Halfedge_handle e_base_naxe_next_o = e_base_naxe_next->opposite();

	e_base_prev_o->face()->set_halfedge(e_base_prev_o->prev());
	e_base_naxe_next_o->face()->set_halfedge(e_base_naxe_next_o->prev());

	_source_mesh->make_hole(e_base);
	v->set_halfedge(e_base_next_o);
	e_base_prev->vertex()->set_halfedge(e_base_o);
	e_base_next->vertex()->set_halfedge(e_base_next);

	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vo->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vo->vertex_begin();
	do 
	{
		assert(vo == edge_begin->vertex());
		edge_begin->set_vertex(v);
		++edge_begin;
	} while (edge_begin != edge_end);

	double x = (v->point().x() + vo->point().x()) / 2.0;
	double y = (v->point().y() + vo->point().y()) / 2.0;
	double z = (v->point().z() + vo->point().z()) / 2.0;
	v->point() = Point_3(x, y, z);

	e_base_prev_o->prev()->set_next(e_base);

	e_base->set_prev(e_base_prev_o->prev());

	e_base_prev_o->next()->set_prev(e_base);

	e_base->set_next(e_base_prev_o->next());

	e_base->set_face(e_base_prev_o->face());


	e_base_naxe_next_o->prev()->set_next(e_base_next);

	e_base_next->set_prev(e_base_naxe_next_o->prev());

	e_base_naxe_next_o->next()->set_prev(e_base_next);

	e_base_next->set_next(e_base_naxe_next_o->next());

	e_base_next->set_face(e_base_naxe_next_o->face());

	decorator.vertices_erase(vo);

	_source_mesh->hds().edges_erase(e_base_prev);
	_source_mesh->hds().edges_erase(e_base_naxe_next);
	return true;
}

bool MeshOptimizer::quad_face_clean_aux1()
{
	bool flag = false;
	std::list<HDS_mesh::Halfedge_handle> edge_list;
	HDS_mesh::Halfedge_iterator it_edge;
	for (it_edge = _source_mesh->halfedges_begin(); it_edge != _source_mesh->halfedges_end(); ++it_edge)
	{
		if (it_edge->is_border() || it_edge->opposite()->is_border())
			continue;
		if (it_edge->testMark(HDS_mesh::Halfedge::VISITED))
			continue;
		if (it_edge->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY)
			|| it_edge->opposite()->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY))
			continue;
		if (it_edge->vertex()->is_trivalent() && it_edge->opposite()->vertex()->is_trivalent())
		{
			it_edge->setMark(HDS_mesh::Halfedge::VISITED);
			it_edge->opposite()->setMark(HDS_mesh::Halfedge::VISITED);
			edge_list.push_back(it_edge);
		}
	}
	std::list<HDS_mesh::Halfedge_handle>::iterator it_e;
	for (it_e = edge_list.begin(); it_e != edge_list.end(); ++it_e)
	{
		if ((*it_e) == HDS_mesh::Halfedge_handle())
			continue;
		if ((*it_e)->testMark(15))
			continue;
		HDS_mesh::Halfedge_handle e_next = (*it_e)->next();
		HDS_mesh::Halfedge_handle oe_next = (*it_e)->opposite()->next();
		if (oe_next->is_border())
			continue;
		if (!e_next->face()->is_quad() || !e_next->opposite()->face()->is_quad()
			|| !oe_next->face()->is_quad() || !oe_next->opposite()->face()->is_quad())
			continue;

		std::vector<HDS_mesh::Vertex_handle> v_vector_temp;
		v_vector_temp.push_back((*it_e)->next()->vertex());
		v_vector_temp.push_back((*it_e)->next()->next()->vertex());
		v_vector_temp.push_back((*it_e)->prev()->opposite()->next()->vertex());
		v_vector_temp.push_back((*it_e)->opposite()->next()->vertex());
		v_vector_temp.push_back((*it_e)->opposite()->next()->next()->vertex());
		v_vector_temp.push_back((*it_e)->opposite()->prev()->opposite()->next()->vertex());
		double max_qua = DBL_MIN;
		int k = -1;
		for (int i = 0; i < 2; ++i)
		{
			double quality1 = quad_quality(v_vector_temp[i], v_vector_temp[i+1], v_vector_temp[i+2], v_vector_temp[i+3]);
			double quality2 = quad_quality(v_vector_temp[i+3], v_vector_temp[(i+4)%6], v_vector_temp[(i+5)%6], v_vector_temp[i]);
			double min_qua_tamp = quality1 < quality2 ? quality1 : quality2;
			if (max_qua < min_qua_tamp)
			{
				max_qua = min_qua_tamp;
				k = i;
			}
		}
		HDS_mesh::Halfedge_handle e_base = _source_mesh->erase_center_vertex(e_next->opposite());
		_source_mesh->erase_center_vertex(oe_next->opposite());
		for (int j = 0; j < k; ++j)
		{
			e_base = e_base->next();
		}
		_source_mesh->split_facet(e_base, e_base->next()->next()->next());
		flag = true;
		//return true;
	}
	return flag;
}

bool MeshOptimizer::quad_bivalent_vertex_clean()
{
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> v_set;
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	HDS_mesh::Vertex_iterator it_v;
	bool flag = false;
	for (it_v = _source_mesh->vertices_begin(); it_v != _source_mesh->vertices_end(); ++it_v)
	{
		if (!it_v->testMark(HDS_mesh::Vertex::BOUNDARY))
		{
			if (it_v->is_bivalent())
			{
				HDS_mesh::Face_handle face_base = it_v->halfedge()->face();
				HDS_mesh::Face_handle oface = it_v->halfedge()->opposite()->face();
				//if(face_base->is_quad() && oface->is_quad())
				{
					v_set.insert(it_v);
				}						
			}
		}
	}
	if (!v_set.empty())
		flag = true;
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_vs;
	for (it_vs = v_set.begin(); it_vs != v_set.end(); ++it_vs)
	{
		HDS_mesh::Halfedge_handle e_base = (*it_vs)->halfedge();
		e_base = _source_mesh->erase_center_vertex(e_base);
		//write_my_mesh_nas("clean.nas");
	}

	return flag;
}

bool MeshOptimizer::quad_face_clean_aux2()
{
	bool flag = false;
	std::list<HDS_mesh::Halfedge_handle> edge_list;///< 入射点内角接近180度的边
	HDS_mesh::Face_iterator it_face;
	for (it_face = _source_mesh->facets_begin(); it_face != _source_mesh->facets_end(); ++it_face)
	{
		if (!it_face->is_quad())
			continue;
		HDS_mesh::Halfedge_handle edge_begin = it_face->halfedge();
		HDS_mesh::Halfedge_handle edge_end = it_face->halfedge();
		do 
		{
			double ang_temp = angle_evaluate(edge_begin->prev()->vertex()->point(),
				edge_begin->vertex()->point(),
				edge_begin->next()->vertex()->point());
			if (ang_temp > 170)
			{
				edge_list.push_back(edge_begin);
				break;
			}
			edge_begin = edge_begin->next();
		} while (edge_begin != edge_end);
	}
	std::list<HDS_mesh::Halfedge_handle>::iterator it_e;
	for (it_e = edge_list.begin(); it_e != edge_list.end(); ++it_e)
	{
		if ((*it_e) == HDS_mesh::Halfedge_handle() || (*it_e)->testMark(15))
			continue;
		if ((*it_e)->prev()->prev()->vertex()->testMark(HDS_mesh::Vertex::LOCKED)
			&& (*it_e)->prev()->vertex()->testMark(HDS_mesh::Vertex::LOCKED))
		{
			Point_3 p_temp = 0.5 * ((*it_e)->prev()->prev()->vertex()->point() + (*it_e)->vertex()->point());
			if ((*it_e)->prev()->is_border_edge())
				continue;
			HDS_mesh::Halfedge_handle g = _source_mesh->join_facet((*it_e)->prev());
			g = _source_mesh->join_facet(g->opposite());
			g = _source_mesh->create_center_vertex(g);
			g->vertex()->resetMark(HDS_mesh::Vertex::LOCKED);
			g->vertex()->point() = p_temp;
			_source_mesh->join_facet(g->next()->opposite()->next());
			_source_mesh->join_facet(g->opposite()->prev()->opposite()->prev());

			_source_mesh->split_vertex(g, g->next()->opposite()->next()->opposite());
			_source_mesh->split_vertex(g->opposite()->prev(), g->next()->opposite()->next()->opposite());
			// 			g = _backgroundMesh.join_facet(g->next());
			// 			_backgroundMesh.join_facet(g->opposite()->prev()->opposite()->prev());
			// 			_backgroundMesh.split_vertex(g->opposite()->prev(), g->next()->opposite());
		}
		else
		{
			_source_mesh->join_facet((*it_e)->prev());
			_source_mesh->split_facet((*it_e), (*it_e)->next()->next()->next());
		}
		// 		std::cout << (*it_e)->prev()->vertex()->point() << std::endl;
		// 		std::cout << (*it_e)->vertex()->point() << std::endl;
		// 		std::cout << "=======================================\n";
		flag = true;
		//return true;
	}
	return flag;
}

bool MeshOptimizer::localReconnect()
{
	for(int iter=0; iter<100; iter++)
	{
		int swapEdgeCnt = 0;
		///< 遍历背景网格的半边链表,收集所有允许翻转的半边
		std::list<HDS_mesh::Halfedge_handle> edge_list_temp;
		HDS_mesh::Halfedge_iterator it_edge;
		for (it_edge=_source_mesh->halfedges_begin(); 
			it_edge!=_source_mesh->halfedges_end(); 
			++it_edge)
		{
			if (it_edge->testMark(HDS_mesh::Halfedge::VISITED))
				continue;
			bool quad_flag = false;
			if (it_edge->testMark(HDS_mesh::Halfedge::QuadEdge))
				quad_flag = true;
			if (!it_edge->testMark(HDS_mesh::Halfedge::BOUNDARY))
			{
				if (it_edge->opposite()->testMark(HDS_mesh::Halfedge::QuadEdge))
					quad_flag = true;
				if (quad_flag)
					continue;
				it_edge->setMark(HDS_mesh::Halfedge::VISITED);
				it_edge->opposite()->setMark(HDS_mesh::Halfedge::VISITED);
				edge_list_temp.push_back(it_edge);
			}
		}
		std::list<HDS_mesh::Halfedge_handle>::iterator it_edge_to_swap;
		for (it_edge_to_swap = edge_list_temp.begin(); 
			it_edge_to_swap != edge_list_temp.end(); 
			++it_edge_to_swap)
		{
			HDS_mesh::Halfedge_handle e_result;
#if 0 //暂时封掉 2016.11.10
			if (flip_edge(*it_edge_to_swap, true))
				++swapEdgeCnt;
			else
			{
				(*it_edge_to_swap)->resetMark(HDS_mesh::Halfedge::VISITED);
				(*it_edge_to_swap)->opposite()->resetMark(HDS_mesh::Halfedge::VISITED);
			}
#endif
		}
		if(swapEdgeCnt == 0) 
		{
			break;
		}
	}
	return true;
}
// bool QMorphMesher::isValid(HDS_mesh::Vertex_handle v0, HDS_mesh::Vertex_handle v1, HDS_mesh::Vertex_handle v2)
// {
// 	Vector3D v01(v0->point(), v1->point());
// 	Vector3D v02(v0->point(), v2->point());
// 	double x01 = v01.getEnd().getX() - v01.getStart().getX();
// 	double y01 = v01.getEnd().getY() - v01.getStart().getY();
// 
// 	double x02 = v02.getEnd().getX() - v02.getStart().getX();
// 	double y02 = v02.getEnd().getY() - v02.getStart().getY();
// 
// 	return x01 * y02 - y01 * x02 > 0;
// }