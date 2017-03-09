// **************************************************************************
// File:     TriMergeImp.cpp
// Author:   Zhang Jun
// Date:     2016-11-07
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 四边形网格生成算法实现(基于三角形合并的算法)
// **************************************************************************
#include "TriMergeGenerator.h"
#include "MeshOptimizer.h"

TriMergeGenerator::TriMergeGenerator()
{
	_source_mesh = nullptr;
	_optimizer = new MeshOptimizer;
}

TriMergeGenerator::~TriMergeGenerator()
{
	delete _optimizer;
}

void TriMergeGenerator::set_SourceMesh(HDS_mesh *source_mesh)
{
	_source_mesh = source_mesh;
	_optimizer->set_SourceMesh(source_mesh);
	//_active_front_list = active_front;
	//_unactive_front_list = unactive_front;

	HDS_mesh::Vertex_iterator it_v;
	HDS_mesh::Halfedge_iterator it_edge;

	for (it_edge = _source_mesh->halfedges_begin(); it_edge != _source_mesh->halfedges_end(); ++it_edge)
	{
		if(it_edge->is_border())
		{
			it_edge->opposite()->vertex()->setMark(HDS_mesh::Vertex::BOUNDARY);///< 标记节点为边界节点
			//it_edge->opposite()->vertex()->setMark(HDS_mesh::Vertex::FrontNode);///< 标记节点为前沿点
			
			it_edge->opposite()->setMark(HDS_mesh::Halfedge::BOUNDARY);///< 标记半边为边界边
			it_edge->opposite()->setMark(HDS_mesh::Halfedge::FRONT);///< 把边界边标记为初始前沿
			it_edge->opposite()->vertex()->set_halfedge(it_edge->opposite());
			//it_edge->set_levelNum(0); ///< 初始前沿层数为0
			//double len_temp = CGAL::squared_distance(it_edge->vertex()->point(), it_edge->opposite()->vertex()->point());
			//it_edge->opposite()->set_SquardLength(len_temp);

			_active_front_list.push_back(it_edge->opposite());
		}
	}
	std::list<HDS_mesh::Halfedge_handle>::iterator it_front;
	for (it_front = _active_front_list.begin(); it_front != _active_front_list.end(); ++it_front)
	{
		HDS_mesh::Halfedge_handle prv_front = (*it_front)->prev()->vertex()->halfedge();
		assert(prv_front->testMark(HDS_mesh::Halfedge::FRONT));
		(*it_front)->set_loopPrev(prv_front);
		prv_front->set_loopNext(*it_front);
	}
// 	assert(_min_len > 0);
// 	_reatio = std::sqrt(_max_len * 1.0 / _min_len);
}

void TriMergeGenerator::deleteFromUnactiveFront(HDS_mesh::Halfedge_handle eh)
{
	std::list<HDS_mesh::Halfedge_iterator>::iterator it_edge_temp;
	for (it_edge_temp = _unactive_front_list.begin(); it_edge_temp != _unactive_front_list.end(); ++it_edge_temp)
	{
		if ((*it_edge_temp)->vertex() == eh->vertex() && (*it_edge_temp)->prev()->vertex() == eh->prev()->vertex())
			break;
	}
	if (it_edge_temp != _unactive_front_list.end())
		_unactive_front_list.erase(it_edge_temp);
}

bool TriMergeGenerator::get_Front_Merge(HDS_mesh::Halfedge_handle& eh)
{
	if (!_active_front_list.empty())
	{
		eh = _active_front_list.front();
		_active_front_list.pop_front();
		return true;
	}
	else if (!_unactive_front_list.empty())
	{
		_active_front_list.clear();
		std::list<HDS_mesh::Halfedge_handle>::iterator it_front_temp;
		for (it_front_temp = _unactive_front_list.begin(); 
			it_front_temp != _unactive_front_list.end(); ++it_front_temp)
		{
			if ((*it_front_temp)->testMark(HDS_mesh::Halfedge::FRONT))
			{
				(*it_front_temp)->resetMark(HDS_mesh::Halfedge::UNACTIVE);
				_active_front_list.push_back((*it_front_temp));
			}
		}
		_unactive_front_list.clear();
		if (!_active_front_list.empty())
		{
			eh = _active_front_list.front();
			_active_front_list.pop_front();
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

void TriMergeGenerator::deleteFromFront(HDS_mesh::Halfedge_handle eh)
{
	std::list<HDS_mesh::Halfedge_handle>::iterator front_temp;
	front_temp = _active_front_list.begin();
	while (front_temp != _active_front_list.end() && (*front_temp) != eh)
	{
		++front_temp;
	}
	if (front_temp != _active_front_list.end())
		_active_front_list.erase(front_temp);	
}

HDS_mesh::Face_handle TriMergeGenerator::merge_Tri(HDS_mesh::Face_handle& tri1, 
											  HDS_mesh::Face_handle& tri2)
{
	assert(tri1->testMark(HDS_mesh::Face::TRI));
	assert(tri2->testMark(HDS_mesh::Face::TRI));

	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());

	HDS_mesh::Halfedge_handle it_edge_begin = tri1->halfedge();
	HDS_mesh::Halfedge_handle it_edge_end = tri1->halfedge();
	do 
	{///< 找到两三角形的公共边
		if (!it_edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY))
		{
			if (it_edge_begin->opposite()->face() == tri2)
				break;
		}
		it_edge_begin = it_edge_begin->next();
	} while (it_edge_begin != it_edge_end);

	HDS_mesh::Halfedge_handle oedge = it_edge_begin->opposite();
	HDS_mesh::Halfedge_handle e11 = it_edge_begin->next();
	HDS_mesh::Halfedge_handle e12 = it_edge_begin->prev();
	HDS_mesh::Halfedge_handle e21 = oedge->next();
	HDS_mesh::Halfedge_handle e22 = oedge->prev();

	HDS_mesh::Halfedge_handle e_down;
	if (e11->testMark(HDS_mesh::Halfedge::ISDOWN))
	{
		e11->resetMark(HDS_mesh::Halfedge::ISDOWN);
		e_down = e11;
	}
	else
	{
		e12->resetMark(HDS_mesh::Halfedge::ISDOWN);
		e_down = e12;
	}
	assert(e_down->testMark(HDS_mesh::Halfedge::FRONT)); ///< 为方便起见，默认新单元的底边为当前前沿

	e11->setMark(HDS_mesh::Halfedge::QuadEdge);
	e11->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e12->setMark(HDS_mesh::Halfedge::QuadEdge);
	e12->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e21->setMark(HDS_mesh::Halfedge::QuadEdge);
	e21->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e22->setMark(HDS_mesh::Halfedge::QuadEdge);
	e22->vertex()->setMark(HDS_mesh::Vertex::QuadNode);

	tri1->resetMark(HDS_mesh::Face::TRI);
	tri1->setMark(HDS_mesh::Face::QUAD);

	decorator.join_face(it_edge_begin);

	HDS_mesh::Halfedge_handle e_right = e_down->next();
	HDS_mesh::Halfedge_handle e_up = e_right->next();
	HDS_mesh::Halfedge_handle e_left = e_up->next();

	///< 更新前沿拓扑关系
	e_down->resetMark(HDS_mesh::Halfedge::FRONT); ///< 取消下侧边的前沿标记

	if (!e_up->testMark(HDS_mesh::Halfedge::FRONT))
	{
		e_up->opposite()->setMark(HDS_mesh::Halfedge::FRONT); ///< 可能会有问题
		e_up->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
		_unactive_front_list.push_back(e_up->opposite());
	}
	else
	{
		if (e_up->testMark(HDS_mesh::Halfedge::UNACTIVE))
		{
			e_up->resetMark(HDS_mesh::Halfedge::FRONT);
			e_up->resetMark(HDS_mesh::Halfedge::UNACTIVE);
			deleteFromUnactiveFront(e_up);
		}
		else
		{
			e_up->resetMark(HDS_mesh::Halfedge::FRONT);
			deleteFromFront(e_up);
		}
	}

	if (e_right->testMark(HDS_mesh::Halfedge::FRONT))
	{
		e_right->resetMark(HDS_mesh::Halfedge::FRONT);
		deleteFromFront(e_right);
	}
	else if(!e_right->testMark(HDS_mesh::Halfedge::BOUNDARY))
	{
		if (e_right->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
		{
			e_right->opposite()->resetMark(HDS_mesh::Halfedge::FRONT);
			deleteFromFront(e_right->opposite());
		}
		else
		{
			e_right->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			e_right->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
			_unactive_front_list.push_back(e_right->opposite());
		}
	}

	if (e_left->testMark(HDS_mesh::Halfedge::FRONT))
	{
		e_left->resetMark(HDS_mesh::Halfedge::FRONT);
		if (e_left->testMark(HDS_mesh::Halfedge::UNACTIVE))
		{
			e_left->resetMark(HDS_mesh::Halfedge::UNACTIVE);
			deleteFromUnactiveFront(e_left);
		}
		else
			deleteFromFront(e_left);
	}
	else if (!e_left->testMark(HDS_mesh::Halfedge::BOUNDARY))
	{
		if (e_left->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
		{
			e_left->opposite()->resetMark(HDS_mesh::Halfedge::FRONT);
			if (e_left->opposite()->testMark(HDS_mesh::Halfedge::UNACTIVE))
			{
				e_left->opposite()->resetMark(HDS_mesh::Halfedge::UNACTIVE);
				deleteFromUnactiveFront(e_left->opposite());
			}
			else
				deleteFromFront(e_left->opposite());
		}
		else
		{
			e_left->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			e_left->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
			_unactive_front_list.push_back(e_left->opposite());
		}
	}

	return tri1;
}

void TriMergeGenerator::mesh()
{
	int n = 0;
	while (1)
	{
		HDS_mesh::Halfedge_handle current_front; ///< 获取前沿边
		if (!get_Front_Merge(current_front))
			break;
		current_front->setMark(HDS_mesh::Halfedge::ISDOWN);
		HDS_mesh::Face_handle tri_face_base = current_front->face();
		assert(tri_face_base->testMark(HDS_mesh::Face::TRI));///< 当前前沿的入射面片必须为三角形
		HDS_mesh::Halfedge_handle it_edge_begin = tri_face_base->halfedge();
		HDS_mesh::Halfedge_handle it_edge_end = tri_face_base->halfedge();
		std::pair<HDS_mesh::Face_handle, double> tri_qual(tri_face_base, -1.0);
		do 
		{
			if (!it_edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY))
			{
				HDS_mesh::Face_handle face_temp = it_edge_begin->opposite()->face();
				if (face_temp->testMark(HDS_mesh::Face::TRI))
				{
					HDS_mesh::Vertex_handle v0 = it_edge_begin->vertex();
					HDS_mesh::Vertex_handle v1 = it_edge_begin->next()->vertex();
					HDS_mesh::Vertex_handle v2 = it_edge_begin->opposite()->vertex();
					HDS_mesh::Vertex_handle v3 = it_edge_begin->opposite()->next()->vertex();
					double quality_temp = _optimizer->quad_quality(v0, v1, v2, v3);
					//assert(quality_temp <= 1.0 && quality_temp >= 0.0);
					if (quality_temp > tri_qual.second)
					{
						tri_qual.first = face_temp;
						tri_qual.second = quality_temp;
					}
				}
			}

			it_edge_begin = it_edge_begin->next();
		} while (it_edge_begin != it_edge_end);
		if (tri_qual.second < TOLERENCE)
		{
			continue;
		}
		assert(tri_qual.second > 0.0);

		HDS_mesh::Face_handle new_quad = merge_Tri(tri_face_base, tri_qual.first);
		++n;
	}
	_optimizer->mesh_OPT();
}