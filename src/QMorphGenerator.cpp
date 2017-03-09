//**************************************************************************
// File:     QMorphGenerator.cpp
// Author:   Zhang Jun
// Date:     2016-11-09
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 四边形网格生成算法实现(基于QMorph的算法)
//**************************************************************************
#include "QMorphGenerator.h"
#include "Vector3D.h"
#include "MeshOptimizer.h"

/** \brief 根据边长对前沿边进行排序
*/
static bool sort_edge(HDS_mesh::Halfedge_handle e1, HDS_mesh::Halfedge_handle e2)
{
	return e1->get_SquardLength() < e2->get_SquardLength();
}

/** \brief 使用编号对节点进行排序
*/
struct sort_edge_set : public std::unary_function<HDS_mesh::Halfedge_handle, bool>
{
	bool operator()(const HDS_mesh::Halfedge_handle &e1, const HDS_mesh::Halfedge_handle &e2)
	{
		int h_halfedge1[2]; ///< 存储半边的两个端点的编号
		h_halfedge1[0] = e1->prev()->vertex()->getIndex();
		h_halfedge1[1] = e1->vertex()->getIndex();
		std::sort(h_halfedge1, h_halfedge1+2);
		int k_halfedge1 = 6151 * h_halfedge1[0] + h_halfedge1[1];

		int h_halfedge2[2]; ///< 存储半边的两个端点的编号
		h_halfedge2[0] = e2->prev()->vertex()->getIndex();
		h_halfedge2[1] = e2->vertex()->getIndex();
		std::sort(h_halfedge2, h_halfedge2+2);
		int k_halfedge2 = 6151 * h_halfedge2[0] + h_halfedge2[1];
		return k_halfedge1 < k_halfedge2;
	}
};

/** \brief 计算向量A与向量B夹角的余弦值
		\param p A向量的起点
		\param q A向量的终点
		\param r B向量的起点
		\param s B向量的终点
*/
// static double angle_evaluate2(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
// {
// 	const double dx1 = q.x() - p.x();
// 	const double dx2 = s.x() - r.x();
// 
// 	const double dy1 = q.y() - p.y();
// 	const double dy2 = s.y() - r.y();
// 
// 	const double dz1 = q.z() - p.z();
// 	const double dz2 = s.z() - r.z();
// 
// 	double ans = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1) * std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2));
// 
// 	if (ans >= 0)
// 		return ans < 1 ? ans : 1;
// 	else
// 		return ans < -1 ? -1 : ans;
// }

/** \brief 计算两共点的向量qp和qr的夹角
	    \return 返回值单位为度
	*/
//static double angle_evaluate(const Point_3 &p, const Point_3 &q, const Point_3 &r)
// {
// 	double cosAngle = angle_evaluate2(q, p, q, r);
// 	double angle_temp = std::acos(cosAngle) * 180.0 / PI;
// 
// 	CGAL::Orientation ori = CGAL::orientation(Point_2(p.x(), p.y()), Point_2(q.x(), q.y()), Point_2(r.x(), r.y()));
// 	if (ori == CGAL::RIGHT_TURN)
// 		angle_temp = 360 - angle_temp;
// 	return angle_temp;
//
static double Epsilon()
{
	double r;

	r = 1.0E+00;

	while ( 1.0E+00 < ( double ) ( 1.0E+00 + r )  )
	{
		r = r / 2.0E+00;
	}

	return ( 2.0E+00 * r );
}

/** \brief 计算直线p1q1与直线p2q2的交点
*/
static bool get_cross(const Point_3 &p1, const Point_3 &q1, 
					  const Point_3 &p2, const Point_3 &q2,
					  Point_3& result)
{
	Line_3 L1(p1, q1);
	Line_3 L2(p2, q2);
	auto v = intersection(L1, L2);
	// 	cpp11::result_of<K::Intersect_3(Line_3, Line_3)>::type
	// 		v = intersection(L1, L1); 
	if(v) 
	{
		if (const Point_3 *p = boost::get<Point_3>(&*v) ) 
		{
			result = *p;
			return true;
		} 
		else 
		{
			return false;
			assert(false);
		}
	} 
	else 
	{
		return false;
		assert(false);
	}
	// 	Point_2 p11(p1.getX(), p1.getY());
	// 	Point_2 p22(q1.getX(), q1.getY());
	// 	Point_2 p33(p2.getX(), p2.getY());
	// 	Point_2 p44(q2.getX(), q2.getY());
	// 	Intersect_Segement_2 segments[] = { Intersect_Segement_2(p11, p22),
	// 		                                Intersect_Segement_2(p33, p44)};
	// 	std::list<Point_2> pts;

	// 	CGAL::compute_intersection_points(segments, segments + 2, std::back_inserter(pts));
}

QMorphGenerator::QMorphGenerator()
{
	//_dimension = 3; ///< 默认维数为三维
	_source_mesh = nullptr;
	_optimizer = new MeshOptimizer;
	_deta = 0.0;
	_max_index = 0;
	_reatio = 0.0;
}

QMorphGenerator::~QMorphGenerator()
{
	delete _optimizer;
}

void QMorphGenerator::set_SourceMesh(HDS_mesh *source_mesh)
{
	_source_mesh = source_mesh;
	_optimizer->set_SourceMesh(source_mesh);
	
	//_active_front_list = active_front;
	//_unactive_front_list = unactive_front;

	double _min_len = DBL_MAX;
	double _max_len = DBL_MIN;
	double min_x = DBL_MAX;
	double min_y = DBL_MAX;
	double min_z = DBL_MAX;

	double max_x = DBL_MIN;
	double max_y = DBL_MIN;
	double max_z = DBL_MIN;

	HDS_mesh::Halfedge_iterator  it_edge;
	for (it_edge = _source_mesh->halfedges_begin(); it_edge != _source_mesh->halfedges_end(); ++it_edge)
	{
		if(it_edge->is_border())
		{
			Point_3 p_temp = it_edge->vertex()->point();
			if (min_x > p_temp.x())
				min_x = p_temp.x();
			if (max_x < p_temp.x())
				max_x = p_temp.x();
			if (min_y > p_temp.y())
				min_y = p_temp.y();
			if (max_y < p_temp.y())
				max_y = p_temp.y();
			if (min_z > p_temp.z())
				min_z = p_temp.z();
			if (max_z < p_temp.z())
				max_z = p_temp.z();
			it_edge->opposite()->vertex()->setMark(HDS_mesh::Vertex::BOUNDARY);///< 标记节点为边界节点
			it_edge->opposite()->vertex()->setMark(HDS_mesh::Vertex::FrontNode);///< 标记节点为前沿点
			set_Front_Node_State(it_edge->opposite()->vertex()); ///< 设置前沿点的状态值
			if (it_edge->opposite()->vertex()->get_angle() > 120 
				&& it_edge->opposite()->vertex()->get_angle() <= 240)
				it_edge->opposite()->vertex()->setMark(HDS_mesh::Vertex::SIDENODE);
			it_edge->opposite()->setMark(HDS_mesh::Halfedge::BOUNDARY);///< 标记半边为边界边
			it_edge->opposite()->setMark(HDS_mesh::Halfedge::FRONT);///< 把边界边标记为初始前沿
			it_edge->opposite()->vertex()->set_halfedge(it_edge->opposite());
			it_edge->set_levelNum(0); ///< 初始前沿层数为0
			double len_temp = CGAL::squared_distance(it_edge->vertex()->point(), it_edge->opposite()->vertex()->point());
			it_edge->opposite()->set_SquardLength(len_temp);
			if (len_temp > _max_len)
				_max_len = len_temp;
			if (len_temp < _min_len)
				_min_len = len_temp;
			_active_front_list.push_back(it_edge->opposite());
		}
	}
	_deta = (1e-5) * std::sqrt((max_x - min_x) * (max_x - min_x)
		                      +(max_y - min_y) * (max_y - min_y)
		                      +(max_z - min_z) * (max_z - min_z));
	std::list<HDS_mesh::Halfedge_handle>::iterator it_front;
	for (it_front = _active_front_list.begin(); it_front != _active_front_list.end(); ++it_front)
	{
		HDS_mesh::Halfedge_handle prv_front = (*it_front)->prev()->vertex()->halfedge();
		assert(prv_front->testMark(HDS_mesh::Halfedge::FRONT));
		(*it_front)->set_loopPrev(prv_front);
		prv_front->set_loopNext(*it_front);
	}
	assert(_min_len > 0);
	_reatio = std::sqrt(_max_len * 1.0 / _min_len);
	_max_index = _source_mesh->size_of_vertices() - 1;
}

void QMorphGenerator::classify_Front_Edge()
{
	std::vector<HDS_mesh::Halfedge_handle> front_edge_11_temp; ///< 状态值为11的前沿边链表
	std::vector<HDS_mesh::Halfedge_handle> front_edge_10_temp; ///< 状态值为10的前沿边链表
	std::vector<HDS_mesh::Halfedge_handle> front_edge_01_temp; ///< 状态值为01的前沿边链表
	std::vector<HDS_mesh::Halfedge_handle> front_edge_00_temp; ///< 状态值为00的前沿边链表

	std::list<HDS_mesh::Halfedge_handle>::iterator it_front;
	for (it_front = _active_front_list.begin(); it_front != _active_front_list.end(); ++it_front)
	{
		(*it_front)->vertex()->set_halfedge((*it_front));
		if ((*it_front)->prev()->vertex()->testMark(HDS_mesh::Vertex::STATE))
		{
			if ((*it_front)->vertex()->testMark(HDS_mesh::Vertex::STATE))
			{
				front_edge_11_temp.push_back(*it_front);
				(*it_front)->set_FrontState(3);
			}
			else
			{
				front_edge_10_temp.push_back(*it_front);
				(*it_front)->set_FrontState(2);
			}
		}
		else
		{
			if ((*it_front)->vertex()->testMark(HDS_mesh::Vertex::STATE))
			{
				front_edge_01_temp.push_back(*it_front);
				(*it_front)->set_FrontState(1);
			}
			else
			{
				front_edge_00_temp.push_back(*it_front);
				(*it_front)->set_FrontState(0);
			}
		}
	}
	///< 对各组前沿边，分别按边长顺序升序排列
	std::sort(front_edge_11_temp.begin(), front_edge_11_temp.end(), sort_edge);
	std::sort(front_edge_10_temp.begin(), front_edge_10_temp.end(), sort_edge);
	std::sort(front_edge_01_temp.begin(), front_edge_01_temp.end(), sort_edge);
	std::sort(front_edge_00_temp.begin(), front_edge_00_temp.end(), sort_edge);

	for (int i = 0; i < front_edge_11_temp.size(); ++i)
		_front_edge_11.push_back(front_edge_11_temp[i]);
	for (int i = 0; i < front_edge_10_temp.size(); ++i)
		_front_edge_10.push_back(front_edge_10_temp[i]);
	for (int i = 0; i < front_edge_01_temp.size(); ++i)
		_front_edge_01.push_back(front_edge_01_temp[i]);
	for (int i = 0; i < front_edge_00_temp.size(); ++i)
		_front_edge_00.push_back(front_edge_00_temp[i]);
}

void QMorphGenerator::set_Front_Node_State(HDS_mesh::Vertex_handle vh)
{
	double angle_temp = 0.0;
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh->vertex_begin();
	do 
	{
		if (!edge_begin->is_border())
		{
			if (edge_begin->face()->testMark(HDS_mesh::Face::TRI))
			{
				double cosAngle = angle_evaluate2(edge_begin->vertex()->point(),
					edge_begin->prev()->vertex()->point(),
					edge_begin->vertex()->point(),
					edge_begin->next()->vertex()->point());
				angle_temp += std::acos(cosAngle) * 180.0 / PI;
			}
		}
		++edge_begin;
	} while (edge_begin != edge_end);
	if (angle_temp < state_tolerence)
		vh->setMark(HDS_mesh::Vertex::STATE);
	else
		vh->resetMark(HDS_mesh::Vertex::STATE);
	vh->set_angle(angle_temp);
}

void QMorphGenerator::updateFront(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list)
{
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> nodeSet; ///< 需要重新计算状态值的节点
	int n_front = 0;
	std::list<HDS_mesh::Halfedge_handle> edge_list_temp; ///< 存储前沿半边
	std::list<HDS_mesh::Halfedge_handle>::iterator it_edge_temp;
	for (it_edge_temp = influenced_front_edge_list.begin(); it_edge_temp != influenced_front_edge_list.end(); ++it_edge_temp)
	{
		HDS_mesh::Halfedge_handle current_edge = *it_edge_temp;
		if (current_edge->testMark(15))
			continue;
		nodeSet.insert(current_edge->vertex());
		nodeSet.insert(current_edge->prev()->vertex());
	}

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_node;
	for (it_node = nodeSet.begin(); it_node != nodeSet.end(); ++ it_node)
	{
		set_Front_Node_State(*it_node);
	}
	for (it_edge_temp = influenced_front_edge_list.begin(); it_edge_temp != influenced_front_edge_list.end(); ++it_edge_temp)
	{ 
		if ((*it_edge_temp)->testMark(HDS_mesh::Halfedge::UNACTIVE))
			continue;
		deleteFromFront(*it_edge_temp); ///< 将半边从相应的前沿半边链表中删除
		if ((*it_edge_temp)->prev()->vertex()->testMark(HDS_mesh::Vertex::STATE))
		{
			if ((*it_edge_temp)->vertex()->testMark(HDS_mesh::Vertex::STATE))
			{
				_front_edge_11.push_back(*it_edge_temp);
				(*it_edge_temp)->set_FrontState(3);
			}
			else
			{
				_front_edge_10.push_back(*it_edge_temp);
				(*it_edge_temp)->set_FrontState(2);
			}
		}
		else
		{
			if ((*it_edge_temp)->vertex()->testMark(HDS_mesh::Vertex::STATE))
			{
				_front_edge_01.push_back(*it_edge_temp);
				(*it_edge_temp)->set_FrontState(1);
			}
			else
			{
				_front_edge_00.push_back(*it_edge_temp);
				(*it_edge_temp)->set_FrontState(0);
			}
		}
	}
}

bool QMorphGenerator::get_Next_Front(HDS_mesh::Halfedge_handle& eh)
{
	bool flag = get_Front_Aux(eh);

	if(flag)
		return true;

	///< 此时需更新活动前沿
	_active_front_list.clear();
	std::list<HDS_mesh::Halfedge_handle>::iterator it_front_temp;
	for (it_front_temp = _unactive_front_list.begin(); 
		it_front_temp != _unactive_front_list.end(); ++it_front_temp)
	{
		if ((*it_front_temp)->testMark(15))
			continue;
		if ((*it_front_temp)->testMark(HDS_mesh::Halfedge::FRONT))
		{
			(*it_front_temp)->resetMark(HDS_mesh::Halfedge::UNACTIVE);
			_active_front_list.push_back((*it_front_temp));
		}
	}
	_unactive_front_list.clear();
	//localReconnect();
	_optimizer->global_smooth(3);
	for (it_front_temp = _active_front_list.begin(); it_front_temp != _active_front_list.end(); ++it_front_temp)
	{
		set_Front_Node_State((*it_front_temp)->vertex());
	}
	classify_Front_Edge();
	if (_active_front_list.size() == 0)
		return false;
	flag = get_Front_Aux(eh);
	return flag;
}

bool QMorphGenerator::get_Front_Aux(HDS_mesh::Halfedge_handle& eh)
{
	if (_front_edge_11.size() != 0)
	{
		eh = _front_edge_11.front();
		_front_edge_11.pop_front();
		return true;
	}
	else if (_front_edge_10.size() != 0)
	{
		eh = _front_edge_10.front();
		_front_edge_10.pop_front();
		return true;
	}
	else if (_front_edge_01.size() != 0)
	{
		eh = _front_edge_01.front();
		_front_edge_01.pop_front();
		return true;
	}
	else if (_front_edge_00.size() != 0)
	{
		eh = _front_edge_00.front();
		_front_edge_00.pop_front();
		return true;
	}

	return false;
}

bool QMorphGenerator::get_Side_Left(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh)
{
	if (front->prev()->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
	{///< 重新计算该点的状态值
		double angle_temp = 0.0;
		HDS_mesh::Halfedge_handle pre_temp = front;
		HDS_mesh::Halfedge_handle next_temp = front->prev();
		do 
		{
			double cosAngle = angle_evaluate2(pre_temp->prev()->vertex()->point(),
				pre_temp->vertex()->point(),
				next_temp->vertex()->point(),
				next_temp->prev()->vertex()->point());
			angle_temp += std::acos(cosAngle) * 180.0 / PI;
			pre_temp = next_temp->opposite();
			next_temp = pre_temp->prev();
		} while (pre_temp->face()->testMark(HDS_mesh::Face::TRI));

		front->prev()->vertex()->set_angle(angle_temp);

		if (angle_temp < state_tolerence)
			front->prev()->vertex()->setMark(HDS_mesh::Vertex::STATE);
		else
			front->prev()->vertex()->resetMark(HDS_mesh::Vertex::STATE);
	}

	if (front->prev()->vertex()->testMark(HDS_mesh::Vertex::STATE))
	{
		eh = front->loopPrev();
		return true;
	}
	HDS_mesh::Vertex_handle v_left = front->prev()->vertex();
	Vector3D vl(v_left->point(), front->loopPrev()->prev()->vertex()->point());
	Vector3D vr(v_left->point(), front->vertex()->point());
	Vector3D vm = bisectVector(vl, vr, front->prev()->vertex()->get_angle());


	std::pair<double, HDS_mesh::Halfedge_handle> edge_candidate(-1.0, front); ///< 候选边
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = v_left->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = v_left->vertex_begin();
	do 
	{///< 从v_left的邻接边表中获取与vm夹角最小的边
		if (!edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY))
		{
			assert(edge_begin->vertex() == v_left);
			double cosAngle = 0.0;
			Vector_3 vtemp(edge_begin->vertex()->point(), edge_begin->prev()->vertex()->point());
			cosAngle = angle_evaluate2(CGAL::ORIGIN, Point_3(vm.x(), vm.y(), vm.z()), 
				CGAL::ORIGIN, Point_3(vtemp.x(), vtemp.y(), vtemp.z()));
			if (cosAngle > edge_candidate.first)
			{
				edge_candidate.first = cosAngle;
				edge_candidate.second = edge_begin;
			}
		}
		++edge_begin;
	} while (edge_begin != edge_end);
	if (edge_candidate.second->testMark(HDS_mesh::Halfedge::BOUNDARY))
		return false;
	double min_angle = std::acos(edge_candidate.first) * 180.0 / PI;

	HDS_mesh::Halfedge_handle edge_swap_split;
	if (min_angle >= side_tolerence)
	{
		assert(!edge_candidate.second->testMark(HDS_mesh::Halfedge::BOUNDARY));
		Point_3 vms_temp = vm.getStart();
		Point_3 vme_temp = vm.getEnd();
		Point_3 p_temp = edge_candidate.second->prev()->vertex()->point();
		CGAL::Orientation ori = CGAL::orientation(Point_2(vms_temp.x(), vms_temp.y()), Point_2(vme_temp.x(), vme_temp.y()), Point_2(p_temp.x(), p_temp.y()));
		switch (ori)
		{
		case CGAL::LEFT_TURN:
			edge_swap_split = edge_candidate.second->prev();
			break;
		case CGAL::RIGHT_TURN:
			edge_swap_split = edge_candidate.second->opposite()->next();
			break;
		default:
			return false;
			break;
		}
	}
	if (min_angle < side_tolerence || edge_swap_split->testMark(HDS_mesh::Halfedge::FRONT))
	{
		eh = edge_candidate.second;
		if (!edge_candidate.second->prev()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
		{
			double leng_eh = CGAL::squared_distance(eh->prev()->vertex()->point(), eh->vertex()->point());
			double leng_left = CGAL::squared_distance(Point_3(0,0,0), Point_3(vl.x(), vl.y(), vl.z()));
			double leng_right = CGAL::squared_distance(Point_3(0,0,0), Point_3(vr.x(), vr.y(), vr.z()));
			if (4 * leng_eh < 3.0 * (leng_left + leng_right + 2 * CGAL::sqrt(leng_left*leng_right)))
				return true;
			else
			{///有错误
				//assert(false);
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->prev()->vertex());
				return true;
			}
		}
		else
		{
			bool new_loop = false;
			int n_sub_loop_left = 0;
			HDS_mesh::Halfedge_handle sub_front_temp = front->loopNext();
			while (sub_front_temp->vertex() != eh->opposite()->vertex())
			{
				++n_sub_loop_left;
				if (sub_front_temp == front)
				{
					eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
					return true;
				}
				sub_front_temp = sub_front_temp->loopNext();
			}
			n_sub_loop_left += 3;
			if (n_sub_loop_left % 2 == 0)
				new_loop = true;

			if (new_loop)
			{
				int n_sub_loop_right = 0;
				sub_front_temp = front->loopPrev();
				while (sub_front_temp->prev()->vertex() != eh->opposite()->vertex())
				{
					++n_sub_loop_right;
					sub_front_temp = sub_front_temp->loopPrev();
				}
				n_sub_loop_right += 2;
				if (n_sub_loop_right % 2 != 0)
					new_loop = false;
			}

			if (new_loop)
			{
				eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
				return true;
			}
			else
			{
				//assert(false);
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->prev()->vertex());
				return true;
			}
		}
		return true;
	}

	Point_3 p_opposite = edge_swap_split->opposite()->next()->vertex()->point();
	double cosAngle_opposite = angle_evaluate2(vm.getStart(), vm.getEnd(), vm.getStart(), p_opposite);
	double length_opposite = CGAL::squared_distance(vm.getStart(), p_opposite);
	double length_left_front = CGAL::squared_distance(vm.getStart(), vl.getEnd());
	double length_right_front = CGAL::squared_distance(vm.getStart(), vr.getEnd());

	if (std::acos(cosAngle_opposite) * 180.0 / PI < state_tolerence
		&& 4 * length_opposite < 3 * (length_left_front + length_right_front + 2 * CGAL::sqrt(length_left_front * length_right_front)))
	{///< 交换对边
		bool flag_temp =  flip_edge(edge_swap_split);
		eh = edge_swap_split;
		if (flag_temp)
		{
			if (eh->vertex() != front->prev()->vertex())
				eh = eh->opposite();
			if (!eh->opposite()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
				return true;
			else
			{
				int n_sub_loop = 0;
				HDS_mesh::Halfedge_handle sub_front_temp = front;
				while (sub_front_temp->vertex() != eh->opposite()->vertex())
				{
					++n_sub_loop;
					sub_front_temp = sub_front_temp->loopNext();
					if (sub_front_temp == front)
					{
						eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
						return true;
					}
				}
				n_sub_loop += 2;
				if (n_sub_loop % 2 == 0)
				{///< 环分裂
					eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
					return true;
				}
				else
				{
					//assert(false);
					Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
					eh = _source_mesh->join_facet(eh);
					eh = _source_mesh->create_center_vertex(eh);
					eh->vertex()->clearMark();
					eh->vertex()->point() = cut_point;
					eh->vertex()->setIndex(_max_index);
					eh = eh->next()->opposite()->next();
					++_max_index;
					assert(eh->vertex() == front->prev()->vertex());
					return true;
				}
			}
		}
	}

	///< 否则分裂
	Point_3 split_point;
	bool isCross = get_cross(vm.getStart(), vm.getEnd(), 
		edge_swap_split->prev()->vertex()->point(), 
		edge_swap_split->vertex()->point(), split_point);
	if (!isCross)
		return false;
	edge_swap_split = _source_mesh->join_facet(edge_swap_split);
	//write_my_mesh_nas("join_facet.nas");
	edge_swap_split = _source_mesh->create_center_vertex(edge_swap_split);
	eh = edge_swap_split;
	eh->vertex()->clearMark();
	eh->vertex()->point() = split_point;
	eh->vertex()->setIndex(_max_index);
	//write_my_mesh_nas("join_facet.nas");
	eh = eh->next();
	++_max_index;
	if (eh->vertex() != front->prev()->vertex())
		return false;
	assert(eh->vertex() == front->prev()->vertex());
	return true;
}

bool QMorphGenerator::get_Side_Right(HDS_mesh::Halfedge_handle& front, HDS_mesh::Halfedge_handle& eh)
{
	if (front->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
	{///< 重新计算该点的状态值
		double angle_temp = 0.0;
		HDS_mesh::Halfedge_handle pre_temp = front;
		HDS_mesh::Halfedge_handle next_temp = front->next();
		do 
		{
			double cosAngle = angle_evaluate2(pre_temp->vertex()->point(),
				pre_temp->prev()->vertex()->point(),
				next_temp->prev()->vertex()->point(),
				next_temp->vertex()->point());
			angle_temp += std::acos(cosAngle) * 180.0 / PI;
			pre_temp = next_temp->opposite();
			next_temp = pre_temp->next();
		} while (pre_temp->face()->testMark(HDS_mesh::Face::TRI));

		if (angle_temp < state_tolerence)
			front->vertex()->setMark(HDS_mesh::Vertex::STATE);
		else
			front->vertex()->resetMark(HDS_mesh::Vertex::STATE);
	}
	if (front->vertex()->testMark(HDS_mesh::Vertex::STATE))
	{
		eh = front->loopNext();
		return true;
	}
	HDS_mesh::Vertex_handle v_right = front->vertex();
	Vector3D vl(v_right->point(), front->prev()->vertex()->point());
	Vector3D vr(v_right->point(), front->loopNext()->vertex()->point());

	Vector3D vm = bisectVector(vl, vr, front->vertex()->get_angle());

	std::pair<double, HDS_mesh::Halfedge_handle> edge_candidate(-1.0, front); ///< 候选边
	std::list<HDS_mesh::Halfedge_handle> edge_valid_list; ///< 有效入射边
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = v_right->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = v_right->vertex_begin();
	do 
	{
		edge_valid_list.push_back(edge_begin);
		if (!edge_begin->testMark(HDS_mesh::Halfedge::FRONT)
			&& !edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY))
		{
			double cosAngle = angle_evaluate2(vm.getStart(), vm.getEnd(), 
				edge_begin->vertex()->point(), edge_begin->prev()->vertex()->point());
			if (cosAngle > edge_candidate.first)
			{
				edge_candidate.first = cosAngle;
				edge_candidate.second = edge_begin;
			}
		}
		++edge_begin;
	} while (edge_begin != edge_end);
	if (edge_valid_list.size() == 2)
	{
		eh = front->loopNext();
		return true;
	}

	double min_angle = std::acos(edge_candidate.first) * 180.0 / PI;
	HDS_mesh::Halfedge_handle edge_swap_split;
	if (min_angle >= side_tolerence)
	{
		assert(!edge_candidate.second->testMark(HDS_mesh::Halfedge::BOUNDARY));
		Point_3 vms_temp = vm.getStart();
		Point_3 vme_temp = vm.getEnd();
		Point_3 p_temp = edge_candidate.second->prev()->vertex()->point();
		CGAL::Orientation ori = CGAL::orientation(Point_2(vms_temp.x(), vms_temp.y()), Point_2(vme_temp.x(), vme_temp.y()), Point_2(p_temp.x(), p_temp.y()));
		switch (ori)
		{
		case CGAL::LEFT_TURN:
			edge_swap_split = edge_candidate.second->prev();
			break;
		case CGAL::RIGHT_TURN:
			edge_swap_split = edge_candidate.second->opposite()->next();
			break;
		default:
			return false;
			break;
		}
		if (edge_swap_split->vertex() == front->prev()->vertex() && edge_swap_split->prev()->vertex()->testMark(HDS_mesh::Vertex::HASCHOICED)
			|| edge_swap_split->prev()->vertex() == front->prev()->vertex() && edge_swap_split->vertex()->testMark(HDS_mesh::Vertex::HASCHOICED))
			return false;
	}
	if ( min_angle < side_tolerence || edge_swap_split->testMark(HDS_mesh::Halfedge::FRONT))
	{
		eh = edge_candidate.second;

		double leng_eh = CGAL::squared_distance(eh->prev()->vertex()->point(), eh->vertex()->point());
		double leng_left = CGAL::squared_distance(vl.getStart(), vl.getEnd());
		double leng_right = CGAL::squared_distance(vr.getStart(), vr.getEnd());

		if (!eh->prev()->vertex()->testMark(HDS_mesh::Vertex::FrontNode) 
			&& !eh->prev()->vertex()->testMark(HDS_mesh::Vertex::HASCHOICED))
		{
			if (eh->testMark(HDS_mesh::Halfedge::QuadEdge))
				return false;
			if (4 * leng_eh < 3.0 * (leng_left + leng_right + 2 * CGAL::sqrt(leng_left * leng_right)))
				return true;
			else
			{
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->vertex());
				return true;
			}
		}
		else
		{
			if (eh->prev()->vertex()->testMark(HDS_mesh::Vertex::HASCHOICED)
				|| (!eh->testMark(HDS_mesh::Halfedge::BOUNDARY) && eh->prev()->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY)))
			{
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->vertex());
				return true;
			}
			else if (4 * leng_eh >= 3.0 * (leng_left + leng_right + 2 * CGAL::sqrt(leng_left * leng_right)))
			{
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->vertex());
				return true;
			}
			else if (eh->opposite()->vertex() == front->loopPrev()->loopPrev()->loopPrev()->vertex())
			{
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->vertex());
				return true;
			}
			else
			{
				bool new_loop = false;
				int n_sub_loop_left = 0;
				HDS_mesh::Halfedge_handle sub_front_temp1 = front->loopNext();
				if (sub_front_temp1->testMark(HDS_mesh::Halfedge::FRONT))
				{
					while (sub_front_temp1->vertex() != eh->opposite()->vertex())
					{
						++n_sub_loop_left;
						if (sub_front_temp1 == front)
						{
							eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
							return true;
						}
						sub_front_temp1 = sub_front_temp1->loopNext();
						//std::cout << sub_front_temp1->vertex()->point() << std::endl;
					}
					n_sub_loop_left += 2;
				}
				else
					n_sub_loop_left = 0;
				if (n_sub_loop_left % 2 == 0)
					new_loop = true;

				if (new_loop)
				{
					int n_sub_loop_right = 0;
					sub_front_temp1 = front->loopPrev();
					while (sub_front_temp1->prev()->vertex() != eh->opposite()->vertex())
					{
						++n_sub_loop_right;
						sub_front_temp1 = sub_front_temp1->loopPrev();
					}
					if (front->get_FrontState() == 2 || front->get_FrontState() == 3)
						n_sub_loop_right += 1;
					else
						n_sub_loop_right += 3;
					if (n_sub_loop_right % 2 != 0)
						new_loop = false;
				}

				if (new_loop)
				{///< 环分裂
					eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
					return true;
				}
				else
				{
					Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
					eh = _source_mesh->join_facet(eh);
					eh = _source_mesh->create_center_vertex(eh);
					eh->vertex()->clearMark();
					eh->vertex()->point() = cut_point;
					eh->vertex()->setIndex(_max_index);
					eh = eh->next()->opposite()->next();
					++_max_index;
					assert(eh->vertex() == front->vertex());
					return true;
				}
			}
		}
	}

	Point_3 p_opposite = edge_swap_split->opposite()->next()->vertex()->point();
	double cosAngle_opposite = angle_evaluate2(vm.getStart(), vm.getEnd(), vm.getStart(), p_opposite);
	double length_opposite = CGAL::squared_distance(vm.getStart(), p_opposite);
	double length_left_front = CGAL::squared_distance(vm.getStart(), vl.getEnd());
	double length_right_front = CGAL::squared_distance(vm.getStart(), vr.getEnd());

	if (std::acos(cosAngle_opposite) * 180.0 / PI < state_tolerence
		&& 4 * length_opposite < 3.0 * (length_left_front + length_right_front + 2 * CGAL::sqrt(length_left_front * length_right_front)))
	{///< 交换对边

		bool flag_temp = flip_edge(edge_swap_split);
		eh = edge_swap_split;
		if (flag_temp)
		{
			if (eh->vertex() != front->vertex())
				eh = eh->opposite();
			if (eh->opposite()->vertex()->testMark(HDS_mesh::Vertex::HASCHOICED))
			{///< 取该边终点将其分裂
				Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
				eh = _source_mesh->join_facet(eh);
				eh = _source_mesh->create_center_vertex(eh);
				eh->vertex()->clearMark();
				eh->vertex()->point() = cut_point;
				eh->vertex()->setIndex(_max_index);
				eh = eh->next()->opposite()->next();
				++_max_index;
				assert(eh->vertex() == front->vertex());
				return true;
			}
			if (!eh->opposite()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
				return true;
			else
			{
				int n_sub_loop = 0;
				HDS_mesh::Halfedge_handle sub_front_temp2 = front->loopNext();
				while (sub_front_temp2->vertex() != eh->opposite()->vertex())
				{
					++n_sub_loop;
					if (sub_front_temp2 == front)
					{
						eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
						return true;
					}
					sub_front_temp2 = sub_front_temp2->loopNext();
				}
				n_sub_loop += 2;
				if (n_sub_loop % 2 == 0)
				{
					eh->opposite()->vertex()->setMark(HDS_mesh::Vertex::LOOPCUT);
					return true;
				}
				else
				{
					Point_3 cut_point = (eh->vertex()->point() + eh->opposite()->vertex()->point()) / 2.0;
					eh = _source_mesh->join_facet(eh);
					eh = _source_mesh->create_center_vertex(eh);
					eh->vertex()->clearMark();
					eh->vertex()->point() = cut_point;
					eh->vertex()->setIndex(_max_index);
					eh = eh->next()->opposite()->next();
					++_max_index;
					assert(eh->vertex() == front->vertex());
					return true;
				}
			}
		}
	}

	///< 否则分裂

	Point_3 split_point;
	bool isCross = get_cross(vm.getStart(), vm.getEnd(), 
		edge_swap_split->prev()->vertex()->point(), 
		edge_swap_split->vertex()->point(),
		split_point);
	if (!isCross)
		return false;
	edge_swap_split = _source_mesh->join_facet(edge_swap_split);
	//write_my_mesh_nas("join_facet.nas");
	edge_swap_split = _source_mesh->create_center_vertex(edge_swap_split);
	eh = edge_swap_split;
	eh->vertex()->clearMark();
	eh->vertex()->point() = split_point;
	//write_my_mesh_nas("join_facet.nas");
	eh->vertex()->setIndex(_max_index);
	eh = eh->next();
	++_max_index;
	assert(eh->vertex() == front->vertex());
	return true;
}

HDS_mesh::Halfedge_handle  QMorphGenerator::get_Top(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr)
{
	HDS_mesh::Halfedge_handle it_edge_temp = HDS_mesh::Halfedge_handle();
	bool flag = recover_edge(vl, vr, it_edge_temp);
	return flag ? it_edge_temp : HDS_mesh::Halfedge_handle();
	//return it_edge_temp;
}

HDS_mesh::Face_handle QMorphGenerator::constructNewFace(HDS_mesh::Halfedge_handle e_down, HDS_mesh::Halfedge_handle e_right, 
													    HDS_mesh::Halfedge_handle e_up,   HDS_mesh::Halfedge_handle e_left)
{
	///< 首先确保四条边绕向为逆时针
	if (e_right->vertex() == e_down->vertex())
		e_right = e_right->opposite();
	if (e_up->vertex() == e_right->vertex())
		e_up = e_up->opposite();
	if (e_left->vertex() != e_down->prev()->vertex())
		e_left = e_left->opposite();
	if (e_down->vertex()->halfedge() != e_down)
		e_down->vertex()->set_halfedge(e_down);

	e_down->setMark(HDS_mesh::Halfedge::QuadEdge);
	e_down->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e_down->vertex()->setMark(HDS_mesh::Vertex::FrontNode);

	e_right->setMark(HDS_mesh::Halfedge::QuadEdge);
	e_right->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e_right->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
	//e_right->vertex()->set_halfedge(e_right);

	e_up->setMark(HDS_mesh::Halfedge::QuadEdge);
	e_up->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e_up->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
	//e_up->vertex()->set_halfedge(e_up);

	e_left->setMark(HDS_mesh::Halfedge::QuadEdge);
	e_left->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
	e_left->vertex()->setMark(HDS_mesh::Vertex::FrontNode);

	HDS_mesh::Face_handle face_base_temp = e_down->face();
	std::list<HDS_mesh::Vertex_handle> old_vertex_list; ///< 待删除节点链表
	std::list<HDS_mesh::Face_handle> old_face_list; ///< 待删除面片链表
	std::list<HDS_mesh::Halfedge_handle> old_edge_list; ///< 待删除边链表
	getFaceToDelete(face_base_temp, old_face_list, old_edge_list, old_vertex_list);

	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());

	face_base_temp->resetMark(HDS_mesh::Face::TRI);
	face_base_temp->setMark(HDS_mesh::Face::QUAD);
	decorator.set_face_halfedge(face_base_temp, e_down);

	std::list<HDS_mesh::Face_handle>::iterator it_old_f;
	for (it_old_f = old_face_list.begin(); it_old_f != old_face_list.end(); ++it_old_f)
	{
		if (it_old_f == old_face_list.begin())
			continue;
		if ((*it_old_f)->testMark(15))
			continue;
		_source_mesh->hds().faces_erase(*it_old_f);
	}
	std::list<HDS_mesh::Halfedge_handle>::iterator it_old_e;
	for (it_old_e = old_edge_list.begin(); it_old_e != old_edge_list.end(); ++it_old_e)
	{///< 删除公共边
		_source_mesh->hds().edges_erase(*it_old_e);
	}
	std::list<HDS_mesh::Vertex_handle>::iterator it_old_v;
	for (it_old_v = old_vertex_list.begin(); it_old_v != old_vertex_list.end(); ++it_old_v)
	{
		_source_mesh->hds().vertices_erase(*it_old_v);			
	}

	///< 设置除e_down边之外的三条边的入射面片
	decorator.set_face(e_right, face_base_temp);
	decorator.set_face(e_up, face_base_temp);
	decorator.set_face(e_left, face_base_temp);

	///< 更新四条边的拓扑关系
	e_down->set_next(e_right);  e_down->set_prev(e_left);
	e_right->set_next(e_up);    e_right->set_prev(e_down);
	e_up->set_next(e_left);     e_up->set_prev(e_right);
	e_left->set_next(e_down);   e_left->set_prev(e_up);

	///< 设置节点入射半边
	decorator.set_vertex_halfedge(e_left);
	decorator.set_vertex_halfedge(e_right);
	decorator.set_vertex_halfedge(e_up);

	///< 更新前沿拓扑关系

	HDS_mesh::Halfedge_handle right_in = HDS_mesh::Halfedge_handle();
	HDS_mesh::Halfedge_handle right_out = HDS_mesh::Halfedge_handle();
	if (e_right->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
	{
		HDS_mesh::Halfedge_handle sub_front_temp = HDS_mesh::Halfedge_handle();
		std::vector<HDS_mesh::Halfedge_handle> sub_front_vector;
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = e_right->vertex()->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = e_right->vertex()->vertex_begin();
		do 
		{
			if (!edge_begin->is_border() && edge_begin->testMark(HDS_mesh::Halfedge::FRONT))
			{
				sub_front_vector.push_back(edge_begin);
			}
			++edge_begin;
		} while (edge_begin != edge_end);
		if (sub_front_vector.size() == 0)
		{
			e_right->vertex()->resetMark(HDS_mesh::Vertex::LOOPCUT);
		}
		else if (sub_front_vector.size() > 1)
		{
			int i = 0;
			for (; i < sub_front_vector.size(); ++i)
			{
				Point_2 ps_temp(sub_front_vector[i]->prev()->vertex()->point().x(), sub_front_vector[i]->prev()->vertex()->point().y());
				Point_2 pt_temp(sub_front_vector[i]->vertex()->point().x(), sub_front_vector[i]->vertex()->point().y());
				Point_2 p_right_temp(e_left->vertex()->point().x(), e_left->vertex()->point().y());
				CGAL::Orientation ori_temp = CGAL::orientation(ps_temp, pt_temp, p_right_temp);
				if (ori_temp == CGAL::LEFT_TURN)
				{///< 位于相应前沿的左侧
					right_in = sub_front_vector[i];
					right_out = sub_front_vector[i]->loopNext();
					break;
				}
			}
			if (i == sub_front_vector.size())
				e_left->prev()->vertex()->resetMark(HDS_mesh::Vertex::LOOPCUT);
		}
		else
		{
			assert(sub_front_vector.size() == 1);
			right_in = sub_front_vector[0];
			right_out = sub_front_vector[0]->loopNext();
		}
	}

	HDS_mesh::Halfedge_handle left_in = HDS_mesh::Halfedge_handle();
	HDS_mesh::Halfedge_handle left_out = HDS_mesh::Halfedge_handle();
	if (e_left->prev()->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
	{
		HDS_mesh::Halfedge_handle sub_front_temp = HDS_mesh::Halfedge_handle();
		std::vector<HDS_mesh::Halfedge_handle> sub_front_vector;
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = e_left->prev()->vertex()->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = e_left->prev()->vertex()->vertex_begin();
		do 
		{
			if (!edge_begin->is_border() && edge_begin->testMark(HDS_mesh::Halfedge::FRONT))
			{
				sub_front_vector.push_back(edge_begin);
			}
			++edge_begin;
		} while (edge_begin != edge_end);
		if (sub_front_vector.size() == 0)
		{
			e_left->prev()->vertex()->resetMark(HDS_mesh::Vertex::LOOPCUT);
		}
		else if (sub_front_vector.size() > 1)
		{
			int i = 0;
			for (; i < sub_front_vector.size(); ++i)
			{
				Point_2 ps_temp(sub_front_vector[i]->prev()->vertex()->point().x(), sub_front_vector[i]->prev()->vertex()->point().y());
				Point_2 pt_temp(sub_front_vector[i]->vertex()->point().x(), sub_front_vector[i]->vertex()->point().y());
				Point_2 p_left_temp(e_left->vertex()->point().x(), e_left->vertex()->point().y());
				CGAL::Orientation ori_temp = CGAL::orientation(ps_temp, pt_temp, p_left_temp);
				if (ori_temp == CGAL::LEFT_TURN)
				{///< 位于相应前沿的左侧
					left_in = sub_front_vector[i];
					left_out = sub_front_vector[i]->loopNext();
					break;
				}
			}
			if (i == sub_front_vector.size())
				e_left->prev()->vertex()->resetMark(HDS_mesh::Vertex::LOOPCUT);
		}
		else
		{
			assert(sub_front_vector.size() == 1);
			left_in = sub_front_vector[0];
			left_out = sub_front_vector[0]->loopNext();
		}
	}

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
		if (e_right->loopNext() != e_up)
			e_right->loopNext()->set_loopPrev(e_up->opposite());
		if (!e_up->testMark(HDS_mesh::Halfedge::BOUNDARY) 
			&& e_up->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
			e_up->opposite()->set_loopNext(e_right->loopNext());

		if (e_right->testMark(HDS_mesh::Halfedge::UNACTIVE))
			deleteFromUnactiveFront(e_right);
		else
			deleteFromFront(e_right);
	}
	else if(!e_right->testMark(HDS_mesh::Halfedge::BOUNDARY))
	{
		if (e_right->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
		{
			e_right->opposite()->resetMark(HDS_mesh::Halfedge::FRONT);
			e_right->opposite()->loopNext()->set_loopPrev(e_up->opposite());
			e_up->opposite()->set_loopNext(e_right->opposite()->loopNext());

			if (e_right->opposite()->testMark(HDS_mesh::Halfedge::UNACTIVE))
				deleteFromUnactiveFront(e_right->opposite());
			else
				deleteFromFront(e_right->opposite());
		}
		else
		{
			e_right->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			e_right->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
			e_right->opposite()->set_loopNext(e_down->loopNext());
			e_down->loopNext()->set_loopPrev(e_right->opposite());

			e_right->opposite()->set_loopPrev(e_up->opposite());
			e_up->opposite()->set_loopNext(e_right->opposite());
			_unactive_front_list.push_back(e_right->opposite());
		}
	}

	if (e_left->testMark(HDS_mesh::Halfedge::FRONT))
	{
		e_left->resetMark(HDS_mesh::Halfedge::FRONT);
		if (e_left->loopPrev() != e_up)
		{
			e_left->loopPrev()->set_loopNext(e_up->opposite());
			e_up->opposite()->set_loopPrev(e_left->loopPrev());
		}
		if (e_left->testMark(HDS_mesh::Halfedge::UNACTIVE))
			deleteFromUnactiveFront(e_left);
		else
			deleteFromFront(e_left);
	}
	else if (!e_left->testMark(HDS_mesh::Halfedge::BOUNDARY))
	{
		if (e_left->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
		{
			e_left->opposite()->resetMark(HDS_mesh::Halfedge::FRONT);
			e_left->opposite()->loopPrev()->set_loopNext(e_up->opposite());
			e_up->opposite()->set_loopPrev(e_left->opposite()->loopPrev());

			if (e_left->opposite()->testMark(HDS_mesh::Halfedge::UNACTIVE))
				deleteFromUnactiveFront(e_left->opposite());
			else
				deleteFromFront(e_left->opposite());
		}
		else
		{
			e_left->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			e_left->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
			e_left->opposite()->set_loopPrev(e_down->loopPrev());
			e_down->loopPrev()->set_loopNext(e_left->opposite());

			e_left->opposite()->set_loopNext(e_up->opposite());
			e_up->opposite()->set_loopPrev(e_left->opposite());

			_unactive_front_list.push_back(e_left->opposite());
		}
	}

	// 	HDS_mesh::Halfedge_handle it_edge_begin = face_base_temp->halfedge();
	// 	HDS_mesh::Halfedge_handle it_edge_end = face_base_temp->halfedge();
	// 	do 
	// 	{
	// 		HDS_mesh::Vertex_handle vh_temp = it_edge_begin->vertex();
	// 		vh_temp->resetMark(HDS_mesh::Vertex::FrontNode);
	// 		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh_temp->vertex_begin();
	// 		HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh_temp->vertex_begin();
	// 		do 
	// 		{
	// 			if (edge_begin->testMark(HDS_mesh::Halfedge::FRONT))
	// 			{
	// 				vh_temp->setMark(HDS_mesh::Vertex::FrontNode);
	// 				break;
	// 			}
	// 			++edge_begin;
	// 		} while (edge_begin != edge_end);
	// 	} while (it_edge_begin != it_edge_end);

	if (e_right->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
	{
		HDS_mesh::Halfedge_handle loop_next1 = right_out;
		HDS_mesh::Halfedge_handle loop_prev1 = e_up->opposite();
		if (right_in ==  e_right)
		{
			right_in->set_loopNext(e_right->loopNext());
		}
		else
		{
			right_in->set_loopNext(e_right->opposite());
			e_right->opposite()->set_loopPrev(right_in);
			//e_down->set_loopNext(e_right);
			e_right->set_loopNext(loop_next1);
		}

		if (loop_next1->testMark(HDS_mesh::Halfedge::FRONT) 
			&& loop_next1 != e_right->opposite())
		{
			loop_next1->set_loopPrev(loop_prev1);
			loop_prev1->set_loopNext(loop_next1);
		}
	}

	if (e_left->prev()->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT) 
		&& e_left->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
	{
		// 		HDS_mesh::Halfedge_handle sub_front_temp = e_down->loopPrev();
		// 		while (sub_front_temp->vertex() != e_left->prev()->vertex())
		// 		{
		// 			sub_front_temp = sub_front_temp->loopPrev();
		// 		}
		HDS_mesh::Halfedge_handle loop_next1 = left_out;
		if (!e_up->testMark(HDS_mesh::Halfedge::BOUNDARY))
		{
			if (e_up->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
			{
				left_in->set_loopNext(e_up->opposite());
				e_up->opposite()->set_loopPrev(left_in);
			}
		}
		if (loop_next1 != e_left->opposite())
		{
			loop_next1->set_loopPrev(e_left->opposite());
			e_left->opposite()->set_loopNext(loop_next1);
		}
	}

	return face_base_temp;
}

void QMorphGenerator::deleteFromFront(HDS_mesh::Halfedge_handle eh)
{
	std::list<HDS_mesh::Halfedge_handle>::iterator front_temp;
	switch (eh->get_FrontState())
	{
	case 0:
		front_temp = _front_edge_00.begin();
		while (front_temp != _front_edge_00.end() && (*front_temp) != eh)
		{
			++front_temp;
		}
		if (front_temp != _front_edge_00.end())
			_front_edge_00.erase(front_temp);
		break;
	case 1:
		front_temp = _front_edge_01.begin();
		while ((*front_temp) != eh)
		{
			++front_temp;
		}
		_front_edge_01.erase(front_temp);
		break;
	case 2:
		front_temp = _front_edge_10.begin();
		while ((*front_temp) != eh)
		{
			++front_temp;
		}
		_front_edge_10.erase(front_temp);
		break;
	case 3:
		if (!_front_edge_11.empty())
		{
			front_temp = _front_edge_11.begin();
			while ((*front_temp) != eh)
			{
				++front_temp;
				if (front_temp == _front_edge_11.end())
					break;
			}
			if (front_temp != _front_edge_11.end())
				_front_edge_11.erase(front_temp);
		}
		break;
	default:
		break;
	}
}

void QMorphGenerator::deleteFromUnactiveFront(HDS_mesh::Halfedge_handle eh)
{
	std::list<HDS_mesh::Halfedge_iterator>::iterator it_edge_temp;
	for (it_edge_temp = _unactive_front_list.begin(); it_edge_temp != _unactive_front_list.end(); ++it_edge_temp)
	{
		if ((*it_edge_temp)->vertex() == eh->vertex() && (*it_edge_temp)->prev()->vertex() == eh->prev()->vertex())
			break;
	}
	_unactive_front_list.erase(it_edge_temp);
}

void QMorphGenerator::getFaceToDelete(HDS_mesh::Face_handle face_base, 
								   std::list<HDS_mesh::Face_handle>& old_face_list,
								   std::list<HDS_mesh::Halfedge_handle>& old_edge_list,
								   std::list<HDS_mesh::Vertex_handle>& old_vertex_list)
{
	old_face_list.clear();
	old_edge_list.clear();
	old_vertex_list.clear();

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> old_vertex_set;

	old_face_list.push_back(face_base);
	face_base->setMark(HDS_mesh::Face::VISITED); ///< 标记面片已被访问
	std::queue<HDS_mesh::Face_handle> face_queue; ///< 待删除面片队列
	HDS_mesh::Halfedge_handle it_edge_begin = face_base->halfedge();
	HDS_mesh::Halfedge_handle it_edge_end = face_base->halfedge();
	do 
	{///< 先把face_base面片的待删除邻接面片压入队列
		if (!it_edge_begin->testMark(HDS_mesh::Halfedge::QuadEdge))
		{
			if (!it_edge_begin->vertex()->testMark(HDS_mesh::Vertex::QuadNode))
				old_vertex_set.insert(it_edge_begin->vertex());
			if (!it_edge_begin->prev()->vertex()->testMark(HDS_mesh::Vertex::QuadNode))
				old_vertex_set.insert(it_edge_begin->prev()->vertex());

			it_edge_begin->opposite()->face()->setMark(HDS_mesh::Face::VISITED);
			face_queue.push(it_edge_begin->opposite()->face());
		}
		it_edge_begin = it_edge_begin->next();
	} while (it_edge_begin != it_edge_end);

	while (!face_queue.empty())
	{
		HDS_mesh::Face_handle face_temp = face_queue.front();
		face_queue.pop();
		old_face_list.push_back(face_temp);
		it_edge_begin = face_temp->halfedge();
		it_edge_end = face_temp->halfedge();
		do 
		{///< 遍历face_temp的邻接单元，将其满足要求的邻接面片压入队列
			if (!it_edge_begin->testMark(HDS_mesh::Halfedge::QuadEdge))
			{
				if (!it_edge_begin->isDead())
				{
					it_edge_begin->setDead();
					it_edge_begin->opposite()->setDead();
					old_edge_list.push_back(it_edge_begin);

					if (!it_edge_begin->vertex()->testMark(HDS_mesh::Vertex::QuadNode))
						old_vertex_set.insert(it_edge_begin->vertex());
					if (!it_edge_begin->prev()->vertex()->testMark(HDS_mesh::Vertex::QuadNode))
						old_vertex_set.insert(it_edge_begin->prev()->vertex());
				}
				if (!it_edge_begin->opposite()->face()->testMark(HDS_mesh::Face::VISITED))
				{
					it_edge_begin->opposite()->face()->setMark(HDS_mesh::Face::VISITED);
					face_queue.push(it_edge_begin->opposite()->face());
				}
			}
			it_edge_begin = it_edge_begin->next();
		} while (it_edge_begin != it_edge_end);
	}

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_v;
	for (it_v = old_vertex_set.begin(); it_v != old_vertex_set.end(); ++it_v)
	{
		old_vertex_list.push_back(*it_v);		
	}
	old_vertex_set.clear();

	std::list<HDS_mesh::Face_handle>::iterator it_face_temp;
	for (it_face_temp = old_face_list.begin(); it_face_temp != old_face_list.end(); ++it_face_temp)
		(*it_face_temp)->resetMark(HDS_mesh::Face::VISITED);
}

void QMorphGenerator::getVertexToSmooth(HDS_mesh::Face_handle quad_face, 
									    std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list)
{
	vertex_smooth_list.clear();

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> v_set;
	HDS_mesh::Halfedge_handle it_edge_begin = quad_face->halfedge();
	HDS_mesh::Halfedge_handle it_edge_end = quad_face->halfedge();
	do 
	{///< 分别获取四边形四个顶点的第一级邻接点
		bool is_front = false;
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = it_edge_begin->vertex()->vertex_begin();	
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = it_edge_begin->vertex()->vertex_begin();
		do 
		{
			if (!is_front)
			{
				if (edge_begin->testMark(HDS_mesh::Halfedge::FRONT)
					|| edge_begin->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
					is_front = true;
			}
			v_set.insert(edge_begin->prev()->vertex());
			++edge_begin;
		} while (edge_begin != edge_end);
		if (!is_front)
			it_edge_begin->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
		it_edge_begin = it_edge_begin->next();
	} while (it_edge_begin != it_edge_end);

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_vs_temp;
	for (it_vs_temp = v_set.begin(); it_vs_temp != v_set.end(); ++it_vs_temp)
	{
		vertex_smooth_list.push_back(*it_vs_temp);
	}
}

void QMorphGenerator::getInfluencedFrontEdge(std::list<HDS_mesh::Vertex_handle>& vertex_smooth_list,
										     std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list)
{
	influenced_front_edge_list.clear();
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> v_set;

	std::list<HDS_mesh::Vertex_handle>::iterator it_v_temp;
	for (it_v_temp = vertex_smooth_list.begin(); it_v_temp != vertex_smooth_list.end(); ++it_v_temp)
	{
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = (*it_v_temp)->vertex_begin();	
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = (*it_v_temp)->vertex_begin();
		do 
		{
			if (edge_begin->prev()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
				v_set.insert(edge_begin->prev()->vertex());
			++edge_begin;
		} while (edge_begin != edge_end);
	}
	std::list<HDS_mesh::Vertex_handle> influenced_vertex_list;
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_vs_temp;
	for (it_vs_temp = v_set.begin(); it_vs_temp != v_set.end(); ++it_vs_temp)
	{
		influenced_vertex_list.push_back(*it_vs_temp);
	}

	std::set<HDS_mesh::Halfedge_handle, sort_edge_set> influenced_edge_set; ///< 收集由于节点光顺而受影响的边(用于更新前沿)

	for (it_v_temp = influenced_vertex_list.begin(); it_v_temp != influenced_vertex_list.end(); ++it_v_temp)
	{///< 从所有被光顺节点的邻接边中获取需要进行缝合检测的前沿边
		if (!(*it_v_temp)->testMark(HDS_mesh::Vertex::FrontNode))
			continue;

		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = (*it_v_temp)->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = (*it_v_temp)->vertex_begin();
		do 
		{
			//assert(!edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY));
			if (!edge_begin->testMark(HDS_mesh::Halfedge::FRONT))
			{
				if (edge_begin->opposite()->testMark(HDS_mesh::Halfedge::FRONT))
					influenced_edge_set.insert(edge_begin->opposite());
				else 
				{
					++edge_begin;
					continue;
				}
			}
			influenced_edge_set.insert(edge_begin);
			++edge_begin;
		} while (edge_begin != edge_end);
	}
	std::set<HDS_mesh::Halfedge_handle, sort_edge_set>::iterator it_es_temp;
	for (it_es_temp = influenced_edge_set.begin(); it_es_temp != influenced_edge_set.end(); ++it_es_temp)
	{
		influenced_front_edge_list.push_back((*it_es_temp));
	}
}

bool QMorphGenerator::MinMaxSwap2D(Point_3& pt0,Point_3& pt1,
								   Point_3& pt2, Point_3& pt3)
{
	double wesp=0.001;
	Vector3D v01(pt0,pt1);
	Vector3D v02(pt0,pt2);
	Vector3D v03(pt0,pt3);

	double area032= v03[0]*v02[1]-v03[1]*v02[0];
	if(area032 < Epsilon()) return false;

	Vector3D v12(pt1,pt2);
	Vector3D v13(pt1,pt3);
	double area123= v12[0]*v13[1]-v12[1]*v13[0];
	if(area123 < Epsilon()) return false;

	v02.normal(); v12.normal();
	v13.normal(); v03.normal();

	double w012 = v02 * v12; ///< 两向量夹角的余弦值
	double w031 = v13 * v03;
	double w0 = w012 < w031 ? w012 : w031; ///< 取余弦值最小(夹角最大)的最为w0

	double w032=v03 * v02;
	double w123=v12 * v13;
	double w1 = w032 < w123 ? w032 : w123;
	double dw=w1-w0;
	return dw>=wesp; ///< 判断边交换前后，最大角是否减小
}

bool QMorphGenerator::flip_edge(HDS_mesh::Halfedge_handle& eh, bool check_angle)
{
	if (eh->testMark(HDS_mesh::Halfedge::BOUNDARY))
		return false;

	HDS_mesh::Halfedge_handle eho = eh->opposite();

	if (!eh->face()->testMark(HDS_mesh::Face::TRI)
		|| !eho->face()->testMark(HDS_mesh::Face::TRI))
		return false;///< 待翻转煸两侧均应是三角形

	HDS_mesh::Halfedge_handle eh_n = eh->next();
	HDS_mesh::Halfedge_handle eh_p = eh->prev();

	HDS_mesh::Halfedge_handle eho_n = eho->next();
	HDS_mesh::Halfedge_handle eho_p = eho->prev();

	HDS_mesh::Vertex_handle   va0 = eh->vertex();
	HDS_mesh::Vertex_handle   va1 = eho->vertex();

	HDS_mesh::Vertex_handle   vb0 = eh_n->vertex();
	HDS_mesh::Vertex_handle   vb1 = eho_n->vertex();

	Point_3 pa0 = va0->point();
	Point_3 pa1 = va1->point();
	Point_3 pb0 = vb0->point();
	Point_3 pb1 = vb1->point();

	if (CGAL::collinear(pa0, pb0, pb1) || CGAL::collinear(pa1, pb1, pb0))
	{
		return false;
	}

	// 	std::cout << CGAL::normal(pa0, pb0, pb1) << std::endl;
	// 	std::cout << CGAL::normal(pa1, pb1, pb0) << std::endl;
	// 	if (CGAL::normal(pa0, pb0, pb1).z() < 1.0e-15)
	// 		return false;
	// 	if (CGAL::normal(pa1, pb1, pb0).z() < 1.0e-15)
	// 		return false;

	if (CGAL::normal(pa0, pb0, pb1) 
		* CGAL::normal(pa1, pb1, pb0) <= 0)
		return false;///< 两面片法向需一致
	// 	CGAL::Orientation ori1 = CGAL::orientation(Point_2(pa0.x(), pa0.y()), Point_2(pb0.x(), pb0.x()), Point_2(pb1.x(), pb1.y()));
	// 	CGAL::Orientation ori2 = CGAL::orientation(Point_2(pa1.x(), pa1.y()), Point_2(pb1.x(), pb1.y()), Point_2(pb0.x(), pb0.y()));
	// 
	// 	if (ori1 == ori2)
	// 		return false;
	if (check_angle)
	{
		if (!MinMaxSwap2D(pa0, pa1, pb1, pb0))
			return false;
	}

	HDS_mesh::Face_handle     f_eh  = eh->face();
	HDS_mesh::Face_handle     f_eho  = eho->face();

	eh->set_vertex(vb0);
	eho->set_vertex(vb1);

	eh->set_next(eh_p);
	eh_p->set_prev(eh);

	eh->set_prev(eho_n);
	eho_n->set_next(eh);

	eho_n->set_prev(eh_p);
	eh_p->set_next(eho_n);

	eho->set_next(eho_p);
	eho_p->set_prev(eho);

	eho->set_prev(eh_n);
	eh_n->set_next(eho);

	eho_p->set_next(eh_n);
	eh_n->set_prev(eho_p);

	eh_n->set_face(f_eho);
	eho_n->set_face(f_eh);

	f_eh->set_halfedge(eh);
	f_eho->set_halfedge(eho);

	if (va0->halfedge() == eh)
		va0->set_halfedge(eho_p);
	if (va1->halfedge() == eho)
		va1->set_halfedge(eh_p);
	return true;
}

bool QMorphGenerator::recover_edge(HDS_mesh::Vertex_handle vl, HDS_mesh::Vertex_handle vr, HDS_mesh::Halfedge_handle& edge)
{
	///< 首先判断待恢复边是否已经存在
	bool exist = false;
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vl->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vl->vertex_begin();
	do 
	{
		if (edge_begin->prev()->vertex() == vr)
		{
			exist = true;
			break;
		}
		++edge_begin;
	} while (edge_begin != edge_end);
	if (exist)
	{
		edge = edge_begin;
		return true;
	}
	else
	{
		int aid = vl->getIndex();
		int bid = vr->getIndex();
		int key = 0;
		if (aid < bid)
			key = 6151 * aid + bid;
		else
			key = 6151 * bid + aid;

		Vector_3 vs(vl->point(), vr->point());
		HDS_mesh::Halfedge_handle ei;
		HDS_mesh::Halfedge_handle e_pre = vl->vertex_begin();
		edge_begin = vl->vertex_begin();
		bool flag_temp = false;
		do 
		{
			Vector3D vk(vl->point(), e_pre->prev()->vertex()->point());
			Vector3D vkk(vl->point(), edge_begin->prev()->vertex()->point());

			///< 把向量vk和vkk绕z轴正向旋转-90度
			vk = vk.rotateZ(-90);
			vkk = vkk.rotateZ(-90);
			Compute_scalar_product_3 scalar_product;
			if (scalar_product(vs, Vector_3(vk.getStart(), vk.getEnd())) > 0
				&& scalar_product(vs, Vector_3(vkk.getStart(), vkk.getEnd())) < 0)
			{
				if (e_pre->testMark(HDS_mesh::Halfedge::QuadEdge))
				{	
					return false;
					HDS_mesh::Halfedge_handle e_temp_nxt = e_pre->next()->next();
					std::list<Point_2> pts;
					Point_3 p_temp;
					bool flag = get_cross(vl->point(), vr->point(),
						e_temp_nxt->prev()->vertex()->point(), 
						e_temp_nxt->vertex()->point(), p_temp);
					if (flag)
					{
						p_temp = e_temp_nxt->vertex()->point() + 0.99 * (p_temp - e_temp_nxt->vertex()->point());
						e_pre->next()->vertex()->point() = p_temp;
						edge_begin = vl->vertex_begin();
						e_pre = edge_begin;
						//write_my_mesh_nas("seam.nas");
						continue;
					}
				}
				ei = e_pre->prev();
				flag_temp = true;
				break;
			}
			else 
				e_pre = edge_begin;
			++edge_begin;
		} while (edge_begin != edge_end);
		if (!flag_temp)
		{
			ei = vl->vertex_begin()->opposite()->next();

			if (ei->testMark(HDS_mesh::Halfedge::QuadEdge) || ei->testMark(HDS_mesh::Halfedge::FRONT))
			{	
				return false;
				HDS_mesh::Halfedge_handle e_temp_nxt = ei;

				Point_3 p_temp;
				bool flag = get_cross(vl->point(), vr->point(), 
					e_temp_nxt->prev()->vertex()->point(), 
					e_temp_nxt->vertex()->point(),
					p_temp);
				if (flag)
				{
					p_temp = e_temp_nxt->vertex()->point() + 0.99 * (p_temp - e_temp_nxt->vertex()->point());
					ei->prev()->vertex()->point() = p_temp;
					e_pre = vl->vertex_begin();
				}
				edge_begin = vl->vertex_begin();
				do 
				{
					Point_3 vk = e_pre->prev()->vertex()->point() + (-vl->point());
					Point_3 vkk = edge_begin->prev()->vertex()->point() + (-vl->point());

					///< 把向量vk和vkk绕z轴正向旋转-90度
					vk = Point_3(-vk.y(), vk.x(), vk.z());
					vkk = Point_3(-vkk.y(), vkk.x(), vkk.z());
					Compute_scalar_product_3 scalar_product;
					if (scalar_product(vs, Vector_3(Point_3(0,0,0), vk)) > 0
						&& scalar_product(vs, Vector_3(Point_3(0,0,0), vkk)) < 0)
					{
						ei = e_pre->prev();
						flag_temp = true;
						break;
					}
					else
						e_pre = edge_begin;
					++edge_begin;
				} while (edge_begin != edge_end);
			}

			Vector3D s_aux(vl->point(), vr->point());
			s_aux = s_aux.normal();
			s_aux = s_aux.rotateZ(90);
			while (1)
			{
				Point_2 p_l(vl->point().x(), vl->point().y());
				Point_2 p_r(vr->point().x(), vr->point().y());
				Point_2 p_m(ei->vertex()->point().x(), ei->vertex()->point().y());
				CGAL::Orientation ori = CGAL::orientation(p_l, p_r, p_m);
				if (ori == CGAL::LEFT_TURN)
					break;
				else
					return false;
				ei->vertex()->point() = ei->vertex()->point() + 0.1 * s_aux;
			}
		}
		if (ei->next()->vertex() != vl)
			return false;
		assert(ei->next()->vertex() == vl);

		if (ei->testMark(HDS_mesh::Halfedge::FRONT))
		{
			return false;
		}

		std::queue<HDS_mesh::Halfedge_handle> edge_queue; ///< 用于存储所有与待恢复边相交的边
		edge_queue.push(ei);
		while (1)
		{///< 收集所有与待恢复边相交的边
			HDS_mesh::Face_handle ti = ei->opposite()->face();
			ei = ei->opposite();
			if (has_node(ti, vr)) ///< 如果当前三角形中包含待恢复边的目标点，则循环结束
				break;
			Vector3D vi(vl->point(), ei->next()->vertex()->point());
			vi = vi.rotateZ(-90);
			Compute_scalar_product_3 scalar_product;
			if (scalar_product(vs, Vector_3(vi.getStart(), vi.getEnd())) < 0)
				ei = ei->prev();				
			else
				ei = ei->next();
			if (!ei->testMark(HDS_mesh::Halfedge::FRONT))
				edge_queue.push(ei);
			else
				return false;
		}

		std::list<HDS_mesh::Halfedge_handle> error_list1; ///< 用于存储恢复失败的边
		std::list<HDS_mesh::Halfedge_handle> error_list2; ///< 用于存储恢复失败的边
		HDS_mesh::Halfedge_handle nedge;
		while (!edge_queue.empty())
		{
			HDS_mesh::Halfedge_handle edge_temp = edge_queue.front();
			edge_queue.pop();

			if (!flip_edge(edge_temp))
			{
				// 				if (edge_queue.size() == 0)
				// 					return false;
				//edge_queue.push(edge_temp);
				error_list1.push_back(edge_temp);
				if (!edge_queue.empty())
					continue;
			}
			else
			{
				nedge = edge_temp;
				int aid = nedge->vertex()->getIndex();
				int bid = nedge->prev()->vertex()->getIndex();
				int key_n = 0;
				if (aid < bid)
					key_n = 6151 * aid + bid;
				else
					key_n = 6151 * bid + aid;

				if (key == key_n)
					break;
				else
				{
					if (nedge->vertex() == vl || nedge->vertex() == vr || nedge->prev()->vertex() == vl || nedge->prev()->vertex() == vr)
					{
						if (!edge_queue.empty())
							continue;
					}
					else
					{
						Point_3 p_temp;
						bool cross_flag = get_cross(nedge->vertex()->point(), 
							nedge->prev()->vertex()->point(),
							vl->point(), vr->point(), p_temp);
						if (cross_flag)
						{
							//edge_queue.push(nedge);
							error_list1.push_back(nedge);
						}
					}
				}///< end if (key == key_n)
			}///< end if (!flip_edge(edge_temp))
			if (edge_queue.empty() && !error_list1.empty())
			{
				std::list<HDS_mesh::Halfedge_handle>::iterator it_error_edge1;
				for (it_error_edge1 = error_list1.begin(); it_error_edge1 != error_list1.end(); ++it_error_edge1)
				{
					edge_queue.push(*it_error_edge1);
				}
				if (error_list1.size() == error_list2.size())
				{
					int n_diff = 0;
					std::list<HDS_mesh::Halfedge_handle>::iterator it_error_edge2;
					for (it_error_edge1 = error_list1.begin(), it_error_edge2 = error_list2.begin(); 
						it_error_edge1 != error_list1.end(); ++it_error_edge1, ++it_error_edge2)
					{
						if (*it_error_edge1 != *it_error_edge2)
							++n_diff;
					}
					if (n_diff > 0)
					{
						error_list2.clear();
						error_list2 = error_list1;
						error_list1.clear();
						continue;
					}
					else
						return false;
				}
				else
				{
					error_list2.clear();
					error_list2 = error_list1;
					error_list1.clear();
				}
			}
		}///< end while
		edge = nedge;
		return true;
	}
}

bool QMorphGenerator::merge_edge(HDS_mesh::Halfedge_handle& edge, double t)
{
	assert(!edge->testMark(HDS_mesh::Halfedge::BOUNDARY)); ///< 仅允许对非边界边执行合并操作
	assert(t >= 0 && t <= 1);

	HDS_mesh::Vertex_handle v_sorce = edge->prev()->vertex();
	HDS_mesh::Vertex_handle v_target = edge->vertex();

	Point_3 p_temp = v_sorce->point() + t * (v_target->point() + (-v_sorce->point())); ///< 合并后的点

	HDS_mesh::Halfedge_handle edgeNext = edge->next();
	HDS_mesh::Halfedge_handle edgePrev = edge->prev();

	HDS_mesh::Halfedge_handle oppEdge = edge->opposite();
	HDS_mesh::Halfedge_handle oppEdgeNext = oppEdge->next();
	HDS_mesh::Halfedge_handle oppEdgePrev = oppEdge->prev();

	if (edgeNext->vertex()->is_trivalent())
	{
		edgeNext = _source_mesh->erase_center_vertex(edgeNext)->next();
		edgePrev = edgeNext->next();
	}

	if (edgeNext->testMark(HDS_mesh::Halfedge::FRONT))
	{
		_source_mesh->flip_edge(oppEdgeNext);
		return true;
	}

	int ids_en = edgeNext->prev()->vertex()->getIndex();
	int idt_en = edgeNext->vertex()->getIndex();
	if (ids_en > idt_en)
	{
		int t = ids_en;
		ids_en = idt_en;
		idt_en = t;
	}

	int ids_eop = oppEdgePrev->prev()->vertex()->getIndex();
	int idt_eop = oppEdgePrev->vertex()->getIndex();
	if (ids_eop > idt_eop)
	{
		int t = ids_eop;
		ids_eop = idt_eop;
		idt_eop = t;
	}

	bool is_adj = false;
	if (ids_en == ids_eop && idt_en == idt_eop)
	{
		is_adj = true;
	}

	HDS_mesh::Face_handle face_base = edge->face();
	HDS_mesh::Face_handle face_opp = oppEdge->face();

	if (is_adj)
	{
		if (edgeNext->vertex()->halfedge() == edgeNext)
			edgeNext->vertex()->set_halfedge(oppEdgeNext);
		if (v_sorce->halfedge() == oppEdge)
			v_sorce->set_halfedge(oppEdgeNext->opposite());

		edgePrev->opposite()->prev()->set_next(oppEdgeNext);

		oppEdgeNext->vertex()->set_halfedge(oppEdgeNext);

		oppEdgeNext->set_prev(edgePrev->opposite()->prev());
		edgePrev->opposite()->next()->set_prev(oppEdgeNext);
		oppEdgeNext->set_next(edgePrev->opposite()->next());
		oppEdgeNext->set_face(edgePrev->opposite()->face());
		//assert(oppEdgeNext->opposite()->face()->testMark(HDS_mesh::Face::QUAD));

		if (edgePrev->testMark(HDS_mesh::Halfedge::FRONT))
		{
			oppEdgeNext->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			oppEdgeNext->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
			oppEdgeNext->opposite()->set_loopPrev(edgePrev->loopPrev());
			edgePrev->loopPrev()->set_loopNext(oppEdgeNext->opposite());
			oppEdgeNext->opposite()->set_loopNext(edgePrev->loopNext());
			edgePrev->loopNext()->set_loopPrev(oppEdgeNext->opposite());
			v_target->setMark(HDS_mesh::Vertex::FrontNode);
			oppEdgeNext->setMark(HDS_mesh::Halfedge::QuadEdge);
			oppEdgeNext->opposite()->set_FrontState(oppEdgeNext->opposite()->get_FrontState());
			if (edgePrev->testMark(HDS_mesh::Halfedge::UNACTIVE))
			{
				deleteFromUnactiveFront(edgePrev);
			}
			else
			{
				deleteFromFront(edgePrev);
			}
			_unactive_front_list.push_back(oppEdgeNext->opposite());
		}

		_source_mesh->hds().edges_erase(edgeNext);
		_source_mesh->hds().edges_erase(edgePrev);
		_source_mesh->hds().edges_erase(edge);
		// 		_backgroundMesh.hds().edges_erase(edgeNext);
		_source_mesh->hds().vertices_erase(v_target);

		_source_mesh->hds().faces_erase(face_base);
		_source_mesh->hds().faces_erase(face_opp);
	}
	else
	{
		if (edgeNext->opposite()->face() == oppEdgePrev->opposite()->face())
		{
			if (oppEdgeNext->testMark(HDS_mesh::Halfedge::FRONT))
			{
				if (edgePrev->testMark(HDS_mesh::Halfedge::FRONT))
				{
					edgeNext->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
					edgeNext->opposite()->setMark(HDS_mesh::Halfedge::UNACTIVE);
					edgeNext->opposite()->set_loopPrev(edgePrev->loopPrev());
					edgePrev->loopPrev()->set_loopNext(edgeNext->opposite());
					edgeNext->opposite()->set_loopNext(edgePrev->loopNext());
					edgePrev->loopNext()->set_loopPrev(edgeNext->opposite());
					v_target->setMark(HDS_mesh::Vertex::FrontNode);
					edgeNext->setMark(HDS_mesh::Halfedge::QuadEdge);
					edgeNext->opposite()->set_FrontState(oppEdgeNext->opposite()->get_FrontState());
					if (edgePrev->testMark(HDS_mesh::Halfedge::UNACTIVE))
					{
						deleteFromUnactiveFront(edgePrev);
					}
					else
					{
						deleteFromFront(edgePrev);
					}
					_unactive_front_list.push_back(edgeNext->opposite());
				}
				_source_mesh->join_facet(oppEdgePrev);
				_source_mesh->join_facet(edgePrev);
				edgeNext->face()->clearMark();
				edgeNext->face()->setMark(HDS_mesh::Face::QUAD);
			}
			else
			{
				_source_mesh->join_facet(oppEdgeNext);
				assert(!edgeNext->testMark(HDS_mesh::Halfedge::FRONT));
				_source_mesh->join_facet(edgeNext);
			}
		}
		else
		{
			_source_mesh->join_facet(oppEdgePrev);
			assert(!edgeNext->testMark(HDS_mesh::Halfedge::FRONT));
			_source_mesh->join_facet(edgeNext);
		}
		//oppEdge->opposite()->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
		_source_mesh->join_vertex(oppEdge)->vertex()->point() = p_temp;
		//v_target->clearMark();
		v_target->setMark(HDS_mesh::Vertex::ERASED);
	}
	/*write_my_mesh_nas("join.nas");*/

	return true;
}

void QMorphGenerator::localSmoothing(std::list<HDS_mesh::Vertex_handle>& smooth_node_list)
{	
	std::list<HDS_mesh::Vertex_handle>::iterator it_vertex;
	std::list<HDS_mesh::Vertex_handle> front_node_list;
	std::list<HDS_mesh::Vertex_handle> interior_node_list;
	for (it_vertex = smooth_node_list.begin(); it_vertex != smooth_node_list.end();  ++it_vertex)
	{
		if ((*it_vertex)->testMark(HDS_mesh::Vertex::BOUNDARY))
			continue; ///< 边界点不执行光顺
		if ((*it_vertex)->testMark(HDS_mesh::Vertex::FrontNode))
			front_node_list.push_back((*it_vertex));
		else
			interior_node_list.push_back((*it_vertex));
	}
	//smooth_node_list.clear();
	for (it_vertex = front_node_list.begin(); it_vertex != front_node_list.end(); ++it_vertex)
	{
		//smooth_node_list.push_back(*it_vertex);
		smoothFrontNode(*it_vertex);
	}
	for (it_vertex = interior_node_list.begin(); it_vertex != interior_node_list.end(); ++it_vertex)
	{
		//smooth_node_list.push_back(*it_vertex);
		_optimizer->smoothInteriorNode(*it_vertex);
	}
	for (it_vertex = front_node_list.begin(); it_vertex != front_node_list.end(); ++it_vertex)
		smoothFrontNode(*it_vertex);
}

void QMorphGenerator::smoothFrontNode(HDS_mesh::Vertex_handle& vh)
{
	if (vh->testMark(HDS_mesh::Vertex::LOCKED))
		return;
	//if (_reatio < RATIO_RANGE_LEFT)
	{
		///< 首先获取点vh的邻接面片链表，根据其中四边形单元的个数采用不同个光顺策略
		std::vector<HDS_mesh::Face_handle> quad_face_list; ///< 收集所有的邻接四边形
		Point_3 new_point(0,0,0);
		HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh->vertex_begin();
		do 
		{
			if (!edge_begin->is_border())
			{
				HDS_mesh::Face_handle f_temp = edge_begin->face();
				if (f_temp->testMark(HDS_mesh::Face::QUAD))
				{
					HDS_mesh::Halfedge_handle it_edge_begin = f_temp->halfedge();
					while (it_edge_begin->vertex() != vh)
					{
						it_edge_begin = it_edge_begin->next();
					}
					Point_3 pj = it_edge_begin->next()->vertex()->point();
					Point_3 pk = it_edge_begin->next()->next()->vertex()->point();
					Point_3 pl = it_edge_begin->prev()->vertex()->point();
					new_point = new_point + pj + pl + (-pk);
					quad_face_list.push_back(f_temp);
				}
			}
			++edge_begin;
		} while (edge_begin != edge_end);

		assert(quad_face_list.size() != 0);
		new_point = new_point / quad_face_list.size();

		if (2 != quad_face_list.size())
		{
			_optimizer->smoothAdjustment(vh, new_point);
			return;
		}
		///< 对于只有两个邻接四边形单元的前沿点(需对上述等参光顺得到的新节点位置做调整)
		HDS_mesh::Halfedge_handle it_edge_begin = quad_face_list[0]->halfedge();
		while (it_edge_begin->vertex() != vh)
		{
			it_edge_begin = it_edge_begin->next();
		}
		if (!it_edge_begin->opposite()->face()->testMark(HDS_mesh::Face::QUAD))
			it_edge_begin = it_edge_begin->next()->opposite();
		assert(it_edge_begin->vertex() == vh);
		Point_3 pj = it_edge_begin->prev()->vertex()->point(); ///< 在两个四边形单元的公共边上与点vh相对的点
		Point_3 detaA = new_point + (-vh->point());
		double lenA = CGAL::sqrt(CGAL::squared_distance(pj, new_point)); ///< pj与等参光顺获得的新节点位置之间的距离
		double lenD = 0.0;
		//assert(it_edge_begin->prev()->vertex()->testMark(HDS_vertex::FrontNode));
		///< 获取pj点处的链条前沿边
		HDS_mesh::Halfedge_handle left_eh = it_edge_begin->prev();
		//assert(left_eh->testMark(HDS_halfedge::FRONT));
		HDS_mesh::Halfedge_handle right_eh = it_edge_begin->opposite()->next();
		assert(left_eh->vertex() == right_eh->prev()->vertex());

		///< 计算距离调整参数
		double d1 = CGAL::sqrt(CGAL::squared_distance(left_eh->prev()->vertex()->point(), 
			left_eh->vertex()->point()));
		double d2 = CGAL::sqrt(CGAL::squared_distance(right_eh->prev()->vertex()->point(), 
			right_eh->vertex()->point()));

		Point_3 ap1 = left_eh->vertex()->point();
		Point_3 ap2 = left_eh->prev()->vertex()->point(); 
		Point_3 ap3 = right_eh->vertex()->point();

		double angle_temp = angle_evaluate(ap2, ap1, ap3);

		if (_reatio < 2.5)
		{
			if (angle_temp < PAVING_SIDE)
			{
				lenD = (d1 + d2) / (2 * sin(angle_temp * PI/ 360.0));
			}
			else if (angle_temp < PAVINE_CORNER)
			{
				lenD = (d1 + d2) / (2 * sin(angle_temp * PI/ 540.0));
				lenD *= std::sqrt(2.0);
			}
			else
				lenD = (d1 + d2) / (2 * sin(angle_temp * PI/ 720.0));
		}
		else if (_reatio < 20)
		{
			double len1 = CGAL::sqrt(CGAL::squared_distance(it_edge_begin->next()->vertex()->point(), 
				it_edge_begin->next()->next()->vertex()->point()));
			double len2 = CGAL::sqrt(CGAL::squared_distance(it_edge_begin->prev()->prev()->vertex()->point(), 
				it_edge_begin->prev()->vertex()->point()));
			HDS_mesh::Halfedge_handle opposite_edge = it_edge_begin->opposite();
			double len3 = CGAL::sqrt(CGAL::squared_distance(opposite_edge->vertex()->point(), 
				opposite_edge->next()->vertex()->point()));
			double len4 = CGAL::sqrt(CGAL::squared_distance(opposite_edge->next()->vertex()->point(), 
				opposite_edge->next()->next()->vertex()->point()));
			double len_sum = len1 + len2 + len3 +len4;
			int ne_temp = 0;
			edge_begin = vh->vertex_begin();
			do 
			{///< 遍历vh的入射边表
				if (edge_begin->face()->testMark(HDS_mesh::Face::QUAD) 
					|| !edge_begin->testMark(HDS_mesh::Halfedge::BOUNDARY) 
					&& edge_begin->opposite()->testMark(HDS_mesh::Face::QUAD))
				{
					++edge_begin;
					continue;
				}
				++ne_temp;
				len_sum += CGAL::sqrt(CGAL::squared_distance(vh->point(), 
					edge_begin->prev()->vertex()->point()));
				++edge_begin;
			} while (edge_begin != edge_end);
			lenD = len_sum * 1.0 / (ne_temp + 4);
		}
		else
		{
			int n_edge = 0;
			double len_sum = 0.0;
			edge_begin = vh->vertex_begin();
			do 
			{
				double len = CGAL::sqrt(CGAL::squared_distance(vh->point(), 
					edge_begin->prev()->vertex()->point()));
				len_sum += len;
				++n_edge;
				++edge_begin;
			} while (edge_begin != edge_end);
			assert(0 != n_edge);
			lenD = len_sum * 1.0 / n_edge;
		}
		assert(lenA != 0);
		Point_3 detaB = pj + (-vh->point()) + (detaA + vh->point() + (-pj)) * lenD / lenA; ///< 距离调整参数

		///< 计算角度调整参数
		Point_3 p_left = it_edge_begin->next()->vertex()->point();
		Point_3 p_right = it_edge_begin->opposite()->next()->next()->vertex()->point();
		Vector3D v_left(pj, p_left);
		Vector3D v_right(pj, p_right);
		Vector3D vm = bisectVector(v_left, v_right);
		Vector3D v_mid(pj, vh->point());
		Vector3D nq = bisectVector(vm, v_mid);
		Point_3 pq;
		bool flag = get_cross(nq.getStart(), nq.getEnd(), p_left, p_right, pq); ///< 计算得到Q点
		assert(flag);
		double lenQ = CGAL::sqrt(CGAL::squared_distance(pj, pq));
		double lenB2 = 0.0;
		if (lenD > lenQ) 
			lenB2 = (lenD + lenQ) / 2.0;
		else
			lenB2 = lenD;

		double t = lenB2 / lenQ; ///< 参数t
		Point_3 B2 = pj + t * (pq + (-pj));
		Point_3 detaC = B2 + (-vh->point()); ///< 角度调整参数
		Point_3 detaI = (detaB + detaC) / 2.0; ///< 综合调整参数
		new_point = detaI + vh->point(); ///< 计算得到的新节点位置
		_optimizer->smoothAdjustment(vh, new_point);
		return;
	}
	// 	else if (_reatio > 20)
	// 	{
	// 		smoothInteriorNode(vh);
	// 		return;
	// 	}
}

bool QMorphGenerator::OptBasedSmooth(HDS_mesh::Vertex_handle& vh)
{
	std::vector<double> miu_old;///< 用于存储扰动前各邻接单元的质量系数
	std::vector<double> miu_new;///< 用于存储扰动后各邻接单元的质量系数
	std::vector<Vector_3> g;
	double min_miu_old = DBL_MAX;
	int min_id = 0;
	int i = 0;
	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = vh->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = vh->vertex_begin();
	do 
	{
		if (edge_begin->face()->is_triangle())
		{
			double mo = _optimizer->tri_quality(edge_begin->vertex(), 
				edge_begin->next()->vertex(), 
				edge_begin->prev()->vertex());
			if (min_miu_old > mo)
			{
				min_miu_old = mo;
				min_id = i;
			}
			miu_old.push_back(mo);

			double mnx = _optimizer->tri_quality(edge_begin->vertex()->point() + Point_3(_deta, 0.0, 0.0), 
				edge_begin->next()->vertex()->point(), 
				edge_begin->prev()->vertex()->point());
			double mny = _optimizer->tri_quality(edge_begin->vertex()->point() + Point_3(0.0, _deta, 0.0), 
				edge_begin->next()->vertex()->point(), 
				edge_begin->prev()->vertex()->point());
			// 			double mnz = tri_quality(edge_begin->vertex()->point() + Point_3(0.0, 0.0, _deta), 
			// 				                     edge_begin->next()->vertex()->point(), 
			// 				                     edge_begin->prev()->vertex()->point());
			g.push_back(Vector_3((mnx - mo) / _deta, (mny - mo) / _deta, 0.0));
		}
		else
		{
			double mo = _optimizer->quad_quality2(edge_begin->vertex(), 
				edge_begin->next()->vertex(), 
				edge_begin->next()->next()->vertex(),
				edge_begin->prev()->vertex());
			if (min_miu_old > mo)
			{
				min_miu_old = mo;
				min_id = i;
			}
			miu_old.push_back(mo);
			double mnx = _optimizer->quad_quality2(edge_begin->vertex()->point() + Point_3(_deta, 0.0, 0.0), 
				edge_begin->next()->vertex()->point(), 
				edge_begin->next()->next()->vertex()->point(),
				edge_begin->prev()->vertex()->point());
			double mny = _optimizer->quad_quality2(edge_begin->vertex()->point() + Point_3(0.0, _deta, 0.0), 
				edge_begin->next()->vertex()->point(), 
				edge_begin->next()->next()->vertex()->point(),
				edge_begin->prev()->vertex()->point());
			g.push_back(Vector_3((mnx - mo) / _deta, (mny - mo) / _deta, 0.0));
		}
		++edge_begin;
		++i;
	} while (edge_begin != edge_end);
	Vector_3 gm = g[min_id];
	double sec_min_miu = DBL_MAX;
	int sec_min_id = 0;
	if (std::abs(gm[0]) < TOLERENCE && std::abs(gm[1]) < TOLERENCE && std::abs(gm[2]) < TOLERENCE)
	{
		for (int im = 0; im < miu_old.size(); ++im)
		{
			if (im == min_id)
				continue;
			if (miu_old[im] < sec_min_miu)
			{
				sec_min_miu = miu_old[im];
				sec_min_id = im;
			}
		}
		gm = g[sec_min_id];
	}
	double garm = DBL_MAX;
	bool flag = false;
	for (int ig = 0; ig < g.size(); ++ig)
	{
		if (g[ig] * gm < 0)
		{
			flag = true;
			double garm_temp = (miu_old[ig] - min_miu_old) / (gm * gm - gm * g[ig]);
			if (garm_temp < garm)
				garm = garm_temp;
		}
	}
	int num = 0;
	while (1 && flag)
	{
		double min_miu_new = DBL_MAX;
		for (int im = 0; im < miu_old.size(); ++im)
		{
			double miu_temp = miu_old[im] + garm * gm * g[im];
			if (min_miu_new > miu_temp)
				min_miu_new = miu_temp;
		}
		if (min_miu_new - min_miu_old >= 1e-5)
		{
			vh->point() = vh->point() + garm * gm;
			return true;
		}
		else
			garm /= 2;
		++num;
		if (num > 2)
			break;
	}
	return false;
}

void QMorphGenerator::rePushFront(HDS_mesh::Halfedge_handle eh)
{
	eh->setMark(HDS_mesh::Halfedge::FRONT);
	eh->setMark(HDS_mesh::Halfedge::UNACTIVE);
	eh->set_FrontState(0);
	_unactive_front_list.push_back(eh);
	// 	switch (eh->get_FrontState())
	// 	{
	// 	case 0:
	// 		_front_edge_00.push_back(eh);
	// 		break;
	// 	case 1:
	// 		_front_edge_01.push_back(eh);
	// 		break;
	// 	case 2:
	// 		_front_edge_10.push_back(eh);
	// 		break;
	// 	case 3:
	// 		_front_edge_11.push_back(eh);
	// 		break;
	// 	default:
	// 		break;
	// 	}
}

bool QMorphGenerator::has_node(HDS_mesh::Face_const_handle tri, HDS_mesh::Vertex_handle vh)
{
	HDS_mesh::Halfedge_const_handle edge_begin = tri->halfedge();
	HDS_mesh::Halfedge_const_handle edge_end = tri->halfedge();
	do 
	{
		if (edge_begin->vertex() == vh)
			return true;
		edge_begin = edge_begin->next();
	} while (edge_begin != edge_end);
	return false;
}

bool QMorphGenerator::canSeam(HDS_mesh::Halfedge_handle& e_left)
{
	int n_front = 0;
	HDS_mesh::Halfedge_handle current_front = e_left;
	while (current_front->vertex() != e_left->prev()->vertex())
	{
		++n_front;
		current_front = current_front->loopNext();
		if (n_front > 3)
			break;
	}
	if (n_front == 3)
		return false;
	return true;
}

HDS_mesh::Halfedge_handle QMorphGenerator::seamFront(HDS_mesh::Halfedge_handle& e_left, 
													 HDS_mesh::Halfedge_handle& e_right, 
													 bool is_top)
{
	if (e_left->opposite()->face() == e_right->opposite()->face())
	{
		return HDS_mesh::Halfedge_handle();
	}
	if (e_left->prev()->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY)
		|| e_right->vertex()->testMark(HDS_mesh::Vertex::BOUNDARY))
	{
		return HDS_mesh::Halfedge_handle();
	}
	double len_left = CGAL::squared_distance(e_left->prev()->vertex()->point(), 
		e_left->vertex()->point());
	double len_right = CGAL::squared_distance(e_right->prev()->vertex()->point(), 
		e_right->vertex()->point());
	double ratio_temp = 0.0;
	HDS_mesh::Halfedge_handle e_biger = e_left;
	HDS_mesh::Halfedge_handle e_shorter = e_right;

	if (len_left > len_right)
		ratio_temp = len_left / len_right;
	else
	{
		ratio_temp = len_right / len_left;
		e_biger = e_right;
		e_shorter = e_left;
	}

	double qua_temp = _optimizer->quad_quality(e_biger->opposite()->vertex(), e_biger->opposite()->next()->vertex(),
		                                       e_biger->opposite()->next()->next()->vertex(), e_biger->vertex());
	if (ratio_temp > 6.25 && qua_temp > 0)
	{
		return seamFrontTrans(e_biger, e_shorter, is_top);
	}

	return seamFrontNormal(e_left, e_right, is_top);
}

HDS_mesh::Halfedge_handle QMorphGenerator::seamFrontNormal(HDS_mesh::Halfedge_handle& e_left,
														HDS_mesh::Halfedge_handle& e_right,
														bool is_top)
{
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	///< 待合并的两点
	HDS_mesh::Vertex_handle v_left = e_left->prev()->vertex();
	HDS_mesh::Vertex_handle v_right = e_right->vertex();
	HDS_mesh::Vertex_handle v_aux = e_left->vertex();

	HDS_mesh::Vertex_handle vk = e_left->vertex();
	HDS_mesh::Halfedge_handle edge_temp;
	bool flag = recover_edge(v_left, v_right, edge_temp); ///< 首先恢复两点之间的边
	if (!flag)
		return HDS_mesh::Halfedge_handle();
	assert(flag);

	if (edge_temp->vertex() != v_right)
		edge_temp = edge_temp->opposite();
	// 	if (!edge_temp->next()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
	// 		smoothInteriorNode(edge_temp->next()->vertex());

	if (edge_temp->next()->vertex()->is_trivalent())
	{
		int n_tri = 0;
		HDS_mesh::Halfedge_around_vertex_circulator e_begin_temp = edge_temp->next()->vertex()->vertex_begin();
		HDS_mesh::Halfedge_around_vertex_circulator e_end_temp = edge_temp->next()->vertex()->vertex_begin();
		do 
		{
			if (e_begin_temp->face()->is_triangle())
				++n_tri;
			++e_begin_temp;
		} while (e_begin_temp != e_end_temp);
		if (3 == n_tri)
			_source_mesh->erase_center_vertex(edge_temp->next());
	}

	if (edge_temp->next()->testMark(HDS_mesh::Halfedge::FRONT) 
		&& edge_temp->prev()->testMark(HDS_mesh::Halfedge::FRONT))
	{
		if (!edge_temp->next()->vertex()->testMark(HDS_mesh::Vertex::LOOPCUT))
			edge_temp->next()->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
	}

	Point_3 ver_d = edge_temp->next()->vertex()->point();
	Point_3 ver_a = edge_temp->prev()->vertex()->point();
	Point_3 ver_c = edge_temp->vertex()->point();///< 备选合并点

	e_left->setMark(HDS_mesh::Halfedge::QuadEdge);
	e_right->setMark(HDS_mesh::Halfedge::QuadEdge);
	edge_temp->next()->setMark(HDS_mesh::Halfedge::QuadEdge);
	edge_temp->prev()->setMark(HDS_mesh::Halfedge::QuadEdge);

	HDS_mesh::Face_handle face_base_temp = edge_temp->face();
	std::list<HDS_mesh::Vertex_handle> old_vertex_list; ///< 待删除节点链表
	std::list<HDS_mesh::Face_handle> old_face_list; ///< 待删除面片链表
	std::list<HDS_mesh::Halfedge_handle> old_edge_list; ///< 待删除边链表
	getFaceToDelete(face_base_temp, old_face_list, old_edge_list, old_vertex_list);

	Point_3 point_temp = (v_left->point() + v_right->point()) / 2.0;///< 合并后的节点位置

	HDS_mesh::Vertex_handle vd = edge_temp->vertex();
	HDS_mesh::Halfedge_handle ed1 = edge_temp->next();
	HDS_mesh::Halfedge_handle ed2 = edge_temp->next()->opposite();
	HDS_mesh::Halfedge_handle ed3 = e_right;
	HDS_mesh::Halfedge_handle ed4 = e_right->opposite();

	HDS_mesh::Halfedge_handle edge_temp_prev = edge_temp->prev()->opposite();

	HDS_mesh::Halfedge_around_vertex_circulator edge_begin = v_right->vertex_begin();
	HDS_mesh::Halfedge_around_vertex_circulator edge_end = v_right->vertex_begin();
	do 
	{
		if (edge_begin->vertex() == v_right)
			edge_begin->set_vertex(v_left);
		++edge_begin;
	} while (edge_begin != edge_end);

	if (is_top)
	{
		e_right->loopNext()->set_loopPrev(e_left->loopPrev());
		e_left->loopPrev()->set_loopNext(e_right->loopNext());

		e_left->opposite()->prev()->set_loopNext(e_left->opposite());
		e_left->opposite()->set_loopPrev(e_left->opposite()->prev());
		e_left->opposite()->set_loopNext(e_right->loopNext());
	}

	if (ed1->testMark(HDS_mesh::Halfedge::FRONT))
	{
		if (edge_temp->prev()->testMark(HDS_mesh::Halfedge::FRONT))
		{
			edge_temp->prev()->resetMark(HDS_mesh::Halfedge::FRONT);
			if (edge_temp->prev()->testMark(HDS_mesh::Halfedge::UNACTIVE))
				deleteFromUnactiveFront(edge_temp->prev());
			else
				deleteFromFront(edge_temp->prev());
			edge_temp->prev()->resetMark(HDS_mesh::Halfedge::UNACTIVE);
		}
		else
		{
			edge_temp->prev()->opposite()->setMark(HDS_mesh::Halfedge::FRONT);
			edge_temp->prev()->opposite()->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			edge_temp->prev()->vertex()->setMark(HDS_mesh::Vertex::FrontNode);

			edge_temp->prev()->opposite()->set_FrontState(0);
			rePushFront(edge_temp->prev()->opposite());
		}

		ed1->opposite()->set_opposite(edge_temp->prev()->opposite());

		e_left->loopPrev()->set_loopNext(edge_temp->prev()->opposite());
		edge_temp->prev()->opposite()->set_loopPrev(e_left->loopPrev());
		ed1->resetMark(HDS_mesh::Halfedge::FRONT);
		ed1->loopNext()->set_loopPrev(edge_temp->prev()->opposite());
		edge_temp->prev()->opposite()->set_loopNext(ed1->loopNext());
		ed1->vertex()->set_halfedge(edge_temp->prev()->opposite());
		e_left->opposite()->prev()->set_loopNext(e_left->opposite());
		e_left->opposite()->set_loopNext(edge_temp->prev()->opposite());
		if (ed1->testMark(HDS_mesh::Halfedge::UNACTIVE))
			deleteFromUnactiveFront(ed1);
		else
			deleteFromFront(ed1);
	}
	else
	{
		e_left->loopPrev()->set_loopNext(e_right->loopNext());
		e_right->loopNext()->set_loopPrev(e_left->loopPrev());
		ed1->vertex()->set_halfedge(edge_temp->prev()->opposite());
		e_right->opposite()->vertex()->set_halfedge(e_left);
	}

	if (e_right->testMark(HDS_mesh::Halfedge::UNACTIVE))
		deleteFromUnactiveFront(e_right);
	else
		deleteFromFront(e_right);

	if (e_left->testMark(HDS_mesh::Halfedge::UNACTIVE))
		deleteFromUnactiveFront(e_left);
	else
		deleteFromFront(e_left);

	e_left->resetMark(HDS_mesh::Halfedge::FRONT);
	e_right->resetMark(HDS_mesh::Halfedge::FRONT);

	edge_temp->next()->resetMark(HDS_mesh::Halfedge::QuadEdge);
	edge_temp->prev()->resetMark(HDS_mesh::Halfedge::QuadEdge);

	e_right->opposite()->prev()->set_next(e_left);
	e_left->set_prev(e_right->opposite()->prev());
	e_right->opposite()->prev()->set_vertex(e_left->opposite()->vertex());
	e_right->opposite()->set_opposite(e_left->opposite());

	e_right->opposite()->next()->set_prev(e_left);
	e_left->set_next(e_right->opposite()->next());

	e_left->set_face(e_right->opposite()->face());

	e_left->vertex()->set_halfedge(e_left);

	e_right->opposite()->face()->set_halfedge(e_left);

	edge_temp->next()->opposite()->prev()->set_next(edge_temp->prev());
	edge_temp->prev()->set_prev(edge_temp->next()->opposite()->prev());

	edge_temp->next()->opposite()->next()->set_prev(edge_temp->prev());
	edge_temp->prev()->set_next(edge_temp->next()->opposite()->next());
	edge_temp->prev()->set_face(edge_temp->next()->opposite()->face());
	edge_temp->next()->opposite()->face()->set_halfedge(edge_temp->prev());

	e_left->prev()->vertex()->point() = point_temp;

	///< 删除多余面片
	std::list<HDS_mesh::Face_handle>::iterator it_face_delete;
	for (it_face_delete = old_face_list.begin(); 
		it_face_delete != old_face_list.end(); ++it_face_delete)
	{
		_source_mesh->hds().faces_erase(*it_face_delete);
	}

	///< 删除多余半边
	std::list<HDS_mesh::Halfedge_handle>::iterator it_edge_delete;
	for (it_edge_delete = old_edge_list.begin(); it_edge_delete != old_edge_list.end(); ++it_edge_delete)
	{
		_source_mesh->hds().edges_erase(*it_edge_delete);
	}

	///< 删除多余节点
	std::list<HDS_mesh::Vertex_handle>::iterator it_vertex_delete;
	for (it_vertex_delete = old_vertex_list.begin(); it_vertex_delete != old_vertex_list.end(); ++it_vertex_delete)
	{
		_source_mesh->hds().vertices_erase(*it_vertex_delete);
	}

	if (vd->testMark(HDS_mesh::Vertex::SEAMCHECK))
		vd->setDead(); ///< 把节点标记为删除
	else
		_source_mesh->hds().vertices_erase(vd);

	_source_mesh->hds().edges_erase(e_right);
	_source_mesh->hds().edges_erase(ed1);

	///< 获取v_left的邻接单元
	decorator.set_vertex_halfedge(e_left->loopPrev());

	std::list<HDS_mesh::Vertex_handle> smoothed_node_list;

	getVertexToSmooth(e_left->face(), smoothed_node_list);

	localSmoothing(smoothed_node_list);

	bool flag_temp = true;
	edge_begin = v_left->vertex_begin();
	edge_end = v_left->vertex_begin();
	do 
	{
		if (!edge_begin->is_border() && edge_begin->face()->testMark(HDS_mesh::Face::TRI))
		{
			Point_3 pa = edge_begin->prev()->vertex()->point();
			Point_3 pb = edge_begin->vertex()->point();
			Point_3 pc = edge_begin->next()->vertex()->point();
			CGAL::Orientation ori = CGAL::orientation(Point_2(pa.x(), pa.y()), Point_2(pb.x(), pb.y()), Point_2(pc.x(), pc.y()));
			if (ori == CGAL::RIGHT_TURN)
			{
				flag_temp = false;
				break;
			}
		}
		++edge_begin;
	} while (edge_begin != edge_end);

	if (!flag_temp)
	{
		flag_temp = true;
		e_left->prev()->vertex()->point() = ver_a;
		edge_begin = v_left->vertex_begin();
		edge_end = v_left->vertex_begin();
		do 
		{
			if (edge_begin->face()->testMark(HDS_mesh::Face::TRI))
			{
				Point_3 pa = edge_begin->prev()->vertex()->point();
				Point_3 pb = edge_begin->vertex()->point();
				Point_3 pc = edge_begin->next()->vertex()->point();
				CGAL::Orientation ori = CGAL::orientation(Point_2(pa.x(), pa.y()), Point_2(pb.x(), pb.y()), Point_2(pc.x(), pc.y()));
				if (ori == CGAL::RIGHT_TURN)
				{
					flag_temp = false;
					break;
				}
			}
			++edge_begin;
		} while (edge_begin != edge_end);
	}
	if (!flag_temp)
	{
		flag_temp = true;
		e_left->prev()->vertex()->point() = ver_c;
		edge_begin = v_left->vertex_begin();
		edge_end = v_left->vertex_begin();
		do 
		{
			if (edge_begin->face()->testMark(HDS_mesh::Face::TRI))
			{
				Point_3 pa = edge_begin->prev()->vertex()->point();
				Point_3 pb = edge_begin->vertex()->point();
				Point_3 pc = edge_begin->next()->vertex()->point();
				CGAL::Orientation ori = CGAL::orientation(Point_2(pa.x(), pa.y()), Point_2(pb.x(), pb.y()), Point_2(pc.x(), pc.y()));
				if (ori == CGAL::RIGHT_TURN)
				{
					flag_temp = false;
					break;
				}
			}
			++edge_begin;
		} while (edge_begin != edge_end);
	}

	if (!flag_temp/* && !edge_temp_prev->testMark(HDS_mesh::Halfedge::QuadEdge)*/)
	{
		flag_temp = true;
		std::cout << "merge_edge()\n";
		e_left->prev()->vertex()->point() = point_temp;		
		//write_my_mesh_nas("join_vertex.nas");
		OptBasedSmooth(e_left->prev()->vertex());
		assert(edge_temp_prev->opposite()->vertex() == v_left);
		//merge_edge(edge_temp_prev, 1);
		//e_left->prev()->vertex()->point() = ver_d;
		//write_my_mesh_nas("join_vertex.nas");
		v_left->set_halfedge(e_left->opposite());
		edge_begin = v_left->vertex_begin();
		edge_end = v_left->vertex_begin();
		do 
		{
			if (edge_begin->face()->testMark(HDS_mesh::Face::TRI))
			{
				Point_3 pa = edge_begin->prev()->vertex()->point();
				Point_3 pb = edge_begin->vertex()->point();
				Point_3 pc = edge_begin->next()->vertex()->point();
				CGAL::Orientation ori = CGAL::orientation(Point_2(pa.x(), pa.y()), Point_2(pb.x(), pb.y()), Point_2(pc.x(), pc.y()));
				if (ori == CGAL::RIGHT_TURN)
				{
					//write_my_mesh_nas("join_vertex.nas");
					OptBasedSmooth(e_left->prev()->vertex());
					//write_my_mesh_nas("join_vertex.nas");
					// 					if (!edge_begin->prev()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
					// 						smoothInteriorNode(edge_begin->prev()->vertex());
					// 					if (!edge_begin->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
					// 						smoothInteriorNode(edge_begin->vertex());
					// 					if (!edge_begin->next()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
					// 						smoothInteriorNode(edge_begin->next()->vertex());

					//write_my_mesh_nas("join_vertex.nas");
					// 					if (edge_begin->next()->testMark(HDS_mesh::Halfedge::FRONT))
					// 						merge_edge(edge_begin->opposite()->next(), 1);
					// 					else
					// 					{
					// 						if (edge_begin->next()->vertex()->testMark(HDS_mesh::Vertex::FrontNode))
					// 							merge_edge(edge_begin->opposite(), 1);
					// 						else
					// 							merge_edge(edge_begin->next(), 1);
					// 					}
					break;
				}
			}
			++edge_begin;
		} while (edge_begin != edge_end);
		//write_my_mesh_nas("join_vertex.nas");
	}
	if (!flag_temp)
		int aaaaaaaaa = 0;
	assert(flag_temp);
	return e_left;
}

HDS_mesh::Halfedge_handle QMorphGenerator::seamFrontTrans(HDS_mesh::Halfedge_handle& e_biger, 
														  HDS_mesh::Halfedge_handle& e_shorter, 
														  bool is_top)
{
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());

	HDS_mesh::Face_handle quad_face = e_biger->opposite()->face();
	HDS_mesh::Face_handle tri_face = e_biger->face();
	assert(quad_face->testMark(HDS_mesh::Face::QUAD));
	assert(tri_face->testMark(HDS_mesh::Face::TRI));

	bool flag_seam = false;
	if (e_shorter->loopNext() == e_biger)
	{
		flag_seam = true;
	}

	///< 计算e_biger的中点
	Point_3 e_mid = (e_biger->prev()->vertex()->point() + e_biger->vertex()->point()) / 2.0;
	///< 计算四边形面片的中心点
	Point_3 quad_mid(0, 0, 0);
	HDS_mesh::Halfedge_handle edge_begin_temp = quad_face->halfedge();
	HDS_mesh::Halfedge_handle edge_end_temp = quad_face->halfedge();
	do 
	{
		quad_mid = quad_mid + edge_begin_temp->vertex()->point();
		edge_begin_temp = edge_begin_temp->next();
	} while (edge_begin_temp != edge_end_temp);
	quad_mid = quad_mid / 4.0;

	// 先把上述两点插入网格
	HDS_mesh::Vertex_handle handle_V_e_mid, handle_V_quad_mid;
	HDS_mesh::Vertex v_e_mid, v_quad_mid;
	int id_e_mid = _max_index + 1;
	int id_quad_mid = _max_index + 2;
	_max_index = id_quad_mid;

	v_e_mid.setIndex(id_e_mid);
	v_e_mid.point() = e_mid;
	handle_V_e_mid = _source_mesh->hds().vertices_push_back(v_e_mid);///< 插入节点

	v_quad_mid.setIndex(id_quad_mid);
	v_quad_mid.point() = quad_mid;
	handle_V_quad_mid = _source_mesh->hds().vertices_push_back(v_quad_mid);///< 插入节点

	///< 创建十二条半边
	HDS_mesh::Halfedge h[12];
	h[0].set_vertex(handle_V_e_mid);  h[1].set_vertex(e_biger->next()->vertex());
	h[2].set_vertex(handle_V_e_mid);  h[3].set_vertex(e_biger->prev()->vertex());
	h[4].set_vertex(handle_V_e_mid);  h[5].set_vertex(handle_V_quad_mid);
	h[6].set_vertex(handle_V_e_mid);  h[7].set_vertex(e_biger->vertex());

	h[8].set_vertex(handle_V_quad_mid); 
	h[10].set_vertex(handle_V_quad_mid); 

	if (flag_seam)
	{
		h[9].set_vertex(e_biger->opposite()->prev()->prev()->vertex());
		h[11].set_vertex(e_biger->opposite()->vertex());
	}
	else
	{
		h[9].set_vertex(e_biger->opposite()->next()->vertex());
		h[11].set_vertex(e_biger->vertex());
	}

	///< 把十二条半边插入网格
	HDS_mesh::Halfedge_handle handle_H[12];
	for (int ih = 0; ih < 11; ih += 2)
	{///< 把新增的两条半边插入背景网格
		handle_H[ih] = _source_mesh->hds().edges_push_back(h[ih], h[ih + 1]);
		handle_H[ih + 1] = handle_H[ih]->opposite();
	}

	if (e_biger->vertex()->halfedge() == e_biger)
		e_biger->vertex()->set_halfedge(handle_H[7]);
	if(e_biger->opposite()->vertex()->halfedge() == e_biger->opposite())
		e_biger->opposite()->vertex()->set_halfedge(handle_H[3]);

	///< 设置各半边之间的拓扑关系
	handle_H[0]->set_next(handle_H[7]); handle_H[0]->set_prev(e_biger->next());
	handle_H[1]->set_next(e_biger->prev());   handle_H[1]->set_prev(handle_H[2]);
	handle_H[2]->set_next(handle_H[1]); handle_H[2]->set_prev(e_biger->prev());
	handle_H[3]->set_prev(handle_H[4]);

	if (flag_seam)
	{
		handle_H[3]->set_next(handle_H[10]);
		handle_H[4]->set_prev(handle_H[10]);
		handle_H[5]->set_next(handle_H[9]);
		handle_H[6]->set_prev(e_biger->opposite()->prev());
		handle_H[8]->set_next(handle_H[11]); 
		handle_H[8]->set_prev(e_biger->opposite()->next()->next());
		handle_H[9]->set_next(e_biger->opposite()->prev());
		handle_H[9]->set_prev(handle_H[5]);
		handle_H[10]->set_next(handle_H[4]);
		handle_H[10]->set_prev(handle_H[3]);
		handle_H[11]->set_next(e_biger->opposite()->next());
		handle_H[11]->set_prev(handle_H[8]);
		e_biger->opposite()->next()->set_prev(handle_H[11]);
		e_biger->opposite()->prev()->set_next(handle_H[6]);
		e_biger->opposite()->prev()->set_prev(handle_H[9]);
		e_biger->opposite()->next()->next()->set_next(handle_H[8]);
	}
	else
	{
		handle_H[3]->set_next(e_biger->opposite()->next());
		handle_H[4]->set_prev(handle_H[8]);
		handle_H[5]->set_next(handle_H[11]);
		handle_H[6]->set_prev(handle_H[11]);
		handle_H[8]->set_next(handle_H[4]);
		handle_H[8]->set_prev(e_biger->opposite()->next());
		handle_H[9]->set_next(e_biger->opposite()->next()->next());
		handle_H[9]->set_prev(handle_H[10]);
		handle_H[10]->set_next(handle_H[9]);
		handle_H[10]->set_prev(e_biger->opposite()->prev());
		handle_H[11]->set_next(handle_H[6]);
		handle_H[11]->set_prev(handle_H[5]);
		e_biger->opposite()->next()->next()->set_prev(handle_H[9]);
		e_biger->opposite()->next()->set_prev(handle_H[3]);
		e_biger->opposite()->next()->set_next(handle_H[8]);
		e_biger->opposite()->prev()->set_next(handle_H[10]); 
	}

	handle_H[4]->set_next(handle_H[3]); 
	handle_H[5]->set_prev(handle_H[6]);
	handle_H[6]->set_next(handle_H[5]); 
	handle_H[7]->set_next(e_biger->next());   handle_H[7]->set_prev(handle_H[0]);

	e_biger->next()->set_prev(handle_H[7]); e_biger->next()->set_next(handle_H[0]);
	e_biger->prev()->set_prev(handle_H[1]); e_biger->prev()->set_next(handle_H[2]);

	///< 创建三个三角形面片
	HDS_mesh::Face tri_face_new[3];
	HDS_mesh::Face_handle handle_tri_face_new[3];
	for (int nf = 0; nf < 3; ++nf)
	{
		tri_face_new[nf].setMark(HDS_mesh::Face::TRI); ///< 标记面片为三角形
		handle_tri_face_new[nf] = _source_mesh->hds().faces_push_back(tri_face_new[nf]); ///< 把三角形面片插入背景网格
	}

	if (flag_seam)
	{
		e_biger->prev()->set_face(handle_tri_face_new[0]);
		e_biger->next()->set_face(handle_tri_face_new[1]);
		handle_H[3]->set_face(handle_tri_face_new[2]);

		handle_tri_face_new[0]->set_halfedge(e_biger->prev());
		handle_tri_face_new[1]->set_halfedge(e_biger->next());
		handle_tri_face_new[2]->set_halfedge(handle_H[3]);
	}
	else
	{
		e_biger->next()->set_face(handle_tri_face_new[0]);
		e_biger->prev()->set_face(handle_tri_face_new[1]);
		handle_H[6]->set_face(handle_tri_face_new[2]);

		handle_tri_face_new[0]->set_halfedge(e_biger->next());
		handle_tri_face_new[1]->set_halfedge(e_biger->prev());
		handle_tri_face_new[2]->set_halfedge(handle_H[6]);
	}

	///< 创建两个四边形
	HDS_mesh::Face quad_face_new[2];
	HDS_mesh::Face_handle handle_quad_face_new[2];
	for (int nf = 0; nf < 2; ++nf)
	{
		quad_face_new[nf].setMark(HDS_mesh::Face::QUAD); ///< 标记面片为三角形
		handle_quad_face_new[nf] = _source_mesh->hds().faces_push_back(quad_face_new[nf]);///< 把四边形面片插入网格
	}

	e_biger->opposite()->next()->set_face(handle_quad_face_new[0]); 
	e_biger->opposite()->prev()->set_face(handle_quad_face_new[1]);

	handle_quad_face_new[0]->set_halfedge(e_biger->opposite()->next()); 
	handle_quad_face_new[1]->set_halfedge(e_biger->opposite()->prev());

	///< 设置各半边的入射面片
	if (flag_seam)
	{
		handle_H[0]->set_face(handle_tri_face_new[1]);
		handle_H[7]->set_face(handle_tri_face_new[1]);
		handle_H[1]->set_face(handle_tri_face_new[0]);
		handle_H[2]->set_face(handle_tri_face_new[0]);

		handle_H[4]->set_face(handle_tri_face_new[2]);
		handle_H[5]->set_face(handle_quad_face_new[1]);
		handle_H[6]->set_face(handle_quad_face_new[1]);
		handle_H[10]->set_face(handle_tri_face_new[2]);
		handle_H[11]->set_face(handle_quad_face_new[0]);

		e_biger->opposite()->next()->next()->set_face(handle_quad_face_new[0]);

		handle_H[10]->setMark(HDS_mesh::Halfedge::FRONT);
		handle_H[4]->setMark(HDS_mesh::Halfedge::FRONT);
		handle_H[7]->setMark(HDS_mesh::Halfedge::FRONT);

		e_biger->loopPrev()->set_loopNext(handle_H[10]);
		handle_H[10]->set_loopPrev(e_biger->loopPrev());
		handle_H[4]->set_loopPrev(handle_H[10]);
		handle_H[10]->set_loopNext(handle_H[4]);
		handle_H[4]->set_loopNext(handle_H[7]);
		handle_H[7]->set_loopPrev(handle_H[4]);
		handle_H[7]->set_loopNext(e_biger->loopNext());
		e_biger->loopNext()->set_loopPrev(handle_H[7]);

		handle_H[10]->set_FrontState(3);

		handle_V_quad_mid->set_halfedge(handle_H[10]);
		handle_V_e_mid->set_halfedge(handle_H[4]); ///< 设置节点的关联半边

		e_shorter->opposite()->prev()->set_loopNext(handle_H[10]);
	}
	else
	{
		handle_H[0]->set_face(handle_tri_face_new[0]);
		handle_H[7]->set_face(handle_tri_face_new[0]);
		handle_H[1]->set_face(handle_tri_face_new[1]);
		handle_H[2]->set_face(handle_tri_face_new[1]);
		handle_H[3]->set_face(handle_quad_face_new[0]);
		handle_H[4]->set_face(handle_quad_face_new[0]);
		handle_H[5]->set_face(handle_tri_face_new[2]);

		handle_H[10]->set_face(handle_quad_face_new[1]);
		handle_H[11]->set_face(handle_tri_face_new[2]);

		handle_H[9]->next()->set_face(handle_quad_face_new[1]);

		handle_H[2]->setMark(HDS_mesh::Halfedge::FRONT);
		handle_H[5]->setMark(HDS_mesh::Halfedge::FRONT);
		handle_H[11]->setMark(HDS_mesh::Halfedge::FRONT);

		e_biger->loopPrev()->set_loopNext(handle_H[2]);
		handle_H[2]->set_loopPrev(e_biger->loopPrev());
		handle_H[5]->set_loopPrev(handle_H[2]);
		handle_H[2]->set_loopNext(handle_H[5]);
		handle_H[5]->set_loopNext(handle_H[11]);
		handle_H[11]->set_loopPrev(handle_H[5]);
		handle_H[11]->set_loopNext(e_shorter);
		e_shorter->set_loopPrev(handle_H[11]);

		handle_H[11]->set_FrontState(3);

		handle_V_quad_mid->set_halfedge(handle_H[5]);
		handle_V_e_mid->set_halfedge(handle_H[2]); ///< 设置节点的关联半边
	}

	handle_H[8]->set_face(handle_quad_face_new[0]); handle_H[9]->set_face(handle_quad_face_new[1]);

	if (e_biger->testMark(HDS_mesh::Halfedge::UNACTIVE))
		deleteFromUnactiveFront(e_biger);
	else
		deleteFromFront(e_biger);

	_source_mesh->hds().faces_erase(tri_face);
	_source_mesh->hds().faces_erase(quad_face);

	_source_mesh->hds().edges_erase(e_biger);

	HDS_mesh::Vertex_handle top_left = handle_V_e_mid;
	HDS_mesh::Vertex_handle top_right = e_shorter->vertex();
	if (flag_seam)
	{
		top_left = e_shorter->prev()->vertex();
		top_right = handle_V_e_mid;
	}

	HDS_mesh::Halfedge_handle top_edge = get_Top(top_left, top_right);
	if (top_edge == HDS_mesh::Halfedge_handle())
	{
		if (flag_seam)
		{
			handle_H[10]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[10]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[10]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[10]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[10]);

			handle_H[4]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[4]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[4]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[4]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[4]);

			handle_H[7]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[7]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[7]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[7]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[7]);
		}
		else
		{
			handle_H[2]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[2]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[2]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[2]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[2]);

			handle_H[5]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[5]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[5]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[5]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[5]);

			handle_H[11]->setMark(HDS_mesh::Halfedge::UNACTIVE);
			handle_H[11]->vertex()->setMark(HDS_mesh::Vertex::QuadNode);
			handle_H[11]->vertex()->setMark(HDS_mesh::Vertex::FrontNode);
			handle_H[11]->set_FrontState(0);
			_unactive_front_list.push_back(handle_H[11]);
		}
		return HDS_mesh::Halfedge_handle();
	}
	//write_my_mesh_nas("edge_recover.nas");
	HDS_mesh::Face_handle new_quad;
	if (flag_seam)
		new_quad = constructNewFace(handle_H[10], handle_H[4], top_edge, e_shorter);
	else
		new_quad = constructNewFace(handle_H[11], e_shorter, top_edge, handle_H[5]);

	std::list<HDS_mesh::Vertex_handle> smoothed_node_list;

	getVertexToSmooth(new_quad, smoothed_node_list);

	localSmoothing(smoothed_node_list);

	e_shorter->setMark(HDS_mesh::Halfedge::ISSPLIT);
	return e_shorter;
}


bool QMorphGenerator::checkFrontToSeam(std::list<HDS_mesh::Halfedge_handle>& influenced_front_edge_list,
									   std::list<HDS_mesh::Vertex_handle>& node_need_smooth,
									   int n)
{
	bool flag = false; ///< 标记是否执行过缝合操作

	std::set<HDS_mesh::Vertex_handle, Vertex_Comp> nodeSet; ///< 无重复地存储前沿点
	std::list<HDS_mesh::Halfedge_handle>::iterator it_edge_temp;
	for (it_edge_temp = influenced_front_edge_list.begin(); it_edge_temp != influenced_front_edge_list.end(); ++it_edge_temp)
	{
		HDS_mesh::Halfedge_handle current_edge = *it_edge_temp;
		nodeSet.insert(current_edge->vertex());
		nodeSet.insert(current_edge->prev()->vertex());
	}

	std::list<HDS_mesh::Vertex_handle> vertex_list_temp;
	std::set<HDS_mesh::Vertex_handle, Vertex_Comp>::iterator it_node;
	for (it_node = nodeSet.begin(); it_node != nodeSet.end(); ++ it_node)
	{
		(*it_node)->setMark(HDS_mesh::Vertex::SEAMCHECK);
		vertex_list_temp.push_back(*it_node);
	}
	std::list<HDS_mesh::Vertex_handle>::iterator it_vertex;
	for (it_vertex = vertex_list_temp.begin(); it_vertex != vertex_list_temp.end(); ++it_vertex)
	{
		HDS_mesh::Vertex_handle v_temp = *it_vertex;
		if (v_temp->isDead())
			continue;
		if (!v_temp->testMark(HDS_mesh::Vertex::FrontNode))
			continue;
		if (v_temp ->testMark(HDS_mesh::Vertex::ERASED))
			continue;

		assert(v_temp->testMark(HDS_mesh::Vertex::FrontNode));
		std::queue<HDS_mesh::Vertex_handle> front_node_queue_temp;

		front_node_queue_temp.push(v_temp);
		while (!front_node_queue_temp.empty())
		{
			HDS_mesh::Vertex_handle v_to_test = front_node_queue_temp.front();
			front_node_queue_temp.pop();
			if (v_to_test ->testMark(HDS_mesh::Vertex::ERASED))
				continue;
			nodeSet.insert(v_to_test);
			if (v_to_test->isDead())
				continue;
			//std::cout << v_to_test->getIndex() << std::endl;
			std::list<HDS_mesh::Halfedge_handle> front_adjecent_edge_list;
			HDS_mesh::Halfedge_around_vertex_circulator edge_begin = v_to_test->vertex_begin();
			HDS_mesh::Halfedge_around_vertex_circulator edge_end = v_to_test->vertex_begin();
			do 
			{
				if (edge_begin->testMark(HDS_mesh::Halfedge::FRONT))
				{
					front_adjecent_edge_list.push_back(edge_begin);
				}
				++edge_begin;
			} while (edge_begin != edge_end);
			if (front_adjecent_edge_list.size() == 0)
			{
				v_to_test->resetMark(HDS_mesh::Vertex::FrontNode);
				continue;
			}
			for (it_edge_temp = front_adjecent_edge_list.begin(); 
				it_edge_temp != front_adjecent_edge_list.end(); 
				++it_edge_temp)
			{
				///< 计算当前前沿与其后继前沿之间的夹角
				HDS_mesh::Halfedge_handle loop_next = (*it_edge_temp)->loopNext();
				if (loop_next->testMark(15))
					continue;
				Point_3 p_p = (*it_edge_temp)->prev()->vertex()->point();
				Point_3 p_m = (*it_edge_temp)->vertex()->point();
				Point_3 p_n = loop_next->vertex()->point();
				double angle_temp = angle_evaluate(p_p, p_m, p_n);
				if (angle_temp < SEAM_ANGLE)
				{
					if ((*it_edge_temp)->loopNext()->loopNext()->loopNext() == (*it_edge_temp))
						continue;
					if ((*it_edge_temp)->testMark(HDS_mesh::Halfedge::BOUNDARY) || loop_next->testMark(HDS_mesh::Halfedge::BOUNDARY))
						continue;
					//front_node_queue_temp.push(loop_next->vertex());
					if (canSeam(*it_edge_temp))
					{
						HDS_mesh::Halfedge_handle new_edge = seamFront((*it_edge_temp), loop_next);
						if (new_edge == HDS_mesh::Halfedge_handle())
							continue;
						// 						if (n == 1046)
						// 							write_my_mesh_nas("seam.nas");
						new_edge->vertex()->setMark(HDS_mesh::Vertex::SEAMCHECK);
						front_node_queue_temp.push(new_edge->vertex());
						new_edge->prev()->vertex()->setMark(HDS_mesh::Vertex::SEAMCHECK);
						front_node_queue_temp.push(new_edge->prev()->vertex());
						flag = true;
					}
				}
			}
		}
	}

	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	for (it_node = nodeSet.begin(); it_node != nodeSet.end(); ++ it_node)
	{///< 将已被标记删除的节点从网格中删除
		if ((*it_node) ->testMark(HDS_mesh::Vertex::ERASED))
			continue;
		(*it_node)->resetMark(HDS_mesh::Vertex::SEAMCHECK);
		if (flag)
		{
			if ((*it_node)->isDead())
				decorator.vertices_erase(*it_node);
			else
				node_need_smooth.push_back(*it_node);
		}
	}
	nodeSet.clear();
	return flag;
}

bool QMorphGenerator::tri_mesh_clean()
{
	CGAL::HalfedgeDS_decorator<HalfedgeDS> decorator(_source_mesh->hds());
	HDS_mesh::Vertex_iterator it_v;
	int n_merge = 0;
	do 
	{
		n_merge = 0;
		for (it_v = _source_mesh->vertices_begin(); it_v != _source_mesh->vertices_end();)
		{
			std::list<HDS_mesh::Face_handle> face_list_temp;
			if (!it_v->testMark(HDS_mesh::Vertex::BOUNDARY))
			{
				if (it_v->degree() == 3)
				{
					HDS_mesh::Halfedge_handle e1 = it_v->halfedge()->prev();
					HDS_mesh::Halfedge_handle e2 = it_v->halfedge()->opposite()->next();
					HDS_mesh::Halfedge_handle e3 = it_v->halfedge()->next()->opposite()->prev();

					HDS_mesh::Halfedge_handle ed1 = e1->next();
					HDS_mesh::Halfedge_handle ed2 = e2->next();
					HDS_mesh::Halfedge_handle ed3 = e3->next();
					decorator.remove_tip(e1);
					decorator.remove_tip(e2);
					decorator.remove_tip(e3);
					_source_mesh->hds().edges_erase(ed1);
					_source_mesh->hds().edges_erase(ed2);
					_source_mesh->hds().edges_erase(ed3);
					if (decorator.get_face(e2) != HDS_mesh::Face_handle())
						decorator.faces_erase(decorator.get_face( e2));
					if (decorator.get_face(e3) != HDS_mesh::Face_handle())
						decorator.faces_erase(decorator.get_face( e3));
					decorator.set_face(e2, decorator.get_face(e1));
					decorator.set_face(e3, decorator.get_face(e1));

					if (decorator.get_face(e1) != HDS_mesh::Face_handle())
						decorator.set_face_halfedge(e1);

					decorator.set_vertex_halfedge(e1);
					decorator.set_vertex_halfedge(e2);
					decorator.set_vertex_halfedge(e3);
					_source_mesh->hds().vertices_erase(it_v++);
					++n_merge;
					continue;
				}			
			}
			++it_v;
		}
	} while (n_merge > 0);
	return true;
}

void QMorphGenerator::mesh()
{
	tri_mesh_clean();

	classify_Front_Edge(); ///< 对活动前沿中的前沿边按状态值进行分类，对每一类前沿边按边长升序排列
	HDS_mesh::Halfedge_handle current_front; ///< 获取前沿边
	get_Next_Front(current_front);
	int n = 0;
	bool top_error = false;
//#if 0

	while (1)
	{
		++n;
		HDS_mesh::Halfedge_handle side_left;
		HDS_mesh::Halfedge_handle side_right;
		HDS_mesh::Vertex_handle top_left;
		HDS_mesh::Vertex_handle top_right;

		if (current_front->loopNext()->loopNext()->loopNext()->loopNext() == current_front)
		{
			side_left = current_front->loopPrev();
			side_right = current_front->loopNext();

			top_left = side_left->prev()->vertex();
			top_right = current_front->loopNext()->vertex();

			top_right->resetMark(HDS_mesh::Vertex::LOOPCUT);
			top_left->resetMark(HDS_mesh::Vertex::LOOPCUT);
		}
		else
		{
			// 			std::cout << "===========================================\n";
			// 			std::cout << current_front->vertex()->point() << std::endl;
			if (!get_Side_Left(current_front, side_left))
			{
				rePushFront(current_front);
				get_Next_Front(current_front);
				--n;
				continue;
			}///< 获取左侧边界边

			if (side_left->vertex() != current_front->prev()->vertex())
				side_left = side_left->opposite();
			top_left = side_left->prev()->vertex();
			top_left->setMark(HDS_mesh::Vertex::HASCHOICED);

			///< 获取右侧边界边
			if (!get_Side_Right(current_front, side_right))
			{
				top_left->resetMark(HDS_mesh::Vertex::HASCHOICED);
				rePushFront(current_front);
				get_Next_Front(current_front);
				--n;
				continue;
			}

			if (side_right->prev()->vertex() != current_front->vertex())
				side_right = side_right->opposite();
			top_right = side_right->vertex();

			top_left->resetMark(HDS_mesh::Vertex::HASCHOICED);

			top_right->resetMark(HDS_mesh::Vertex::HASCHOICED);
		}

		HDS_mesh::Halfedge_handle top_edge = get_Top(top_left, top_right);

		if (top_edge == HDS_mesh::Halfedge_handle())
		{
			rePushFront(current_front);
			get_Next_Front(current_front);
			--n;	
			continue;
		}

		if (top_edge == side_left)
		{
			rePushFront(current_front);
			get_Next_Front(current_front);
			--n;	
			continue;
		}

		HDS_mesh::Face_handle new_quad = constructNewFace(current_front, side_right, top_edge, side_left);

		std::list<HDS_mesh::Vertex_handle> smoothed_node_list;
		getVertexToSmooth(new_quad, smoothed_node_list);

		localSmoothing(smoothed_node_list);

		std::list<HDS_mesh::Halfedge_handle> influenced_front_edge_list;
		std::list<HDS_mesh::Vertex_handle> node_need_smooth;
		getInfluencedFrontEdge(smoothed_node_list, influenced_front_edge_list);
		if (checkFrontToSeam(influenced_front_edge_list, node_need_smooth, n))
		{
			localSmoothing(node_need_smooth);
			getInfluencedFrontEdge(node_need_smooth, influenced_front_edge_list);
		}
		updateFront(influenced_front_edge_list);

		if (3995 == n)
		{
			//return;
			//write_my_mesh_nas("0811.nas");
			int xxxx = 0;
			std::cout << "fsdfs\n";
			int cca = 1;
		}
		if (!get_Next_Front(current_front))
			break;
		while (current_front->loopNext()->loopNext()->loopNext() == current_front)
		{
			current_front->resetMark(HDS_mesh::Halfedge::FRONT);
			current_front->loopNext()->resetMark(HDS_mesh::Halfedge::FRONT);
			current_front->loopNext()->loopNext()->resetMark(HDS_mesh::Halfedge::FRONT);
			current_front->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
			current_front->loopNext()->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
			current_front->loopNext()->loopNext()->vertex()->resetMark(HDS_mesh::Vertex::FrontNode);
			if (!get_Next_Front(current_front))
				break;
		}
		if (!current_front->testMark(HDS_mesh::Halfedge::FRONT))
		{
			if (!get_Next_Front(current_front))
				break;
		}
	}
//#endif
	_optimizer->mesh_OPT();
}