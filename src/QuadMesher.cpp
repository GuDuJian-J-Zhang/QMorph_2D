// **************************************************************************
// File:     QuadMesher.cpp
// Author:   Zhang Jun
// Date:     2016-05-05
// Email:    zhangjun_dg@mail.dlut.edu.cn
// **************************************************************************
#include "QuadMesher.h"
#include "Vector3D.h"
#include "sldelaunay2d.h"
#include "sldelaunay2d_std.h"

#include "TriMergeGenerator.h"
#include "QMorphGenerator.h"

//#include <process.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <boost/algorithm/string.hpp>

// #define TOLERENCE 1.0e-6
// 
// #define MINLENGTH 1.0e30
// 
// #define MAXLENGTH -1.0e30
// 
// #define RATIO_RANGE_LEFT 2.5 ///< 用于前沿点光顺
// 
// #define RATIO_RANGE_RIGHT 20 ///< 用于前沿点光顺
// 
// #define PAVING_SIDE 240 ///< Paving算法中边节点处前沿边夹角的上限
// 
// #define PAVINE_CORNER 300 ///< Paving算法中角节点处前沿边夹角的上限
// 
// #define PAVING_REVERSAL 360///< Paving算法中转节点处前沿边夹角的上限
// 
// #define SEAM_ANGLE 45 ///< 缝合角(当相邻前沿夹角小于该值时执行前沿缝合操作)
// 
// #define ASPECT_RATIO 0.2 ///< 四边形单元纵横比下限


template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> 
{
public:
	Build_triangle(DelaunayAFT2DInterface& tri_mesh): _source_tri_mesh(tri_mesh) {}
	void operator()( HDS& hds)
	{
		// Postcondition: hds is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true);
		builder.begin_surface(1000, 1000, 3000);
		typedef typename HDS::Vertex          Vertex;
		typedef typename HDS::Vertex_handle   Vertex_handle;
		typedef typename Vertex::Point        Point;
		int id; 
		double uvxyz[5];
		DelaunayAFT2DInterface::ResultIterator mesherIter(&_source_tri_mesh);
		std::map<int, int> id_to_index;
		mesherIter.node_begin();
		int index = 0;
		while(1)
		{
			if(!mesherIter.node_next(id,2,uvxyz)) break;
			if(id<0) 
			{
				continue;
			}
			double x = uvxyz[0];
			double y = uvxyz[1];
			double z = 0.0;
			Vertex_handle v_temp = builder.add_vertex(Point(x, y, z));
			v_temp->setIndex(id);
			id_to_index[id] = index;
			++index;
		}
		mesherIter.triangle_begin();

		int nodes[3];
		mesherIter.triangle_begin();
		while(1)
		{
			if(!mesherIter.triangle_next(nodes)) break;
			int ncount = 0;
			int jfield = 0;
			std::list<int> tri_temp;
			for (jfield = 4; jfield <= 9; jfield++)
			{
				ncount++;
				if (ncount>3) break;
				int idi = nodes[ncount-1];
				if(idi<0) idi = -idi;
				tri_temp.push_back(id_to_index.find(idi)->second);
			} 
			builder.add_facet(tri_temp.begin(), tri_temp.end());
		}
		builder.end_surface();
	}
private:
	DelaunayAFT2DInterface& _source_tri_mesh; 
};

template <class HDS>
class Build_triangle_INP : public CGAL::Modifier_base<HDS> 
{
public:
	Build_triangle_INP(std::string& tri_mesh_file): _tri_mesh_file(tri_mesh_file) {}
	void operator()( HDS& hds)
	{
		// Postcondition: hds is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true);
		builder.begin_surface(1000, 1000, 3000);
		typedef typename HDS::Vertex          Vertex;
		typedef typename HDS::Vertex_handle   Vertex_handle;
		typedef typename Vertex::Point        Point;
		
		std::map<int, int> id_to_index;
		int index = 0;
		std::ifstream infile_point_tri;
		infile_point_tri.open(_tri_mesh_file.c_str(),std::ios::in);
		std::string buffer;
		while (std::getline(infile_point_tri, buffer) && buffer != "*Node")
		{
			;
		}
		//读取节点坐标
		std::stringstream sstr;
		while (std::getline(infile_point_tri, buffer) && !boost::contains(buffer, "*Element"))
		{
			std::vector<std::string> result;
			boost::split(result, buffer, boost::is_any_of(","), boost::token_compress_on);
			sstr.str("");
			sstr.clear();
			for (int i = 0; i < result.size(); ++i)
				sstr << result[i];
			Point_3 p_temp;
			int id;
			sstr >> id;
			double x = 0, y = 0, z = 0;
			sstr >> x >> y;
			if (3 == result.size())
				sstr >> z;
			Vertex_handle v_temp = builder.add_vertex(Point(x, y, z));
			v_temp->setIndex(index);
			id_to_index[id] = index;
			++index;
		}

		if (index == 0)
		{
			std::cerr << "There are no points! File error...\n ";
			exit(-1);
		}

		int n_face = 0;
		//读取三角面片
		while (std::getline(infile_point_tri, buffer) && buffer != "*End Part")
		{
			std::vector<std::string> result;
			boost::split(result, buffer, boost::is_any_of(","), boost::token_compress_on);
			assert(result.size() == 4);
			sstr.str("");
			sstr.clear();
			for (int i = 0; i < 4; ++i)
			{
				sstr << result[i];
			}
			int id;
			sstr >> id;
			int n0 = 0;
			int n1 = 0;
			int n2 = 0;
			sstr >> n0 >> n1 >> n2;

			n0 = id_to_index.find(n0)->second;
			n1 = id_to_index.find(n1)->second;
			n2 = id_to_index.find(n2)->second;

			std::list<int> tri_temp;
			tri_temp.push_back(n0); tri_temp.push_back(n1); tri_temp.push_back(n2);
			builder.add_facet(tri_temp.begin(), tri_temp.end());
			++n_face;
		}
		builder.end_surface();

		if (n_face == 0)
		{
			std::cerr << "There are no faces! File error...\n ";
			exit(-1);
		}
		infile_point_tri.close();
	}
private:
	std::string& _tri_mesh_file; 
};

QuadMesher::QuadMesher(HDS_mesh &backgroundMesh) : _backgroundMesh(backgroundMesh)
{}

bool QuadMesher::mesh(const char* method)
{
	if (!method)
		return false;
	//std::cout << method << std::endl;
	if (strcmp(method, "TriMerge") == 0)
		_meshGenerator = new TriMergeGenerator;
	else
		_meshGenerator = new QMorphGenerator;
	_meshGenerator->set_SourceMesh(&_backgroundMesh);
	_meshGenerator->mesh();
	return true;
}

bool QuadMesher::read_file(std::string file_name)
{
	_backgroundMesh.delegate(Build_triangle_INP<HalfedgeDS>(file_name));
	HDS_mesh::Face_handle it_face;
	for (it_face = _backgroundMesh.facets_begin(); it_face != _backgroundMesh.facets_end(); ++it_face)
	{
		it_face->setMark(HDS_mesh::Face::TRI);
	}
	std::cout << "Number of points = " << _backgroundMesh.size_of_vertices() << std::endl;
	std::cout << "Number of Tris = " << _backgroundMesh.size_of_facets() << std::endl;
	std::cout << "Read file over.\n";
	return true;
}

bool QuadMesher::read_tri_mesh(DelaunayAFT2DInterface& tri_mesh)
{
	//Build_triangle<HalfedgeDS> source_mesh(tri_mesh);
	_backgroundMesh.delegate(Build_triangle<HalfedgeDS>(tri_mesh));
	HDS_mesh::Face_handle it_face;
	for (it_face = _backgroundMesh.facets_begin(); it_face != _backgroundMesh.facets_end(); ++it_face)
	{
		it_face->setMark(HDS_mesh::Face::TRI);
	}
	std::cout << "Number of points = " << _backgroundMesh.size_of_vertices() << std::endl;
	std::cout << "Number of Tris = " << _backgroundMesh.size_of_facets() << std::endl;
	std::cout << "Read file over.\n";
	return true;
}

void QuadMesher::print_node_long_format_nas(int id,double x,double y,double z,char* strBuffer)
{
	char bufx[32],bufy[32],bufz[32];
	sprintf(bufx,"%16.16lf",x);bufx[15]='\0';
	sprintf(bufy,"%16.16lf",y);bufy[15]='\0';
	sprintf(bufz,"%16.16lf",z);bufz[15]='\0';
	sprintf(strBuffer,"%-16s%8d%-16s%-16s%-16s%-1s\n%-8s%-16s\n", "GRID*", id, "", bufx,bufy, "*" ,"*",bufz);
}

void QuadMesher::write_my_mesh_nas(std::string file_name)
{
	std::ofstream out(file_name);

	char buffer[256];
	out << "$$\n$$  GRID Data\n$$\n";
	int id = 0;
	int m = 0;
	HDS_mesh::Vertex_iterator it_vertex;
	for(it_vertex=_backgroundMesh.vertices_begin(); 
		it_vertex!=_backgroundMesh.vertices_end();
		++it_vertex)
	{
// 		if (id == 564)
// 			int aaaaaa = 0;
 		it_vertex->setIndex(id++);
		print_node_long_format_nas(it_vertex->getIndex() + 1,it_vertex->point().x(),
			                          it_vertex->point().y(),
									  it_vertex->point().z(),buffer);
		out << buffer;
	}

	out<<"$\n";

	sprintf(buffer,"CTRIA3  ");
	out<<"$  CTRIA3 Data\n$\n";

	HDS_mesh::Face_iterator it_face;
	int nf = 0; ///< 单元域数
	int nTri = 0; ///< 三角形单元总数
	int nQuad = 0; ///< 四边形单元总数
	int prop=309;
	for (it_face = _backgroundMesh.facets_begin(); it_face != _backgroundMesh.facets_end(); ++it_face)
	{
		HDS_mesh::Halfedge_iterator edge_begin = it_face->halfedge();
		HDS_mesh::Halfedge_iterator edge_end   = it_face->halfedge();
		int ne = 0;
		//if (it_face->testMark(HDS_mesh::Facet::TRI))
		if (it_face->is_triangle())
		{
			++nTri;
			out << "CTRIA3  " << std::setw(8) << ++nf << std::setw(8) << prop;
		}
		else if (it_face->is_quad())
		{
			++nQuad;
			out << "CQUAD4  " << std::setw(8) << ++nf << std::setw(8) << prop + 1;
		}
		else
		{
			++nQuad;
			out << "CQUAD4  " << std::setw(8) << ++nf << std::setw(8) << prop + 1;
		}
		do 
		{
			out << std::setw(8) << edge_begin->vertex()->getIndex() + 1;
			edge_begin = edge_begin->next();
		} while (edge_begin != edge_end);
		out << "\n";
	}
	out.close();
	std::cout << "Number of Tri = " << nTri << std::endl;
	std::cout << "Number of Quad = " << nQuad << std::endl;
}

QMorphInput::QMorphInput(unsigned nnode,const double* nodes,
						 unsigned  nedge,const int* edges)
						 :_nnode(nnode),_nedge(nedge),_nodes(nodes),_edges(edges)
{
	_ei=_ni=0;
}

bool QMorphInput::node_next(int& id,double xy[2])
{
	if(_ni>=_nnode)
		return false;
	xy[0]=_nodes[_ni*2];
	xy[1]=_nodes[_ni*2+1];
	++_ni;
	id=_ni;
	return true;
}

bool QMorphInput::edge_next(int nodes[2])
{
	if(_ei>=_nedge)
		return false;
	nodes[0]=_edges[_ei*2]+1;
	nodes[1]=_edges[_ei*2+1]+1;
	++_ei;
	return true;
}