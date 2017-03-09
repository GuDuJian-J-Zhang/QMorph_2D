#include "QuadMesher.h"

#include "quadMeshAPI.h"
#include "mg_def.h"
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

class PropertyTree : public boost::property_tree::ptree 
{
public:
	typedef boost::property_tree::ptree Base;
	PropertyTree(){}
	PropertyTree(const Base& atree)
		:Base(atree){}

	bool loadJsonStream(std::istream& stream)
	{
		try
		{
			boost::property_tree::json_parser::read_json<Base>(stream, *this);
		}catch(boost::property_tree::ptree_error&)
		{
			return false;
		}
		return true;
	}

	bool get(const std::string& str,PropertyTree& subTree)const
	{
		try
		{
			subTree=this->get_child(str);
		}catch(...)
		{
			return false;
		}
		return true;
	}

	bool get(const std::string& path,std::string& val)const
	{
		try
		{
			val=this->get_child(path).data();
		}catch(...)
		{
			return false;
		}
		return true;
	}

	bool get(const std::string& path,int& val)const
	{
		std::string valStr;
		if(!get(path,valStr))
			return false;
		val = atoi(valStr.c_str());
		return true;
	}

	bool get(const std::string& path,double& val)const
	{
		std::string valStr;
		if(!get(path,valStr))
			return false;
		val = atof(valStr.c_str());
		return true;
	}
};

ErrorCodes QMorph2D(MR_NodeArray* nodes, MR_ElemArray* elems,const char* ctrls)
{
	PropertyTree json_tree;
	std::stringstream jstr(ctrls);
	json_tree.loadJsonStream(jstr);
	std::string method = "";
	json_tree.get("data.options.surfaceMesh.method", method);
	std::cout << method;

	double* myNodes = new double[2 * nodes->length];
	int* myEdges = new int[2 * elems->length];
	for (int ni = 0; ni < nodes->length; ++ni)
	{
		myNodes[2 * ni] = nodes->array[ni].xyz[0];
		myNodes[2 * ni + 1] = nodes->array[ni].xyz[1];
	}
	for (int ei = 0; ei < elems->length; ++ei)
	{
		myEdges[2 * ei] = elems->array[ei].nodes[0] - 1;
		myEdges[2 * ei + 1] = elems->array[ei].nodes[1] - 1;
	}
	QMorphInput input(nodes->length, myNodes,elems->length, myEdges);
	DelaunayAFT2DInterface mesher2d;
	int rt=mesher2d.doMesh(&input);
	if(rt!=0)
		return rt;
	HDS_mesh backoundMesh;
	QuadMesher myMesher(backoundMesh);
	myMesher.read_tri_mesh(mesher2d);
	myMesher.mesh(method.c_str());
// 	int rt = mesher2d.doMeshInp(file_name);
// 	if(rt!=0)
// 		return rt;
// 	HDS_mesh backoundMesh;
// 	QuadMesher myMesher(backoundMesh);
// 	myMesher.read_tri_mesh(mesher2d);
// 
// 	myMesher.mesh();
// 	myMesher.tri_mesh_clean();
// 	myMesher.quad_generation_QMorph();
// 	myMesher.quad_mesh_OPT();

// 	MG_NS::ArrayProxy<double> nodeProxy(nodes.array,nodes.length);
// 
// 	size_t nnode = backoundMesh.size_of_vertices();
// 	nodeProxy.resize(nnode * 2);
// 	HDS_mesh::Vertex_iterator it_vertex;
// 	int idx = 0;
// 	for(it_vertex=backoundMesh.vertices_begin(); 
// 		it_vertex!=backoundMesh.vertices_end();
// 		++it_vertex)
// 	{
// 		it_vertex->setIndex(idx);
// 		nodeProxy[idx*2] = it_vertex->point().x(); 
// 		nodeProxy[idx*2+1] = it_vertex->point().y(); 
// 		++idx;
// 	}
// 
// 	MG_NS::ArrayProxy<int> quad_face(quads.array,quads.length);
// 	MG_NS::ArrayProxy<int> tri_face(tris.array,tris.length);
// 
// 	HDS_mesh::Face_iterator it_face;
// 	for (it_face = backoundMesh.facets_begin(); it_face != backoundMesh.facets_end(); ++it_face)
// 	{
// 		HDS_mesh::Halfedge_iterator edge_begin = it_face->halfedge();
// 		HDS_mesh::Halfedge_iterator edge_end   = it_face->halfedge();
// 		if (it_face->is_triangle())
// 		{
// 			do 
// 			{
// 				tri_face.push_back(edge_begin->vertex()->getIndex());
// 				edge_begin = edge_begin->next();
// 			} while (edge_begin != edge_end);
// 		}
// 		else if (it_face->is_quad())
// 		{
// 			do 
// 			{
// 				quad_face.push_back(edge_begin->vertex()->getIndex());
// 				edge_begin = edge_begin->next();
// 			} while (edge_begin != edge_end);
// 		}
// 	}
	myMesher.write_my_mesh_nas("0811.nas");
	return 0;
}