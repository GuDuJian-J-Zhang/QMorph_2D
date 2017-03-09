#include "QuadMesher.h"

#include "qmorphAPI.h"
#include "mg_def.h"

ErrorCodes QMorph2D(MR_DoubleArray& nodes,
					const MR_IntArray& edges,
					MR_IntArray& quads,
					MR_IntArray& tris)
{
	QMorphInput input(nodes.length/2,nodes.array,edges.length/2,edges.array);
	DelaunayAFT2DInterface mesher2d;
	int rt=mesher2d.doMesh(&input);
	if(rt!=0)
		return rt;
	HDS_mesh backoundMesh;
	QMorphMesher myMesher(backoundMesh);
	myMesher.read_tri_mesh(mesher2d);
	//myMesher.write_my_mesh_nas("0811.nas");
// 	myMesher.tri_mesh_clean();
// 	myMesher.quad_generation_QMorph();
// 	myMesher.quad_mesh_OPT();
	myMesher.mesh();

	MG_NS::ArrayProxy<double> nodeProxy(nodes.array,nodes.length);

	size_t nnode = backoundMesh.size_of_vertices();
	nodeProxy.resize(nnode * 2);
	HDS_mesh::Vertex_iterator it_vertex;
	int idx = 0;
	for(it_vertex=backoundMesh.vertices_begin(); 
		it_vertex!=backoundMesh.vertices_end();
		++it_vertex)
	{
		it_vertex->setIndex(idx);
		nodeProxy[idx*2] = it_vertex->point().x(); 
		nodeProxy[idx*2+1] = it_vertex->point().y(); 
		++idx;
	}

	MG_NS::ArrayProxy<int> quad_face(quads.array,quads.length);
	MG_NS::ArrayProxy<int> tri_face(tris.array,tris.length);

	HDS_mesh::Face_iterator it_face;
	for (it_face = backoundMesh.facets_begin(); it_face != backoundMesh.facets_end(); ++it_face)
	{
		HDS_mesh::Halfedge_iterator edge_begin = it_face->halfedge();
		HDS_mesh::Halfedge_iterator edge_end   = it_face->halfedge();
		if (it_face->is_triangle())
		{
			do 
			{
				tri_face.push_back(edge_begin->vertex()->getIndex());
				edge_begin = edge_begin->next();
			} while (edge_begin != edge_end);
		}
		else if (it_face->is_quad())
		{
			do 
			{
				quad_face.push_back(edge_begin->vertex()->getIndex());
				edge_begin = edge_begin->next();
			} while (edge_begin != edge_end);
		}
	}
	myMesher.write_my_mesh_nas("0811.nas");
	return 0;
}

ErrorCodes QMorph2D(const char* file_name)
{
	DelaunayAFT2DInterface mesher2d;

	int rt = mesher2d.doMeshInp(file_name);
	if(rt!=0)
		return rt;
	HDS_mesh backoundMesh;
	QuadhMesher myMesher(backoundMesh);
	myMesher.read_tri_mesh(mesher2d);

	myMesher.mesh();
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