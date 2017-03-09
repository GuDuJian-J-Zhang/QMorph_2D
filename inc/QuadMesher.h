//**************************************************************************
// File:     QuadMesher.h
// Author:   Zhang Jun
// Date:     2016-05-05
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 半边结构下三角形面片编号规则(相对的顶点和边的编号一致  逆时针方向)
//**************************************************************************

#ifndef _QMORPHMESHER_H
#define _QMORPHMESHER_H

#include "USE_CGAL.h"
#include "toolsdefine.h"
#include "sldelaunay2d.h"

#define EPS_QMorph 1.0e-4

/*class Vector3D;*/
class QuadMeshGenerator;
struct DelaunayAFT2DInterface;
class QuadMesher
{
public:
	QuadMesher(HDS_mesh &backgroundMesh);
	~QuadMesher() {};

	/** \brief 从inp文件中读取三角形网格
	*    根据节点坐标判断网格维数
	   \param file_name  文件名(.inp格式)
	*/
	bool read_file(std::string file_name);

	bool read_tri_mesh(DelaunayAFT2DInterface& tri_mesh);

	bool mesh(const char* method);

	/** \brief 输出网格
	*
		\param file_name 文件名
	*/
	void write_my_mesh(std::string file_name);

	/** \brief nas格式的网格文件
	*
		\param file_name 文件名
	*/
	void write_my_mesh_nas(std::string file_name);

private:

	/** \brief 通过对满足条件的三角形边进行交换来提高参与三角形网格的质量，以利于后续四边形网格的生成
	*/
	bool localReconnect();

	//----------------------------------------------------------------------------//

	void  print_node_long_format_nas(int id,double x,double y,double z,char* strBuffer);

private:
	std::string _file_name; ///< inp文件名

	HDS_mesh &_backgroundMesh; ///< 背景网格

	///< 四边形网格生成算法实现
	QuadMeshGenerator* _meshGenerator;
};

class QMorphInput : public DelaunayAFT2DInterface::Input
{
public:
	QMorphInput(unsigned nnode,const double* nodes,
		unsigned  nedge,const int* edges);

	void      node_begin(){_ni=0;}
	bool      node_next(int& id,double xy[2]);
	void      edge_begin(){_ei=0;}
	bool      edge_next(int nodes[2]);
private:
	unsigned  _nnode,_nedge,_ni,_ei;
	const double* _nodes;
	const int*    _edges;
};
#endif