//**************************************************************************
// File:     main.cpp
// Author:   Zhang Jun
// Date:     2016.06.08
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 测试QMorph算法
//**************************************************************************
#include <iostream>
#include "QuadMesher.h"
//#include "USE_OpenGL.h"
#include "sldelaunay2d.h"
#include <boost/filesystem.hpp>
#include <time.h>

//#define Tri_Merge ///< 启用三角形合并算法

#ifndef Tri_Merge // !Tri_Merge

    #define QMorph ///< 启用QMorph算法

#endif // Tri_Merge
int main()
{
	HDS_mesh backoundMesh;

	QuadMesher myMesher(backoundMesh);

#ifdef QMorph
	std::string file_name = "D:\\VS2010_CODE\\CGAL_QMorph_2D\\CGAL_QMorph\\test_0520_b.inp";

	boost::filesystem::path file_path(file_name);

	if (file_path.extension() != ".inp")
	{
		std::cerr << "File format error!\n";
		std::cerr << "Please input an inp file.\n";
		return -1;
	}
	else if (!boost::filesystem::exists(file_name))
	{
		std::cerr << "The file is not exist!\n";
		std::cerr << "Please input an inp file.\n";
		return -1;
	}

// 	if (!myMesher.read_file(file_name)) //读取背景网格(三角形网格)
// 	{
// 		std::cin.get();
// 		return -1;
// 	}
	DelaunayAFT2DInterface mesher2d;

	mesher2d.doMeshInp(file_name.c_str());

	myMesher.read_tri_mesh(mesher2d);

	double dur;
	clock_t start,end;
	start = clock();

	myMesher.mesh("");

#endif
 	
#ifdef Tri_Merge
	#define INP
	std::string file_name = "";
	#ifdef INP
		std::string file_for_tri_mesh = "test_5_2.inp";
	#else
		std::string file_for_tri_mesh = "circle_2.msh";
	#endif

		file_name = file_for_tri_mesh;
	 	DelaunayAFT2DInterface mesher2d;
	#ifdef INP
		mesher2d.doMeshInp(file_for_tri_mesh.c_str());
	#else
		mesher2d.doMesh(file_for_tri_mesh.c_str());
	#endif

	myMesher.read_tri_mesh(mesher2d);

	double dur;
	clock_t start,end;
	start = clock();
	myMesher.mesh();
#endif

	end = clock();
	dur = (double)(end - start);

	std::cout << "Use Time: " << dur/CLOCKS_PER_SEC << "s\n";

 	myMesher.write_my_mesh_nas(file_name + ".nas");

	std::cout << "QMorph Mesher\n";
	std::cin.get();
	return 0;
}