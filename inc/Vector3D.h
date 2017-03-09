//**************************************************************************
// File:     USE_CGAL.h
// Author:   Zhang Jun
// Date:     2016-06-21
// Email:    zhangjun_dg@mail.dlut.edu.cn
//**************************************************************************
#include "USE_CGAL.h"
class Vector3D : public Vector_3
{
	typedef MeshKernel::Point_3 Point_3;
public:
	Vector3D(Point_3& a, Point_3& b);
	Vector3D(double x, double y, double z);

	void setStart(Point_3 ps) { _p_start = ps; }
	Point_3 getStart()        { return _p_start; }

	void setEnd(Point_3 pe)   { _p_end = pe; }
	Point_3 getEnd()          { return _p_end; }

	Vector3D normal(); ///< 返回单位向量

	Vector3D     translate(Point_3& p_target);

	Vector3D     rotateZ(double t);///< 将矢量绕Z轴正向旋转t度

	friend Vector3D     operator- (Vector3D& w);

	friend Vector3D     operator+ (Vector3D& v1, Vector3D& v2);
	/** \brief 计算从同一点出发的两个向量的角平分线向量
	*   \param 返回从两向量公共点出发的单位向量
	*/
	friend Vector3D bisectVector(Vector3D& v1, Vector3D& v2, double ang = 0.0);

	friend Point_3 operator+ (Point_3&p, Point_3& q);

	friend Point_3 operator- (Point_3&p);

	friend Point_3 operator/ (Point_3&p, double t);
private:
	Point_3 _p_start; ///< 起点
	Point_3 _p_end; ///< 终点 
};