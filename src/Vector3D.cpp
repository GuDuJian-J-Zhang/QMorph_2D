//**************************************************************************
// File:     Vector3D.cpp
// Author:   Zhang Jun
// Date:     2016-05-08
// Email:    zhangjun_dg@mail.dlut.edu.cn
// Description: 三维向量
//**************************************************************************
#include "USE_CGAL.h"
#include "Vector3D.h"
#include <cmath>

Vector3D::Vector3D(MeshKernel::Point_3& a, MeshKernel::Point_3& b)
	:Vector_3(a, b)
{
	_p_start = a;
	_p_end = b;
}

Vector3D::Vector3D(double x, double y, double z)
	:Vector_3(x, y, z)
{
	_p_start = Point_3(0, 0, 0);
	_p_end = Point_3(x, y, z);
}

Vector3D  operator- (Vector3D& w)
{
	return Vector3D(w.getStart(), -w.getEnd());
}

Vector3D operator+ (Vector3D& v1, Vector3D& v2)
{
	return Vector3D(Point_3(0,0,0), Point_3(v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z()));
}

Vector3D Vector3D::translate(Point_3& p_target)
{
	return Vector3D(p_target, _p_end + p_target + (-_p_start));
}

Vector3D Vector3D::rotateZ(double t)
{
	Point_3 old_start = _p_start; ///< 矢量的原始起点
	Vector3D v_temp = this->translate(Point_3(0, 0, 0)); ///< 先把矢量平移至原点
	MeshKernel::RT x = v_temp.x() * std::cos(t * PI / 180) - v_temp.y() * std::sin(t * PI / 180);
	MeshKernel::RT y = v_temp.x() * std::sin(t * PI / 180) + v_temp.y() * std::cos(t * PI / 180);
	MeshKernel::RT z = v_temp.z();
	if (std::abs(x) < TOLERENCE)     x = 0.0;
	if (std::abs(y) < TOLERENCE)     y = 0.0;
	if (std::abs(this->z()) < TOLERENCE) z = 0.0;
	Point_3 p_e(x, y, z);
	v_temp = Vector3D(v_temp.getStart(), p_e);
	return v_temp.translate(old_start);
}

Vector3D Vector3D::normal()
{
	double length_temp = CGAL::sqrt(this->squared_length());
	return Vector3D(MeshKernel::Point_3(0,0,0), 
		            MeshKernel::Point_3(this->x()/length_temp, this->y()/length_temp, this->z()/length_temp));
}

Vector3D bisectVector(Vector3D& v1, Vector3D& v2, double ang)
{
	Vector3D vl_normal = v1.normal();
	Vector3D vr_normal = v2.normal();

	Vector3D vm = vl_normal + vr_normal;

	if (std::abs(vm.x()) < TOLERENCE && std::abs(vm.y()) < TOLERENCE && std::abs(vm.z()) < TOLERENCE)
	{
		vm = vr_normal;
		double t = vm.x();
		double x = -vm.y();
		double y = t;
		vm = Vector3D(x, y, vm.z());///< 绕z轴负向旋转90度作为角平分线向量
	}	
	else
	{
		if (ang > 180)
		{
			vm = -vm;
		}
	}
	return vm.translate(v1.getStart()); ///< 角平分线所在的单位向量
}

// double Vector3D::dot(Vector3D& v)
// {
// // 	double x1 = _p_end.getX() - _p_start.getX();
// // 	double y1 = _p_end.getY() - _p_start.getY();
// // 	double z1 = _p_end.getZ() - _p_start.getZ();
// // 
// // 	double x2 = v.getEnd().getX() - v.getStart().getX();
// // 	double y2 = v.getEnd().getY() - v.getStart().getY();
// // 	double z2 = v.getEnd().getZ() - v.getStart().getZ();
// // 
// // 	return x1 * x2 + y1 * y2 + z1 * z2;
// 
// 	return _data[0] * v[0] + _data[1] * v[1] + _data[2] * v[2];
// }

/*
Vector3D Vector3D::cross(Vector3D& v)
{
	double x_temp = _data[1] * v[2] - _data[2] * v[1];
	double y_temp = _data[2] * v[0] - _data[0] * v[2];
	double z_temp = _data[0] * v[1] - _data[1] * v[0];
	Point3D p_temp(x_temp, y_temp, z_temp);
	Vector3D resultv(Point3D(0, 0, 0), p_temp);
	return resultv;
}*/