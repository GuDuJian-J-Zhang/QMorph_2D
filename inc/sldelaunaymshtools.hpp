//sldelaunaymshtools.hpp
#ifndef __SLDELAUNAYMSHTOOLS_H_
#define __SLDELAUNAYMSHTOOLS_H_
#include "sldelaunay2d_std.h"
SL_MSH_NS_BEGIN
class  RandomGen
{
public:
    Int64 u,v,w;
    RandomGen(Int64 j=23);
    Int64 int64();
    UInt randInt(UInt num);
};
inline Real ToRadian(Real degree){return (Real)(0.01745329252*degree);}
Real Incircle_Math(const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d);
Real Orient2d_Math(const Point2D &a, const Point2D &b, const Point2D &c);
Real Insphere_Math(const Point3D &pa, const Point3D &pb, const Point3D &pc, const Point3D &pd, const Point3D &pe);
Real Orient3d_Math(const Point3D &pa, const Point3D &pb, const Point3D &pc, const Point3D &pd);
template <typename T>
inline void Swap(T& a, T& b) 
{
    T tmp = a;
    a = b;
    b = tmp;
}

class UIntBits
{
public:
    UIntBits() : _flag(0) {}

    void set(int i) 
    { 
        unsigned int bi=1<<i;
        _flag|=bi;
    }
    void reset() {_flag = 0;}
    void clear(int i)
    {
        unsigned int bi=~(1<<i);
        _flag&=bi;
    }

    bool test(int i) const 
    { 
        unsigned int bi=1<<i;
        bi&=_flag;
        return bi>0; 
    }
    unsigned int get() const { return _flag; }
private:
    unsigned int  _flag;
};


SL_MSH_NS_END
#endif
