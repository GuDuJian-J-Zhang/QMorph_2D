//sldelaunay2d_std.h
#ifndef __SLDELAUNAY2D_STD_H_
#define __SLDELAUNAY2D_STD_H_
#   pragma warning (disable : 4267)
#include <time.h>
#include <vector>
#include <fstream>
#include <cassert>
#include <math.h>
#include <stack>
#include <iomanip>
#include <float.h>
#include "predicates.h"
#include <map>
#include <bitset>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

# define SL_MSH_NS DELAUNAY_MSHGEN
# define SL_MSH_NS_BEGIN namespace SL_MSH_NS {
# define SL_MSH_NS_END }

SL_MSH_NS_BEGIN

#define Array1D     std::vector
#define Stack       std::stack
#define StdPair     std::pair
#define MakePair    std::make_pair
#define Assert      assert
#define StdBitSet    std::bitset
#define SL_NEW      new
typedef  double              Real;
typedef  long  long          Int64;
typedef  size_t              UInt;
typedef  unsigned long       ULong;
typedef  unsigned long long  HashValue;
typedef std::map<int,int>    IdMap;
#define InValidEntity        UInt(-1)
#define PI                   3.141592653589793238462643

/*********************************************
三角形单元节点/边编号示意图
                 n2
			  |      |
			e1         e0
		 |               |
		n0 - - - e2 - - - n1
*********************************************/
static int TRI_EDGE_VERTEX[3][2]      =   {{1,2},{2,0},{0,1}}; ///< 三角形单元各边的顶点编号
static  int TRI_NEXT[3]         =  {1,2,0};
static  int TRI_PRE[3]          =  {2,0,1};
static  int QUAD_NEXT[4]        =  {1,2,3,0};
static  int QUAD_PRE[4]         =  {3,0,1,2};

static int TRI_VERTEX_EDGE[3][2]       =   {{TRI_PRE[0],TRI_NEXT[0]},{TRI_PRE[1],TRI_NEXT[1]},{TRI_PRE[2],TRI_NEXT[2]}}; ///< 三角形各边在本单元中的前驱边及后继边的编号(三角形各顶点在本单元中的邻接边)
struct Tolerance
{
    static Real EPS_COMPUTER;
    static Real EPS_RELATIVE;
    static Real EPS;
    static Real MAX_REAL;
};

template<class Key>
class HashCompare
{
public:
    size_t operator( )( const Key& key ) const
    {
        return (size_t)(key.hash());
    }
    bool operator( )( const Key& _Key1,const Key& _Key2) const
    {
        return _Key1==_Key2;
    }
};

template<class T>
class HashSetDefault : public boost::unordered_set<T,HashCompare<T>,HashCompare<T> >
{
};

typedef  boost::unordered_set<unsigned int>    HashUIntSet;

template <class T, int N>
class Tuple
{
public:
    Tuple ();
    Tuple (const T* v);
    Tuple (const Tuple& rhs);
    Tuple& operator= (const Tuple& rhs);
    T& operator[] (int i);
    T  operator[] (int i) const;
    const T*        memory()const;
    T*              memory();
    bool            operator == (const Tuple& vec)const;
    bool            operator != (const Tuple& vec)const;

    T data[N];
};

struct  CoupleReal : public Tuple<Real,2>
{
    CoupleReal(Real x=0.0,Real y=0.0){data[0]=x;data[1]=y;}
    CoupleReal(const CoupleReal&,const CoupleReal&);
    const CoupleReal& setToMin (const CoupleReal & p2);
    const CoupleReal& setToMax (const CoupleReal & p2);
    void             set(const CoupleReal&,const CoupleReal&);
    Real             operator [] (int i) const{return data[i];}
    Real&            operator [] (int i){return data[i];}
    CoupleReal       operator * (Real val) const;
    CoupleReal       operator - ()const;
    CoupleReal&      operator *= (Real val);
    CoupleReal       operator + (const CoupleReal & p2) const;
    CoupleReal&      operator += (const CoupleReal & p2);
    CoupleReal&      operator -= (const CoupleReal & p2);

    Real             distance2(const CoupleReal & p2)const;
    Real             length()const{return sqrt(data[0]*data[0]+data[1]*data[1]);}

    Real             distanceTo(const CoupleReal & p2)const{return sqrt(distance2(p2));}
    Real             dot(const CoupleReal& vec) const;
    CoupleReal&      normalize();
    CoupleReal&      set(Real x,Real y){data[0]=x;data[1]=y;return *this;}
    const Real*      memory()const{return &(data[0]);}
    CoupleReal       perpVector()const{return CoupleReal(-data[1],data[0]);}
	CoupleReal&      rotateZ(double ang);
};

struct  TripleReal : public Tuple<Real,3>
{
    TripleReal(Real x=0.0,Real y=0.0,Real z=0.0){data[0]=x;data[1]=y;data[2]=z;}
    TripleReal(const TripleReal&,const TripleReal&);
    TripleReal&      set(const TripleReal&,const TripleReal&);
    Real  distance2(const TripleReal&)const;
    Real  distanceTo(const TripleReal& src)const{return sqrt(distance2(src));}
    TripleReal&      set(Real x,Real y,Real z){data[0]=x;data[1]=y;data[2]=z;return *this;}
    Real  len2()const{return dot(*this);}
    Real  length()const{return sqrt(len2());}
    TripleReal& normalize();
    Real  dot(const TripleReal&)const;
    TripleReal cross(const TripleReal& vec)const;
    TripleReal       operator + (const TripleReal & p2) const;
    TripleReal&      operator += (const TripleReal & p2);
    TripleReal&      operator -= (const TripleReal & p2);
    TripleReal       operator * (Real val) const;
    TripleReal&      operator *= (Real val);
    TripleReal       operator - ()const;
};



inline CoupleReal   CoupleReal::operator * (Real val) const
{
    return CoupleReal(data[0]*val, data[1]*val);
}

inline CoupleReal&  CoupleReal::operator *= (Real val)
{
    data[0]*=val;data[1]*=val;
    return *this;
}

inline CoupleReal   CoupleReal::operator + (const CoupleReal & p2) const
{
    return CoupleReal(data[0]+p2.data[0],data[1]+p2.data[1]);
}

inline CoupleReal&  CoupleReal::operator += (const CoupleReal & p2)
{
    data[0]+=p2.data[0];data[1]+=p2.data[1];
    return *this;
}

inline CoupleReal&  CoupleReal::operator -= (const CoupleReal & p2)
{
    data[0]-=p2.data[0];data[1]-=p2.data[1];
    return *this;
}

inline TripleReal&      TripleReal::operator += (const TripleReal & p2)
{
    data[0]+=p2.data[0];data[1]+=p2.data[1]; data[2]+=p2.data[2];
    return *this;
}

inline TripleReal&      TripleReal::operator -= (const TripleReal & p2)
{
    data[0]-=p2.data[0];data[1]-=p2.data[1]; data[2]-=p2.data[2];
    return *this;
}


inline const CoupleReal & CoupleReal::setToMin (const CoupleReal & p2)
{
    if (p2[0] < data[0]) data[0] = p2[0];
    if (p2[1] < data[1]) data[1] = p2[1];
    return *this;
}

inline Real Area2(const CoupleReal& point1, const CoupleReal& point2,const CoupleReal& point3)
{
    return (point1[0]*(point2[1]-point3[1])+
        point2[0]*(point3[1]-point1[1])+
        point3[0]*(point1[1]-point2[1]));
}

inline  const CoupleReal & CoupleReal::setToMax (const CoupleReal & p2)
{
    if (p2[0] > data[0]) data[0] = p2[0];
    if (p2[1] > data[1]) data[1] = p2[1];
    return *this;
}

inline void CoupleReal::set(const CoupleReal& a,const CoupleReal& b)
{
    data[0]=b.data[0]-a.data[0];
    data[1]=b.data[1]-a.data[1];
}

inline CoupleReal CoupleReal::operator - () const
{
    return CoupleReal (-data[0], -data[1]);
}


inline CoupleReal::CoupleReal(const CoupleReal& a,const CoupleReal& b)
{
    data[0]=b.data[0]-a.data[0];
    data[1]=b.data[1]-a.data[1];
}

inline TripleReal::TripleReal(const TripleReal& a,const TripleReal& b)
{
    data[0]=b.data[0]-a.data[0];
    data[1]=b.data[1]-a.data[1];
    data[2]=b.data[2]-a.data[2];
}

inline TripleReal& TripleReal::set(const TripleReal& a,const TripleReal& b)
{
    data[0]=b.data[0]-a.data[0];
    data[1]=b.data[1]-a.data[1];
    data[2]=b.data[2]-a.data[2];
}

inline Real CoupleReal::distance2(const CoupleReal & p2)const
{
    Real dx=data[0]-p2.data[0];
    Real dy=data[1]-p2.data[1];
    return dx*dx+dy*dy;
}

inline Real  TripleReal::distance2(const TripleReal& p2)const
{
    Real dx=data[0]-p2.data[0];
    Real dy=data[1]-p2.data[1];
    Real dz=data[2]-p2.data[2];
    return dx*dx+dy*dy+dz*dz;
}

inline TripleReal      TripleReal::operator * (Real val) const
{
    return TripleReal(data[0]*val, data[1]*val,data[2]*val);
}

inline TripleReal&      TripleReal::operator *= (Real val)
{
    data[0]*=val;
    data[1]*=val;
    data[2]*=val;
    return *this;
}

inline TripleReal TripleReal::operator - ()const
{
    return TripleReal (-data[0], -data[1],-data[2]);
}

inline TripleReal      TripleReal:: operator + (const TripleReal & p2) const
{
    return TripleReal(data[0]+p2.data[0],data[1]+p2.data[1],data[2]+p2.data[2]);
}

inline CoupleReal& CoupleReal::normalize()
{
    Real len=length();
    if(len<Tolerance::EPS_COMPUTER)
    {
        data[0]=data[1]=0.0;
    }else
    {
        len=1.0/len;
        data[0]*=len;data[1]*=len;
    }
    return *this;
}

/** \brief 将矢量绕Z轴逆时针旋转角度ang
*/
inline CoupleReal& CoupleReal::rotateZ(double ang)
{
	double x_new = data[0] * std::cos(ang * PI / 180) - data[1] * std::sin(ang * PI / 180);
	double y_new = data[0] * std::sin(ang * PI / 180) + data[1] * std::cos(ang * PI / 180);
	data[0] = x_new;
	data[1] = y_new;
	return *this;
}

inline Real CoupleReal::dot(const CoupleReal& vec) const
{
    return data[0] * vec.data[0] + data[1] * vec.data[1] ;
}

inline Real  TripleReal::dot(const TripleReal& vec)const
{
    return data[0] * vec.data[0] + data[1] * vec.data[1] + data[2] * vec.data[2] ;
}

inline  TripleReal TripleReal::cross(const TripleReal& vec)const
{
    return TripleReal(data[1]*vec.data[2]-data[2]*vec.data[1],
        data[2]*vec.data[0]-data[0]*vec.data[2],
        data[0]*vec.data[1]-data[1]*vec.data[0]);
}

inline TripleReal& TripleReal::normalize()
{
    Real len=length();
    if(len<Tolerance::EPS_COMPUTER)
    {
        data[0]=data[1]=data[2]=0.0;
    }else
    {
        len=1.0/len;
        data[0]*=len;data[1]*=len;data[2]*=len;
    }
    return *this;
}

inline TripleReal NormalTri(const TripleReal& pt1,
                            const TripleReal& pt2,
                            const TripleReal& pt3)
{
    TripleReal vector1(pt1,pt2);
    TripleReal vector2(pt1,pt3);
    return vector1.cross(vector2) ;
}

inline CoupleReal  operator- (const CoupleReal& vec1,const CoupleReal& vec2)
{
    return CoupleReal(vec1[0]-vec2[0],vec1[1]-vec2[1]);
}

template <class T, int N>
Tuple<T,N>::Tuple()
{
    for(int i=0;i<N;i++)
    {
        data[i]=T();
    }
}

template <class T, int N>
T& Tuple<T,N>::operator[] (int i)
{
    return data[i];
}

template <class T, int N>
T Tuple<T,N>::operator[] (int i)const
{
    return data[i];
}


template <class T, int N>
Tuple<T,N>::Tuple (const T* v)
{
    for (int i = 0; i < N; ++i)
        data[i] = v[i];
}

template <class T, int N>
Tuple<T,N>::Tuple (const Tuple<T,N>& rhs)
{
    for (int i = 0; i < N; ++i)
        data[i] = rhs.data[i];
}

template <class T, int N>
Tuple<T,N>&
Tuple<T,N>::operator= (const Tuple<T,N>& rhs)
{
    if(this!=&rhs)
    {
        for (int i = 0; i < N; ++i)
            data[i] = rhs.data[i];
    }
    return *this;
}

template <class T, int N>
const T*  Tuple<T,N>::memory()const{return data;}

template <class T, int N>
T*  Tuple<T,N>::memory(){return data;}

template <class T, int N>
bool Tuple<T,N>::operator == (const Tuple<T,N>& vec)const
{
    return (memcmp( data, vec.data, N*sizeof(T) )==0);
}

template <class T, int N>
bool  Tuple<T,N>::operator != (const Tuple<T,N>& vec)const
{
    return (memcmp( data, vec.data, N*sizeof(T) )!=0);
}

template<class T>
class  KeyIndex2
{
public:
    KeyIndex2(const T& k0=T(),const T& k1=T())
    {
        reset(k0,k1);
    }

    KeyIndex2(const T ks[2])
    {
        reset(ks[0],ks[1]);
    }

    void reset(const T& k0,const T& k1)
    {
        if(k0<k1)
        {
            _key[0]=k0;
            _key[1]=k1;
        }else
        {
            _key[0]=k1;
            _key[1]=k0;
        }
    }


    bool operator<(const KeyIndex2& src)const
    {
        if(_key[0]<src._key[0]) return true;
        if(src._key[0]<_key[0]) return false;
        return _key[1]<src._key[1];
    }


    bool  operator==(const KeyIndex2& src)const
    {
        if(_key[0]==src._key[0])
        {
            return _key[1]==src._key[1];
        }
        return false;
    }


    bool operator!=(const KeyIndex2& src)const
    {
        return !(*this==src);
    }

    T    operator[](int index)const{return _key[index];}
private:
    T   _key[2];
};

template<class T>
class  KeyIndex3
{
public:
    KeyIndex3(const T& k0=T(),const T& k1=T(),const T& k2=T())
    {
        set(k0,k1,k2);
    }

    void  set(const T& k0,const T& k1,const T& k2)
    {
        if(k0<k1)
        {
            if(k1<k2)
            {
                _key[0]=k0;
                _key[1]=k1;
                _key[2]=k2;
            }else
            {
                _key[2]=k1;
                if(k0<k2)
                {
                    _key[0]=k0;
                    _key[1]=k2;
                }else
                {
                    _key[0]=k2;
                    _key[1]=k0;
                }
            }
        }else
        {
            if(k2<k1)
            {
                _key[0]=k2;
                _key[1]=k1;
                _key[2]=k0;
            }else
            {
                _key[0]=k1;
                if(k2<k0)
                {
                    _key[1]=k2;
                    _key[2]=k0;
                }else
                {
                    _key[1]=k0;
                    _key[2]=k2;
                }
            }

        }
    }

    T  operator[](int index)const{return _key[index];}

    bool operator<(const KeyIndex3& src)const
    {
        if(_key[0]<src._key[0]) return true;
        if(src._key[0]<_key[0]) return false;
        if(_key[1]<src._key[1]) return true;
        if(src._key[1]<_key[1]) return false;
        return _key[2]<src._key[2];
    }

    bool  operator==(const KeyIndex3& src)const
    {
        if(_key[0]==src._key[0])
        {
            if(_key[1]==src._key[1])
            {
                return  _key[2]==src._key[2];
            }
        }
        return false;
    }

    bool operator!=(const KeyIndex3& src)const
    {
        return !(*this==src);
    }
private:
    T   _key[3];
};

template<class T>
class  KeyIndex4
{
public:
    KeyIndex4(const T& k0=T(),const T& k1=T(),const T& k2=T(),const T& k3=T())
    {

        set(k0,k1,k2,k3);
    }

    void  set(const T& k0,const T& k1,const T& k2,const T& k3)
    {
        _key[0]=k0;
        _key[1]=k1;
        _key[2]=k2;
        _key[3]=k3;
        sort();
    }

    T  operator[](int index)const{return _key[index];}

    bool operator<(const KeyIndex4& src)const
    {
        if(_key[0]<src._key[0]) return true;
        if(src._key[0]<_key[0]) return false;
        if(_key[1]<src._key[1]) return true;
        if(src._key[1]<_key[1]) return false;
        if(_key[2]<src._key[2]) return true;
        if(src._key[2]<_key[2]) return false;
        return _key[3]<src._key[3];
    }

    bool  operator==(const KeyIndex4& src)const
    {
        if(_key[0]==src._key[0])
        {
            if(_key[1]==src._key[1])
            {
                if(_key[2]==src._key[2])
                {
                    return  _key[3]==src._key[3];
                }
            }
        }
        return false;
    }

    bool operator!=(const KeyIndex4& src)const
    {
        return !(*this==src);
    }

private:
    void  sort()
    {
        if (_key[1] < _key[0])std::swap (_key[0], _key[1]);
        if (_key[2] < _key[1])
        {
            std::swap (_key[1], _key[2]);
            if (_key[1]<_key[0]) std::swap (_key[0], _key[1]);
        }
        /* step2: isort4p */
        if (_key[3] < _key[2])
        {
            if (_key[3] < _key[0])
            {
                std::swap (_key[2], _key[3]);
                std::swap (_key[1], _key[2]);
                std::swap (_key[0], _key[1]);
            }
            else if (_key[3] < _key[1])
            {
                std::swap (_key[2], _key[3]);
                std::swap (_key[1], _key[2]);
            }
            else std::swap (_key[2], _key[3]);
        }
    }
private:
    T   _key[4];
};



template <class T>
struct HashCouple  : public  Tuple<T,2>
{
    HashCouple(const T& k0=T(),const T& k1=T())
    {
        this->data[0]=k0;
        this->data[1]=k1;
        _hash=0;
    }

    bool  operator==(const HashCouple<T>& src)const
    {
        if(this->hash()!=src.hash()) return false;
        KeyIndex2<T> s1(this->data[0],this->data[1]);
        KeyIndex2<T> s2(src.data[0],src.data[1]);
        return s1==s2;
    }

    HashValue  hash()const
    {
        if(_hash==0)
        {
            KeyIndex2<T> s1(this->data[0],this->data[1]);
            _hash= Hash(s1[0],s1[1]);
        }
        return _hash;
    }

private:
    mutable  HashValue _hash;
};

template <class T>
class TripleT  : public Tuple<T,3>
{
public:
    TripleT():Tuple<T,3>(){}
    TripleT(const T*);
    TripleT(const T& a, const T& b, const T& c);
    TripleT<T>&   set(const T& a, const T& b, const T& c);
	double get_min_angle() { return min_angle; }
	void set_min_angle(double t) { min_angle = t; }
private:
	double min_angle;
	double max_angle;
};

template <class T>
TripleT<T>::TripleT(const T& a, const T& b, const T& c)
{
    set(a,b,c);
}

template <class T>
TripleT<T>::TripleT(const T* a3)
{
    set(a3[0],a3[1],a3[2]);
}

template <class T>
TripleT<T>&  TripleT<T>::set (const T& a , const T& b, const T& c)
{
    this->data[0]=a;
    this->data[1]=b;
    this->data[2]=c;
    return *this;
}

template <class T>
class QuaternionT  : public Tuple<T,4>
{
public:
    QuaternionT(const T& a, const T& b, const T& c,const T& d);
    QuaternionT<T>&   set(const T& a, const T& b, const T& c,const T& d);
    QuaternionT():Tuple<T,4>(){}
};

template <class T>
QuaternionT<T>::QuaternionT(const T& a, const T& b, const T& c,const T& d)
{
    set(a,b,c,d);
}

template <class T>
QuaternionT<T>&  QuaternionT<T>::set (const T& a , const T& b, const T& c,const T& d)
{
    this->data[0]=a;
    this->data[1]=b;
    this->data[2]=c;
    this->data[3]=d;
    return *this;
}

inline HashValue Hash(HashValue h1,HashValue h2)   {return HashValue(((h1 << 16) | (h1 >> 16)) ^ h2);}

typedef CoupleReal Vector2D;
typedef CoupleReal Point2D;
typedef TripleReal Vector3D;
typedef TripleReal Point3D;

/////////////////////////////////////////////////////////////////////

#define  USE_POSITIONFORSEARCH 0

namespace SlObjectPoolDetails
{
#if  USE_POSITIONFORSEARCH
    template <class T>
    struct PageData_ : public T
    {
        PageData_(const T& data=T(),UInt p=UInt(-1))
            :T(data),pos(p)
        {
        }

        PageData_<T>&  operator=(const T& data)
        {
            T::operator=(data);
            return *this;
        }
        UInt   pos;
    };
#endif

    template <class T,int PAGE_SIZE>
    struct PageEx_
    {
#if    USE_POSITIONFORSEARCH
        typedef PageData_<T>   PageData;
#else
        typedef T              PageData;
#endif

        PageEx_(UInt pn,const T& val)
        {
            for(UInt i=0;i<PAGE_SIZE;i++)
            {
#if    USE_POSITIONFORSEARCH
                data[i]=PageData(val,pn*PAGE_SIZE+i);
#else
                data[i]=val;
#endif
                flag[i]=true;
            }
        }

        PageEx_(const PageEx_<T,PAGE_SIZE>& src)
        {
            operator=(src);
        }

        PageEx_<T,PAGE_SIZE>& operator=(const PageEx_<T,PAGE_SIZE>& src)
        {
            if(this!=&src)
            {
                for(UInt i=0;i<PAGE_SIZE;i++)
                {
                    data[i]=src.data[i];
                    flag[i]=src.flag[i];
                }
            }
            return *this;
        }

        UInt erase(UInt pos)
        {
            if (pos>=PAGE_SIZE) return  UInt(-1);
            if(!flag[pos])
            {
                flag[pos]=true;
                return pos;
            }
            return  UInt(-1);
        }

        bool  isValid(int i)const{return !flag[i];}

        PageData             data[PAGE_SIZE];
        StdBitSet<PAGE_SIZE>   flag;// erased flag;
    };
};

//=========================================SlObjectPool===========================
template <class T,int PAGE_SIZE=1024>
class  SlObjectPool
{
public:
    typedef SlObjectPoolDetails::PageEx_<T,PAGE_SIZE>  Page;
    typedef T                           type;
    typedef T*                          pointer;
    typedef const T*                    const_pointer;
    typedef typename Page::PageData     PageData;
    typedef PageData*                   PageDataPtr;
    typedef const PageData*             ConstPageDataPtr;
    typedef Page*                       PagePtr;

    SlObjectPool()
        :size_(0)
#if  !USE_POSITIONFORSEARCH
        ,sortedpagesPtr_(0)
#endif
    {
#if  !USE_POSITIONFORSEARCH
        flag_[INIT_PAGES_SORTED]=true;
        flag_[SORTEDPAGEOK]=false;
#endif
    }

    ~SlObjectPool()
    {
        clear();
    }

#if    USE_POSITIONFORSEARCH
    bool  erase(const T* data0)
    {
        PageData* data=(PageData*)(data0);
        return erase(data);
    }
#else
    bool compare(const PageData* dat0,const PageData* dat1)const
    {
        return dat0<dat1;
    }
    bool compare(const PageData* dat0,const PagePtr dat1)const
    {
        return dat0<dat1->data;
    }

    bool compare(const PagePtr dat1,const PageData* dat0)const
    {
        return ((dat0-dat1->data)>=PAGE_SIZE);
    }


#endif

    T*  add(const T& obj)
    {
        UInt pi,off;

        if(!erased_.empty())
        {
            UInt pos=erased_[erased_.size()-1];
            erased_.pop_back();
            pi= getPage(pos,off);
        }else
        {
            pi= getPage(size_,off);

            if(pi==UInt(-1))
            {
                PagePtr p = SL_NEW Page(pages_.size(),T());
                pi=pages_.size();
                off=0;
                pages_.push_back(p);

#if  !USE_POSITIONFORSEARCH
                if(flag_[INIT_PAGES_SORTED] && pi>1 &&
                    compare(pages_[pi]->data,pages_[pi-1]->data))
                {
                    flag_[INIT_PAGES_SORTED]=false;
                }
                flag_[SORTEDPAGEOK]=false;
#endif
            }
        }
        pages_[pi]->data[off]=obj;
        pages_[pi]->flag[off]=false;
        size_++;
        return &(pages_[pi]->data[off]);
    }

    bool  isContains(const PageData* data)const
    {
        if(data==0)
        {
            return false;
        }
        UInt pidx,pn;
#if    USE_POSITIONFORSEARCH
        pn=getPage(data->pos,pidx);
#else
        pn=getPage(data,pidx);
#endif
        if(pn==UInt(-1))
        {
            return false;
        }
        return pages_[pn]->isValid(pidx);
    }

    bool  erase(const PageData* data)
    {
        if(data==0)
        {
            return false;
        }
        UInt pidx,pn,pos;
#if    USE_POSITIONFORSEARCH
        pos=data->pos;
        pn=getPage(pos,pidx);
#else
        pn=getPage(data,pidx);
        pos=pn*PAGE_SIZE+pidx;
#endif
        if(pn==UInt(-1))
        {
            return false;
        }
        UInt rt=pages_[pn]->erase(pidx);
        if(rt==UInt(-1))
        {
            return false;
        }
        size_-=1;
        erased_.push_back(pos);
        return true;
    }

    void clear()
    {
        if(size()==0) return;
        UInt i,sz=pages_.size();
        for(i=0;i<sz;i++)
        {
            delete pages_[i];
        }
        pages_.clear();
        erased_.clear();
#if  !USE_POSITIONFORSEARCH
        flag_[INIT_PAGES_SORTED]=true;
        flag_[SORTEDPAGEOK]=false;
        delete[] sortedpagesPtr_;
        sortedpagesPtr_=0;
#endif
    }

    UInt size()const{return size_;}
    bool   empty()const{return size_==0;}
    void   begin()const;
    bool   end()const{return iter_==UInt(-1);}

    void   next()const;
    const_pointer current()const;
    pointer current();

    const_pointer nextValidLoop(UInt pos)const;
private:
    UInt  getPage(UInt pos,UInt& pidx)const
    {
        UInt pn=pos/PAGE_SIZE;
        if(pn>=pages_.size()) return UInt(-1);
        pidx=pos%PAGE_SIZE;
        return pn;
    }

#if  !USE_POSITIONFORSEARCH
    struct Sorter
    {
        bool operator()(UInt i,UInt j)const
        {
            return pages_->at(i)<pages_->at(j);
        }

        Sorter(const Array1D<PagePtr>* pages)
            :pages_(pages)
        {
        }
        const Array1D<PagePtr>*        pages_;
    };

    UInt  getPage(const PageData* data,UInt& pidx)const
    {
        UInt ps=pages_.size();
        if(ps==0) return UInt(-1);
        if(flag_[INIT_PAGES_SORTED])
        {
            int begin=0,mid=0,end=pages_.size();
            while(begin<=end)
            {
                mid=(end+begin)>>1;
                if(ps<=mid) return UInt(-1);
                if(compare(data,pages_[mid]))end=mid-1;
                else if(compare(pages_[mid],data)) begin=mid+1;
                else
                {
                   // _SLAssert((data >= pages_[mid]->data));
                    pidx=data-pages_[mid]->data;
                    //_SLAssert(pidx<PAGE_SIZE);
                    return mid;
                }
            }
        }else
        {
            if(!flag_[SORTEDPAGEOK])
            {
                delete[] sortedpagesPtr_;
                UInt sz=pages_.size();
                sortedpagesPtr_=SL_NEW UInt[sz];
                for(UInt i=0;i<sz;i++)
                {
                    sortedpagesPtr_[i]=i;
                }
                Sorter sorter(&pages_);
                std::sort(sortedpagesPtr_,sortedpagesPtr_+sz,sorter);
                flag_[SORTEDPAGEOK]=true;
            }

            int pos,begin=0,mid=0,end=pages_.size();
            while(begin<=end)
            {
                mid=(end+begin)>>1;
                pos=sortedpagesPtr_[mid];
                if(ps<=pos) return UInt(-1);
                if(compare(data,pages_[pos]))end=mid-1;
                else if(compare(pages_[pos],data)) begin=mid+1;
                else
                {
                    //_SLAssert(data>=pages_[pos]->data);
                    pidx=data-pages_[pos]->data;
                    //_SLAssert(pidx<PAGE_SIZE);
                    return pos;
                }
            }
        }
        return UInt(-1);
    }
#endif
    void   next(UInt& pos)const
    {
        pos+=1;
        UInt pn=pos/PAGE_SIZE;
        UInt off=pos%PAGE_SIZE;
        UInt szp=pages_.size();
        while(1)
        {
            if(pn>=szp)
            {
                pos=UInt(-1);
                return;
            }
            for(UInt i=off;i<PAGE_SIZE;i++)
            {
                if(!(pages_[pn]->flag[i]))
                {
                    pos=pn*PAGE_SIZE+i;
                    return;
                }
            }
            pn++;
            off=0;
        }
    }

private:
    SlObjectPool(const SlObjectPool<T,PAGE_SIZE>& src){};
    SlObjectPool<T,PAGE_SIZE>& operator=(const SlObjectPool<T,PAGE_SIZE>& src);
private:
    Array1D<PagePtr>        pages_;
    Array1D<UInt>         erased_;
    UInt                    size_;
#if  !USE_POSITIONFORSEARCH
    enum {INIT_PAGES_SORTED=0,SORTEDPAGEOK,UNDEFINED_STATE};
    mutable UInt  *        sortedpagesPtr_;
    mutable StdBitSet<UNDEFINED_STATE> flag_;
#endif

    mutable UInt       iter_;

};


template <class T,int PAGE_SIZE>
void   SlObjectPool<T,PAGE_SIZE>::begin()const
{
    iter_=0;
    if(pages_[0]->flag[0])
    {
        next(iter_);
    }
}

template <class T,int PAGE_SIZE>
void   SlObjectPool<T,PAGE_SIZE>::next()const
{
    next(iter_);
}

template <class T,int PAGE_SIZE>
typename SlObjectPool<T,PAGE_SIZE>::pointer
SlObjectPool<T,PAGE_SIZE>::current()
{
    UInt off=0;
    UInt sz= getPage(iter_,off);
    //_SLAssert(sz!=UInt(-1) && off<PAGE_SIZE);
    return &(pages_[sz]->data[off]);
}

template <class T,int PAGE_SIZE>
typename  SlObjectPool<T,PAGE_SIZE>::const_pointer
SlObjectPool<T,PAGE_SIZE>::current()const
{
    UInt off=0;
    UInt sz= getPage(iter_,off);
    //_SLAssert(sz!=UInt(-1) && off<PAGE_SIZE);
    return &(pages_[sz]->data[off]);
}

template <class T,int PAGE_SIZE>
typename  SlObjectPool<T,PAGE_SIZE>::const_pointer
SlObjectPool<T,PAGE_SIZE>::nextValidLoop(UInt pos)const
{
    int times=0;
    while(1)
    {
        next(pos);
        while(pos==UInt(-1))
        {
            if(times>0) return 0;
            times++;
            pos=0;
            next(pos);
        }
        UInt off=0;
        UInt sz= getPage(pos,off);
        //_SLAssert(sz!=UInt(-1) && off<PAGE_SIZE);
        return &(pages_[sz]->data[off]);
    }

}

//###################################################################

template <class T>
class GeneralCriterion
{
public:
    virtual  bool operator()(const T& data)=0;
};

template <typename Real>
bool  IsInD1(int Dim,const Real cmin[],const Real cmax[],const Real pt[])
{
    for(int i=0;i<Dim;i++)
    {
        if(pt[i]>cmax[i]) return false;
        if(pt[i]<cmin[i]) return false;
    }
    return true;
}

template<class Real,class Data>
struct AdTreeExNode;
#define   EPS_RALATIVE  1.e-5

template <int Dim,class Real,class Data>
class AdTreeEx
{
public:
    typedef AdTreeExNode<Real,Data> Node;
    AdTreeEx(Real cmin[],Real cmax[]);
    void   insert (const Data& dt);
    Data   uniqueInsert (const Data& dt,bool& isSuccess,Real tol=Tolerance::EPS);
    void   getIntersecting (Real bmin[],Real bmax[],Array1D<Data>& pis, GeneralCriterion<Data>* pFun=0);
    bool   search(Real pt[],Data&,const GeneralCriterion<Data>&,Real eps=Tolerance::EPS);
private:
    void   reset ();
    enum{STACKSIZE=1000};
    Node*  _root;
    Real  _cmin[Dim], _cmax[Dim],_eps;
    Array1D<Node*> _stack;
    Array1D<int> _stackdir;
    int            _stackindex;
    SlObjectPool<Node > _dataPool;
};

template<class Real,class Data>
struct AdTreeExNode
{
    AdTreeExNode()
        :left(0),
        right(0),
        fulled(false)
    {}



    AdTreeExNode<Real,Data>  *left,*right;//*father;
    Real sep;
    Data data;
    bool fulled;
};
///////////////////////////////////////////////////////////////////////

template <int Dim,class Real,class Data>
AdTreeEx<Dim,Real,Data>::AdTreeEx(Real acmin[],Real acmax[])
:_stack(STACKSIZE),
_stackdir(STACKSIZE)
{
    memcpy (_cmin, acmin, Dim * sizeof(Real));
    memcpy (_cmax, acmax, Dim * sizeof(Real));
    //_root = _dataPool.add(AdTreeExNode<Real,Data>);
    _root->sep = (_cmin[0] + _cmax[0])*0.5;
    _eps=_root->sep*EPS_RALATIVE;
}

template <int Dim,class Real,class Data>
void  AdTreeEx<Dim,Real,Data>::insert (const Data& pi)
{
    int lr;
    Real bmin[Dim];
    Real bmax[Dim];
    memcpy (bmin, _cmin, Dim * sizeof(Real));
    memcpy (bmax, _cmax, Dim * sizeof(Real));

    Node *node=0,*next = _root;
    int dir = 0;
    while (next)
    {
        node = next;
        if (!node->fulled)
        {
            node->data= pi;
            node->fulled=true;
            return;
        }

        if (node->sep > pi[dir])
        {
            next = node->left;
            bmax[dir] = node->sep;
            lr = 0;
        } else
        {
            next = node->right;
            bmin[dir] = node->sep;
            lr = 1;
        }

        dir++;
        if (dir == Dim)
            dir = 0;
    }


    //next = _dataPool.add(AdTreeExNode<Real,Data>);
    next->data= pi;
    next->fulled=true;
    next->sep = (bmin[dir] + bmax[dir])*0.5;
    Real tmp=next->sep*EPS_RALATIVE;
    if(tmp<_eps) _eps=tmp;
    if (lr)
        node->right = next;
    else
        node->left = next;
    // next -> father = node;
}

template <int Dim,class Real,class Data>
Data   AdTreeEx<Dim,Real,Data>::uniqueInsert (const Data& dt,bool& isSuccess,Real tol)
{
    Real bmin[Dim], bmax[Dim];
    for(int i=0;i<Dim;i++)
    {
        bmin[i]=dt[i]-tol;
        bmax[i]=dt[i]+tol;
    }

    Array1D<Data> pis;
    isSuccess=false;
    getIntersecting (bmin,bmax,pis);
    if(!pis.empty()) return pis[0];
    isSuccess=true;
    insert(dt);
    return dt;
}

template <int Dim,class Real,class Data>
void   AdTreeEx<Dim,Real,Data>:: reset ()
{
    _stack[1] = _root;
    _stackdir[1] = 0;
    _stackindex = 1;
}


template <int Dim,class Real,class Data>
void AdTreeEx<Dim,Real,Data>::getIntersecting (Real bmin[],
                                               Real bmax[],
                                               Array1D<Data>& pis,
                                               GeneralCriterion<Data>* pFun)
{

    Node * node;
    int dir;
    reset();

    while (_stackindex)
    {
        node = _stack[_stackindex];
        dir =  _stackdir[_stackindex];
        _stackindex--;

        if (node->fulled)
        {
            if (IsInD2(Dim,bmin,bmax,node->data))
            {
                if(pFun)
                {
                    if((*pFun)(node->data))
                    {
                        pis.push_back (node->data);
                    }
                }else
                {
                    pis.push_back (node->data);
                }
            }
        }


        int ndir = dir+1;
        if (ndir == Dim)
            ndir = 0;

        if (node->left && bmin[dir] <= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->left;
            _stackdir[_stackindex] = ndir;
        }
        if (node->right && bmax[dir] >= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->right;
            _stackdir[_stackindex] = ndir;
        }
    }
}

template <int Dim,class Real,class Data>
bool   AdTreeEx<Dim,Real,Data>::search(Real p[],
                                       Data& pos,
                                       const GeneralCriterion<Data>& fun,
                                       Real eps)
{
    Real bmin[Dim], bmax[Dim];
    for(int i=0;i<Dim;i++)
    {
        bmin[i]=p[i]-eps;
        bmax[i]=p[i]+eps;
    }

    Array1D<Data> pis;
    getIntersecting (bmin,bmax,pis);
    if(pis.empty())return false;
    if(pis.size()==1)
    {
        if(fun(pis[0])){pos= pis[0];return true;}
        return false;
    }

    for(int i=1;i<pis.size();i++)
    {
        if(fun(pis[i])){pos= pis[i];return true;}
    }
    return false;
}


///////////////////////////////////////////////////////////////*/



template<int Dim,class Real=Real>
struct AdTreeExNodeRef
{
    AdTreeExNodeRef():left(0),right(0),pi(0){}


    AdTreeExNodeRef<Dim,Real>  *left,*right;
    int pi;
    Real sep,data[Dim];
};


template <int Dim,class Real=double>
class AdTreeExRef
{
public:
    typedef AdTreeExNodeRef<Dim,Real> Node;
    struct AdTreeExRefDistance
    {
        virtual ~AdTreeExRefDistance(){}
        virtual Real distance(const Real* p, int pi)const=0;
    };
    enum{STACKSIZE=1000};
    AdTreeExRef(const Real* cmin=NULL,const Real* cmax=NULL);
    void   insert (const Real* p, int pi);
    int    find(const Real* p,Real tol)const;
    void   deleteElement (int pi);
    void   getIntersecting (const Real* bmin,const Real* bmax,Array1D<int>& pis)const;
    void   getIntersecting (const Real* bmin,const Real* bmax, GeneralCriterion<int>& pFun)const;
    bool   search(const Real* pt,int&, const GeneralCriterion<int>&);
    void   reset ()const;
    bool   isIn(const Real* p,Real eps)const;

private:
    Real  _cmin[Dim], _cmax[Dim],_eps;
    Node*            _root;
    Array1D<Node*>         _ela;
    mutable Array1D<Node*> _stack;
    mutable Array1D<int>   _stackdir;
    mutable int            _stackindex;
    SlObjectPool<Node>     _dataPool;
};


template <int Dim,class Real>
AdTreeExRef<Dim,Real>::AdTreeExRef(const Real* acmin,const Real* acmax)
:_stack(STACKSIZE),
_stackdir(STACKSIZE)
{
    if(acmin!=NULL && acmax!=NULL)
    {
        memcpy (_cmin, acmin, Dim * sizeof(Real));
        memcpy (_cmax, acmax, Dim * sizeof(Real));
    }else
    {
        for(int i=0;i<Dim;i++)
        {
            _cmin[i]= -Tolerance::MAX_REAL;
            _cmax[i]= -Tolerance::MAX_REAL;
        }

    }

    _root =_dataPool.add(Node());
    _root->sep = (_cmin[0] + _cmax[0])*0.5;
    _eps=_root->sep*EPS_RALATIVE;
}

template <int Dim,class Real>
void   AdTreeExRef<Dim,Real>:: reset ()const
{
    _stack[1] = _root;
    _stackdir[1] = 0;
    _stackindex = 1;
}

template <int Dim,class Real>
bool   AdTreeExRef<Dim,Real>::isIn(const Real* p,Real eps)const
{
    for(int i=0;i<Dim;i++)
    {
        if(p[i]>=(_cmax[i]+eps)) return false;
        if(p[i]<=(_cmin[i]-eps)) return false;
    }
    return true;
}

template <int Dim,class Real>
void  AdTreeExRef<Dim,Real>::insert (const Real* p, int pi)
{
    int lr;
    Real bmin[Dim];
    Real bmax[Dim];
    memcpy (bmin, _cmin, Dim * sizeof(Real));
    memcpy (bmax, _cmax, Dim * sizeof(Real));

    Node *node=0,*next = _root;
    int dir = 0;
    while (next)
    {
        node = next;
        if (node->pi == 0)
        {
            memcpy (node->data, p, Dim * sizeof(Real));
            node->pi = pi+1;
            int cap=_ela.capacity();
            if (cap <(int) pi+1)_ela.resize(2*(pi+1));
            _ela[pi] = node;
            return;
        }

        if (node->sep > p[dir])
        {
            next = node->left;
            bmax[dir] = node->sep;
            lr = 0;
        } else
        {
            next = node->right;
            bmin[dir] = node->sep;
            lr = 1;
        }

        dir++;
        if (dir == Dim)
            dir = 0;
    }


    next = _dataPool.add(Node());
    memcpy (next->data, p, Dim * sizeof(Real));
    next->pi = pi+1;
    next->sep = (bmin[dir] + bmax[dir])*0.5;
    int cap=_ela.capacity();
    if (cap <(int) pi+1)_ela.resize(2*(pi+1));
    _ela[pi] = next;


    if (lr)
        node->right = next;
    else
        node->left = next;
}

template <int Dim,class Real>
void  AdTreeExRef<Dim,Real>::deleteElement (int pi)
{
    Node * node = _ela[pi];
    node->pi = 0;

    node = node->father;
    while (node)
    {
        node->nchilds--;
        node = node->father;
    }
}

template <int Dim,class Real>
void AdTreeExRef<Dim,Real> :: getIntersecting (const Real* bmin,
                                               const Real*  bmax,
                                               Array1D<int>& pis)const
{
    Node * node;
    int dir;
    reset ();

    while (_stackindex)
    {
        node = _stack[_stackindex];
        dir =  _stackdir[_stackindex];
        _stackindex--;

        if (node->pi != 0)
        {
            if (IsInD1(Dim,bmin,bmax,node->data))
            {
                pis.push_back (node->pi-1);
            }
        }


        int ndir = dir+1;
        if (ndir == Dim)
            ndir = 0;

        if (node->left && bmin[dir] <= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->left;
            _stackdir[_stackindex] = ndir;
        }
        if (node->right && bmax[dir] >= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->right;
            _stackdir[_stackindex] = ndir;
        }
    }
}



template <int Dim,class Real>
void AdTreeExRef<Dim,Real> ::getIntersecting (const Real*  bmin,
                                              const Real*  bmax,
                                              GeneralCriterion<int>& pFun)const
{

    Node * node;
    int dir;
    reset ();

    while (_stackindex)
    {
        node = _stack[_stackindex];
        dir =  _stackdir[_stackindex];
        _stackindex--;

        if (node->pi != 0)
        {
            if (IsInD1(Dim,bmin,bmax,node->data))
            {
                pFun(node->pi-1);
            }
        }


        int ndir = dir+1;
        if (ndir == Dim)
            ndir = 0;

        if (node->left && bmin[dir] <= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->left;
            _stackdir[_stackindex] = ndir;
        }
        if (node->right && bmax[dir] >= node->sep)
        {
            _stackindex++;
            _stack[_stackindex] = node->right;
            _stackdir[_stackindex] = ndir;
        }
    }
}

template <int Dim,class Real>
int AdTreeExRef<Dim,Real>::find(const Real* p,Real eps)const
{
    int pos=-1;
    Real bmin[Dim], bmax[Dim];
    for(int j=0;j<Dim;j++)
    {
        bmin[j]=p[j]-eps;
        bmax[j]=p[j]+eps;
    }

    Array1D<int> pis;
    getIntersecting (bmin,bmax,pis);
    if(pis.empty())return pos;
    if(pis.size()==1)
    {
        pos= pis[0];
        return pos;
    }

    Real dis=DistanceSqrt(p,_ela[0]->data);
    int index=0;
    for(UInt i=1;i<pis.size();i++)
    {
        Real tmp=DistanceSqrt(p,_ela[i]->data);
        if(tmp<dis)
        {
            dis=tmp;
            index=i;
        }
    }

    pos= pis[index];
    return pos;
}

template <int Dim,class Real>
bool   AdTreeExRef<Dim,Real>::search(const Real*  p,int& pos,const GeneralCriterion<int>& fun)
{
    Real bmin[Dim], bmax[Dim];
    for(int j=0;j<Dim;j++)
    {
        bmin[j]=p[j]-_eps;
        bmax[j]=p[j]+_eps;
    }

    Array1D<int> pis;
    getIntersecting (bmin,bmax,pis);
    if(pis.empty())return false;
    if(pis.size()==1)
    {
        if(fun(0))
		{
			pos= pis[0];
			return true;
		}
        return false;
    }

    for(int i=1;i<pis.size();i++)
    {
        if(fun(i)){pos= pis[i];return true;}
    }
    return false;
}
SL_MSH_NS_END

#endif

