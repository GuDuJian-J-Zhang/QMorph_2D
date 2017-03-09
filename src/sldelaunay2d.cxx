//#include <sldelaunay2d_std.h>
#include "sldelaunaymshtools.hpp"
//#include <slmesh25config.h>
#include <sstream>///< added by ZhangJun 2016.06.13
#include <boost/algorithm/string.hpp>///< added by ZhangJun 2016.06.13
#include <iostream>///< added by ZhangJun 2016.06.13

#ifdef __SLMESH25CFG_SLNS_H_
SL_NS_BEGIN
#include <sldelaunay2d.h>
namespace MESH25 {
#else
#include "sldelaunay2d.h"
SL_MSH_NS_BEGIN
#endif
#ifdef __SLMESH25CFG_SLNS_H_
#define  D_NEW  SL_NEW
inline HashValue Hash(HashValue h1,HashValue h2)   {return HashValue(((h1 << 16) | (h1 >> 16)) ^ h2);}
#endif
//######################################################################
#ifdef __SLDELAUNAY2D_STD_H_
#define  D_NEW  new
Real Epsilon ()
{
    Real r;

    r = 1.0E+00;

    while ( 1.0E+00 < ( Real ) ( 1.0E+00 + r )  )
    {
        r = r / 2.0E+00;
    }

    return ( 2.0E+00 * r );
}

Real Tolerance::EPS_COMPUTER = Epsilon ();
Real Tolerance::EPS          = 1.e-10;
Real Tolerance::EPS_RELATIVE = 1.e-4;
Real Tolerance::MAX_REAL     = 1.e20;


static int GeoPredicateInited=0;

Real Incircle_Math(const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d)
{
    if(GeoPredicateInited==0)
    {
        exactinit();
        GeoPredicateInited=1;
    }
    return incircle((Real*)(a.memory()),
        (Real*)(b.memory()),
        (Real*)(c.memory()),
        (Real*)(d.memory()));
}

Real Insphere_Math(const Point3D& a,const Point3D& b,const Point3D& c,const Point3D& d,const Point3D& e)
{
    if(GeoPredicateInited==0)
    {
        exactinit();
        GeoPredicateInited=1;
    }

    return insphere((Real*)(a.memory()),
        (Real*)(c.memory()),
        (Real*)(b.memory()),
        (Real*)(d.memory()),
        (Real*)(e.memory()));
}


/** \brief 判断三点的位置关系
    \return + 点c位于有向线段a--b的左侧
	\return - 点c位于有向线段a--b的右侧
	\return 0 三点共线
*/
Real Orient2d_Math(const Point2D &a, const Point2D &b, const Point2D &c)
{
    if(GeoPredicateInited==0)
    {
        exactinit();
        GeoPredicateInited=1;
    }
    return orient2d((Real*)(a.memory()),
        (Real*)(b.memory()),
        (Real*)(c.memory()));
}

Real Orient3d_Math(const Point3D& a,const Point3D& b,const Point3D& c,const Point3D& d)
{
    if(GeoPredicateInited==0)
    {
        exactinit();
        GeoPredicateInited=1;
    }
    return orient3d((Real*)(a.memory()),
        (Real*)(c.memory()),
        (Real*)(b.memory()),
        (Real*)(d.memory()));
}
#endif


//########################################################################
typedef int (*Fun_Message)(int type,char*msg);
int MessageDefault(int type,char*msg) {printf(msg); return 0;}
#define ErrorCodes  DelaunayAFT2DInterface
Real  QuaTriangle(Real a,Real b,Real c);
///
bool     MinMaxSwap2D(const Point2D&, const Point2D&,const Point2D&,const Point2D&);
bool     MinMaxSwap3D(const Point3D&, const Point3D&,const Point3D&,const Point3D&);
int      Circumcircle2D(const Point2D&,const Point2D& pt2,const Point2D& pt3,Point2D&,Real&);



RandomGen::RandomGen(Int64 j) 
: v(4101842887655102017),
w(1) 
{
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
}

Int64 RandomGen::int64() 
{
    u = u * 2862933555777941757 + 7046029254386353087;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    Int64 x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
}

UInt RandomGen::randInt(UInt num)
{
    UInt rand=(UInt)int64();
    return rand%num;
}
//#######################################################

//###############################################################
struct SwapCheck
{
    virtual bool  pass(const UInt nds[4],const UInt tris[2])=0;
};

typedef TripleT<UInt>           Triangle;
typedef Tuple<UInt,4>           Quad;
typedef Array1D<Tuple<UInt,2> > EdgeArray; ///< 存储边的两个端点编号(节点编号从0开始)
typedef Array1D<Point3D>        Point3DArray;
typedef Array1D<Triangle>       TriangleArray;
typedef Array1D<Quad >          QuadArray;
TriangleArray* gTriArray=NULL;
typedef StdPair<UInt,char>      Link;
template<typename PT>
class MesherBasic
{
public:
    enum {RESERVE0=0,ERASED,IN_DOMAIN,OUT_DOMAIN,FROZEN,QUAD_TRI,USER0,MAX_FLAG, FRONT_NODE_USED, BOUNDARY_NODE, FRONT_NODE, UNACTIVE};
    typedef Array1D<PT>   PointArray;
    typedef Tuple<Link,3>      TriangleLink; ///< 三角形单元三条边的链接信息(每个Link中存储了邻接面片的编号，及该边在邻接单元中的局部编号)
    typedef Tuple<Real,4>      TriangleCircle; ///< 存储三角形单元的外接圆圆心及半径(二维情况下共存储三个数据，三维情况下存储四个数据)
    struct LinkVisitor
    {
        virtual bool visit(const Link&)=0;
    };

    class Edge 
    {
    public:
        Edge(UInt tri=0,int whi=-1,bool reverse=false)
            :_tri(tri),_whi(whi),_reversed(reverse),_hash(0)
        {

        }
        Edge(const Link& lnk,bool reverse)
            :_tri(lnk.first),_whi(lnk.second),_reversed(reverse),_hash(0)
        {

        }
        HashValue   hash()const;
        bool        operator==(const Edge&)const;
        void        link(UInt& tri,int& whi,bool& reversed)const
        {
            tri=_tri;whi=(int)(_whi);
            reversed=(_reversed==1);
        }

        Link       link()const
        {
            return Link(_tri,(char)_whi);
        }

    protected:
        UInt                _tri; ///< 边的入射三角形单元的编号
        short               _whi; ///< 该边在当前三角形中的局部编号
        bool                _reversed;
        mutable  HashValue  _hash;
    };

    typedef   HashCouple<UInt>                           Edge2;
    typedef   HashSetDefault<Edge>                       EdgeSet;
    typedef   HashSetDefault<Edge2>                      EdgeSet2;
    typedef   Array1D<TriangleLink>   TriangleAdjArray;

    MesherBasic()
    {
        gTriArray=&_workingEleArray;
    }
    ~MesherBasic()
    {
        gTriArray=NULL;
    }
public:
    void   getTriangles(Array1D<Triangle>&,bool all=false);
    UInt   triangleNum()const{return _workingEleArray.size();}
    UInt   nodeNum()const{return _nodesArray.size();}
    const  Point2D& nodeCoord(UInt i)const{return _nodesArray[i];}
    int    getNodeId(UInt id);
    const  Triangle& getTriangle(UInt i)const{return _workingEleArray[i];}
    bool  isCheckedTri(UInt tri,int i)const{return _triMark[tri].test(i);}
    void  checkinTri(UInt tri,int i){_triMark[tri].set(i);}
    void  checkoutTri(UInt tri,int i){_triMark[tri].clear(i);}
    virtual void getNodeCoord(UInt ni,int /*dim*/,Real coord[3])const;
    Link  adjacentTri(const Link& lnk)const;
    void  setNodeLink(UInt,const Link&);
    void  setTriLink(UInt,int,const Link&);
    void  exportBlk(std::ostream& out,int dim,bool all);
    bool  visitNodeBall(UInt ndi,LinkVisitor&);
    bool  visitNodeBall(UInt ndi,Array1D<Link>&); ///< 获取节点的邻接三角形数组

	bool  getAdjacentEdge(UInt ndi,LinkVisitor&); ///< 获取节点的邻接边数组(added by ZhangJun 2016.06.07)
	bool  getAdjacentEdge(UInt ndi,Array1D<Link>&); ///< 获取节点的邻接边数组(added by ZhangJun 2016.06.07)
	Real  length(UInt id1, UInt id2); ///< 计算两点之间的距离(added by ZhangJun 2016.06.14)

    inline bool  isSuperNode(UInt ni)const;
    virtual int getTriGroupId(UInt)const{return -1;}
protected:
    PointArray        _nodesArray; ///< 网格节点数组
    EdgeArray         _edges;
    EdgeSet2          _boundaryEdgeSet;
    TriangleArray     _workingEleArray; ///< 网格单元数组
    TriangleAdjArray  _adjElementsArray; ///< 存储三角形单元的邻接单元(first: 邻接单元的编号；second: 公共边在first面片中的局部编号)
    Array1D<Link>     _linkNodeArray; ///< 每个节点存储其邻接面片组中的任一面片(first: 邻接面片的局部编号；second: 点在该面片中的局部编号)
    Array1D<UIntBits> _triMark,_nodeMark; ///< 用于对三角形单元和节点做标记

    Array1D<UInt>     _inValidArray; ///< 无效单元数据，用于存储被标记删除的三角形单元的编号
    IdMap             _idMap,_idMapReverse;
    UInt              _inputNodeSize;
    int                _maxNdId,_lastError;
    Fun_Message       _msgFun;
};

class Delaunay2D : public MesherBasic<Point2D>
{
public:
    typedef MesherBasic<Point2D>::PointArray PointArray;
    typedef DelaunayAFT2DInterface::Input Input;
    Delaunay2D(const  Array1D<Point2D>& pts,const EdgeArray& edges,Fun_Message fun=NULL);
    Delaunay2D(Input* inputPtr,Fun_Message fun=NULL);
    virtual ~Delaunay2D();
    int    delaunay();

    inline Point2D orginalCoord(const Point2D&)const;
public:
    int   delaunayKernal();
    int   recovery();
    int   recovery(UInt,UInt);
    enum 
    {
        SAMPLEFACTOR=11,
        OUTSIDE,INSIDE,ONEDGE,ONVERTEX,OUTDOMAIN,
        FINISHED,CONTINUE,SWAPED,NEXT_ADJ,ERRO_INPUT,
    };

    typedef Array1D<TriangleCircle> TriCircleArray;
    void  createSuper(const  Array1D<Point2D>& pts);
    int   locate(const Point2D& pt,Link& whi);
    Real  distance2(UInt whi,const Point2D&)const;
    int   orientTri(UInt whi,const Point2D&,int&);
    int   nextTri(const Point2D&,UInt& whiTri,int&);
    bool  calCircumcircle(UInt whi);
    bool  isInCircumcircle(UInt whi,const Point2D& pt);
    bool  allNodesInserted(Array1D<UInt>& lostNodes);
    bool  formCave(UInt ni,UInt whiTri,bool correct,EdgeSet&,Array1D<UInt>&);
    bool  getBoundaryEdge(Array1D<UInt>&,EdgeSet&,UInt,bool);
    int   validTri(UInt,UInt,UInt)const;
    int   validTri(UInt nd0,UInt nd1,const Point2D&)const;
    int   insert(UInt tri,UInt ndHd,bool force);
    int   insert(UInt nodeHd,bool force);
    int   optConnection();
    void  localReconnect(SwapCheck* fun);
    void  exportBoundaryBlk(std::ostream& out,int dim);
    void  getEdgeNodes(UInt tri,int whi,UInt nodes[2])const;
    int   getPiple(UInt edge[2],Array1D<Edge>& edgestack);

    int   recoverySwap(UInt edge[2],const Link& lnk);
    bool  getSwapData(const Link&,UInt tris[2],UInt nd[4],Link lnk[4],bool&)const;
    int   swapIt(const Link&,bool,SwapCheck*,Link*,Link*);
    UInt  addTri(UInt,UInt,UInt);
    void  eraseTri(UInt);
    //###############################################################################

    //#################################################################################
    void   markout();
    void   markSuperNode(const Link& lnk);
    void   markAdjacent(UInt tri, int mk,UInt&);
    UInt   markedTriNum()const;
    bool   isInDomain(const Point2D&)const;

	Real  getFScale() { return _fScale; } ///< added by ZhangJun 2016.06.04
	
	Point2D getMinPoint() { return _min; } ///< added by ZhangJun 2016.06.04

protected:
    EdgeArray*        _edgesPtr;
    TriCircleArray    _triCircleArray; ///< 存储三角形单元的外心坐标及外接圆半径
    Point2D           _min,_max,_minSc,_maxSc;
    UInt              _samples,_lastTri;
    Real              _fScale,_fScaleInv; ///< 缩放因子
};

Delaunay2D::~Delaunay2D()
{
    gTriArray=NULL;
}
struct SwapCheckerRecover : public SwapCheck
{
public:
    SwapCheckerRecover(Delaunay2D* msh)
        :_meshPtr(msh)
    {
    }

    bool pass(const UInt nds[4],const UInt tris[2]);
protected:
    Delaunay2D*  _meshPtr;
};

struct MinMaxChecker2D : public SwapCheckerRecover
{
    MinMaxChecker2D(Delaunay2D* msh)
        :SwapCheckerRecover(msh){}
    bool  pass(const UInt nds[4],const UInt tris[2]);
};



Delaunay2D::Delaunay2D(const  Array1D<Point2D>& pts,const EdgeArray& edges,Fun_Message fun)
:_edgesPtr((EdgeArray*)(&edges))
{
    _lastError=ErrorCodes::eOk;
    _msgFun=MessageDefault;
    if(fun!=NULL) _msgFun=fun;
    _samples=0;
    createSuper(pts);
}

/** \brief 边界恢复操作
*/
int   Delaunay2D::recovery()
{
    if(0)
    {
        std::ofstream out("boundary.nas");
        exportBoundaryBlk(out,2);
    }
    UInt i,sz=_edgesPtr->size();
    for(i=0;i<sz;i++)
    {
        int rt= recovery((*_edgesPtr)[i][0],(*_edgesPtr)[i][1]);
        if(rt!=ErrorCodes::eOk) return rt;
        _boundaryEdgeSet.insert(Edge2((*_edgesPtr)[i][0],(*_edgesPtr)[i][1]));
    }
    return ErrorCodes::eOk;
}

int   Delaunay2D::delaunay()
{
    if(0)
    {
        std::ofstream out("boundary.nas");
        exportBoundaryBlk(out,2);
    }
    _msgFun(1,"2-d delaunay meshing...\n");
    if(_lastError!=ErrorCodes::eOk) return _lastError;
    _lastError = delaunayKernal(); ///< 插点操作
    //_msgFun(2,"*");
    if(_lastError!=ErrorCodes::eOk) return _lastError;
    if(_edgesPtr==NULL) return _lastError;
    _lastError=recovery(); ///< 边界恢复
    //_msgFun(2,"*");
    if(_lastError!=ErrorCodes::eOk) 
    {
        return _lastError;
    }
    markout(); ///< 标记边界内外单元
    _msgFun(1,"2-d delaunay mesh improving...\n");
    return optConnection(); ///< 拓扑优化
}

/** \brief Delaunay插点操作
*/
int   Delaunay2D::delaunayKernal()
{
    _samples=0;
    UInt i,sz = _nodesArray.size() - 3;
    RandomGen rand(sz);
    Array1D<UInt> indexArray(sz);
    for(i = 0; i < sz; i++) 
        indexArray[i]=i;
    for(i = 0;i < sz; i++) 
        Swap(indexArray[i],indexArray[rand.randInt(sz)]);

    for(i = 0; i < sz; i++)
    {
        if(0)
        {
            std::ofstream out("delaunay_debug.nas");
            exportBlk(out,2,true);
        }
        if(ErrorCodes::eInValidInput==insert(indexArray[i],false))
        {
            return ErrorCodes::eInValidInput;
        }
    }

    if(0)
    {
        std::ofstream out("delaunay_debug.nas");
        exportBlk(out,2,true);
    }
    bool allNodeInserted=false;
    int count=0;
    Array1D<UInt> lostNodes;
    allNodesInserted(lostNodes);
    if(lostNodes.empty()) 
    {
        return ErrorCodes::eOk;
    }
    while(!allNodeInserted && count++<10)
    {
        sz=lostNodes.size();
        for(i=0;i<sz;i++)
        {
            if(ErrorCodes::eInValidInput==insert(lostNodes[i],false))
            {
                return ErrorCodes::eInValidInput;
            }
        } 

        Array1D<UInt> lostNodes2;
        allNodesInserted(lostNodes2);
        UInt inserted=lostNodes2.size()-sz;
        if(inserted==0)
        {
            sz=lostNodes2.size();
            for(i=0;i<sz;i++)
            {
                if(ErrorCodes::eInValidInput==insert(lostNodes2[i],true))
                {
                    return ErrorCodes::eInValidInput;
                }
            }
            lostNodes.clear();
            allNodesInserted(lostNodes);
        }else
        {
            lostNodes.swap(lostNodes2);
        }
        if(lostNodes.empty()) 
        {
            return ErrorCodes::eOk;
        }
    }

    return ErrorCodes::eMeshingFailed;
}

Delaunay2D::Delaunay2D(Input* inputPtr,Fun_Message fun)
{
    _maxNdId=-INT_MAX;
    if(inputPtr->maxNodeId()>_maxNdId) _maxNdId=inputPtr->maxNodeId();
    _lastError=ErrorCodes::eOk;
    _msgFun=MessageDefault;
    if(fun!=NULL) _msgFun=fun;
    int id;double xy[2];
    inputPtr->node_begin();
    _min=Point2D(DBL_MAX,DBL_MAX);
    _max=Point2D(-DBL_MAX,-DBL_MAX);
    _inputNodeSize=0;
    while(1)
    {
        if(!inputPtr->node_next(id,xy)) break;
        Point2D apt(xy[0],xy[1]);
        _nodesArray.push_back(apt);
        _min.setToMin(apt);
        _max.setToMax(apt);
        _idMap.insert(MakePair(id,_inputNodeSize++));
        if(id>_maxNdId) _maxNdId=id;
    }
    int nodes[2];
    inputPtr->edge_begin();
    Tuple<UInt,2> aEdge;
    while(1)
    {
        if(!inputPtr->edge_next(nodes)) break;
        IdMap::iterator iter=_idMap.find(nodes[0]);
        if(iter==_idMap.end())
        {
            _lastError=ErrorCodes::eInValidInput;
            char buffer[1024];
            sprintf(buffer,"There is no node as:%d",nodes[0]);
            _msgFun(0,buffer);
            return;
        }
        aEdge[0]=iter->second;
        iter=_idMap.find(nodes[1]);
        if(iter==_idMap.end())
        {
            char buffer[1024];
            sprintf(buffer,"There is no node as:%d",nodes[1]);
            _msgFun(0,buffer);
            _lastError=ErrorCodes::eInValidInput;
            return;
        }
        aEdge[1]=iter->second;
        _edges.push_back(aEdge);
    }
    _edgesPtr=&_edges;

    Real width = _max[0] - _min[0];
    Real tmp=_max[1] - _min[1];
    if (tmp > width) 
    {
        width = tmp;
    }
    if (width == 0.0) 
    {
        width = 1.0;
    }
    _fScaleInv=width;
    _fScale=1.0/width;
    _linkNodeArray.resize(_inputNodeSize+3);
    _minSc.set(-Tolerance::EPS,-Tolerance::EPS);
    _maxSc= (_max - _min)*_fScale;
    _maxSc[0]+=Tolerance::EPS;_maxSc[1]+=Tolerance::EPS;
    for (UInt i = 0; i < _inputNodeSize; i++)
    {
        _nodesArray[i] = (_nodesArray[i] - _min)*_fScale;
        _linkNodeArray[i]=Link(InValidEntity,-1);
    }

    _nodesArray.push_back(Point2D(-1.0,-1.0));
    _nodesArray.push_back(Point2D(4.0,-1.0));
    _nodesArray.push_back(Point2D(-1.0,4.0));
    UInt tri=addTri(_inputNodeSize,_inputNodeSize+1,_inputNodeSize+2);
    _linkNodeArray[_inputNodeSize]=Link(tri,0);
    _linkNodeArray[_inputNodeSize+1]=Link(tri,1);
    _linkNodeArray[_inputNodeSize+2]=Link(tri,2);
}

inline Point2D Delaunay2D::orginalCoord(const Point2D& pt)const
{
    return pt*_fScaleInv+_min;
}

template<typename PT>
int    MesherBasic<PT>::getNodeId(UInt id)
{
    if(_idMap.empty()) return (int)id;
    if(_idMapReverse.empty())
    {
        IdMap::iterator iter=_idMap.begin();
        for(;iter!=_idMap.end();++iter)
        {
            _idMapReverse.insert(MakePair(iter->second,iter->first));
        }
        UInt i=0,sz=_nodesArray.size();
        for(i=_inputNodeSize+3;i<sz;i++)
        {
            _idMapReverse.insert(MakePair(i,++_maxNdId));
        }
        for(i=_inputNodeSize;i<_inputNodeSize+3;i++)
        {
            _idMapReverse.insert(MakePair(i,++_maxNdId));
        }
    }
    IdMap::iterator iter=_idMapReverse.find(id);
    if(iter!=_idMapReverse.end())
    {
        return iter->second;
    }
    Assert(0);
    return _maxNdId;
}

/** \brief 对待剖分区域进行归一化并创建超级三角形
*/
void  Delaunay2D::createSuper(const  Array1D<Point2D>& pts)
{
    Real EPS=Tolerance::EPS_COMPUTER;
    UInt i,sz=pts.size();
    _inputNodeSize=sz;
    _min=_max=pts[0];
    for(i=1;i<sz;i++)
    {
        _min.setToMin(pts[i]);
        _max.setToMax(pts[i]);
    }

    Real width = _max[0] - _min[0];
    Real tmp=_max[1] - _min[1];
    if (tmp > width) {
        width = tmp;
    }
    if (width == 0.0) {
        width = 1.0;
    }
    _fScaleInv=width;
    _fScale=1.0/width;
    _nodesArray.resize(sz+3);
    _linkNodeArray.resize(sz+3);
    for (i = 0; i < sz; i++)
    {
        _nodesArray[i] = (pts[i] - _min)*_fScale;
        _linkNodeArray[i]=Link(InValidEntity,-1);
    }
    _minSc.set(-Tolerance::EPS,-Tolerance::EPS);
    _maxSc= (_max - _min)*_fScale;
    _maxSc[0]+=Tolerance::EPS;_maxSc[1]+=Tolerance::EPS;
    _nodesArray[sz].set(-1.0,-1.0);
    _nodesArray[sz+1].set(4.0,-1.0);
    _nodesArray[sz+2].set(-1.0,4.0);
    UInt tri=addTri(sz,sz+1,sz+2);
    _linkNodeArray[sz]=Link(tri,0);///< 设置超级三角形三个顶点的链接三角形
    _linkNodeArray[sz+1]=Link(tri,1);
    _linkNodeArray[sz+2]=Link(tri,2);
}

/** \brief 使用三点构造三角形
    \pram nd1~nd3 三个点的编号
    \return 新生成的三角形单元的编号
*/
UInt  Delaunay2D::addTri(UInt nd1,UInt nd2,UInt nd3)
{
    UInt ti=_workingEleArray.size();
    if(!_inValidArray.empty())
    {
        ti=_inValidArray[_inValidArray.size()-1];
        _inValidArray.pop_back();
        _workingEleArray[ti][0]=nd1;
        _workingEleArray[ti][1]=nd2;
        _workingEleArray[ti][2]=nd3;
        checkoutTri(ti,ERASED);
    }else
    {
        _workingEleArray.push_back(Triangle(nd1,nd2,nd3));
        _adjElementsArray.push_back(TriangleLink());
        _triCircleArray.push_back(TriangleCircle());
        _triMark.push_back(UIntBits());
    }
    _triCircleArray[ti][3]=-1.0;
    _adjElementsArray[ti][0]=_adjElementsArray[ti][1]=_adjElementsArray[ti][2]=Link(InValidEntity,-1);
    _lastTri=ti;
    return ti;
}

template<typename PT>
void MesherBasic<PT>::getNodeCoord(UInt ni,int dim,Real coord[3])const
{
    coord[0]=_nodesArray[ni][0];
    coord[1]=_nodesArray[ni][1];
    coord[2]=0.0;
    if(dim==3)
        coord[2]=_nodesArray[ni][2];
}

/** \brief 把编号为ti的三角形标记删除
*/
void  Delaunay2D::eraseTri(UInt ti)
{
    _inValidArray.push_back(ti);
    checkinTri(ti,ERASED);
    _triCircleArray[ti][3]=-1.0;
    for(int i=0;i<3;i++)
    {
        _adjElementsArray[ti][i]=Link(InValidEntity,-1);
    }
}

template<typename PT>
void  MesherBasic<PT>::setNodeLink(UInt whi,const Link& lnk)
{
    _linkNodeArray[whi]=lnk;
}

template<typename PT>
void  MesherBasic<PT>::setTriLink(UInt tri,int whi,const Link& lnk)
{
    _adjElementsArray[tri][whi]=lnk;
    if(lnk.first!=InValidEntity)
    {
        _adjElementsArray[lnk.first][lnk.second]=Link(tri,whi);
    }
}

/** \brief 判断点pt是否在待剖分域内
*/
bool   Delaunay2D::isInDomain(const Point2D& pt)const
{
    if(pt[0]<_minSc[0]) return false;
    if(pt[1]<_minSc[1]) return false;
    if(pt[0]>_maxSc[0]) return false;
    if(pt[1]>_maxSc[1]) return false;
    return true;
}

/** \brief 判断点pt与三角形whi的位置关系
    \param whi 三角形的编号
	\param pt 二维节点
	\param weg 用于存储节点在三角形中的位置 \n
	     case1 在边上 weg = 所在边的局部编号 \n
		 case2 在顶点上 weg = 所在顶点的局部编号 \n
		 case3 在三角形内 weg = 3
*/
int   Delaunay2D::orientTri(UInt whi,const Point2D& pt,int& weg)
{
    int inCnt = 0, onEdgeCnt = 0;
    int onEdge[3] = {-1,-1,-1};

    for(int i = 0; i < 3; i++)
    {
        int isValidi = validTri(_workingEleArray[whi][TRI_EDGE_VERTEX[i][0]],_workingEleArray[whi][TRI_EDGE_VERTEX[i][1]],pt);
        if(isValidi==-1)
        {
            weg = i;
            return OUTSIDE;
        }else if(isValidi == 1)
        {
            inCnt++;
        }else if(isValidi == 0)
        {
            onEdge[onEdgeCnt++] = i;
        }
    }

    if(inCnt == 3) return INSIDE;
    if(onEdgeCnt >= 3)
    {
        if(0)
        {
            std::ofstream out("debugit.nas");
            exportBlk(out,2,true);
        }
        Assert(onEdgeCnt < 3);
    }
    if(onEdgeCnt == 1)
    {
        weg = onEdge[0]; ///< 此时weg为pt所在边的局部编号
        return ONEDGE;
    }else
    {
        weg = 3 - onEdge[0] - onEdge[1]; ///< 此时weg为pt所在顶点的局部编号
        return ONVERTEX;
    }

}

/** \brief 通过三角形单元之间的邻接关系来确定点pt所在的三角形
    \param pt 待定位点
	\param whiTri 作为输出时表示上一步判断失败的三角形, 作为输出时表示pt所在的三角形
	\param whi 作为输入时表示三角形单元whiTri中局部编号为whi的边，点pt落在该边的外侧，作为输出时，表示点与三角形whiTri的位置关系
	\return loc 点pt与whiTri的具体位置关系
*/
int   Delaunay2D::nextTri(const Point2D& pt,UInt& whiTri,int& whi)
{
    if(whiTri == InValidEntity) return OUTSIDE;
    UInt ti = _adjElementsArray[whiTri][whi].first; ///< 获取邻接单元
    int ei= _adjElementsArray[whiTri][whi].second;  ///< 获取边whi在邻接单元中的编号
    whiTri=ti;
    int inCnt=0,onEdgeCnt=0;
    int onEdge[2]={-1,-1};

    for(int i=0;i<2;i++)
    {
        int ej=TRI_VERTEX_EDGE[ei][i];
        int isValidi=validTri(_workingEleArray[whiTri][TRI_EDGE_VERTEX[ej][0]],
            _workingEleArray[whiTri][TRI_EDGE_VERTEX[ej][1]],pt);
        if(isValidi==-1)
        {
            whi=ej;
            return OUTSIDE;
        }else if(isValidi==1)
        {
            inCnt++;
        }else if(isValidi==0)
        {
            onEdge[onEdgeCnt++]=ej;
        }
    }

    if(inCnt==2) return INSIDE;
    if(onEdgeCnt>=3)
    {
        if(0)
        {
            std::ofstream out("debugit.nas");
            exportBlk(out,2,true);
        }
        Assert(onEdgeCnt<3);
    }
    if(onEdgeCnt==1)
    {
        whi=onEdge[0];
        return ONEDGE;
    }else
    {
        whi=3-onEdge[0]-onEdge[1];
        return ONVERTEX;
    }
}

/** \brief 确定点pt所在的三角形
    \param pt 待定位点
	\pram whiTri 从该三角形开始查找，且查找结果存储在whiTri中\n
	      whiTri.first = 点pt所在的三角形
		  whiTri.second = 3 表示点在三角形内
		  whiTri.second = 0~2 表示点在三角形的边上或与某一点重合
    \return loc 13:在三角形内；14:在边上；15:与某一点重合
*/
int Delaunay2D::locate(const Point2D& pt,Link& whiTri)
{
    if(!isInDomain(pt)) return  OUTDOMAIN; 
    static RandomGen rand;
    UInt trinum = _workingEleArray.size();
    if(whiTri.first >= trinum) whiTri.first = _lastTri;
    if(isCheckedTri(whiTri.first,ERASED)) whiTri.first = _lastTri;
    Real searchdist = distance2(whiTri.first,pt); ///< 计算点与三角形外心的距离


    while (SAMPLEFACTOR * _samples * _samples * _samples* _samples < trinum)
    {
        _samples++;
    }

    UInt i =1;

    while(i < _samples) ///< 随机找几个三角形进行判断，取其中外心与待插点最近的三角形
    {
        UInt whiTrii = rand.randInt(trinum);
        if(isCheckedTri(whiTrii,ERASED)) continue;
        Real dis = distance2(whiTrii,pt);
        if(dis < searchdist)
        {
            searchdist = dis;
            whiTri.first = whiTrii;
        }
        i++;
    }

    int whi = 3;
    int loc = orientTri(whiTri.first,pt,whi); ///< 判断点与三角形的位置关系
    if(loc==INSIDE || loc==ONEDGE || loc==ONVERTEX)
    {
        whiTri.second=(char)whi;
        return loc;
    }

    for (i = 0; i < trinum; i++)
    {
        loc = nextTri(pt,whiTri.first,whi); ///< 通过邻接关系取下一单元进行判断
        if(loc==INSIDE || loc==ONEDGE || loc==ONVERTEX)
        {
            whiTri.second=(char)whi;
            return loc;
        }
    }

    // when call here, there may be a problem
    for (i = 0; i < trinum; i++)
    {
        if(isCheckedTri(i,ERASED)) continue;
        whiTri.first=i;

        loc=orientTri(whiTri.first,pt,whi);
        if(loc==INSIDE || 
            loc==ONEDGE ||
            loc==ONVERTEX)
        {
            whiTri.second=(char)whi;
            return loc;
        }

    }
    return OUTDOMAIN;
}

/** \brief 判断点pt是否落在三角形whi的外接圆内
*/
bool  Delaunay2D::isInCircumcircle(UInt whi,const Point2D& pt)
{
    if(calCircumcircle(whi))
    {
        Real dx=pt[0]-_triCircleArray[whi][0];
        Real dy=pt[1]-_triCircleArray[whi][1];
        Real diff=_triCircleArray[whi][3]-sqrt(dx*dx+dy*dy);
        if(diff<=0.0) return false;
        return diff>Tolerance::EPS*_triCircleArray[whi][3];
    }else
    {
        Real det= Incircle_Math(_nodesArray[_workingEleArray[whi][0]],
            _nodesArray[_workingEleArray[whi][1]],
            _nodesArray[_workingEleArray[whi][2]],pt);
        return det>0.0;
    }

    return false;
}

/** \breif 计算点pt与三角形whi外心距离的平方值
    \param whi 三角形编号
	\param pt 二维节点
*/
Real  Delaunay2D::distance2(UInt whi,const Point2D& pt)const
{
    Point2D ptc=_nodesArray[_workingEleArray[whi][0]]+_nodesArray[_workingEleArray[whi][1]]+_nodesArray[_workingEleArray[whi][2]];
    ptc*=0.3333333333333333333;
    return ptc.distance2(pt);
}

/** \brief 判断由三点组成的三角形单元是否有效
*   对有效单元，点nd2应位于有向线段nd0--nd1的左侧
*/
int   Delaunay2D::validTri(UInt nd0,UInt nd1,UInt nd2)const
{
    Real det= Orient2d_Math(_nodesArray[nd0],_nodesArray[nd1],_nodesArray[nd2]);
    if(det>0) return 1;
    if (det<0) return -1;
    return 0;

    //return (det > Tolerance::EPS_COMPUTER ? +1 : (det < -Tolerance::EPS_COMPUTER ? -1 : 0));
}

/** \breif 判断由nd0、nd1、nd2三点组成的三角形是否合法
    三点成逆时针方向，则合法
*/
int   Delaunay2D::validTri(UInt nd0,UInt nd1,const Point2D& nd2)const
{
    Real det= Orient2d_Math(_nodesArray[nd0],_nodesArray[nd1],nd2);
    if(det>0) return 1;
    if (det<0) return -1;
    return 0;
    //return (det > Tolerance::EPS_COMPUTER ? +1 : (det < -Tolerance::EPS_COMPUTER ? -1 : 0));
}

void  Delaunay2D::getEdgeNodes(UInt tri,int whi,UInt nodes[2])const
{
    nodes[0]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]];
    nodes[1]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];
}

template<typename PT>
bool  MesherBasic<PT>::Edge::operator==(const Edge& edge)const
{
    if(hash()!=edge.hash()) return false;
    UInt nodes[2],nodes2[2];
    nodes[0]=(*gTriArray)[_tri][TRI_EDGE_VERTEX[_whi][0]];
    nodes[1]=(*gTriArray)[_tri][TRI_EDGE_VERTEX[_whi][1]];
    nodes2[0]=(*gTriArray)[edge._tri][TRI_EDGE_VERTEX[edge._whi][0]];
    nodes2[1]=(*gTriArray)[edge._tri][TRI_EDGE_VERTEX[edge._whi][1]];

    if(nodes[0]==nodes2[0]) return nodes[1]==nodes2[1];
    if(nodes[0]==nodes2[1]) return nodes[1]==nodes2[0];
    return false;
}

template<typename PT>
HashValue MesherBasic<PT>::Edge::hash()const
{
    if(_hash==0)
    {
        UInt nodes[2];
        nodes[0]=(*gTriArray)[_tri][TRI_EDGE_VERTEX[_whi][0]];
        nodes[1]=(*gTriArray)[_tri][TRI_EDGE_VERTEX[_whi][1]];
        if(nodes[0]>nodes[1])
        {
            _hash= Hash(nodes[1],nodes[0]);
        }else
        {
            _hash= Hash(nodes[0],nodes[1]);
        }
    }
    return _hash;

}

/** \brief 计算并保存三角形的外接圆圆心坐标及半径
    \return true 计算成功
	\return false 计算失败
*/
bool  Delaunay2D::calCircumcircle(UInt whi)
{
    if(_triCircleArray[whi][3]>0.0) return true;

    Point2D cen;Real rad;
    if(Circumcircle2D(_nodesArray[_workingEleArray[whi][0]],
        _nodesArray[_workingEleArray[whi][1]],
        _nodesArray[_workingEleArray[whi][2]],
        cen,rad)==ErrorCodes::eOk)
    {
        _triCircleArray[whi][0]=cen[0];
        _triCircleArray[whi][1]=cen[1];
        _triCircleArray[whi][3]=rad;
        return true;
    }
    return false;
}

/** \brief 收集所有未被插入的节点
    \return false 不存在未被插入的节点
*/
bool Delaunay2D::allNodesInserted(Array1D<UInt>& lostNodes)
{
    UInt i,num,sz=_nodesArray.size();
    Array1D<int> flag(sz,0);

    num=_workingEleArray.size();
    for(i=0;i<num;i++)
    {
        if(isCheckedTri(i,ERASED)) continue;
        flag[_workingEleArray[i][0]]=1;
        flag[_workingEleArray[i][1]]=1;
        flag[_workingEleArray[i][2]]=1;
    }

    for(i=0;i<sz;i++)
    {
        if(flag[i]==0)
        {
            lostNodes.push_back(i);
        }
    }

    return lostNodes.empty();
}

/** \brief 收集内核面片及内核边
    \param[in] point 待插点的编号
	\param[in] whiTri 待插点所在的三角形单元的编号
	\param[out] caveEdges 用于存储收集到的内核边(待插点将于这些边组成新的三角形单元)
	\param[out] caveElements 用于存储收集到的内核单元(将被标记删除)
*/
bool  Delaunay2D::formCave(UInt point,UInt whiTri,bool correct,
                           EdgeSet& caveEdges,
                           Array1D<UInt>& caveElements)
{
    HashUIntSet elementsFilter;
    Stack<UInt> faces;
    faces.push(whiTri);
    elementsFilter.insert(whiTri);
    caveElements.push_back(whiTri);
    StdPair<EdgeSet::iterator,bool> inserter;
    while(!faces.empty())
    {///< 从whiTri出发，利用邻接关系收集所有外接圆圆心包含待插点point的三角形单元
        UInt ele=faces.top();
        faces.pop();
        for(int i=0;i<3;i++)
        {
            UInt elei=_adjElementsArray[ele][i].first;
            if(elei==InValidEntity) continue;

            if(elementsFilter.insert(elei).second)
            {
                if(isInCircumcircle(elei,_nodesArray[point]))
                {
                    faces.push(elei);
                    caveElements.push_back(elei);
                }
            }
        }
    }
    return getBoundaryEdge(caveElements,caveEdges,point,correct);
}

/** \brief 使用内核单元收集内核边
    \param[in] caveElements 内核单元，即外接圆圆心包含待插点pt的三角形单元
	\param[in] pt 待插点
	\param[out] caveEdges 内核边(即内核区域的边界边)
*/
bool Delaunay2D::getBoundaryEdge(Array1D<UInt>& caveElements,EdgeSet& caveEdges,
                                 UInt pt,bool correct)
{
    StdPair<EdgeSet::iterator,bool> inserter;
    if(!correct)
    {
        UInt i,j,sz=caveElements.size();
        for(i=0; i<sz; i++)
        {
            UInt aTri=caveElements[i];
            for(j=0; j<3; j++)
            {
                UInt elei=_adjElementsArray[aTri][j].first;///< 获取邻接三角形单元的编号
                if(elei!=InValidEntity)
                {
                    inserter=caveEdges.insert(Edge(_adjElementsArray[aTri][j],true));
                    if(!inserter.second)
                    {
                        caveEdges.erase(inserter.first);
                    }

                }else
                {
                    inserter=caveEdges.insert(Edge(aTri,j,false));
                    Assert(inserter.second);
                }
            }
        }

        for(EdgeSet::iterator iter=caveEdges.begin();iter!=caveEdges.end();++iter)
        {
            UInt tri;int whi;bool rev;
            iter->link(tri,whi,rev);
            int rt=-1;
            if(rev)
            {
                rt=validTri(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]],_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]],pt);
            }else
            {
                rt=validTri(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]],pt);
            }
            if(rt<=0) return false;
        }
    }else
    {
        UInt i,j,sz=caveElements.size();
        if(sz<1) return false;
        while(1)
        {
            sz=caveElements.size();
            HashUIntSet eleSet;
            caveEdges.clear();
            for(i=0;i<sz;i++)
            {
                UInt aTri=caveElements[i];
                for(j=0;j<3;j++)
                {
                    UInt elei=_adjElementsArray[aTri][j].first;
                    inserter=caveEdges.insert(Edge(aTri,j,false));
                    if(!inserter.second)
                    {
                        caveEdges.erase(inserter.first);
                    }
                }
                eleSet.insert(aTri);
            }


            for(EdgeSet::iterator iter=caveEdges.begin();iter!=caveEdges.end();++iter)
            {
                UInt tri;int whi;bool rev;
                iter->link(tri,whi,rev);
                int rt=validTri(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]],pt);
                if(rt<=0)
                {
                    eleSet.erase(tri);
                }
            }
            if(eleSet.size()==sz)
            {
                sz=caveElements.size();
                caveEdges.clear();
                for(i=0;i<sz;i++)
                {
                    UInt aTri=caveElements[i];
                    for(j=0;j<3;j++)
                    {
                        UInt elei=_adjElementsArray[aTri][j].first;
                        if(elei!=InValidEntity)
                        {
                            inserter=caveEdges.insert(Edge(_adjElementsArray[aTri][j],true));
                            if(!inserter.second)
                            {
                                caveEdges.erase(inserter.first);
                            }

                        }else
                        {
                            inserter=caveEdges.insert(Edge(aTri,j,false));
                            Assert(inserter.second);
                        }
                    }
                }
                break;
            }else
            {
                caveElements.clear();
                for(HashUIntSet::iterator eiter=eleSet.begin();eiter!=eleSet.end();++eiter)
                {
                    caveElements.push_back(*eiter);
                }
            }
        }
    }
    return caveElements.size()>0;
}

/** \breif Delaunay插点操作
* 首先判断待插点与三角形的位置关系\n
* 然后根据判断结果执行插点操作
*/
int Delaunay2D::insert(UInt nodeHd,bool force)
{
    Link whiTri;
    whiTri.first=_lastTri;
    int rt = locate(_nodesArray[nodeHd],whiTri); ///< 获取包含待插点的三角形及该点与三角形的位置关系
    if(rt == OUTDOMAIN)
    {
        return ErrorCodes::eInValidInput;
    }
    return insert(whiTri.first,nodeHd,force);
}

/** \brief Delaunay插点操作
    \param triSt 待插点所在的三角形单元的编号
	\param ndHd 待插点的编号
*/
int Delaunay2D::insert(UInt triSt,UInt ndHd,bool force)
{
    EdgeSet caveEdges;
    Array1D<UInt> caveElements;
    if(!formCave(ndHd,triSt,force,caveEdges,caveElements)) 
        return ErrorCodes::eMeshingFailed;

    UInt tri;int whi;bool rev;

    EdgeSet adjEdges;
    EdgeSet::iterator iter=caveEdges.begin();
    UInt nodes[3];
    nodes[2]=ndHd; ///< 待插点作为三角形单元的第三点

    StdPair<EdgeSet::iterator,bool> inserter;
    for(;iter!=caveEdges.end();++iter)
    {
        iter->link(tri,whi,rev);
        int rt=-1;
        if(rev)
        {
            nodes[0]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];
            nodes[1]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]];
        }else
        {
            nodes[0]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]];
            nodes[1]=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];
        }
        UInt triNew=addTri(nodes[0],nodes[1],nodes[2]);
        setNodeLink(nodes[0],Link(triNew,0)); ///< 设置各点的连接三角形单元(first: 三角形单元的局部编号；second: 点在该三角形中的局部编号)
        setNodeLink(nodes[1],Link(triNew,1));
        setNodeLink(nodes[2],Link(triNew,2));
        if(rev)
        {
            setTriLink(triNew,2,Link(tri,whi)); ///< 待插点所对的边(及内核区域的边界边) 如果它不是几何边界边，则需要重建拓扑关系
        }

        for(int i=0;i<2;i++)
        {
            inserter=adjEdges.insert(Edge(triNew,i,false));
            if(!inserter.second)
            {
                inserter.first->link(tri,whi,rev);
                setTriLink(triNew,i,Link(tri,whi));
            }
        }
    }

    UInt num=caveElements.size();
    for(UInt i=0;i<num;i++)
    {
        eraseTri(caveElements[i]);
    }
    return ErrorCodes::eOk;
}

struct ElementVisitorRec : public Delaunay2D::LinkVisitor
{
    ElementVisitorRec(Delaunay2D* meshPtr,UInt* edge)
        :_meshPtr(meshPtr),_edge(edge)
    {
    }

    bool visit(const Link& lnk)
    {
        _rt=_meshPtr->recoverySwap(_edge,lnk);//try to swap first intersected edge
        if(_rt==Delaunay2D::FINISHED) return false;
        if(_rt==Delaunay2D::CONTINUE) return true;
        if(_rt==Delaunay2D::SWAPED) return false;
        if(_rt==Delaunay2D::ERRO_INPUT) return false;
        Link cur=lnk;

        int cnt=0;
        while(_rt==Delaunay2D::NEXT_ADJ && cnt++<1000)
        {
            cur=_meshPtr->adjacentTri(cur);
            for(int i=1;i<3;i++)
            {
                Link lnki(cur.first,(cur.second+i)%3);
                _rt=_meshPtr->recoverySwap(_edge,lnki);
                if(_rt==Delaunay2D::FINISHED) return false;
                if(_rt==Delaunay2D::SWAPED) return false;
                if(_rt==Delaunay2D::NEXT_ADJ)
                {
                    cur=lnki; 
                    break;
                }
            }
        }
        Assert(0);
        return true;
    }

    Delaunay2D*  _meshPtr;
    UInt*        _edge;
    int          _rt;
};

/** \brief 边界恢复操作
*/
int  Delaunay2D::recovery(UInt nd1,UInt nd2)
{
    if(0)
    {
        std::ofstream out("delaunay_debug_rec.nas");
        exportBlk(out,2,true);
    }
    UInt i,sz=_workingEleArray.size()+100;
    {
        UInt edge[2]={nd1,nd2};
        ElementVisitorRec vi(this,edge);
        for(i=0;i<sz;i++)
        {
            visitNodeBall(edge[0],vi);
            if(vi._rt==FINISHED)
            {
                return ErrorCodes::eOk;
            }else if(vi._rt==ERRO_INPUT)
            {
                return ErrorCodes::eInValidInput;
            }
        }
    }

    {
        UInt edge[2]={nd2,nd1};
        ElementVisitorRec vi(this,edge);
        for(i=0;i<sz;i++)
        {
            visitNodeBall(edge[0],vi);
            if(vi._rt==FINISHED)
            {
                return ErrorCodes::eOk;
            }else if(vi._rt==ERRO_INPUT)
            {
                return ErrorCodes::eInValidInput;
            }
        }
    }
    return ErrorCodes::eMeshingFailed;
}

template<typename PT>
Link MesherBasic<PT>::adjacentTri(const Link& lnk)const
{
    UInt ti=lnk.first;
    int whi=lnk.second;
    return _adjElementsArray[ti][whi];
}

int Delaunay2D::recoverySwap(UInt edge[2],const Link& lnk)
{
    UInt ti=lnk.first;
    int whi=lnk.second;
    UInt nd0=_workingEleArray[ti][TRI_EDGE_VERTEX[whi][0]];
    if(nd0==edge[1]) return FINISHED;
    if(nd0==edge[0]) return CONTINUE;
    UInt nd1=_workingEleArray[ti][TRI_EDGE_VERTEX[whi][1]];
    if(nd1==edge[1]) return FINISHED;
    if(nd1==edge[0]) return CONTINUE;

    int rt0=validTri(edge[0],edge[1],nd0);
    int rt1=validTri(edge[0],edge[1],nd1);
    if(rt0*rt1>0) return CONTINUE;

    int rt2=validTri(nd0,nd1,edge[0]);
    int rt3=validTri(nd0,nd1,edge[1]);
    if(rt2*rt3>0) return CONTINUE;
    SwapCheckerRecover swaper(this);
    if(0==swapIt(lnk,true,&swaper,NULL,NULL))
    {
        return SWAPED;
    }
    if(_boundaryEdgeSet.find(Edge2(nd0,nd1))!=_boundaryEdgeSet.end())
    {
        char buffer[1024];
        sprintf(buffer,"Edge(%d,%d) is intersected with edge(%d,%d)",nd0+1,nd1+1,edge[0]+1,edge[1]+1);
        _msgFun(0,buffer);
        return ERRO_INPUT;
    }
    return NEXT_ADJ;

}

template<typename PT>
bool   MesherBasic<PT>::isSuperNode(UInt ni)const
{
    if(ni>=_inputNodeSize)
    {
        return ni<_inputNodeSize+3;
    }
    return false;
}

/** \brief 获取待翻转边的两个邻接三角形及两三角形的四个顶点及待翻转边在两邻接三角形中的四条邻接边
*/
bool   Delaunay2D::getSwapData(const Link& lnk,UInt tris[2],
                               UInt nd[4],
                               Link lnks[4],
                               bool& hasSuper)const
{
    tris[0]=lnk.first;
    int whi=lnk.second;
    if(isCheckedTri(tris[0],ERASED)) return false;

    hasSuper=false;

    const Link& adj=_adjElementsArray[tris[0]][whi];

    tris[1]=adj.first;
    int whiadj=adj.second;
    if(isCheckedTri(tris[1],ERASED)) return false;

    nd[0]=_workingEleArray[tris[0]][TRI_EDGE_VERTEX[whi][0]];

    if(!hasSuper) hasSuper=isSuperNode(nd[0]);

    lnks[0]=_adjElementsArray[tris[0]][TRI_EDGE_VERTEX[whi][0]];

    nd[1]=_workingEleArray[tris[0]][TRI_EDGE_VERTEX[whi][1]];
    if(!hasSuper) hasSuper=isSuperNode(nd[1]);
    lnks[1]=_adjElementsArray[tris[0]][TRI_EDGE_VERTEX[whi][1]];

    nd[2]=_workingEleArray[tris[0]][whi];

    if(!hasSuper) hasSuper=isSuperNode(nd[2]);
    lnks[2]=_adjElementsArray[tris[1]][TRI_EDGE_VERTEX[whiadj][0]];

    nd[3]=_workingEleArray[tris[1]][whiadj];
    if(!hasSuper) hasSuper=isSuperNode(nd[3]);
    lnks[3]=_adjElementsArray[tris[1]][TRI_EDGE_VERTEX[whiadj][1]];
    return true;
}


int  Delaunay2D:: swapIt(const Link& elnk,bool doRecover,SwapCheck* checkPtr,Link* newEdge,Link* outlink)
{
    UInt tris[2];
    UInt nd[4];
    Link lnk[4];
    bool hasSuper;
    if(!getSwapData(elnk,tris,nd,lnk,hasSuper)) return -1;
	if (isCheckedTri(tris[0],FROZEN) || isCheckedTri(tris[1],FROZEN))
		return -1;
    if((!doRecover) && hasSuper)return 1;
    if(checkPtr)
    {
        if(!checkPtr->pass(nd,tris))///< 检查边翻转后生成的两个三角形是否有效，如无效，则不能翻转
        {
            return 2;
        }
    }

    _workingEleArray[tris[0]][0]=nd[0];
    _workingEleArray[tris[0]][1]=nd[3];
    _workingEleArray[tris[0]][2]=nd[2];

    _workingEleArray[tris[1]][0]=nd[1];
    _workingEleArray[tris[1]][1]=nd[2];
    _workingEleArray[tris[1]][2]=nd[3];
    //setNodeLink(nd[0],Link(tris[0],0));

    setNodeLink(nd[0],Link(tris[0],0));
    setNodeLink(nd[3],Link(tris[0],1));
    setNodeLink(nd[2],Link(tris[0],2));

    setNodeLink(nd[1],Link(tris[1],0));
    setNodeLink(nd[2],Link(tris[1],1));
    setNodeLink(nd[3],Link(tris[1],2));

    //===
    setTriLink(tris[0],0,Link(tris[1],0));
    setTriLink(tris[0],1,lnk[1]);
    setTriLink(tris[0],2,lnk[2]);

    setTriLink(tris[1],1,lnk[3]);
    setTriLink(tris[1],2,lnk[0]);
    if(outlink)
    {
        outlink[0]=_adjElementsArray[tris[0]][1];
        outlink[1]=_adjElementsArray[tris[0]][2];
        outlink[2]=_adjElementsArray[tris[1]][1];
        outlink[3]=_adjElementsArray[tris[1]][2];
    }
    if(newEdge)*newEdge=Link(tris[0],0);
    return 0;
}

bool  SwapCheckerRecover::pass(const UInt nds[4],const UInt tris[2])
{
    int det=_meshPtr->validTri(nds[0],nds[3],nds[2]);
    if(det<=0) return false;
    det=_meshPtr->validTri(nds[1],nds[2],nds[3]);
    if(det<=0.0) return false;
    return true;
}

/** \brief 首先判断边是否允许翻转(即翻转后不出现负单元), 
*   如果允许翻转，则继续判断翻转操作是否可以使三角形最大角减小
    \return true 边允许翻转且翻转后三角形最大角减小
	\return false 边不允许翻转或翻转后三角形最大角不会减小
*/
bool  MinMaxChecker2D::pass(const UInt nds[4],const UInt tris[2])
{
    if(!SwapCheckerRecover::pass(nds,tris))
    {
        return false;
    }

    return  MinMaxSwap2D(_meshPtr->nodeCoord(nds[0]),
        _meshPtr->nodeCoord(nds[1]),
        _meshPtr->nodeCoord(nds[2]),
        _meshPtr->nodeCoord(nds[3]));
}

template<typename PT>
bool MesherBasic<PT>::visitNodeBall(UInt ndi,LinkVisitor& vi)
{
    Stack<Link>     frontStack;
    Array1D<UInt>   visitedLnks;

    frontStack.push(_linkNodeArray[ndi]);

    while(!frontStack.empty())
    {
        Link lnk=frontStack.top();
        frontStack.pop();
        UInt fhandle=lnk.first;
        if(isCheckedTri(fhandle,ERASED)) continue;
        if(isCheckedTri(fhandle,RESERVE0))continue;
        int   whi=lnk.second;
        if(!vi.visit(lnk))
        {
            break;
        }
        checkinTri(fhandle,RESERVE0);
        visitedLnks.push_back(fhandle);
        for(int i=0;i<2;i++)
        {
            int ei=TRI_VERTEX_EDGE[whi][i];
            UInt fhandlei= _adjElementsArray[fhandle][ei].first;
            if(fhandlei==InValidEntity) continue;
            if(isCheckedTri(fhandlei,ERASED)) continue;
            if(isCheckedTri(fhandlei,RESERVE0))continue;
            ei=_adjElementsArray[fhandle][ei].second;
            if(ndi==_workingEleArray[fhandlei][TRI_EDGE_VERTEX[ei][0]])
            {
                frontStack.push(Link(fhandlei,(char)(TRI_EDGE_VERTEX[ei][0])));
            }else if(ndi==_workingEleArray[fhandlei][TRI_EDGE_VERTEX[ei][1]])
            {
                frontStack.push(Link(fhandlei,(char)(TRI_EDGE_VERTEX[ei][1])));
            }else
            {
                Assert(0);
            }
        }
    }

    UInt i,sz=visitedLnks.size();
    for(i=0;i<sz;i++)
    {
        checkoutTri(visitedLnks[i],RESERVE0);
    }
    return true;
}

template<typename PT>
bool MesherBasic<PT>::visitNodeBall(UInt ndi,Array1D<Link>& lnks)
{
    struct ElementVisitor : public LinkVisitor
    {
        bool visit(const Link& lnk)
        {
            links->push_back(lnk);
            return true;
        }
        Array1D<Link>*  links;
    };

    ElementVisitor vi;
    vi.links=&lnks;
    return visitNodeBall(ndi,vi);
}

template<typename PT>
bool MesherBasic<PT>::getAdjacentEdge(UInt ndi,LinkVisitor& vi)
{
	Stack<Link>     frontStack;
	Array1D<UInt>   visitedLnks;

	frontStack.push(_linkNodeArray[ndi]); ///< 把节点的链接面片压入栈

	while(!frontStack.empty())
	{
		Link lnk=frontStack.top();///< 存储了邻接面片的局部编号及点在该面片中的局部编号
		frontStack.pop();
		UInt fhandle=lnk.first;
		if(isCheckedTri(fhandle,ERASED)) continue;
		if(isCheckedTri(fhandle,RESERVE0))continue;
		int whi=lnk.second; ///< 点ndi在当前面片中的局部编号
// 		if(!vi.visit(lnk))
// 		{
// 			break;
// 		}
		checkinTri(fhandle,RESERVE0);
		visitedLnks.push_back(fhandle);
		for(int i=0;i<2;i++) ///< 遍历该点在当前单元中的两条邻接边
		{
			int ei=TRI_VERTEX_EDGE[whi][i];
			Link edge_link; ///< 顶点的邻接边
			edge_link.first = fhandle;
			edge_link.second = ei;

			UInt fhandlei= _adjElementsArray[fhandle][ei].first;
			if(fhandlei==InValidEntity) continue;
			if(isCheckedTri(fhandlei,ERASED)) continue;
			if(isCheckedTri(fhandlei,RESERVE0))continue;
			if(!vi.visit(edge_link))
			{
				break;
			}
			ei=_adjElementsArray[fhandle][ei].second;
			if(ndi==_workingEleArray[fhandlei][TRI_EDGE_VERTEX[ei][0]])
			{
				frontStack.push(Link(fhandlei,(char)(TRI_EDGE_VERTEX[ei][0])));
			}else if(ndi==_workingEleArray[fhandlei][TRI_EDGE_VERTEX[ei][1]])
			{
				frontStack.push(Link(fhandlei,(char)(TRI_EDGE_VERTEX[ei][1])));
			}else
			{
				Assert(0);
			}
		}
	}

	UInt i,sz=visitedLnks.size();
	for(i=0;i<sz;i++)
	{
		checkoutTri(visitedLnks[i],RESERVE0);
	}
	return true;
}

template<typename PT>
bool MesherBasic<PT>::getAdjacentEdge(UInt ndi,Array1D<Link>& lnks)
{
	struct EdgeVisitor : public LinkVisitor
	{
		bool visit(const Link& lnk)
		{
			links->push_back(lnk);
			return true;
		}
		Array1D<Link>*  links;
	};

	EdgeVisitor vi;
	vi.links=&lnks;
	return getAdjacentEdge(ndi,vi);
}

struct ElementVisitorMark : public Delaunay2D::LinkVisitor
{
    ElementVisitorMark(Delaunay2D* meshPtr)
        :_meshPtr(meshPtr)
    {
    }

    bool visit(const Link& lnk)
    {
        _meshPtr->markSuperNode(lnk); 
        return true;
    }
    Delaunay2D*  _meshPtr;
};

void Delaunay2D::markSuperNode(const Link& lnk)
{
    UInt tri=lnk.first;
    checkinTri(tri,OUT_DOMAIN);

    int whi=lnk.second;
    Edge2 k0(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]]);
    if(_boundaryEdgeSet.find(k0)!=_boundaryEdgeSet.end())
    {
        checkinTri(_adjElementsArray[tri][whi].first,IN_DOMAIN);
    }
}

void Delaunay2D::markAdjacent(UInt tri, int mk,UInt& markedEleNum)
{
    Stack<UInt> eleStack;
    eleStack.push(tri);
    while(!eleStack.empty())
    {
        UInt thi=eleStack.top();
        eleStack.pop();
        for(int i=0;i<3;i++)
        {
            UInt lnkj=_adjElementsArray[thi][i].first;
            if(lnkj==InValidEntity) continue;
            if(isCheckedTri(lnkj,ERASED)) continue;
            if(isCheckedTri(lnkj,OUT_DOMAIN)) continue;
            if(isCheckedTri(lnkj,IN_DOMAIN)) continue;
            Edge2 k0(_workingEleArray[thi][TRI_EDGE_VERTEX[i][0]],_workingEleArray[thi][TRI_EDGE_VERTEX[i][1]]);
            if(_boundaryEdgeSet.find(k0)!=_boundaryEdgeSet.end()) continue;
            checkinTri(lnkj,mk);
            eleStack.push(lnkj);
            markedEleNum++;
        }
    }
}

UInt  Delaunay2D:: markedTriNum()const
{
    UInt i,sz=_workingEleArray.size();
    UInt triNumMarked=0;
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,ERASED)) continue;
        if(isCheckedTri(i,OUT_DOMAIN) || isCheckedTri(i,IN_DOMAIN))
        {
            triNumMarked++;
        }
    }
    return triNumMarked;
}

void  Delaunay2D::markout()
{
    ElementVisitorMark mkSu(this);
    visitNodeBall(_inputNodeSize,mkSu);
    visitNodeBall(_inputNodeSize+1,mkSu);
    visitNodeBall(_inputNodeSize+2,mkSu);

    if(0)
    {
        std::ofstream out("delaunay_debug_mark.nas");
        exportBlk(out,2,true);
    }

    UInt i,j,triNum=0,sz=_workingEleArray.size();
    UInt triNumMarked=0;
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,ERASED)) continue;
        triNum++;
        if(isCheckedTri(i,OUT_DOMAIN)) 
        {
            markAdjacent(i,OUT_DOMAIN,triNumMarked);
        }else if(isCheckedTri(i,IN_DOMAIN)) 
        {
            markAdjacent(i,IN_DOMAIN,triNumMarked);
        }
    }

    if(0)
    {
        std::ofstream out("delaunay_debug_mark.nas");
        exportBlk(out,2,true);
    }
    triNumMarked=markedTriNum();
    if(triNumMarked>=triNum) return;


    for(i=0;i<sz;i++)
    {
        Array1D<UInt> newOuts;
        for(j=0;j<sz;j++)
        {
            if(isCheckedTri(j,ERASED)) continue;
            if(isCheckedTri(j,OUT_DOMAIN)) continue;
            if(isCheckedTri(j,IN_DOMAIN)) continue;
            for(int k=0;k<3;k++)
            {
                UInt lnkj=_adjElementsArray[j][k].first;
                if(lnkj==InValidEntity) continue;
                if(isCheckedTri(lnkj,ERASED)) continue;
                if(isCheckedTri(lnkj,IN_DOMAIN)) 
                {
                    triNumMarked++;
                    checkinTri(j,OUT_DOMAIN);
                    newOuts.push_back(j);
                    break;
                }
            }
        }
        if(0)
        {
            std::ofstream out("delaunay_debug_mark.nas");
            exportBlk(out,2,true);
        }
        if(triNumMarked>=triNum) 
        {
            triNumMarked=markedTriNum();
            if(triNumMarked>=triNum) return;
        }
        UInt sz2=newOuts.size();
        for(j=0;j<sz2;j++)
        {
            markAdjacent(newOuts[j],OUT_DOMAIN,triNumMarked);
        }
        if(0)
        {
            std::ofstream out("delaunay_debug_mark.nas");
            exportBlk(out,2,true);
        }
        if(triNumMarked>=triNum) 
        {
            triNumMarked=markedTriNum();
            if(triNumMarked>=triNum) return;
        }

        Array1D<UInt> newIn;
        for(j=0;j<sz;j++)
        {
            if(isCheckedTri(j,ERASED)) continue;
            if(isCheckedTri(j,OUT_DOMAIN)) continue;
            if(isCheckedTri(j,IN_DOMAIN)) continue;
            for(int k=0;k<3;k++)
            {
                UInt lnkj=_adjElementsArray[j][k].first;
                if(lnkj==InValidEntity) continue;
                if(isCheckedTri(lnkj,ERASED)) continue;
                if(isCheckedTri(lnkj,OUT_DOMAIN)) 
                {
                    triNumMarked++;
                    checkinTri(j,IN_DOMAIN);
                    newIn.push_back(j);
                    break;
                }
            }
        }
        if(0)
        {
            std::ofstream out("delaunay_debug_mark.nas");
            exportBlk(out,2,true);
        }
        if(triNumMarked>=triNum) 
        {
            triNumMarked=markedTriNum();
            if(triNumMarked>=triNum) return;
        }
        sz2=newOuts.size();
        for(j=0;j<sz2;j++)
        {
            markAdjacent(newOuts[j],IN_DOMAIN,triNumMarked);
        }
        if(0)
        {
            std::ofstream out("delaunay_debug_mark.nas");
            exportBlk(out,2,true);
        }
        if(triNumMarked>=triNum) 
        {
            triNumMarked=markedTriNum();
            if(triNumMarked>=triNum) return;
        }
    }
}

int   Delaunay2D::optConnection()
{
    MinMaxChecker2D fun(this);
    localReconnect(&fun);
    if(0)
    {
        std::ofstream out("delaunay_debug_opt.nas");
        exportBlk(out,2,true);
    } 
    return ErrorCodes::eOk;
}

void  Delaunay2D::localReconnect(SwapCheck* fun)
{
    for(int iter=0;iter<100;iter++)
    {
        if(0)
        {
            std::ofstream out("delaunay_debug_opt.nas");
            exportBlk(out,2,true);
        }
        UInt i=0,sz=_workingEleArray.size();
        int swapEdgeCnt=0;
        for(i=0;i<sz;i++)
        {
            if(isCheckedTri(i,ERASED)) continue;
            if(isCheckedTri(i,OUT_DOMAIN)) continue;
            if(isCheckedTri(i,FROZEN)) continue;
            for(int j=0;j<3;j++)
            {
                UInt adji=_adjElementsArray[i][j].first;
                if(!isCheckedTri(adji,IN_DOMAIN)) continue;
                Link ei(i,j);
                if(0==swapIt(ei,false,fun,NULL,NULL))
                {
                    checkoutTri(adji,FROZEN);
                    checkoutTri(i,FROZEN);
                    swapEdgeCnt++;
                    break;
                }
            }
        }
        if(swapEdgeCnt==0) 
        {
            break;
        }
    }
}

void  print_node_long_format_nas(int id,Real x,Real y,Real z,char* strBuffer)
{
    char bufx[32],bufy[32],bufz[32];
    sprintf(bufx,"%16.16lf",x);bufx[15]='\0';
    sprintf(bufy,"%16.16lf",y);bufy[15]='\0';
    sprintf(bufz,"%16.16lf",z);bufz[15]='\0';
    sprintf(strBuffer,"%-16s%8d%-16s%-16s%-16s%-1s\n%-8s%-16s\n", "GRID*", id, "", bufx,bufy, "*" ,"*",bufz);
}

template<typename PT>
void  MesherBasic<PT>::exportBlk(std::ostream& out,int dim,bool all)
{
    char buffer[256];
    UInt i, sz=_nodesArray.size();

    out << "$$\n$$  GRID Data\n$$\n";
    Real coord[3];
    for(i = 0; i < sz; i++)
    {
        if(dim == 3 && isSuperNode(i)) continue;
        getNodeCoord(i,dim,coord);
        print_node_long_format_nas(i+1,coord[0],coord[1],coord[2],buffer);
        out << buffer;
    }

    out<<"$\n";

    UInt cnt = 1;
    sz = _workingEleArray.size();

    sprintf(buffer,"CTRIA3  ");
    out<<"$  CTRIA3 Data\n$\n";

    for(i=0;i<sz;i++)
    {
        int prop=309;
        if(isCheckedTri(i,ERASED)) continue;
        if(!all && isCheckedTri(i,OUT_DOMAIN)) continue;
        if(isCheckedTri(i,OUT_DOMAIN)) 
        {
            prop=300;
        }else if(isCheckedTri(i,IN_DOMAIN)) 
        {
            prop=301;
            if(isCheckedTri(i,FROZEN)) 
            {
                prop=308;
            }
        }

        int gp=getTriGroupId(i);
        if(gp>0)prop=gp;

        out<<buffer<<std::setw(8)<<cnt++<<std::setw(8)<<prop;
        int ncount = 0;
        int jfield = 0;
        // first line
        for (jfield = 4; jfield <=9; jfield++)
        {
            ncount++;
            if (ncount>3) break;
            int idi=_workingEleArray[i][ncount-1]+1;
            out <<std::setw(8) << idi;
        } 
        out<<"\n";
    }
}

void  Delaunay2D::exportBoundaryBlk(std::ostream& out,int dim)
{
    char buffer[256];
    UInt i, sz = _inputNodeSize;

    out << "$$\n$$  GRID Data\n$$\n";
    Real coord[3];
    for( i = 0; i < sz; i++)
    {
        getNodeCoord(i,dim,coord);
        print_node_long_format_nas(i+1, coord[0], coord[1], coord[2], buffer);
        out << buffer;
    }

    out<<"$\n";

    UInt cnt=1;


    sprintf(buffer,"PLOTEL  ");
    out<<"$  PLOTEL  Data\n$\n";
    sz = _edgesPtr->size();
    int prop = 201;
    for(i = 0; i < sz; i++)
    {
        out << buffer << std::setw(8) << cnt++;
        int ncount = 0;
        int jfield = 0;
        // first line
        for (jfield = 4; jfield <= 9; jfield++)
        {
            ncount++;
            if (ncount > 2) break;
            int idi=(*_edgesPtr)[i][ncount-1]+1;
            out <<std::setw(8) << idi;
        } 
        out<<"\n";
    }
}

template<typename PT>
void  MesherBasic<PT>::getTriangles(Array1D<Triangle>& tris,bool all)
{
    UInt i,sz=_workingEleArray.size();
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,ERASED)) continue;
        if(!all && isCheckedTri(i,OUT_DOMAIN)) continue;
        tris.push_back(_workingEleArray[i]);
    }
}


//############################################################################
int ConstrianedDealunay2D(const  Array1D<Point2D>& pts, 
                          const  Array1D<Tuple<UInt,2> >& edges,
                          Array1D<TripleT<UInt> >& tris)
{
    Delaunay2D delauy(pts,edges);
    int rt= delauy.delaunay();
    if(rt!=ErrorCodes::eOk) return rt;
    delauy.getTriangles(tris);
    return ErrorCodes::eOk;
}

int PolygonTriangulate2D(const  Array1D<Point2D>& pts, 
                         Array1D<TripleT<UInt> >& tris)
{
    UInt i,sz=pts.size();
    if(sz==3)
    {
        tris.push_back(TripleT<UInt>(0,1,2));
        return ErrorCodes::eOk;
    }

    if(sz==4)
    {
        if(Orient2d_Math(pts[0],pts[1],pts[2])<=0 ||
            Orient2d_Math(pts[0],pts[2],pts[3])<=0)
        {
            if(Orient2d_Math(pts[0],pts[1],pts[3])<=0 ||
                Orient2d_Math(pts[1],pts[2],pts[3])<=0)
            {
                return ErrorCodes::eInValidInput;
            }
            tris.push_back(TripleT<UInt>(0,1,3));
            tris.push_back(TripleT<UInt>(1,2,3));
            return ErrorCodes::eOk;
        }

        if(Orient2d_Math(pts[0],pts[1],pts[3])<=0 ||
            Orient2d_Math(pts[1],pts[2],pts[3])<=0)
        {
            tris.push_back(TripleT<UInt>(0,1,2));
            tris.push_back(TripleT<UInt>(0,2,3));
            return ErrorCodes::eOk;
        }

        if(MinMaxSwap2D(pts[2],pts[0],pts[1],pts[3]))
        {
            tris.push_back(TripleT<UInt>(0,1,3));
            tris.push_back(TripleT<UInt>(1,2,3));
        }else
        {
            tris.push_back(TripleT<UInt>(0,1,2));
            tris.push_back(TripleT<UInt>(0,2,3));
        }

        return ErrorCodes::eOk;
    }



    Array1D<Tuple<UInt,2> > edges(sz);
    Tuple<UInt,2>  aedge;
    sz-=1;
    for(i=0;i<sz;i++)
    {
        edges[i][0]=i;
        edges[i][1]=i+1;
    }
    edges[sz][0]=sz;
    edges[sz][1]=0;
    return ConstrianedDealunay2D(pts,edges,tris);
}

/** \brief 判断边交换前后，最大角是否减小
*/
bool MinMaxSwap2D(const Point2D& pt0,
                  const Point2D& pt1,
                  const Point2D& pt2,
                  const Point2D& pt3)
{
    Real wesp=0.001;
    Vector2D v01(pt0,pt1);
    Vector2D v02(pt0,pt2);
    Vector2D v03(pt0,pt3);
//     Real area012= v01[0]*v02[1]-v01[1]*v02[0];
//     if(area012<Tolerance::EPS_COMPUTER) return false;
//     Real area031= v03[0]*v01[1]-v03[1]*v01[0];
//     if(area031<Tolerance::EPS_COMPUTER) return false;

    Real area032= v03[0]*v02[1]-v03[1]*v02[0];
    if(area032<Tolerance::EPS_COMPUTER) return false;

    Vector2D v12(pt1,pt2);
    Vector2D v13(pt1,pt3);
    Real area123= v12[0]*v13[1]-v12[1]*v13[0];
    if(area123<Tolerance::EPS_COMPUTER) return false;
    v02.normalize(); v12.normalize();
    v13.normalize();v03.normalize();

    Real w012=v02.dot(v12); ///< 两向量夹角的余弦值
    Real w031=v13.dot(v03);
    Real w0=w012<w031?w012:w031; ///< 取余弦值最小(夹角最大)的最为w0

    Real w032=v03.dot(v02);
    Real w123=v12.dot(v13);
    Real w1=w032<w123?w032:w123;
    Real dw=w1-w0;
    return dw>=wesp; ///< 判断边交换前后，最大角是否减小
}

bool Normalize(Vector3D& vec)
{
    Real len=vec.length();
    if(len<Tolerance::EPS_COMPUTER) return false;
    len=1.0/len;
    vec*=len;
    return true;
}

bool     MinMaxSwap3D(const Point3D& pt0,
                      const Point3D& pt1,
                      const Point3D& pt2,
                      const Point3D& pt3)
{
    Real wesp=0.001;
    Vector3D v01(pt0,pt1);
    Vector3D v02(pt0,pt2);
    Vector3D v03(pt0,pt3);
    Vector3D v12(pt1,pt2);
    Vector3D v13(pt1,pt3);

    Real dismax012=v01.len2();
    Real dismax031=dismax012;
    Real tmp=v02.len2(); 
    if(tmp>dismax012)dismax012=tmp;
    tmp=v12.len2(); 
    if(tmp>dismax012)dismax012=tmp;

    tmp=v03.len2(); 
    if(tmp>dismax031)dismax031=tmp;
    tmp=v13.len2(); 
    if(tmp>dismax031)dismax031=tmp;
    Real sc=dismax012/dismax031;
    if(sc<1.0) sc=1.0/sc;
    if(sc>2.0) return false;

    if(!Normalize(v02)) return false;
    if(!Normalize(v12)) return false;
    if(!Normalize(v13)) return false;
    if(!Normalize(v03)) return false;
    if(!Normalize(v12)) return false;
    if(!Normalize(v13)) return false;

    Real w012=v02.dot(v12);
    Real w031=v13.dot(v03);
    Real w0=w012<w031?w012:w031;

    Real w032=v03.dot(v02);
    Real w123=v12.dot(v13);
    Real w1=w032<w123?w032:w123;
    Real dw=w1-w0;
    return dw>=wesp;
}

/** \brief 计算三角形的外接圆圆心坐标及半径
    \param[in] pt1~pt3 三角形的三个顶点
	\param[out] cen 三角形的外接圆圆心
	\param[out] rad 三角形外接圆的半径
*/
int Circumcircle2D(const Point2D& pt1,
                   const Point2D& pt2,
                   const Point2D& pt3,
                   Point2D& cen,
                   Real& rad)
{
    Vector2D vec1(pt1,pt2);
    Vector2D vec2(pt2,pt3);
    Point2D pta=(pt1+pt2)*.5;
    Point2D ptb=pta+vec1.perpVector();
    Point2D ptc=(pt2+pt3)*.5;
    Point2D ptd=ptc+vec2.perpVector();
    Real  up=(pta[1]-ptc[1])*(ptd[0]-ptc[0])-(pta[0]-ptc[0])*(ptd[1]-ptc[1]);
    Real  down=(ptb[0]-pta[0])*(ptd[1]-ptc[1])-(ptb[1]-pta[1])*(ptd[0]-ptc[0]);
    if(fabs(down)<Tolerance::EPS) return ErrorCodes::eDegeneratedGeo;

    Real  r=up/down;
    cen= (pta+Vector2D(pta,ptb)*r);
    rad=cen.distanceTo(pt1);
    return ErrorCodes::eOk;
}

/** \brief 三角形面积与边长平方和之比
* 正三角形的质量系数为1
* 等腰直角三角形的质量系数为 2 * sqrt(3)
*/
Real  QuaTriangle(Real a,Real b,Real c)
{
    Real min=a<b?a:b;
    min=min<c?min:c;
    min*=0.1;
    min=1.0/min;
    a*=min;b*=min;c*=min;
    //Real ditor=3.4641016151377545870548926830117;///< 2 * sqrt(3)
	Real ditor = 2 * sqrt(3.0);///< modified by ZhangJun 2016.06.07
    Real bottom=a*a+b*b+c*c;
    Real det=(a+b+c)*(a+b-c)*(a+c-b)*(b+c-a);
    if(det<Tolerance::EPS) return -Tolerance::MAX_REAL;
    Real up=sqrt(det)*0.5;
    return (ditor*up/bottom);
}

/*#########################################################################################

DelaunayAFT2D

//###########################################################################################*/

/** \brief 计算线段在黎曼度量下的长度
*/
Real  NormalLength(Real len,Real h0,Real h1)
{
    if(fabs(h0-h1)<Tolerance::EPS)
    {
        return len/h0;
    }
    return len*log(h0/h1)/(h0-h1);
}
static const Real LENGTH_DET=sqrt(2.0);
static const Real LENGTHINTERVAL[4]={LENGTH_DET*0.5,LENGTH_DET,LENGTH_DET*0.25,0.01};



class DelaunayAFT2D : public Delaunay2D
{
public:
    enum{E_TRI=0,R_TRI};
    DelaunayAFT2D(const  Array1D<Point2D>& pts,const EdgeArray& edges,Fun_Message fun=NULL)
        :Delaunay2D(pts,edges,fun)
    {
        //_gr=1.05;
		_gr = 1.0;
        _maxTriSize=-1.0;
        _gdsc=0.618;
    }

    DelaunayAFT2D(Input* inputPtr,Fun_Message fun=NULL);
    virtual int  doMesh();

protected:
    virtual Real distance(const Point2D&,const Point2D&)const;
    virtual Real distance(UInt nd1,UInt nd2)const;
    virtual UInt createFieldPoint(const Point2D& pt);
    virtual void localReconnect();
    virtual void getSmoothData(UInt ti,Real& qua,Real &wi,Point2D& ptc);
    virtual void setNodePos(UInt ni,const Point2D& newPos);
    virtual void getTriQua(UInt ti,Real& qua);
    virtual void outputDebugMesh(const char*);
    virtual Real distance(UInt nd1,const Point2D& pt2)const;
    virtual Real normalDistanceEdge(UInt e0,UInt e1,const Point2D& pt2,Real sz)const;
    virtual int  preAFT() ;
    virtual Real anglecosine(UInt, UInt,UInt,UInt)const;
	virtual Real angle_evaluate(UInt, UInt,UInt,UInt)const;
    struct FrontNode
    {
        UInt ni;
        mutable Link edge[2];
        HashValue hash()const{return ni;}
        bool operator==(const FrontNode& src)const{return ni==src.ni;}
        FrontNode(UInt nii=InValidEntity)
            :ni(nii)
        {
            edge[0]=edge[1]=Link(InValidEntity,-1);
        }
    };
    typedef HashSetDefault<FrontNode>                       FrontNodeSet;
    typedef StdPair<Point2D,Real> NodeSizePair;
    virtual int     getRTriPoint(const FrontNode&,NodeSizePair nodert[2]);
    virtual void    discardClosetNodes(Array1D<NodeSizePair>&){};

public:
    bool  isCheckedNode(UInt ni,int i)const{return _nodeMark[ni].test(i);}
    void  checkinNode(UInt ni,int i){_nodeMark[ni].set(i);}
    void  checkoutNode(UInt ni,int i){_nodeMark[ni].clear(i);}
    void  setupNodeSize();
    int   advancing();
    UInt  insertDirect(const Link& triHandle,const Point2D& pt,Link newTri[4]);
    void  markTriangle(Array1D<UInt>& badEles);
    void  markETriangle(Array1D<UInt>& badEles);
    void  markRTriangle(Array1D<UInt>& badEles);
    Real  normalLength(UInt nd1,UInt nd2)const;
    Real  normalLength2(UInt nd1,const Point2D& pt,Real sz1,Real sz2)const;
    Real  normalLength(const Point2D&,const Point2D&,Real)const;
    Real  normalLength(const Point2D&,const Point2D&,Real,Real)const;
    bool  get3rdPoint(UInt tri,Point2D& ptrt,Real& nodesizert);
    bool  get3rdPoint_mid(UInt,Point2D&,Real&);
    bool  canBeInserted(const Point2D&,Real,const Link&,Array1D<UInt>&);
    bool  canBeInserted(const Point2D& pt,Link& whiTri,int& loc,Real szNode);
    void  doSmooth(int maxIter);
    bool  smooth(UInt ni);
    void  meshOpt();
    void  erase34();
    void  erase3(UInt ndi,Array1D<Link>& lnks);
    void  erase4(UInt ndi,Array1D<Link>& lnks);
    Real  sizeInterpolate(UInt tri,const Point2D& pos)const;
	bool  isLongest(const Link& edge); ///< 判断edge在其所在三角形中是否为最长边
	Real  length(const Point2D& pos1, const Point2D& pos2);

    void    generateNodes(const Array1D<UInt>&,Array1D<NodeSizePair>&,bool at_middle=false);
    void    trans2quad();
    Real    qualityQuad(const Link&)const;
    void    exportQuadBlk(std::ostream& out,int dim);
protected:
    Array1D<Real>     _nodeSize; ///< 节点尺寸数组
    QuadArray         _quadArray;
    Real              _gr,_gdsc,_maxTriSize,_minTriSize;
    Real              _avgLen; ///< 节点尺寸平均值
    int               _triType;
};

DelaunayAFT2D::DelaunayAFT2D(Input* inputPtr,Fun_Message fun)
:Delaunay2D(inputPtr,fun)
{
    _triType=E_TRI;
    //_triType=R_TRI;
    _gr=1.05;
    _maxTriSize=-1.0;
    _minTriSize=DBL_MAX;
    _gdsc=0.618;
    if(inputPtr!=NULL)
    {
        Real gr=inputPtr->growth_ratio();
        if(gr>=1.0) _gr=gr;

        Real maxS,minS;
        inputPtr->max_min_size(maxS,minS);
        if(maxS>0)  _maxTriSize=maxS;
        if(minS<_minTriSize)_minTriSize=minS;
    }
}

/** \brief AFT预处理
* 设置网格各节点处的尺寸值\n
* 实际上是各边界节点处的尺寸值
* 节点处的尺寸之取该点最短邻接边的长度
*/
int  DelaunayAFT2D::preAFT() 
{
    setupNodeSize();
    return ErrorCodes::eOk;
}

int  DelaunayAFT2D::doMesh()
{
    int  rt= delaunay(); ///< 完成边界点集的约束Delaunay三角化
    if(rt!=ErrorCodes::eOk) return rt;
    rt=preAFT(); ///< 设置各节点的尺寸值为该点最短邻接边的长度
    if(rt!=ErrorCodes::eOk) return rt;
    rt=advancing(); ///< 前沿推进
    if(rt!=ErrorCodes::eOk) return rt;
    if(R_TRI==_triType)
    {
        _triType=E_TRI;
        rt=advancing();
        _triType=R_TRI;
        if(rt!=ErrorCodes::eOk) return rt;
        trans2quad();
    }

    //meshOpt(); ///< 网格优化
    _msgFun(1,"Mesh done.");
    return rt;
}



void  DelaunayAFT2D::generateNodes(const Array1D<UInt>& badEle,
                                   Array1D<NodeSizePair>& newPoints,bool middle)
{
    if(_triType==R_TRI)
    {
        UInt i,j,sz = badEle.size();
        EdgeSet edgeSet; ///< 用于存储所有坏单元中被冻结的边或边界边
        for(i = 0; i < sz; i++) ///< 遍历所有坏单元，收集所有边界边和被冻结的边存储在edgeSet中
        {
            UInt tri = badEle[i]; ///< 当前单元
            for(j = 0; j < 3; j++)
            {
                UInt adj=_adjElementsArray[tri][j].first; ///< 遍历该单元的三个邻接单元
                if(adj!=InValidEntity && (!isCheckedTri(adj,ERASED))) ///< 邻接单元有效且未被标记删除
                {
                    if(isCheckedTri(adj,FROZEN)) ///< 如果邻接单元已被冻结
                    {
                        edgeSet.insert(Edge(tri,j)); ///< 存储当前单元中相应的边
                    }else
                    {
                        Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[j][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[j][1]]);
                        if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei)) ///< 如果该边是边界边
                        {
                            edgeSet.insert(Edge(tri,j)); ///< 存储当前单元中相应的边
                        }
                    }
                }else ///< 邻接单元无效或已被标记删除
                {
                    Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[j][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[j][1]]);
                    if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei))///< 如果该边是边界边
                    {
                        edgeSet.insert(Edge(tri,j));///< 存储当前单元中相应的边
                    }
                }
            }
        }

        FrontNodeSet nodeSet;
        StdPair<FrontNodeSet::iterator,bool> inserter;
        EdgeSet::iterator eiter=edgeSet.begin();
        for(;eiter!=edgeSet.end();++eiter) ///< 遍历上一步收集到的边
        {
            UInt tri;int whi;bool rv;
            eiter->link(tri,whi,rv);
            for(int i = 0; i < 2; i++)
            {
                inserter = nodeSet.insert(FrontNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][i]]));
                if((inserter.first->edge[0]).first==InValidEntity)
                {
                    inserter.first->edge[0]=Link(tri,whi);
                }else if((inserter.first->edge[1]).first==InValidEntity)
                {
                    inserter.first->edge[1]=Link(tri,whi);
                }
            }
        }
        newPoints.reserve(nodeSet.size());
        FrontNodeSet::iterator niter=nodeSet.begin();
        NodeSizePair node_sz[2];
        for(;niter!=nodeSet.end();++niter)
        {
            int num=getRTriPoint(*niter,node_sz);
            if(num<=0) continue;
            newPoints.push_back(node_sz[0]);
            if(num==2)newPoints.push_back(node_sz[1]);
        }
    }else
    {
        Point2D pt;
        Real szNode;
        UInt i,sz=badEle.size();
        newPoints.reserve(sz);
        if(middle)
        {
            for(i=0;i<sz;i++)
            {
                if(!get3rdPoint_mid(badEle[i],pt,szNode)) continue;
                newPoints.push_back(NodeSizePair(pt,szNode));
            }
        }else
        {
            for(i=0;i<sz;i++)
            {
                if(!get3rdPoint(badEle[i],pt,szNode)) continue;
                newPoints.push_back(NodeSizePair(pt,szNode));
            }
        }

    }

}

/** \brief 前沿推进
*/
int   DelaunayAFT2D::advancing()
{
    int layerNum=0;
    _msgFun(1,"Advancing front...\n");
    bool insert_middle=false;
	Array1D<UInt> new_node_id;
    while(1)
    {
        //_msgFun(2,"*");
        Array1D<UInt> badEle;
        markTriangle(badEle); ///< 收集在黎曼度量下不满足要求的单元(三角形单元的最长边不应大于sqrt(2.0))
        if(badEle.empty())
        {
            break;
        }
        if(0)
        {
            outputDebugMesh("marked");
        }
		for (int in = 0; in < new_node_id.size(); ++in)
		{
			checkoutNode(new_node_id[in], UNACTIVE);
		}
		new_node_id.clear();
        Array1D<NodeSizePair> newPoints;///< 节点--尺寸数据对
        generateNodes(badEle,newPoints,insert_middle);
        //discardClosetNodes(newPoints);
        if(newPoints.empty())break;
        UInt i,sz=newPoints.size();
        Link whiTrilink;
        int loc;
        int insertCnt=0;
        for(i=0;i<sz;i++)
        {
            if(0)
            {
                outputDebugMesh("insertdebug");
            }
            if(!canBeInserted(newPoints[i].first,whiTrilink,loc,newPoints[i].second)) 
				continue ;
            if(loc!=ONEDGE)
            {
                whiTrilink.second=-1;
            }
            Link newTri[4];
            UInt ndnew=insertDirect(whiTrilink,newPoints[i].first,newTri);
            _nodeSize[ndnew]=newPoints[i].second;
			checkinNode(ndnew, FRONT_NODE);///< 把心cherub的点标记为非活动前沿点
			checkinNode(ndnew, UNACTIVE);
			new_node_id.push_back(ndnew);
            insertCnt++;
// 			if (insertCnt == 2)
// 			{
// 			    int a = 0;
// 				break;
// 			}
        }
        if(insertCnt==0) 
        {
            if(insert_middle) 
                break;

            insert_middle=true;
            continue;
        }

        insert_middle=false;
        if(0)
        {
            outputDebugMesh("inserted");
        }
        //_msgFun(2,"*");
		//if (layerNum < 2)
        localReconnect();
        //_msgFun(2,"*");
        if(0)
        {
            outputDebugMesh("reconnected");
        }

        layerNum++;
// 		if (layerNum == 5)
// 			break;
    }

    if(0)
    {
        outputDebugMesh("final_0");
    }

    return ErrorCodes::eOk;
}

void  DelaunayAFT2D::meshOpt()
{
    _msgFun(1,"Meshing improving...\n");
    doSmooth(3);
    localReconnect();
    doSmooth(3);
    if(0)
    {
        outputDebugMesh("final");
    }
}

void DelaunayAFT2D::localReconnect()
{
    MinMaxChecker2D fun(this);
    Delaunay2D::localReconnect(&fun);
}

/** \brief 设置各节点处的尺寸值为其邻接边中的最短值
* 同时计算所有网格边长的平均值_avgLen
*/
void  DelaunayAFT2D::setupNodeSize()
{
    UInt i,sz=_nodesArray.size();
    _nodeSize.resize(sz);
    _nodeMark.resize(sz);

    for(i=0;i<sz;i++) ///< 初始化节点尺寸值
    {
		checkinNode(i, FRONT_NODE);///< 把边界点标记为前沿点
        _nodeSize[i]=-1.0;
    }
    _avgLen=0;
    sz=_edgesPtr->size();
    for(i=0;i<sz;i++)
    {
        UInt nd1=(*_edgesPtr)[i][0];
        UInt nd2=(*_edgesPtr)[i][1];
        Real dis=distance(nd1,nd2);
        if(_nodeSize[nd1]<0 || _nodeSize[nd1]>dis)
        {
            _nodeSize[nd1]=dis;
        }

        if(_nodeSize[nd2]<0 || _nodeSize[nd2]>dis)
        {
            _nodeSize[nd2]=dis;
        }
        _avgLen+=dis;
    }
    if(sz>0) _avgLen=_avgLen/sz;

}

/** \brief 计算两点距离
*/
Real DelaunayAFT2D::distance(const Point2D& pt1,const Point2D& pt2)const
{
    return pt1.distanceTo(pt2);
}

/** \brief 计算两点距离
* nd1 和 nd2为两点编号
*/
Real DelaunayAFT2D::distance(UInt nd1,UInt nd2)const
{
    return distance(_nodesArray[nd1],_nodesArray[nd2]);
}

bool  DelaunayAFT2D::isLongest(const Link& edge)
{
	if (-1 == edge.first)
		return false;
	UInt current_tri = edge.first;
	UInt edge_id = edge.second;
	Real currend_length = distance(_workingEleArray[current_tri][TRI_EDGE_VERTEX[edge_id][0]],
		_workingEleArray[current_tri][TRI_EDGE_VERTEX[edge_id][1]]);

	Real length_prev = distance(_workingEleArray[current_tri][TRI_EDGE_VERTEX[(edge_id + 2) % 3][0]],
		_workingEleArray[current_tri][TRI_EDGE_VERTEX[(edge_id + 2) % 3][1]]);

	Real length_next = distance(_workingEleArray[current_tri][TRI_EDGE_VERTEX[(edge_id + 1 ) % 3][0]],
		                        _workingEleArray[current_tri][TRI_EDGE_VERTEX[(edge_id + 1 ) % 3][1]]);

	if (currend_length < length_prev || currend_length < length_next)
		return false;
	else
		return true;
}

Real  DelaunayAFT2D::length(const Point2D& pos1, const Point2D& pos2)
{
	return distance(pos1,pos2);
}

/** \brief 计算线段在黎曼度量下的长度
*/
Real  DelaunayAFT2D::normalLength(UInt nd1,UInt nd2)const
{
    Real dis=distance(nd1,nd2);
    dis=NormalLength(dis,_nodeSize[nd1],_nodeSize[nd2]);
    return dis;
}

Real  DelaunayAFT2D::normalLength(const Point2D& nd1,const Point2D& pt,Real sz)const
{
    return normalLength(nd1,pt,sz,sz);
}

Real  DelaunayAFT2D::normalLength2(UInt nd1,const Point2D& pt,Real sz1,Real sz2)const
{
    Real dis=distance(nd1,pt);
    dis=NormalLength(dis,sz1,sz2);
    return dis;
}

Real DelaunayAFT2D::distance(UInt nd1,const Point2D& pt2)const
{
    return distance(_nodesArray[nd1],pt2);
}

Real DelaunayAFT2D::normalDistanceEdge(UInt e0,UInt e1,const Point2D& pt2,Real sz)const
{
    Vector2D vedge(_nodesArray[e0],_nodesArray[e1]);
    vedge.normalize();
    Vector2D vdiff(_nodesArray[e0],pt2);
    Real len=vdiff.dot(vedge);
    if(len<=0.0) return normalLength(_nodesArray[e0],pt2,_nodeSize[e0],sz);
    Real vlen=vedge.length();
    if(len>=vlen) return normalLength(_nodesArray[e1],pt2,_nodeSize[e1],sz);
    Real t=len/vlen;
    Real szt=_nodeSize[e0]*(1.0-t)+_nodeSize[e1]*t;
    Point2D pt=_nodesArray[e0]*(1.0-t)+_nodesArray[e1]*t;
    return normalLength(pt,pt2,szt,sz);
}

Real  DelaunayAFT2D::normalLength(const Point2D& nd1,const Point2D& pt,Real sz1,Real sz2)const
{
    Real dis=distance(nd1,pt);
	if (!(sz1>0.0 && sz2>0.0))
	{
		int aaaaaaaaaaaaa = 0;
	}
    Assert(sz1>0.0 && sz2>0.0);
    dis=NormalLength(dis,sz1,sz2);
    return dis;
}

UInt DelaunayAFT2D::createFieldPoint(const Point2D& pt)
{
    _nodeMark.push_back(UIntBits());
    _nodesArray.push_back(pt);
    _linkNodeArray.push_back(Link(InValidEntity,-1));
    _nodeSize.push_back(-1.0);
    return _nodesArray.size()-1;
}

void Center(const Point2D& pta,const Point2D& ptb,const Point2D& ptc,const Point2D& pt,Real  coord[3])
{
    Real areaMe=1.0/Area2(pta,ptb,ptc);
    coord[0]=Area2(ptb,ptc,pt)*areaMe;
    coord[1]=Area2(ptc,pta,pt)*areaMe;
    coord[2]=Area2(pta,ptb,pt)*areaMe;
}

Real  DelaunayAFT2D::sizeInterpolate(UInt tri,const Point2D& pos)const
{
    Real  ptc[3];
    Center(_nodesArray[_workingEleArray[tri][0]],
        _nodesArray[_workingEleArray[tri][1]],
        _nodesArray[_workingEleArray[tri][2]],
        pos,ptc);
    Real sz=_nodeSize[_workingEleArray[tri][0]]*ptc[0]+_nodeSize[_workingEleArray[tri][1]]*ptc[1]+_nodeSize[_workingEleArray[tri][2]]*ptc[2];
    sz*=_gr;
    if(_maxTriSize>0.0)
    {
        if(sz>_maxTriSize) sz= _maxTriSize;
    }
    return sz;
}

UInt DelaunayAFT2D::insertDirect(const Link& triHandle,const Point2D& pt,Link newTri[4])
{
    if(triHandle.second<0)
    {
        UInt newNodeHd(createFieldPoint(pt));
        UInt nodeH[3];
        nodeH[0]= _workingEleArray[triHandle.first][0];
        nodeH[1]= _workingEleArray[triHandle.first][1];
        nodeH[2]= _workingEleArray[triHandle.first][2];
        Link lnk[3];
        lnk[0]= _adjElementsArray[triHandle.first][0];
        lnk[1]= _adjElementsArray[triHandle.first][1];
        lnk[2]= _adjElementsArray[triHandle.first][2];
        _workingEleArray[triHandle.first][0]=nodeH[0];
        _workingEleArray[triHandle.first][1]=nodeH[1];
        _workingEleArray[triHandle.first][2]=newNodeHd;

        UInt newTri0(addTri(nodeH[1],nodeH[2],newNodeHd));
        UInt newTri1(addTri(nodeH[2],nodeH[0],newNodeHd));
        checkinTri(newTri0,IN_DOMAIN);
        checkinTri(newTri1,IN_DOMAIN);
        checkoutTri(triHandle.first,FROZEN);
        checkoutTri(newTri0,FROZEN);
        checkoutTri(newTri1,FROZEN);
        newTri[0]=Link(triHandle.first,2);
        newTri[1]=Link(newTri0,2);
        newTri[2]=Link(newTri1,2);
        //=======================================
        setNodeLink(nodeH[0],Link(triHandle.first,0));
        setNodeLink(nodeH[1],Link(triHandle.first,1));
        setNodeLink(nodeH[2],Link(newTri1,0));
        setNodeLink(newNodeHd,Link(triHandle.first,2));

        //=========================================
        setTriLink(triHandle.first,0,Link(newTri0,1));
        setTriLink(triHandle.first,1,Link(newTri1,0));
        setTriLink(triHandle.first,2,lnk[2]);
        setTriLink(newTri0,0,Link(newTri1,1));
        setTriLink(newTri0,2,lnk[0]);
        setTriLink(newTri1,2,lnk[1]);
        return newNodeHd;
    }else
    {

        if(0)
        {
            std::ofstream out("debugit.nas");
            exportBlk(out,2,true);
        }
        UInt tris[2],nd[4];
        Link lnk[4];bool super;
        if(getSwapData(triHandle,tris,nd,lnk,super))
        {
            UInt newNodeHd(createFieldPoint(pt));
            _workingEleArray[tris[0]][0]=nd[0];
            _workingEleArray[tris[0]][1]=newNodeHd;
            _workingEleArray[tris[0]][2]=nd[2];

            _workingEleArray[tris[1]][0]=nd[1];
            _workingEleArray[tris[1]][1]=newNodeHd;
            _workingEleArray[tris[1]][2]=nd[3];

            UInt newTri0(addTri(newNodeHd,nd[1],nd[2]));
            UInt newTri1(addTri(newNodeHd,nd[0],nd[3]));

            checkinTri(newTri0,IN_DOMAIN);
            checkinTri(newTri1,IN_DOMAIN);
            checkoutTri(tris[0],FROZEN);
            checkoutTri(tris[1],FROZEN);
            checkoutTri(newTri0,FROZEN);
            checkoutTri(newTri1,FROZEN);
            newTri[0]=Link(tris[0],1);
            newTri[1]=Link(tris[1],1);
            newTri[2]=Link(newTri0,0);
            newTri[3]=Link(newTri1,0);
            //=======================================
            setNodeLink(nd[0],Link(tris[0],0));
            setNodeLink(nd[1],Link(tris[1],0));
            setNodeLink(nd[2],Link(tris[0],2));
            setNodeLink(nd[3],Link(tris[1],2));

            setNodeLink(newNodeHd,Link(tris[0],1));

            //=========================================
            setTriLink(tris[0],0,Link(newTri0,1));
            setTriLink(tris[0],1,lnk[1]);
            setTriLink(tris[0],2,Link(newTri1,2));

            setTriLink(tris[1],0,Link(newTri1,1));
            setTriLink(tris[1],1,lnk[3]);
            setTriLink(tris[1],2,Link(newTri0,2));

            setTriLink(newTri0,0,lnk[0]);
            setTriLink(newTri1,0,lnk[2]);
            if(0)
            {
                std::ofstream out("debugit2.nas");
                exportBlk(out,2,true);
            }
            return newNodeHd;
        }
    }
    Assert(0);
    return InValidEntity;
}

/** \brief 产生新节点ptrt并设置新节点处的尺寸值nodesizert
*/
bool DelaunayAFT2D::get3rdPoint(UInt tri,
                                Point2D& ptrt,
                                Real& nodesizert)
{
    Link trifound;
    Real nodesize;
    const Point2D* pt2ds[3]; ///< 用于存储三角形tri的三个顶点
    Real  h[3]; ///< 存储三角形三顶点处的尺寸值

    int whis[3],num=0;
    int i=0,j=0,kk;
	int n_front_node = 0;
    for(i=0;i<3;i++)
    {
		if (isCheckedNode(_workingEleArray[tri][i], FRONT_NODE))
			++n_front_node;
        pt2ds[i]=&(_nodesArray[_workingEleArray[tri][i]]); ///< 节点坐标
        h[i]=_nodeSize[_workingEleArray[tri][i]]; ///< 节点处的尺寸值

        UInt adj=_adjElementsArray[tri][i].first; ///< 邻接单元
        if(adj!=InValidEntity && (!isCheckedTri(adj,ERASED)))
        {
            if(isCheckedTri(adj,FROZEN))
            {
                whis[num++]=i;
            }
			else
            {
                Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[i][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[i][1]]);
                if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei))
                {
                    whis[num++]=i;
                }
            }
        }else
        {
            Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[i][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[i][1]]);
            if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei))
            {
                whis[num++]=i;
            }
        }
    }
    if(num>=2 && 3 != n_front_node)
    {
// 		Real quai = 0.0;
// 		getTriQua(tri, quai);
// 		if (std::abs(quai - std::sqrt(3.0) / 2.0) <= 0.2)
        checkinTri(tri,FROZEN);
        return false;
    }else if(num==0)
    {
        return false;
    }

	int i_shortest = 0; ///< 最短前沿边的编号
	Real min_l = DBL_MAX;
	if (num >= 2)///< 计算获得最短前沿边
	{
		for(i=0;i<num;i++)
		{
			int whi=whis[i];
			const Point2D& pta=*(pt2ds[TRI_EDGE_VERTEX[whi][0]]);
			const Point2D& ptb=*(pt2ds[TRI_EDGE_VERTEX[whi][1]]);
			Real l_temp = length(pta, ptb);
			if (min_l > l_temp)
			{
				i_shortest = i;
				min_l = l_temp;
			}
		}
	}

    Real diffmin=DBL_MAX,diffMax=0;
    bool rt=false;
    int cnt=0;
    //for(i=0;i<num;i++)
    {
        int whi=whis[i_shortest];
        const Point2D& pta=*(pt2ds[TRI_EDGE_VERTEX[whi][0]]);
        const Point2D& ptb=*(pt2ds[TRI_EDGE_VERTEX[whi][1]]);

		UInt nia = _workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]];
		UInt nib = _workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];

        Point2D midPos=(ptb+pta)*0.5; ///< 计算当前前沿的中点
        Real szmid=(h[TRI_EDGE_VERTEX[whi][0]]+h[TRI_EDGE_VERTEX[whi][1]])*0.5; ///< 计算该边两端点尺寸值的平均值(即边中点处的尺寸值)
        Real gr_temp = 1.0;///< 节点尺寸缩放因子
		Vector2D dir(pta,ptb);
        Real off=dir.length(); ///< 计算当前边的边长
        Point2D newPos,pos;
//         Real tmp=dir[0];
//         dir[0]=-dir[1];
//         dir[1]=tmp; ///< 上述三步的交换操作相当于将向量dir绕逆时针方向旋转90度
//         dir*=0.866;

		Point2D p_left, p_mid, p_right;
		int left_id, mid_id, right_id;
		if (isCheckedNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]], FRONT_NODE_USED)
			&& isCheckedNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]], FRONT_NODE_USED))
		{
			checkoutNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]], FRONT_NODE_USED);
			checkoutNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]], FRONT_NODE_USED);
		}
		bool right_big = false;
		bool left_big = false;
		if (!isCheckedNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]], FRONT_NODE_USED))
		{
			p_mid = pta;
			midPos = pta; ///< added by ZhangJun 2016.06.05
			mid_id = nia;

			Array1D<Link> edge_links;
			getAdjacentEdge(nia, edge_links); ///< 获取节点的邻接边组
			UInt i,sz = edge_links.size();

			int n_front = 0; ///< 用于统计点pta的邻接边中前沿边的个数
			//Array1D<Link> front_edge;
			Array1D<std::pair<Link, bool>> front_edge; ///< first:边, second: 该边在已冻结三角形中是否为最长边
			for(i=0; i<sz; i++) ///< 遍历节点的邻接边组
			{
				UInt elei = edge_links[i].first;
				UInt adj = _adjElementsArray[edge_links[i].first][edge_links[i].second].first;
				if (isCheckedTri(edge_links[i].first,OUT_DOMAIN))
				{
					if (adj != InValidEntity && isCheckedTri(adj,IN_DOMAIN))
					{
						if (!isCheckedTri(adj,FROZEN))
						{
							++n_front;
							front_edge.push_back( std::make_pair(_adjElementsArray[edge_links[i].first][edge_links[i].second], false));
						}
					}
				}
				else
				{
					if (isCheckedTri(edge_links[i].first,FROZEN) 
						&& isCheckedTri(adj,IN_DOMAIN) 
						&& !isCheckedTri(adj,FROZEN))
					{
						++n_front;
						bool flag_temp = false;
						if (isLongest(edge_links[i]))
							flag_temp = true;
						front_edge.push_back(std::make_pair(_adjElementsArray[edge_links[i].first][edge_links[i].second], flag_temp));
					}
					else if (!isCheckedTri(edge_links[i].first,FROZEN) 
						     && (isCheckedTri(adj,OUT_DOMAIN) || isCheckedTri(adj,FROZEN)))
					{
						++n_front;
						bool flag_temp = false;
						if (isCheckedTri(adj,FROZEN) && isLongest(_adjElementsArray[edge_links[i].first][edge_links[i].second]))
							flag_temp = true;
						front_edge.push_back(std::make_pair(edge_links[i], flag_temp));
					}
				}
			}
			if (2 != n_front)
				return false;
			assert(n_front == 2);
			if (_workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][0]] == nia)
			{
				left_id = _workingEleArray[front_edge[1].first.first][TRI_EDGE_VERTEX[front_edge[1].first.second][0]];
				p_left = _nodesArray[left_id];
				right_id = _workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][1]];
				p_right = _nodesArray[right_id];

				if (front_edge[1].second)
					left_big = true;
				if (front_edge[0].second)
					right_big = true;
			} 
			else
			{
				left_id = _workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][0]];
				p_left = _nodesArray[left_id];
				right_id = _workingEleArray[front_edge[1].first.first][TRI_EDGE_VERTEX[front_edge[1].first.second][1]];
				p_right = _nodesArray[right_id];
				if (front_edge[0].second)
					left_big = true;
				if (front_edge[1].second)
					right_big = true;
			}
			if (!left_big && !right_big)
				checkinNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]], FRONT_NODE_USED);
		}
		else if (!isCheckedNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]], FRONT_NODE_USED))
		{
			p_mid = ptb;
			midPos = ptb; ///< added by ZhangJun 2016.06.05
			mid_id = nib;

			Array1D<Link> edge_links;
			getAdjacentEdge(nib, edge_links); ///< 获取节点的邻接边组
			UInt i,sz=edge_links.size();

			int n_front = 0; ///< 用于统计邻接单元中被冻结或位于待剖分域外部的单元
			Array1D<std::pair<Link, bool>> front_edge; ///< first:边, second: 该边在已冻结三角形中是否为最长边
			for(i=0; i<sz; i++) ///< 遍历邻接面片组
			{
				UInt elei = edge_links[i].first;
				UInt adj = _adjElementsArray[edge_links[i].first][edge_links[i].second].first;
				if (isCheckedTri(edge_links[i].first,OUT_DOMAIN))
				{
					if (adj != InValidEntity && isCheckedTri(adj,IN_DOMAIN))
					{
						if (!isCheckedTri(adj,FROZEN))
						{
							++n_front;
							front_edge.push_back(std::make_pair(_adjElementsArray[edge_links[i].first][edge_links[i].second], false));
						}
					}
				}
				else
				{
					if (isCheckedTri(edge_links[i].first,FROZEN) 
						&& isCheckedTri(adj,IN_DOMAIN) 
						&& !isCheckedTri(adj,FROZEN))
					{
						++n_front;
						bool flag_temp = false;
						if (isLongest(edge_links[i]))
							flag_temp = true;
						front_edge.push_back(std::make_pair(_adjElementsArray[edge_links[i].first][edge_links[i].second], flag_temp));
					}
					else if (!isCheckedTri(edge_links[i].first,FROZEN) 
						&& (isCheckedTri(adj,OUT_DOMAIN) || isCheckedTri(adj,FROZEN)))
					{
						++n_front;
						bool flag_temp = false;
						if (isCheckedTri(adj,FROZEN) && isLongest(_adjElementsArray[edge_links[i].first][edge_links[i].second]))
							flag_temp = true;
						front_edge.push_back(std::make_pair(edge_links[i], flag_temp));
					}
				}
			}
			if (2 != n_front)
				return false;
			assert(2 == n_front);
			if (_workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][0]] == nib)
			{
				left_id = _workingEleArray[front_edge[1].first.first][TRI_EDGE_VERTEX[front_edge[1].first.second][0]];
				p_left = _nodesArray[left_id];
				right_id = _workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][1]];
				p_right = _nodesArray[right_id];

				if (front_edge[1].second)
					left_big = true;
				if (front_edge[0].second)
					right_big = true;
			} 
			else
			{
				left_id = _workingEleArray[front_edge[0].first.first][TRI_EDGE_VERTEX[front_edge[0].first.second][0]];
				p_left = _nodesArray[left_id];
				right_id = _workingEleArray[front_edge[1].first.first][TRI_EDGE_VERTEX[front_edge[1].first.second][1]];
				p_right = _nodesArray[right_id];

				if (front_edge[0].second)
					left_big = true;
				if (front_edge[1].second)
					right_big = true;
			}
			if (!left_big && !right_big)
				checkinNode(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]], FRONT_NODE_USED);
		}

		double angle_left_mid_right = angle_evaluate(mid_id, left_id, mid_id, right_id);

		if (left_big && mid_id == nia)
			angle_left_mid_right /= 1.5;
		else if (left_big && mid_id == nib)
			angle_left_mid_right /= 3.0;
		else if (right_big && mid_id == nia)
			angle_left_mid_right /= 3.0;
		else if (right_big && mid_id == nib)
			angle_left_mid_right /= 1.5;
		else
		    angle_left_mid_right /= 1.999999999;///< 矢量dir需要旋转的角度

		dir.rotateZ(angle_left_mid_right);

		szmid = h[TRI_EDGE_VERTEX[whi][0]]; ///< added by ZhangJun 2016.06.05

        Point2D ptIdel = midPos+dir; ///< 计算得到的理想点(该点与当前边组成正三角形)
        newPos = ptIdel;

        bool finded=false;
        UInt trifoundTemp=InValidEntity;
        Real minDiff=DBL_MAX;
        for(j=0; j<20; j++)
        {
            if(locate(newPos,trifound) != OUTDOMAIN) ///< 找到点所在的三角形
            {
                if(!isCheckedTri(trifound.first,OUT_DOMAIN))
                {
                    Real nh0= 0.866/normalLength(midPos,newPos,szmid);

                    if(fabs(1.0-nh0)<LENGTHINTERVAL[3])
                    {
                        trifoundTemp=trifound.first;
                        ptIdel=newPos;
                        pos=newPos;
                        break;
                    }
					else 
                    {
                        if(fabs(1.0-nh0)<minDiff) 
                        {
                            trifoundTemp=trifound.first;
                            ptIdel=newPos;
                            minDiff=fabs(1.0-nh0);
                        }
                        Vector2D vec1(midPos,newPos);
                        newPos+=vec1*(nh0-1.0);
                    }

                }
				else  
					newPos=(midPos+newPos)*0.5;
            }
			else 
				newPos=(midPos+newPos)*0.5;
        }
        if(trifoundTemp==InValidEntity || isCheckedTri(trifoundTemp,ERASED)) return false;
        if(j>=10)
        {
            trifound.first=trifoundTemp;
            pos=ptIdel;
        }


        nodesize=sizeInterpolate(trifound.first,pos);

        Real hi=(h[TRI_EDGE_VERTEX[whi][0]]+h[TRI_EDGE_VERTEX[whi][1]])*0.5;
        hi*=_gr;
        nodesize=hi*(1.0-_gdsc)+nodesize*_gdsc;
        if(_maxTriSize>0.0)
        {
            if(nodesize>_maxTriSize) nodesize= _maxTriSize;
        }

		nodesize = szmid * gr_temp;

        for(kk=0;kk<20;kk++)
        {
            bool breakthisloop=true;
            Vector2D vec1(pta,pos);
            //Vector2D vec2(ptb,pos);
			Vector2D vec2(pta, ptb);
            Point2D pt11(pos),pt12(pos);
            Real nh0= 1.0/normalLength2(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]],pos,nodesize,nodesize);
            //Real nh1= 1.0/normalLength2(_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]],pos,nodesize,nodesize);
            if(fabs(nh0-1.0)>LENGTHINTERVAL[3])
            {
                pt11+=vec1*(nh0-1.0);
                breakthisloop=false;
            }
//             if(fabs(nh1-1.0)>LENGTHINTERVAL[3])
//             {
//                 breakthisloop=false;
//                 pt12+=vec2*(nh1-1.0);
//             }
            if(breakthisloop) break;
            //Real all=fabs(nh0-1)+fabs(nh1-1.0);
           // all=1.0/all;
            Point2D cen=pt11/**fabs(nh0-1)*all+pt12*fabs(nh1-1.0)*all*/;
            pos+=cen;
            pos*=0.5;
            if(!isInDomain(pos)) 
            {
                return false;
            }
        }


        if(num==1 || i==0)
        {
            cnt++;
            ptrt=pos;
            nodesizert=nodesize;
            if(num==1) return true;
        }else
        {
            if(cnt==0)
            {
                ptrt=pos;
				nodesizert=nodesize;
                return true;
            }else
            {
                pos+=ptrt;

                pos*=0.5;
                ptrt=pos;
				nodesizert=nodesize;
//                 nodesizert+=nodesize;
//                 nodesizert*=0.5;
                return true;
            }
        }  
    }
    return false;
}

bool DelaunayAFT2D::get3rdPoint_mid(UInt tri,
                                    Point2D& mid_pt,
                                    Real& nodesizert)
{
    Link trifound;
    const Point2D* pt2ds[3];
    Real  h[3];
    UInt tri_nodes[3];
    int whis[3],num=0;
    int i=0,j=0;
    for(i=0;i<3;i++)
    {
        tri_nodes[i]=_workingEleArray[tri][i];
        pt2ds[i]=&(_nodesArray[tri_nodes[i]]);
        h[i]=_nodeSize[tri_nodes[i]];

        UInt adj=_adjElementsArray[tri][i].first;
        if(adj!=InValidEntity && (!isCheckedTri(adj,ERASED)))
        {
            if(isCheckedTri(adj,FROZEN))
            {
                whis[num++]=i;
            }else
            {
                Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[i][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[i][1]]);
                if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei))
                {
                    whis[num++]=i;
                }
            }
        }else
        {
            Edge2 ei(_workingEleArray[tri][TRI_EDGE_VERTEX[i][0]],_workingEleArray[tri][TRI_EDGE_VERTEX[i][1]]);
            if(_boundaryEdgeSet.end()!= _boundaryEdgeSet.find(ei))
            {
                whis[num++]=i;
            }
        }
    }

    if(num!=2) 
        return false;
    int mid_whi=3-whis[0]-whis[1];
    UInt tri_adj=_adjElementsArray[tri][mid_whi].first;
    if(isCheckedTri(tri_adj,FROZEN) || isCheckedTri(tri_adj,OUT_DOMAIN))
    {
        return false;
    }

    Real dis=normalLength(tri_nodes[TRI_EDGE_VERTEX[mid_whi][0]],tri_nodes[TRI_EDGE_VERTEX[mid_whi][1]]);
    if(dis<LENGTHINTERVAL[1])
        return false;

    Real eps=dis*Tolerance::EPS_RELATIVE;
    int cnt=0;

    Point2D pta=*(pt2ds[TRI_EDGE_VERTEX[mid_whi][0]]);
    Point2D ptb=*(pt2ds[TRI_EDGE_VERTEX[mid_whi][1]]);

    h[0]=_nodeSize[tri_nodes[TRI_EDGE_VERTEX[mid_whi][0]]];
    h[1]=_nodeSize[tri_nodes[TRI_EDGE_VERTEX[mid_whi][1]]];
    h[2]=(h[0]+h[1])*0.5;
    while(cnt<50)
    {
        mid_pt=(pta+ptb)*0.5;
        Real len1=normalLength2(tri_nodes[TRI_EDGE_VERTEX[mid_whi][0]],mid_pt,h[0],h[2]);
        Real len2=normalLength2(tri_nodes[TRI_EDGE_VERTEX[mid_whi][1]],mid_pt,h[1],h[2]);
        if(fabs(len1-len2)<eps)
        {
            break;
        }else
        {
            if(len1>len2) ptb=mid_pt;
            else pta=mid_pt;
        }
        cnt++;
    }
    Real dis0=normalLength2(tri_nodes[TRI_EDGE_VERTEX[mid_whi][0]],mid_pt,h[0],h[2]);
    Real dis1=normalLength2(tri_nodes[TRI_EDGE_VERTEX[mid_whi][1]],mid_pt,h[1],h[2]);
    Real t=dis0/(dis0+dis1);
    nodesizert=_nodeSize[tri_nodes[TRI_EDGE_VERTEX[mid_whi][0]]]*(1.0-t)+_nodeSize[tri_nodes[TRI_EDGE_VERTEX[mid_whi][1]]]*t;
    return true;
}



Real DelaunayAFT2D::anglecosine(UInt nd0, UInt nd1,UInt nd2,UInt nd3)const
{
    Vector2D v0(_nodesArray[nd0],_nodesArray[nd1]);
    Vector2D v1(_nodesArray[nd2],_nodesArray[nd3]);
    v0.normalize();v1.normalize();
    return v0.dot(v1);
}

/** \brief 计算矢量nd0--nd1 和 矢量nd2--nd3的夹角
    \param nd0 第一个矢量的起点编号
	\param nd1 第一个矢量的终点编号
	\param nd2 第二个矢量的起点编号
	\param nd3 第二个质量的终点编号
	\return 矢量夹角(单位：度)
*/
Real DelaunayAFT2D::angle_evaluate(UInt nd0, UInt nd1,UInt nd2,UInt nd3)const
{
	Real cosAngle = anglecosine(nd0, nd1, nd2, nd3);
	Real angle_temp = std::acos(cosAngle) * 180.0 / PI;

	Real det= Orient2d_Math(_nodesArray[nd1],_nodesArray[nd0],_nodesArray[nd3]);
	if (det<0)
	{
		angle_temp = 360 - angle_temp;
	}
		
	return angle_temp;
}

int   DelaunayAFT2D:: getRTriPoint(const FrontNode& nd,NodeSizePair  nodesizert[2])
{
    static UInt tmp=849;
    UInt nd0=nd.ni,nd1(InValidEntity),nd2(InValidEntity);
    if(nd0==tmp)
    {
        int xxxx=0;
    }
    UInt tri;int whi;
    Vector2D v1,v2;
    Real len0=-1;
    if(nd.edge[0].first!=InValidEntity)
    {
        tri=nd.edge[0].first; ///< 边所属的三角形编号
        whi=nd.edge[0].second;///< 边在三角形中的局部编号
        nd1=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]]; ///< 边的左端点的编号
        if(nd1 == nd0)
        {
            nd1 = _workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];
        }
        Vector2D v1tmp(_nodesArray[_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]]],
            _nodesArray[_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]]]);
        v1.set(-v1tmp[1],v1tmp[0]);
        len0=v1.length();
    }

    if(nd.edge[1].first!=InValidEntity)
    {
        tri=nd.edge[1].first;
        whi=nd.edge[1].second;
        nd2=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]];
        if(nd2==nd0)
        {
            nd2=_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]];
        }
        Vector2D v2tmp(_nodesArray[_workingEleArray[tri][TRI_EDGE_VERTEX[whi][0]]],
            _nodesArray[_workingEleArray[tri][TRI_EDGE_VERTEX[whi][1]]]);
        v2.set(-v2tmp[1],v2tmp[0]);
        len0+=v2.length();
        len0*=0.5;
    }

    int dircnt=0;
    Vector2D dir[2];
    if(InValidEntity!=nd1 && InValidEntity!=nd2)
    {
        Vector2D v01(_nodesArray[nd0],_nodesArray[nd1]);
        Vector2D v02(_nodesArray[nd0],_nodesArray[nd2]);
        v01.normalize();v02.normalize();
        Real det0=v01.dot(v02)+1;
        if(det0<0.001)
        {
            // near 180 degree
            v01=v1;v02=v2;
            v01.normalize();v02.normalize();
            Point2D ptm=(_nodesArray[nd0]+v01+_nodesArray[nd0]+v02)*0.5;
            dir[0].set(_nodesArray[nd0],ptm);
            dir[0].normalize();
            dircnt=1;
        }else
        {
            dircnt=1;
            Point2D ptm=(_nodesArray[nd0]+v01+_nodesArray[nd0]+v02)*0.5;
            dir[0].set(_nodesArray[nd0],ptm);
            dir[0].normalize();
            if(dir[0].dot(v1)<0.0)//96
            {
                // >180
                dir[0]=-dir[0];
                Real tdet=det0-1.0;
                if(tdet>0.25882) return 0;
                if(tdet>-0.258812)
                {
                    dircnt=2;
                    dir[0]=-v01;
                    dir[1]=-v02;
                    dir[0].normalize();
                    dir[1].normalize();
                    dir[0]*=len0;
                    dir[1]*=len0;
                } 
            }
            if(dircnt==1)
            {
                if(dir[0].dot(v2)<-0.10453) return 0;
                v01.normalize();
                Real det=dir[0].dot(v01);
                if(det>0.4226)// 65 degree
                    return 0;
            }
        }
        if(dircnt==1) 
        {
            dir[0]*=len0;
        }
    }else 
    {
        dircnt=1;
        if(InValidEntity==nd1)
        {
            dir[0]=v2;
        }
        else if(InValidEntity==nd2)dir[0]=v1;
    }
    int nd_cnt=0;
    for(int it=0;it<dircnt;it++)
    {
        UInt j,trifoundTemp(InValidEntity);
        Link trifound;
        Real minDiff=DBL_MAX;
        Point2D newPos=_nodesArray[nd0]+dir[it];
        Point2D ptIdel(newPos),pos(newPos);
        Point2D midPos(_nodesArray[nd0]);
        Real szmid=_nodeSize[nd0];
        for(j=0;j<20;j++)
        {
            if(locate(newPos,trifound)!=OUTDOMAIN)
            {
                if(!isCheckedTri(trifound.first,OUT_DOMAIN))
                {
                    Real nh0= normalLength(midPos,newPos,szmid);

                    if(fabs(1.0-nh0)<LENGTHINTERVAL[3])
                    {
                        trifoundTemp=trifound.first;
                        ptIdel=newPos;
                        pos=newPos;
                        break;
                    }else 
                    {
                        if(fabs(1.0-nh0)<minDiff) 
                        {
                            trifoundTemp=trifound.first;
                            ptIdel=newPos;
                            minDiff=fabs(1.0-nh0);
                        }
                        Vector2D vec1(midPos,newPos);
                        newPos+=vec1*(nh0-1.0);
                    }

                }else  newPos=(midPos+newPos)*0.5;
            }else  newPos=(midPos+newPos)*0.5;
        }
        if(trifoundTemp==InValidEntity || isCheckedTri(trifoundTemp,ERASED)) continue;
        if(j>=10)
        {
            trifound.first=trifoundTemp;
            pos=ptIdel;
        }

        nodesizert[nd_cnt].first=pos;
        nodesizert[nd_cnt].second=sizeInterpolate(trifound.first,pos);
        nd_cnt++;
    }
    return nd_cnt;
}


struct  ElementLnkViAFT2 : public Delaunay2D::LinkVisitor
{
    bool visit(const Link& lnk)
    {
        _canInsert=_parent->canBeInserted(_pt,_nodesz,lnk,_nodes);
        if(!_canInsert) return false;

        return true;
    }
    bool canbeinsert(){return _canInsert;}

    ElementLnkViAFT2(DelaunayAFT2D* parent,const Point2D& pt,Real sz)
        :_parent(parent),_pt(pt),_nodesz(sz)
    {
        _canInsert=true;
    }
    ~ElementLnkViAFT2()
    {
        UInt i,sz= _nodes.size();
        for(i=0;i<sz;i++)
        {
            _parent->checkoutNode(_nodes[i],DelaunayAFT2D::MAX_FLAG);
        }
    }
protected:
    DelaunayAFT2D*        _parent;
    Point2D              _pt;
    Real                 _nodesz;
    Array1D<UInt>         _nodes;
    bool                 _canInsert;
};

bool DelaunayAFT2D::canBeInserted(const Point2D& pt,Real nodesz,
                                  const Link& lnk,
                                  Array1D<UInt>& checkedNodes)
{
    UInt triHandle=lnk.first;
    int whi=lnk.second;

    UInt ndi[2];
    ndi[0]=_workingEleArray[triHandle][TRI_NEXT[whi]]; 
    ndi[1]=_workingEleArray[triHandle][TRI_PRE[whi]];
    for(int i=0;i<2;i++)
    {
        if(isCheckedNode(ndi[i],MAX_FLAG)) continue ;
        if(isSuperNode(ndi[i])) continue;
        Real nhi=normalLength(_nodesArray[ndi[i]],pt,_nodeSize[ndi[i]],nodesz); 
        //if(nhi<LENGTHINTERVAL[0]) 
		if (nhi < nodesz * 9)
        {
            return false;
        }else
        {
            checkinNode(ndi[i],MAX_FLAG);
            checkedNodes.push_back(ndi[i]);
        }
    }
    return true;
}

void DelaunayAFT2D:: outputDebugMesh(const char* filename)
{
    char buffer[1024];
    sprintf(buffer,"%s.nas",filename);
    std::ofstream out(buffer);
    exportBlk(out,2,true);
}

bool DelaunayAFT2D::canBeInserted(const Point2D& pt,Link& whiTri,int& loc,Real szNode)
{
    loc=locate(pt,whiTri);
    if(loc==OUTDOMAIN) return false;
    if(isCheckedTri(whiTri.first,OUT_DOMAIN)) return false ;
    if(isCheckedTri(whiTri.first,FROZEN)) return false ;
    ElementLnkViAFT2 vi(this,pt,szNode);
    int j=0;
    for(j=0;j<3;j++)
    {
        UInt ndj=_workingEleArray[whiTri.first][j];
        UInt en0=_workingEleArray[whiTri.first][TRI_EDGE_VERTEX[j][0]];
        UInt en1=_workingEleArray[whiTri.first][TRI_EDGE_VERTEX[j][1]];
        Real nhi=normalLength(_nodesArray[ndj],pt,szNode);
        // Real nhi=normalDistanceEdge(en0,en1,pt,szNode);
        if(nhi<LENGTHINTERVAL[2]) 
			break;
        visitNodeBall(ndj,vi);
        if(!vi.canbeinsert())
        {
            break;
        }
    }
    if(j<3) return false;
    return true;
}

/** \brief 收集在黎曼度量下不满足边长要求的单元，并将满足条件的三角形标记冻结
* 三角形单元最长边的长度不应大于sqrt(2.0)
*/
void  DelaunayAFT2D::markETriangle(Array1D<UInt>& badEles)
{
    Array1D<UInt> badEle;
    UInt i,j,sz=_workingEleArray.size();
    int newFrozenNum=0;
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,OUT_DOMAIN))
            continue;
        if(isCheckedTri(i,ERASED))
            continue;
        if(isCheckedTri(i,FROZEN))
            continue;

        for(j=0;j<3;j++)
        {///< 查看当前单元的三个邻接单元
            UInt adj=_adjElementsArray[i][j].first;
            if(adj!=InValidEntity && !isCheckedTri(adj,ERASED))
            {///< 如果邻接单元有效且未被标记删除
                if((!isCheckedTri(adj,FROZEN)) && (!isCheckedTri(adj,OUT_DOMAIN)))
                {///< 如果该邻接单元未被冻结且在待剖分域内部
                    break;
                }
            }
        }
        if(j>=3)
        {
            checkinTri(i,FROZEN); ///< 将当前单元标记冻结
            newFrozenNum++;
            continue;
        }

        bool frozenit=true;
		int n_edge_correct = 0;
		Real min_angle = 180.0;
		Real max_angle = 0.0;
		Real max_l = DBL_MIN;
		REAL min_l = DBL_MAX;
		int min_edge_id = 0;
		Real quai = 0.0;
		getTriQua(i,quai);
		int n_front_node_active = 0;
		int n_front_node_unactive = 0;
		//double deta = std::abs(quai - sqrt(3.0) / 4.0); ///< 计算当前单元的质量系数与等腰直角三角形质量系数之差
        for(j=0;j<3;j++)///< 计算当前单元的三条边长，如果某一边长值>1.414， 则将当前单元标记为坏单元
        {
            Real sz= normalLength(_workingEleArray[i][TRI_EDGE_VERTEX[j][0]],
                _workingEleArray[i][TRI_EDGE_VERTEX[j][1]]); ///< 计算线段在黎曼度量下的长度
			if (isCheckedNode(_workingEleArray[i][j], FRONT_NODE))
			{
				if (isCheckedNode(_workingEleArray[i][j], UNACTIVE))
					++n_front_node_unactive;
				else
					++n_front_node_active;
			}
//             if(sz>LENGTHINTERVAL[1] )
//             {
//                 break;
//             }
			if (sz > max_l)
				max_l = sz;
			if (sz < min_l)
			{
				min_l = sz;
				min_edge_id = j;
			}
			double angle_temp = angle_evaluate(_workingEleArray[i][j], _workingEleArray[i][(j+2)%3],
			                                   _workingEleArray[i][j], _workingEleArray[i][(j+1)%3]);
		    if (max_angle < angle_temp)
				max_angle = angle_temp;
			if (min_angle > angle_temp)
				min_angle = angle_temp;
// 			if (sz <= LENGTHINTERVAL[1])
// 			{
// 				++n_edge_correct;
// 			}
        }
// 		if (isCheckedNode(_workingEleArray[i][TRI_EDGE_VERTEX[min_edge_id][0]],FRONT_NODE)
// 			&& isCheckedNode(_workingEleArray[i][TRI_EDGE_VERTEX[min_edge_id][1]],FRONT_NODE))
		{
			if(isLongest(_adjElementsArray[i][min_edge_id]))
			{///< 如果最短边是斜边
				badEle.push_back(i);
				_workingEleArray[i].set_min_angle(min_angle);
				continue;
			}
		}

// 		if (min_angle / max_angle < 0.45 || min_angle / max_angle < 0.55)
// 		{
// 		}

        //if(j<3) 
		//if (min_angle / max_angle < 0.45 || min_angle / max_angle > 0.55)
		//if (min_angle < 30 || max_l > 1.5 * LENGTHINTERVAL[1])
		if (std::abs(quai - std::sqrt(3.0) / 2.0) > 0.2 || max_l > 1.5 * LENGTHINTERVAL[1])
        {
            badEle.push_back(i);
			_workingEleArray[i].set_min_angle(min_angle);
            continue;
        }
// 		else if (3 == n_front_node_active || 3 == n_front_node_unactive)
// 		{
// 			badEle.push_back(i);
// 			_workingEleArray[i].set_min_angle(min_angle);
// 			continue;
// 		}
        checkinTri(i,FROZEN); ///< 如果两条边长值都小于1.414，则将当前单元标记冻结
        newFrozenNum++;
    }

    Array1D<UInt> newFrozen;
    while(newFrozenNum>0)
    {
        newFrozenNum=0;
        sz = badEle.size();
        for(i=0; i<sz; i++)
        {
            UInt curTri=badEle[i];
            if(isCheckedTri(curTri,FROZEN))  continue;
            int cnt=0;
            for(j=0;j<3;j++)
            {
                UInt adj=_adjElementsArray[curTri][j].first;
                if(adj!=InValidEntity && !isCheckedTri(adj,ERASED))
                {
                    //if(!isCheckedTri(adj,FROZEN))  break;
					if (isCheckedTri(adj,FROZEN))
						++cnt;
                }
            }
            //if(j<3) continue;
			if(cnt < 3 /*|| _workingEleArray[i].get_min_angle() < 30*/) continue; ///< modified by ZhangJun 2016.06.06
            checkinTri(curTri,FROZEN); ///< 如果当前单元的三个邻接单元均被标记为冻结，则应将当前单元标记为冻结
            newFrozenNum++;
        }
    }

    sz=badEle.size();
    for(i=0;i<sz;i++)
    {///< 从首次收集的坏单元中把被标记为冻结的单元删除
        UInt curTri=badEle[i];
        if(isCheckedTri(curTri,FROZEN))  continue;
        badEles.push_back(curTri);
    }
}

Real DelaunayAFT2D::qualityQuad(const Link& lnk)const
{
    UInt tris[2],nds[4];
    Link lnks[4];bool sp;
    getSwapData(lnk,tris,nds,lnks,sp);
    int det=validTri(nds[0],nds[3],nds[2]);
    if(det<=0) return -1;
    det=validTri(nds[1],nds[2],nds[3]);
    if(det<=0.0) return -1;
    Real d01=distance(nds[0],nds[1]);
    Real d02=distance(nds[0],nds[2]);
    Real d03=distance(nds[0],nds[3]);
    Real d12=distance(nds[1],nds[2]);
    Real d13=distance(nds[1],nds[3]);
    Real d23=distance(nds[2],nds[3]);
    KeyIndex4<Real> keyQua(QuaTriangle(d01,d02,d12),
        QuaTriangle(d01,d03,d13),
        QuaTriangle(d03,d02,d23),
        QuaTriangle(d13,d12,d23));

    return (keyQua[0]*keyQua[1])/(keyQua[2]*keyQua[3]);

}

void  DelaunayAFT2D::exportQuadBlk(std::ostream& out,int dim)
{
    char buffer[256];
    UInt i, sz=_nodesArray.size();

    out<<"$$\n$$  GRID Data\n$$\n";
    Real coord[3];
    for(i=0;i<sz;i++)
    {
        if( isSuperNode(i)) continue;
        getNodeCoord(i,dim,coord);
        print_node_long_format_nas(i+1,coord[0],coord[1],coord[2],buffer);
        out<<buffer;
    }

    out<<"$\n";

    UInt cnt=1;
    sz=_quadArray.size();
    Array1D<UInt> tris;
    for(i=0;i<sz;i++)
    {
        if(_quadArray[i][2]==_quadArray[i][3]) 
            tris.push_back(i);
    }

    if(!tris.empty())
    {
        sprintf(buffer,"CTRIA3  ");
        out<<"$  CTRIA3 Data\n$\n";
        sz=tris.size();
        for(i=0;i<sz;i++)
        {
            int prop=301;

            out<<buffer<<std::setw(8)<<cnt++<<std::setw(8)<<prop;
            int ncount = 0;
            int jfield = 0;
            // first line
            for (jfield = 4; jfield <=9; jfield++)
            {
                ncount++;
                if (ncount>3) break;
                int idi=_quadArray[tris[i]][ncount-1]+1;
                out <<std::setw(8) << idi;
            } 
            out<<"\n";
        }
    }

    sz=_quadArray.size();
    if(tris.size()!=sz)
    {
        out<<"$\n";
        sprintf(buffer,"CQUAD4  ");
        out<<"$  CQUAD4 Elements\n";

        for(i=0;i<sz;i++)
        {
            int prop=401;
            if(_quadArray[i][2]==_quadArray[i][3]) continue;
            out<<buffer<<std::setw(8)<<cnt++<<std::setw(8)<<prop;
            int ncount = 0;
            int jfield = 0;
            // first line
            for (jfield = 4; jfield <=9; jfield++)
            {
                ncount++;
                if (ncount>4) break;
                int idi=_quadArray[i][ncount-1]+1;
                out <<std::setw(8) << idi;
            } 
            out<<"\n";
        }
    }
}

void DelaunayAFT2D::erase4(UInt ndi,Array1D<Link>& lnks)
{

}

void DelaunayAFT2D::erase3(UInt nderase,Array1D<Link>& lnks)
{
    /*UInt chgNodes[3];
    UInt tri0=lnks[0].first;
    UInt tri1=lnks[1].first;
    UInt tri2=lnks[2].first;
    _adjElementsArray[tri0][lnks[0].second];

    Link adj0=_adjElementsArray[tri0][lnks[0].second];
    Link adj1=_adjElementsArray[tri1][lnks[1].second];
    Link adj2=_adjElementsArray[tri2][lnks[2].second];

    chgNodes[0]=_workingEleArray[tri0][TRI_EDGE_VERTEX[lnks[0].second][0]];
    chgNodes[1]=_workingEleArray[tri0][TRI_EDGE_VERTEX[lnks[0].second][0]];
    chgNodes[2]=_workingEleArray[tri0][TRI_EDGE_VERTEX[lnks[0].second][0]];
    int flag=0;
    if(chgNodes[2]==chgNodes[0])
    {
    chgNodes[2]=tri1->getEdgeNode(lk1.which(),0);
    flag=1;
    }
    #if _DEBUG_SLMESH25_
    {
    Assert(chgNodes[2]!=chgNodes[1] && chgNodes[2]!=chgNodes[0]);
    }
    #endif
    tri1->flagCheckin(Triangle::ERASED);
    tri2->flagCheckin(Triangle::ERASED);

    nderase->flagCheckin(Node::ERASED);
    tri0->set(chgNodes[0],chgNodes[1],chgNodes[2]);
    tri0->setLink(2,adj0);
    if(flag==0)
    {
    tri0->setLink(0,adj1);
    tri0->setLink(1,adj2);
    }else if(flag==1)
    {
    tri0->setLink(1,adj1);
    tri0->setLink(0,adj2);
    }
    return true;*/
}

void DelaunayAFT2D::erase34()
{
    UInt i, sz=_nodesArray.size();
    for(i=0;i<sz;i++)
    {
        if(i<_inputNodeSize) continue;
        if(isSuperNode(i)) continue;
        Array1D<Link> links;
        visitNodeBall(i,links);
        if(links.size()==4)
        {
            erase4(i,links);
        }else if(links.size()==3)
        {
            erase3(i,links);
        }
    }

}
void DelaunayAFT2D::trans2quad()
{
    Quad aQuad;
    Real Size_Max=LENGTHINTERVAL[1]*LENGTHINTERVAL[1];
    Real Size_Max_Min=Size_Max*0.5;
    UInt i,j,k,sz=_workingEleArray.size(),newFrozenNum=0;
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,OUT_DOMAIN))
            continue;
        if(isCheckedTri(i,ERASED))
            continue;
        if(isCheckedTri(i,QUAD_TRI)) continue;
        bool frozenit=true;
        Real szs[3],s_max=0;
        int whi=-1;
        for(j=0;j<3;j++)
        {
            szs[j]= normalLength(_workingEleArray[i][TRI_EDGE_VERTEX[j][0]],
                _workingEleArray[i][TRI_EDGE_VERTEX[j][1]]);
            if(szs[j]>Size_Max )
            {
                break;
            }
            if(szs[j]>s_max)
            {
                s_max=szs[j];
                whi=j;
            }
        }

        if(j<3) 
        {
            continue;
        }else
        {
            if(szs[TRI_NEXT[whi]]>LENGTHINTERVAL[1] || 
                szs[TRI_PRE[whi]]>LENGTHINTERVAL[1])
            {
                continue;
            }

            UInt adj= _adjElementsArray[i][whi].first;
            if(adj==InValidEntity)
            {
                continue;
            }
            if(isCheckedTri(adj,QUAD_TRI)) continue;
            int whiadj=_adjElementsArray[i][whi].second;
            UInt nd0=_workingEleArray[adj][whiadj];
            UInt nd1=_workingEleArray[adj][(whiadj+1)%3];
            UInt nd2=_workingEleArray[adj][(whiadj+2)%3];
            Real detAngle=anglecosine(nd0,nd1,nd0,nd2);
            if(fabs(detAngle)<0.1736482)
            {
                aQuad[0]=nd0;
                aQuad[1]=nd1;
                aQuad[2]=_workingEleArray[i][whi];
                aQuad[3]=nd2;
                checkinTri(i,QUAD_TRI);
                checkinTri(adj,QUAD_TRI);
                _quadArray.push_back(aQuad);
                newFrozenNum+=2;
                continue;
            }
            whiadj=TRI_NEXT[_adjElementsArray[i][whi].second];
            Real szadj= normalLength(_workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][0]],
                _workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][1]]);
            if(szadj>LENGTHINTERVAL[1])
            {
                continue;
            }
            whiadj=TRI_PRE[_adjElementsArray[i][whi].second];
            szadj= normalLength(_workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][0]],
                _workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][1]]);
            if(szadj>LENGTHINTERVAL[1])
            {
                continue;
            }
            aQuad[0]=nd0;
            aQuad[1]=nd1;
            aQuad[2]=_workingEleArray[i][whi];
            aQuad[3]=nd2;
            checkinTri(i,QUAD_TRI);
            checkinTri(adj,QUAD_TRI);
            _quadArray.push_back(aQuad);
            newFrozenNum+=2;
        }
    }

    if(0)
    {
        std::ofstream out("quad_out0.nas");
        exportQuadBlk(out,2);
    }

    //step 2 for a node linked 4 triangle
    sz=_nodesArray.size();
    for(i=0;i<sz;i++)
    {
        if(i<_inputNodeSize) continue;
        if(isSuperNode(i)) continue;
        Array1D<Link> links;
        visitNodeBall(i,links);
        if(links.size()!=4)  continue;
        for(j=0;j<4;j++)
        {
            if(isCheckedTri(links[j].first,QUAD_TRI)) break;
            if(isCheckedTri(links[j].first,OUT_DOMAIN)) break;
        }
        if(j<4) continue;
        UInt tri=links[0].first;
        int whi=links[0].second;
        aQuad[0]=_workingEleArray[tri][(whi+1)%3];
        aQuad[1]=_workingEleArray[tri][(whi+2)%3];
        int cnt=2;
        for(k=1;k<4;k++)
        {
            for(j=1;j<4;j++)
            {
                tri=links[j].first;
                whi=links[j].second;
                UInt nj=_workingEleArray[tri][(whi+1)%3];
                if(nj==aQuad[cnt-1])
                {
                    aQuad[cnt]=_workingEleArray[tri][(whi+2)%3];
                    cnt++;
                    if(cnt==4) break;
                }
            }
            if(cnt==4) break;
        }
        checkinTri(links[0].first,QUAD_TRI);
        checkinTri(links[1].first,QUAD_TRI);
        checkinTri(links[2].first,QUAD_TRI);
        checkinTri(links[3].first,QUAD_TRI);
        _quadArray.push_back(aQuad);
        newFrozenNum+=4;
    }

    if(0)
    {
        std::ofstream out("quad_out1.nas");
        exportQuadBlk(out,2);
    }
    sz=_workingEleArray.size();
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,OUT_DOMAIN)) continue;
        if(isCheckedTri(i,ERASED)) continue;
        if(isCheckedTri(i,QUAD_TRI)) continue;
        int whi=-1;
        Real quamax=0;
        for(j=0;j<3;j++)
        {
            UInt adj= _adjElementsArray[i][j].first;
            if(adj==InValidEntity)continue;
            if(isCheckedTri(adj,QUAD_TRI) || 
                isCheckedTri(adj,OUT_DOMAIN))
            {
                continue;
            }else
            {
                Real tmp=qualityQuad(Link(i,j));
                if(tmp>quamax)
                {
                    quamax=tmp;
                    whi=j;
                }
            }
        }
        if(whi<0) continue;
        UInt adj= _adjElementsArray[i][whi].first;
        int whiadj=_adjElementsArray[i][whi].second;
        UInt nd0=_workingEleArray[adj][whiadj];
        UInt nd1=_workingEleArray[adj][(whiadj+1)%3];
        UInt nd2=_workingEleArray[adj][(whiadj+2)%3];
        aQuad[0]=nd0;
        aQuad[1]=nd1;
        aQuad[2]=_workingEleArray[i][whi];
        aQuad[3]=nd2;
        checkinTri(i,QUAD_TRI);
        checkinTri(adj,QUAD_TRI);
        _quadArray.push_back(aQuad);
        newFrozenNum+=2;
    }

    if(0)
    {
        std::ofstream out("quad_out2.nas");
        exportQuadBlk(out,2);
    }
    // 
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,OUT_DOMAIN)) continue;
        if(isCheckedTri(i,ERASED)) continue;
        if(isCheckedTri(i,QUAD_TRI)) continue;
        checkinTri(i,QUAD_TRI);
        aQuad[0]=_workingEleArray[i][0];
        aQuad[1]=_workingEleArray[i][1];
        aQuad[2]=_workingEleArray[i][2];
        aQuad[3]=aQuad[2];
        _quadArray.push_back(aQuad);
    }

}

void  DelaunayAFT2D::markRTriangle(Array1D<UInt>& badEles)
{
    Array1D<UInt> badEle;
    UInt i,j,sz=_workingEleArray.size();
    int newFrozenNum=0;
    Real Size_Max=LENGTHINTERVAL[1]*LENGTHINTERVAL[1];
    Real Size_Max_Min=Size_Max*0.5;
    for(i=0;i<sz;i++)
    {
        if(isCheckedTri(i,OUT_DOMAIN))
            continue;
        if(isCheckedTri(i,ERASED))
            continue;
        if(isCheckedTri(i,FROZEN))
            continue;

        for(j=0;j<3;j++)
        {
            UInt adj=_adjElementsArray[i][j].first;
            if(adj!=InValidEntity && !isCheckedTri(adj,ERASED))
            {
                if((!isCheckedTri(adj,FROZEN)) && (!isCheckedTri(adj,OUT_DOMAIN)))
                {
                    break;
                }
            }
        }
        if(j>=3)
        {
            checkinTri(i,FROZEN);
            newFrozenNum++;
            continue;
        }

        bool frozenit=true;
        Real szs[3],s_max=0;
        int whi=-1;
        for(j=0;j<3;j++)
        {
            szs[j]= normalLength(_workingEleArray[i][TRI_EDGE_VERTEX[j][0]],
                _workingEleArray[i][TRI_EDGE_VERTEX[j][1]]);
            if(szs[j]>Size_Max )
            {
                break;
            }
            if(szs[j]>s_max)
            {
                s_max=szs[j];
                whi=j;
            }
        }

        if(j<3) 
        {
            badEle.push_back(i);
            continue;
        }else
        {
            if(szs[TRI_NEXT[whi]]>LENGTHINTERVAL[1] || 
                szs[TRI_PRE[whi]]>LENGTHINTERVAL[1])
            {
                badEle.push_back(i);
                continue;
            }

            UInt adj= _adjElementsArray[i][whi].first;
            if(adj==InValidEntity)
            {
                badEle.push_back(i);
                continue;
            }
            int whiadj=_adjElementsArray[i][whi].second;
            UInt nd0=_workingEleArray[adj][whiadj];
            UInt nd1=_workingEleArray[adj][(whiadj+1)%3];
            UInt nd2=_workingEleArray[adj][(whiadj+2)%3];
            Real detAngle=anglecosine(nd0,nd1,nd0,nd2);
            if(fabs(detAngle)<0.1736482)
            {
                checkinTri(i,FROZEN);
                checkinTri(adj,FROZEN);
                newFrozenNum+=2;
                continue;
            }
            whiadj=TRI_NEXT[_adjElementsArray[i][whi].second];
            Real szadj= normalLength(_workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][0]],
                _workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][1]]);
            if(szadj>LENGTHINTERVAL[1])
            {
                badEle.push_back(i);
                continue;
            }
            whiadj=TRI_PRE[_adjElementsArray[i][whi].second];
            szadj= normalLength(_workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][0]],
                _workingEleArray[adj][TRI_EDGE_VERTEX[whiadj][1]]);
            if(szadj>LENGTHINTERVAL[1])
            {
                badEle.push_back(i);
                continue;
            }
            checkinTri(i,FROZEN);
            checkinTri(adj,FROZEN);
            newFrozenNum+=2;
        }
    }

    Array1D<UInt> newFrozen;
    while(newFrozenNum>0)
    {
        newFrozenNum=0;
        sz=badEle.size();
        for(i=0;i<sz;i++)
        {
            UInt curTri=badEle[i];
            if(isCheckedTri(curTri,FROZEN))  continue;
            int cnt=0;
            for(j=0;j<3;j++)
            {
                UInt adj=_adjElementsArray[curTri][j].first;
                if(adj!=InValidEntity && !isCheckedTri(adj,ERASED))
                {
                    if(!isCheckedTri(adj,FROZEN))  break;
                }
            }
            if(j<3) continue;
            checkinTri(curTri,FROZEN);
            newFrozenNum++;
        }
    }

    sz=badEle.size();
    for(i=0;i<sz;i++)
    {
        UInt curTri=badEle[i];
        if(isCheckedTri(curTri,FROZEN))  continue;
        badEles.push_back(curTri);
    }
}

void   DelaunayAFT2D::markTriangle(Array1D<UInt>& badEles)
{
    if(_triType==R_TRI)
    {
        markRTriangle(badEles);
    }else
    {
        markETriangle(badEles);
    }
    /////////////////////////////
}

/** \brief 计算三角形ti的质量、面积、外心坐标
    \param[in] ti 三角形单元的编号
	\param[out] qua 三角形的质量系数
	\param[out] wi 三角形的面积
	\param[out] ptc 三角形外心坐标
*/
void DelaunayAFT2D::getSmoothData(UInt ti,Real& qua,Real &wi,Point2D& ptc)
{
    Vector2D v01(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][1]]);
    Vector2D v02(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][2]]);
    Vector2D v12(_nodesArray[_workingEleArray[ti][1]],_nodesArray[_workingEleArray[ti][2]]);
    qua=QuaTriangle(v12.length(),v02.length(),v01.length());///< 三角形质量
    ptc=_nodesArray[_workingEleArray[ti][0]]+_nodesArray[_workingEleArray[ti][1]]+_nodesArray[_workingEleArray[ti][2]];
    ptc*=0.333333333333333;///< 外心坐标
    wi=Area2(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][1]],_nodesArray[_workingEleArray[ti][2]]);
    wi*=0.5;///< 三角形面积
}

/** \brief 计算三角形单元ti的质量系数
    \param[in] ti 三角形单元的编号
	\param[out] qua 三角形单元的质量系数
*/
void DelaunayAFT2D::getTriQua(UInt ti,Real& qua)
{
    Vector2D v01(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][1]]);
    //Vector2D v02(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][1]]);
	Vector2D v02(_nodesArray[_workingEleArray[ti][0]],_nodesArray[_workingEleArray[ti][2]]); ///< modified by ZhangJun 2016.06.07
    Vector2D v12(_nodesArray[_workingEleArray[ti][1]],_nodesArray[_workingEleArray[ti][2]]);
    qua=QuaTriangle(v12.length(),v02.length(),v01.length());
}

void DelaunayAFT2D::setNodePos(UInt ni,const Point2D& newPos)
{
    _nodesArray[ni]=newPos;
}

void  DelaunayAFT2D::doSmooth(int maxIter)
{
    UInt i,sz=_nodesArray.size();
    for(int it=0;it<maxIter;it++)
    {
        int cnt=0;
        for(i=_inputNodeSize+3;i<sz;i++)
        {
            if(smooth(i))
            {
                cnt++;
            }
        }
        if(cnt==0) break;
    }
}

/** \brief 对节点ni进行光顺
* 光顺标准：应使三角形单元尽可能接近等腰直角三角形 \n
* 等腰直角三角形的质量系数为 2 * sqrt(3)
    \param ni 待光顺点的编号
*/
bool DelaunayAFT2D::smooth(UInt ni)
{
    if(ni<_inputNodeSize) return false;
    if(isSuperNode(ni)) return false;
    if(_nodesArray.size()<=ni) return false;
    Array1D<Link> links;
    visitNodeBall(ni,links);///< 获取待光顺点的邻接面片组
    UInt i,sz=links.size();
    //Real quaMin=DBL_MAX;
	Real quaMax = DBL_MIN;
    Real quai,wiall=0;
    Array1D<Real> weights(sz);
    Array1D<Point2D> centers(sz);
    for(i=0;i<sz;i++)///< 求出邻接面片组的各面片质量系数与等腰直角三角形质量系数偏差的最大值(偏差量越小说明质量越好)，及各面片面积之和
    {
        getSmoothData(links[i].first,quai,weights[i],centers[i]);
		double deta = std::abs(quai - 2 * sqrt(3.0)); ///< 计算当前单元的质量系数与等腰直角三角形质量系数之差
        if(quaMax < deta) quaMax = deta; ///< modified by ZhangJun 2016.06.07
        wiall+=weights[i];
    }

    Point2D oldPos=_nodesArray[ni];
    if(wiall>0.0)
    {
        wiall=1.0/wiall;
        Point2D newPos(0.0,0.0);
        for(i=0;i<sz;i++)
        {
            newPos+=(centers[i]*(weights[i]*wiall));
        }

        setNodePos(ni,newPos);
    }

    for(i=0;i<sz;i++)
    {
        getTriQua(links[i].first,quai);

		double deta = std::abs(quai - sqrt(3.0) / 2.0); ///< 计算当前单元的质量系数与等腰直角三角形质量系数之差

        //if(quai<=quaMin) 
		if(deta > quaMax)
        {
            setNodePos(ni,oldPos); ///< 如果光顺后偏差量增大，则放弃此次光顺
            return false;
        }
    }

    return true;
}



//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//############################################################################

//###########################################################################
int Mesh2D_DelaunayAFT(const  Array1D<Point2D>& pts, 
                       const  Array1D<Tuple<UInt,2> >& edges,
                       Array1D<TripleT<UInt> >& tris)
{
    DelaunayAFT2D delauy(pts,edges);
    int rt= delauy.doMesh();
    if(rt!=ErrorCodes::eOk) return rt;
    delauy.getTriangles(tris);
    return ErrorCodes::eOk;
}
}

#ifdef __SLDELAUNAY2D_STD_H_
using namespace DELAUNAY_MSHGEN;
#endif

#ifdef __SLMESH25CFG_SLNS_H_
using namespace MESH25;
#endif
struct DelaunayAFT2DInterface::Data
{
    Delaunay2D*  mesherPtr;
    int          type; ///< 输入数据类型
    UInt         ni,trii,tri_num,node_num;
    Data():mesherPtr(NULL),type(-1)
    {
    }
    static    void Clear(Data*&);
};

DelaunayAFT2DInterface::DelaunayAFT2DInterface()
:_data(NULL)
{

}


DelaunayAFT2DInterface::~DelaunayAFT2DInterface()
{
    Data::Clear(_data);
    _data=NULL;
}

void DelaunayAFT2DInterface::Data::Clear(Data*& data) 
{
    if(data==NULL) return;
    if(Input::MESH_POLYGON==data->type)
    {
        delete data->mesherPtr;
    }else if(Input::MESH_2D==data->type)
    {
        DelaunayAFT2D* meshPtr=(DelaunayAFT2D*)(data->mesherPtr);
        delete meshPtr;
    } 
    delete data;
}

DelaunayAFT2DInterface::Input* gInputPtr=NULL;
int MessageInput(int type,char*msg) 
{
    if(gInputPtr!=NULL)
    {
        gInputPtr->message(type,msg);
    }
    return 0;
}

int DelaunayAFT2DInterface::doMesh(Input* inputPtr)
{
    gInputPtr=inputPtr;
    Data::Clear(_data);
    _data=D_NEW Data();
    _data->type=inputPtr->mesh_type();
    int rt=eOk;
    if(Input::MESH_POLYGON==_data->type)
    {
        _data->mesherPtr=D_NEW Delaunay2D(inputPtr,MessageInput);
        rt=_data->mesherPtr->delaunay();
    }else if(Input::MESH_2D==_data->type)
    {
        DelaunayAFT2D* meshPtr=D_NEW DelaunayAFT2D(inputPtr,MessageInput);
        _data->mesherPtr=meshPtr;
        rt=meshPtr->doMesh();
    }

    gInputPtr=NULL;
    return rt;
}

/*#######################################################################
node_num edge_num  growth_ratio  max_size
x y 
...
...
e0 e1
######################################################################*/
int DelaunayAFT2DInterface::doMesh(const char* inputFile)
{
    int rt=eOk;
    std::ifstream input(inputFile);
    if(!input) return eFileOpen;
    int i,node_num, edge_num;
    Real growth_ratio,  max_size;
    input >> node_num; 
	input >> edge_num;
    input >> growth_ratio; ///<  没保存
	input >> max_size;
    Array1D<Point2D> pts(node_num);
    EdgeArray edges(edge_num);
    Real x,y;
    for(i = 0;i < node_num; i++)
    {///< 读入边界节点坐标
        input >> x; input >> y;
        pts[i][0] = x; pts[i][1] = y;
    }

    int nd1,nd2;
    for(i = 0;i < edge_num; i++)
    {///< 读入离散边界
        input>>nd1; input>>nd2;
        edges[i][0] = nd1-1; edges[i][1] = nd2-1;
    }

    Data::Clear(_data);
    _data = D_NEW Data();
    _data->type = Input::MESH_2D;
    DelaunayAFT2D* meshPtr = D_NEW DelaunayAFT2D(pts, edges);
    _data->mesherPtr = meshPtr;
    rt = meshPtr->doMesh();
    return rt;
}

int DelaunayAFT2DInterface::doMeshInp(const char* inputFile)
{
	int rt=eOk;
	std::ifstream input(inputFile);
	if(!input) return eFileOpen;

	std::string buffer;
	while (std::getline(input, buffer) && buffer != "*Node")
	{
		;
	}

	//读取节点坐标
	Array1D<Point2D> pts;
	pts.reserve(1000);
	std::stringstream sstr;
	while (std::getline(input, buffer) && buffer != "*Element, type=B21")
	{
		boost::replace_all(buffer, ",", " ");
		sstr.str("");
		sstr.clear();
		sstr << buffer;
		Point2D p_temp;
		int id;
		sstr >> id;
		sstr >> p_temp[0] >> p_temp[1];
		pts.push_back(p_temp);
	}

	if (pts.size() == 0)
	{
		std::cerr << "There are no points! File error...\n ";
		return false;
	}

	//读取约束边
	EdgeArray edges(pts.size());
	int i = 0;
	while (std::getline(input, buffer) && buffer != "*End Part")
	{
		boost::replace_all(buffer, ",", " ");
		sstr.str("");
		sstr.clear();
		sstr << buffer;
		int id;
		sstr >> id;
		int nd1;
		int nd2;
		sstr >> nd1 >> nd2;
		edges[i][0] = nd1-1; edges[i][1] = nd2-1;
		++i;
	}

	Data::Clear(_data);
	_data = D_NEW Data();
	_data->type = Input::MESH_2D;
	DelaunayAFT2D* meshPtr = D_NEW DelaunayAFT2D(pts, edges);
	_data->mesherPtr = meshPtr;
	rt = meshPtr->doMesh();
	return rt;
}

void      DelaunayAFT2DInterface::ResultIterator::node_begin()
{
    _mesherPtr->_data->ni=0;
    _mesherPtr->_data->node_num=_mesherPtr->_data->mesherPtr->nodeNum();
    _mesherPtr->_data->mesherPtr->getNodeId(_mesherPtr->_data->ni);
}

bool      DelaunayAFT2DInterface::ResultIterator::node_next(int& id, int dim, double uvxyz[5])
{
    while(1)
    {
        if(_mesherPtr->_data->ni>=_mesherPtr->_data->node_num) return false;
        if(_mesherPtr->_data->mesherPtr->isSuperNode(_mesherPtr->_data->ni)) 
        {
            _mesherPtr->_data->ni++;
            continue;
        }
        break;
    }
    const Point2D& nd2d=_mesherPtr->_data->mesherPtr->nodeCoord(_mesherPtr->_data->ni);
    id=_mesherPtr->_data->mesherPtr->getNodeId(_mesherPtr->_data->ni);
    uvxyz[0]=nd2d[0];
    uvxyz[1]=nd2d[1];
    _mesherPtr->_data->ni++;
    return true;
}

void  DelaunayAFT2DInterface::ResultIterator::triangle_begin()
{
    _mesherPtr->_data->trii=0;
    _mesherPtr->_data->tri_num=_mesherPtr->_data->mesherPtr->triangleNum();
}

bool  DelaunayAFT2DInterface::ResultIterator::triangle_next(int nodes[3])
{
    while(1)
    {
        if(_mesherPtr->_data->tri_num<=_mesherPtr->_data->trii) return false;
        if(_mesherPtr->_data->mesherPtr->isCheckedTri(_mesherPtr->_data->trii,Delaunay2D::ERASED))
        {
            _mesherPtr->_data->trii++;
        }else if(_mesherPtr->_data->mesherPtr->isCheckedTri(_mesherPtr->_data->trii,Delaunay2D::OUT_DOMAIN))
        {
            _mesherPtr->_data->trii++;
        }else break;
    }

    const Triangle& ti= _mesherPtr->_data->mesherPtr->getTriangle(_mesherPtr->_data->trii);
    nodes[0]=_mesherPtr->_data->mesherPtr->getNodeId(ti[0]);
    nodes[1]=_mesherPtr->_data->mesherPtr->getNodeId(ti[1]);
    nodes[2]=_mesherPtr->_data->mesherPtr->getNodeId(ti[2]);
    _mesherPtr->_data->trii++;
    return true;
}

void DelaunayAFT2DInterface::writeGMesh(const char* filename)
{
	std::ofstream mshFile(filename);
	mshFile << "$MeshFormat\n";
	mshFile << "2.2 0 8\n";
	mshFile << "$EndMeshFormat\n";

	int id; Real uvxyz[5];
	DelaunayAFT2DInterface::ResultIterator mesherIter(this);
	int node_number = 0; ///< 统计节点个数
	mesherIter.node_begin();
	while(1)
	{
		if(!mesherIter.node_next(id, 2, uvxyz)) break;
		if(id<0) 
		{
			continue;
		}
		++node_number;
	}
	mshFile << "$Nodes\n" << node_number << std::endl;
	mesherIter.node_begin();
	while(1)
	{
		if(!mesherIter.node_next(id, 2, uvxyz)) break;
		if(id<0) 
		{
			continue;
		}
		mshFile << id + 1
			    << "	" << uvxyz[0]/* / _data->mesherPtr->getFScale() + _data->mesherPtr->getMinPoint()[0]*/
		        << "	" << uvxyz[1]/* / _data->mesherPtr->getFScale() + _data->mesherPtr->getMinPoint()[1]*/ 
				<< "	" << 0.0
			    << std::endl;
	}
	mshFile << "$EndNodes\n";
	int ele_number = 0; ///< 统计单元数
	mesherIter.triangle_begin();

	int nodes[3];
	while(1)
	{
		if(!mesherIter.triangle_next(nodes)) break;
		++ele_number;
	}
	mshFile << "$Elements\n" << ele_number << std::endl;
	ele_number = 0;
	mesherIter.triangle_begin();
	while(1)
	{
		if(!mesherIter.triangle_next(nodes)) break;
		int ncount = 0;
		int jfield = 0;
		++ele_number;
		mshFile << ele_number << "	2	2	" << 1 << "	" <<  ele_number;
		for (jfield = 4; jfield <= 9; jfield++)
		{
			ncount++;
			if (ncount>3) break;
			int idi = nodes[ncount-1];
			if(idi<0) idi = -idi;
			mshFile << "	" << idi + 1;
		} 
		mshFile << "\n";
	}
	mshFile << "$EndElements\n";
	mshFile.close();
}

/** \brief 输出剖分结果
*/
void DelaunayAFT2DInterface::writeResult(const char* filename)
{
    std::ofstream out(filename);
    if(!out) return;
    out<<"$$\n$$  GRID Data\n$$\n";
    char buffer[1024];
    int id; Real uvxyz[5];
    DelaunayAFT2DInterface::ResultIterator mesherIter(this);
    mesherIter.node_begin();
    while(1)
    {
        if(!mesherIter.node_next(id, 2, uvxyz)) break;
        if(id<0) 
        {
            continue;
        }
        //print_node_long_format_nas(id,uvxyz[2],uvxyz[3],uvxyz[4],buffer);
		out << id << "	" << uvxyz[0] / _data->mesherPtr->getFScale() + _data->mesherPtr->getMinPoint()[0]
			      << "	" << uvxyz[1] / _data->mesherPtr->getFScale() + _data->mesherPtr->getMinPoint()[1] << std::endl;
        //out<<buffer;
    }

    out<<"$\n";
    sprintf(buffer,"CTRIA3  ");
    out<<"$  CTRIA3 Data\n$\n";
    mesherIter.triangle_begin();

    int nodes[3],cnt=1;
    while(1)
    {
        if(!mesherIter.triangle_next(nodes)) break;
        out<<buffer<<std::setw(8)<<cnt++<<std::setw(8)<<300;
        int ncount = 0;
        int jfield = 0;
        // first line
        for (jfield = 4; jfield <=9; jfield++)
        {
            ncount++;
            if (ncount>3) break;
            int idi=nodes[ncount-1];
            if(idi<0) idi=-idi;
            out <<std::setw(8) << idi;
        } 
        out<<"\n";
    }

SL_MSH_NS_END

#ifdef __SLMESH25CFG_SLNS_H_
        SL_NS_END
#endif
