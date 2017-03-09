#ifndef __SLDELAUNAY2D_H_
#define __SLDELAUNAY2D_H_

struct DelaunayAFT2DInterface
{
    enum 
    {
        eOk=0,
        eFailedFound,
        eDuplicateEntity,
        eWatertight,
        ePointerNULL,
        eInValidInput,
        eErrorNun,
        eFileOpen, ///< 文件打开错误
        eMeshingFailed,
        eDegeneratedGeo
    };

    struct Input
    {
        enum{MESH_POLYGON=0,MESH_2D,MESH_SURF};
        virtual void      node_begin()=0;
        virtual bool      node_next(int& id,double xy[2])=0;
        virtual void      edge_begin()=0;
        virtual bool      edge_next(int nodes[2])=0;
        virtual int       maxNodeId(){return -1;}

        virtual int       mesh_type(){return MESH_2D;}
        //##################################################
        virtual void      max_min_size(double& max,double& min){max=min -1.0;}
        virtual double    max_angle_span(){return -1.0;}
        virtual double    growth_ratio(){return 1.0;}
        virtual int       get_xyz_uv(const double* uv,double* xyz){return -1;}
        //#########################################################################
        // type:0 error message
        //      1 user  message
        virtual void      message(int type,char*){}
    };

    

    struct ResultIterator
    {
        ResultIterator(DelaunayAFT2DInterface*meshptr)
            :_mesherPtr(meshptr){}
        void      node_begin();
        bool      node_next(int& id, int dim, double xyzuv[5]);
        void      triangle_begin();
        bool      triangle_next(int nodes[3]);

        DelaunayAFT2DInterface*  _mesherPtr;
    };
    
    DelaunayAFT2DInterface();
    ~DelaunayAFT2DInterface();
    int doMesh(Input*);
    int doMesh(const char* inputFile);
	int doMeshInp(const char* inputFile); ///< 从inp文件中读取边界数据生成网格
    void writeResult(const char*);
	void writeGMesh(const char*); ///< added by Zhangjun 2016.06.04
public:
    struct  Data;
    Data*   _data;
};


#endif
