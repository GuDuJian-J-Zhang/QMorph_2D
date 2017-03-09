#pragma once
#ifndef QUADMESH_API
    #ifdef WIN32
        #ifdef quadMesh_EXPORTS
            #define QUADMESH_API __declspec(dllexport)
        #else
            #define QUADMESH_API __declspec(dllimport)
        #endif // qmorph_EXPORTS
    #else
        #define QUADMESH_API
    #endif // WIN32

#endif // QMORPH_API
extern "C" {

	QUADMESH_API int QuadMesh(const char* fileName, const char* ctrls);

};
