#ifndef __PREDICATES_MSH_GEO_H
#define __PREDICATES_MSH_GEO_H
typedef double REAL ;
REAL exactinit();
REAL incircle(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
REAL orient2d(REAL *pa, REAL *pb, REAL *pc);
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);

#endif
