//
//  dlamov.c
//
//  Written by Lee Killough 04/19/2012
//  

#define TYPE  double
#define FUNC  "DLAMOV"
#if ( defined Add_ )
#define LAMOV dlamov_
#define LACPY dlacpy_
#elif ( defined UpCase )
#define LAMOV DLAMOV
#define LACPY DLACPY
#elif ( defined NoChange )
#define LAMOV dlamov
#define LACPY dlacpy
#endif
#include "lamov.h"
