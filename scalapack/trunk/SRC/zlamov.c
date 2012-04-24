//
//  zlamov.c
//
//  Written by Lee Killough 04/19/2012
//  

#define TYPE  complex16
#define FUNC  "ZLAMOV"
#if ( defined Add_ )
#define LAMOV zlamov_
#define LACPY zlacpy_
#elif ( defined UpCase )
#define LAMOV ZLAMOV
#define LACPY ZLACPY
#elif ( defined NoChange )
#define LAMOV zlamov
#define LACPY zlacpy
#endif
#include "lamov.h"
