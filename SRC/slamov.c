//
//  slamov.c
//
//  Written by Lee Killough 04/19/2012
//  

#define TYPE  float
#define FUNC  "SLAMOV"
#if ( defined Add_ )
#define LAMOV slamov_
#define LACPY slacpy_
#elif ( defined UpCase )
#define LAMOV SLAMOV
#define LACPY SLACPY
#elif ( defined NoChange )
#define LAMOV slamov
#define LACPY slacpy
#endif
#include "lamov.h"
