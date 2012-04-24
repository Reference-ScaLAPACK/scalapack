//
//  clamov.c
//
//  Written by Lee Killough 04/19/2012
//  

#define TYPE  complex
#define FUNC  "CLAMOV"
#if ( defined Add_ )
#define LAMOV clamov_
#define LACPY clacpy_
#elif ( defined UpCase )
#define LAMOV CLAMOV
#define LACPY CLACPY
#elif ( defined NoChange )
#define LAMOV clamov
#define LACPY clacpy
#endif
#include "lamov.h"
