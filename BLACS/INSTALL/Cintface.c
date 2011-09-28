#include <stdio.h>

void c_intface_(int *i)
{
   fprintf(stderr, "For this platform, set INTFACE = -DAdd_\n");
}

void c_intface(int *i)
{
   fprintf(stderr, "For this platform, set INTFACE = -DNoChange\n");
}

void c_intface__(int *i)
{
   fprintf(stderr, "For this platform, set INTFACE = -Df77IsF2C\n");
}

void C_INTFACE(int *i)
{
   fprintf(stderr, "For this platform, set INTFACE = -DUpCase\n");
}
