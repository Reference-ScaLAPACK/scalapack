#include <stdio.h>

void c_intface_(int *i)
{
   fprintf(stdout, "Add_\n");
}

void c_intface(int *i)
{
   fprintf(stdout, "NoChange\n");
}

void c_intface__(int *i)
{
   fprintf(stdout, "f77IsF2C\n");
}

void C_INTFACE(int *i)
{
   fprintf(stdout, "UpCase\n");
}
