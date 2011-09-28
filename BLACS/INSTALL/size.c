#include <stdio.h>
main()
{
   printf("ISIZE=%d\nSSIZE=%d\nDSIZE=%d\nCSIZE=%d\nZSIZE=%d\n",
          sizeof(int), sizeof(float), sizeof(double), 
          2*sizeof(float), 2*sizeof(double));
}
