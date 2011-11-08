set (ALLAUX  
   pgemraux.c)

set (IMRSRC  
   pigemr.c pigemr2.c pitrmr.c pitrmr2.c)

set (SMRSRC  
   psgemr.c psgemr2.c pstrmr.c pstrmr2.c)

set (CMRSRC  
   pcgemr.c pcgemr2.c pctrmr.c pctrmr2.c)

set (DMRSRC  
   pdgemr.c pdgemr2.c pdtrmr.c pdtrmr2.c)

set (ZMRSRC  
   pzgemr.c pzgemr2.c pztrmr.c pztrmr2.c)
   
set(redist ${ALLAUX} ${IMRSRC} ${SMRSRC} ${CMRSRC} ${DMRSRC} ${ZMRSRC})

