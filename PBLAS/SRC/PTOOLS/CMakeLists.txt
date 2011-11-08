set( ALLCTOOLS  
    PB_CGatherV.c    PB_CInV.c        PB_CInV2.c       PB_CInOutV.c     
    PB_CInOutV2.c    PB_COutV.c       PB_CScatterV.c   PB_CVMinit.c     
    PB_CVMloc.c      PB_CVMnpq.c      PB_CVMpack.c     PB_CVMswp.c      
    PB_CVMupdate.c   PB_CVMcontig.c   PB_Cabort.c      PB_Cainfog2l.c   
    PB_Cbinfo.c      PB_Cchkmat.c     PB_Cchkvec.c     PB_Cconjg.c      
    PB_Cgetbuf.c     PB_Cinfog2l.c    PB_Citypeset.c   PB_Cgcd.c        
    PB_Clcm.c        PB_Cmalloc.c     PB_Cnumroc.c     PB_Cg2lrem.c     
    PB_Cindxg2p.c    PB_Cnnxtroc.c    PB_Cnpreroc.c    PB_CpswapNN.c    
    PB_CpswapND.c    PB_Cpdot11.c     PB_CpdotNN.c     PB_CpdotND.c     
    PB_CpaxpbyNN.c   PB_CpaxpbyND.c   PB_CpaxpbyDN.c   PB_Cpaxpby.c     
    PB_CpgemmBC.c    PB_CpgemmAC.c    PB_CpgemmAB.c    PB_Cplaprnt.c    
    PB_Cplapad.c     PB_Cplapd2.c     PB_Cplascal.c    PB_Cplasca2.c    
    PB_Cplacnjg.c    PB_Cpsym.c       PB_CpsymmAB.c    PB_CpsymmBC.c    
    PB_Cpsyr.c       PB_CpsyrkA.c     PB_CpsyrkAC.c    PB_Cpsyr2.c      
    PB_Cpsyr2kA.c    PB_Cpsyr2kAC.c   PB_Cptrm.c       PB_Cpgeadd.c     
    PB_Cptradd.c     PB_Cptran.c      PB_CptrmmAB.c    PB_CptrmmB.c     
    PB_Cptrsm.c      PB_CptrsmAB.c    PB_CptrsmAB0.c   PB_CptrsmAB1.c   
    PB_CptrsmB.c     PB_Cptrsv.c      PB_Ctop.c        PB_Ctzahemv.c    
    PB_Ctzasymv.c    PB_Ctzatrmv.c    PB_Ctzhemm.c     PB_Ctzhemv.c     
    PB_Ctzher.c      PB_Ctzherk.c     PB_Ctzher2.c     PB_Ctzher2k.c    
    PB_Ctzsymm.c     PB_Ctzsymv.c     PB_Ctzsyr.c      PB_Ctzsyrk.c     
    PB_Ctzsyr2.c     PB_Ctzsyr2k.c    PB_Ctztrmm.c     PB_Ctztrmv.c     
    PB_Cwarn.c       PB_freebuf_.c    PB_topget_.c     PB_topset_.c     
    PB_Cdescset.c    PB_Cdescribe.c   PB_CargFtoC.c    PB_Cfirstnb.c    
    PB_Clastnb.c     PB_Cspan.c)

set( SCTOOLS PB_Cstypeset.c)

set( DCTOOLS PB_Cdtypeset.c)

set( CCTOOLS PB_Cctypeset.c)

set( ZCTOOLS PB_Cztypeset.c)

set(ptools ${ALLCTOOLS} ${SCTOOLS} ${DCTOOLS} ${CCTOOLS} ${ZCTOOLS})
