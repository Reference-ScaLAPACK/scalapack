set (APBTZ  
    pxerbla.f)

set (IPBTZ  
    immadd.f       immdda.f       immtadd.f      immddat.f)

set (SPBTZ  
    svasum.f       sset.f         scshft.f       srshft.f     
    svvdot.f       smmadd.f       smmcadd.f      smmtadd.f    
    smmtcadd.f     smmdda.f       smmddac.f      smmddat.f    
    smmddact.f     stzpad.f       stzpadcpy.f    stzscal.f    
    sagemv.f       sasymv.f       satrmv.f)

set (SCPBTZ  
    sasqrtb.f      sascal.f)

set (CPBTZ  
    scvasum.f      cset.f         ccshft.f       crshft.f     
    cvvdotu.f      cvvdotc.f      cmmadd.f       cmmcadd.f    
    cmmtadd.f      cmmtcadd.f     cmmdda.f       cmmddac.f    
    cmmddat.f      cmmddact.f     ctzpad.f       ctzpadcpy.f  
    chescal.f      ctzscal.f      ctzcnjg.f      cagemv.f     
    cahemv.f       catrmv.f       casymv.f       csymv.f      
    csyr.f         csyr2.f)

set (DPBTZ  
    dvasum.f       dset.f         dcshft.f       drshft.f     
    dvvdot.f       dmmadd.f       dmmcadd.f      dmmtadd.f    
    dmmtcadd.f     dmmdda.f       dmmddac.f      dmmddat.f    
    dmmddact.f     dtzpad.f       dtzpadcpy.f    dtzscal.f    
    dagemv.f       dasymv.f       datrmv.f)

set (DZPBTZ  
    dasqrtb.f      dascal.f)

set (ZPBTZ  
    dzvasum.f      zset.f         zcshft.f       zrshft.f     
    zvvdotu.f      zvvdotc.f      zmmadd.f       zmmcadd.f    
    zmmtadd.f      zmmtcadd.f     zmmdda.f       zmmddac.f    
    zmmddat.f      zmmddact.f     ztzpad.f       ztzpadcpy.f  
    zhescal.f      ztzscal.f      ztzcnjg.f      zagemv.f     
    zahemv.f       zatrmv.f       zasymv.f       zsymv.f      
    zsyr.f         zsyr2.f)


set(ptzblas 
      ${APBTZ} ${IPBTZ} ${SPBTZ} ${SCPBTZ}
      ${CPBTZ} ${DPBTZ} ${DZPBTZ} ${ZPBTZ} )
