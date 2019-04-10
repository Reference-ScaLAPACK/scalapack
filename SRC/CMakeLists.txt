set (ALLAUX  pjlaenv.f  pilaenvx.f piparmq.f pilaver.f pmpim2.f pmpcol.f)
set (ALLAUX-C  pbchkvect.c getpbbuf.c pcrot.c pslaiect.c pdlaiect.c pzrot.c slamov.c clamov.c dlamov.c zlamov.c)

set (SCLAUX  
   pslabad.f pslaed0.f pslaed1.f pslaed2.f pslaed3.f pslaedz.f 
   pslamch.f pslared1d.f pslasrt.f psstebz.f psstedc.f slapst.f slasrt2.f 
   sstein2.f
   slar1va.f slarrb2.f slarrd2.f slarre2.f slarre2a.f slarrf2.f
   slarrv2.f sstegr2.f sstegr2a.f sstegr2b.f)

set (DZLAUX  
   pdlabad.f pdlaed0.f pdlaed1.f pdlaed2.f pdlaed3.f pdlaedz.f 
   pdlamch.f pdlared1d.f pdlasrt.f pdstebz.f pdstedc.f dlapst.f dlasrt2.f 
   dstein2.f
   dlar1va.f dlarrb2.f dlarrd2.f dlarre2.f dlarre2a.f dlarrf2.f
   dlarrv2.f dstegr2.f dstegr2a.f dstegr2b.f)

set (SLASRC  
   psdbsv.f  psdbtrf.f psdbtrs.f psdbtrsv.f 
   psdtsv.f  psdttrf.f psdttrs.f psdttrsv.f 
   psgbsv.f  psgbtrf.f psgbtrs.f 
   psgebd2.f psgebrd.f psgecon.f pstrcon.f psgeequ.f psgehd2.f psgehrd.f
   psgelq2.f psgelqf.f psgels.f  psgeql2.f psgeqlf.f psgeqpf.f psgeqr2.f 
   psgeqrf.f psgerfs.f psgerq2.f psgerqf.f                               
   psgesv.f  psgesvd.f psgesvx.f psgetf2.f psgetrf.f psgetri.f psgetrs.f 
   psggqrf.f psggrqf.f pslabrd.f 
   pslacon.f pslacp2.f pslacpy.f pslahrd.f pslange.f pslanhs.f pslansy.f 
   pslantr.f pslapiv.f pslapv2.f pslaqge.f pslaqsy.f pslarf.f  pslarfb.f 
   pslarfg.f pslarft.f pslase2.f pslaset.f pslascl.f pslassq.f pslaswp.f 
   pslatra.f pslatrd.f pslatrs.f pslauu2.f pslauum.f psorg2l.f psorg2r.f 
                       psorgl2.f psorglq.f psorgql.f psorgqr.f psorgr2.f 
   psorgrq.f           psorm2l.f psorm2r.f psormbr.f psormhr.f psorml2.f 
   psormlq.f psormql.f psormqr.f psormr2.f psormrq.f psormtr.f pspocon.f 
   pspbsv.f  pspbtrf.f pspbtrs.f pspbtrsv.f 
   psptsv.f  pspttrf.f pspttrs.f pspttrsv.f 
   pspoequ.f psporfs.f psposv.f  psposvx.f pspotf2.f pspotrf.f pspotri.f 
   pspotrs.f psrscl.f  psstein.f pssyev.f  pssyevd.f pssyevx.f 
   pssygs2.f pssygst.f 
   pssygvx.f pssyngst.f pssyntrd.f  pssyttrd.f pssytd2.f pssytrd.f 
   pstrti2.f pstrtri.f pstrtrs.f 
   pslaevswp.f 
   pslarzb.f pslarzt.f pslarz.f pslatrz.f pstzrzf.f psormr3.f psormrz.f 
   pslahqr.f pslaconsb.f pslacp3.f pslawil.f 
   pslasmsub.f pslared2d.f pslamr1d.f slaref.f slamsh.f slasorte.f ssteqr2.f 
   sdbtf2.f  sdbtrf.f  sdttrf.f sdttrsv.f spttrsv.f strmvt.f pssyevr.f
   bslaapp.f bslaexc.f bstrexc.f pstrord.f pstrsen.f psgebal.f pshseqr.f
   pslamve.f pslaqr0.f pslaqr1.f pslaqr2.f pslaqr3.f pslaqr4.f pslaqr5.f
   psrot.f slaqr6.f)
   
set (CLASRC  
   pcdbsv.f  pcdbtrf.f pcdbtrs.f pcdbtrsv.f 
   pcdtsv.f  pcdttrf.f pcdttrs.f pcdttrsv.f 
   pcgbsv.f  pcgbtrf.f pcgbtrs.f 
   pcgebd2.f pcgebrd.f pcgecon.f pctrcon.f pcgeequ.f pcgehd2.f pcgehrd.f
   pcgelq2.f pcgelqf.f  pcgels.f pcgeql2.f pcgeqlf.f pcgeqpf.f pcgeqr2.f 
   pcgeqrf.f pcgerfs.f pcgerq2.f pcgerqf.f                               
   pcgesv.f  pcgesvd.f pcgesvx.f pcgetf2.f pcgetrf.f pcgetri.f pcgetrs.f 
   pcggqrf.f 
   pcggrqf.f pcheev.f  pcheevd.f pcheevx.f pchegs2.f pchegst.f pchegvx.f 
   pchengst.f pchentrd.f pchettrd.f pchetd2.f 
   pchetrd.f pclabrd.f pclacon.f pclacgv.f pclacp2.f pclacpy.f pclahrd.f 
   pclahqr.f pclaconsb.f pclasmsub.f pclacp3.f pclawil.f
   pclange.f pclanhe.f pclanhs.f pclansy.f pclantr.f pclapiv.f pclapv2.f 
   pclaqge.f pclaqsy.f pclarf.f  pclarfb.f pclarfc.f pclarfg.f pclarft.f 
   pclascl.f pclase2.f pclaset.f pclassq.f pclaswp.f pclatra.f pclatrd.f 
   pclatrs.f pclauu2.f pclauum.f pcpocon.f pcpoequ.f pcporfs.f pcposv.f  
   pcpbsv.f  pcpbtrf.f pcpbtrs.f pcpbtrsv.f 
   pcptsv.f  pcpttrf.f pcpttrs.f pcpttrsv.f 
   pcposvx.f pcpotf2.f pcpotrf.f pcpotri.f pcpotrs.f pcsrscl.f pcstein.f 
   pctrevc.f pctrti2.f pctrtri.f pctrtrs.f pcung2l.f pcung2r.f 
             pcungl2.f pcunglq.f pcungql.f pcungqr.f pcungr2.f pcungrq.f 
             pcunm2l.f pcunm2r.f pcunmbr.f pcunmhr.f pcunml2.f pcunmlq.f 
   pcunmql.f pcunmqr.f pcunmr2.f pcunmrq.f pcunmtr.f 
   pclaevswp.f 
   pclarzb.f pclarzt.f pclarz.f  pclarzc.f pclatrz.f pctzrzf.f 
   pclattrs.f 
   pcunmr3.f pcunmrz.f pcmax1.f pscsum1.f pclamr1d.f 
   cdbtf2.f  cdbtrf.f  cdttrf.f cdttrsv.f  cpttrsv.f csteqr2.f ctrmvt.f 
   clamsh.f  claref.f  clanv2.f clahqr2.f pcheevr.f)

set (DLASRC  
   pddbsv.f  pddbtrf.f pddbtrs.f pddbtrsv.f 
   pddtsv.f  pddttrf.f pddttrs.f pddttrsv.f 
   pdgbsv.f  pdgbtrf.f pdgbtrs.f 
   pdgebd2.f pdgebrd.f pdgecon.f pdtrcon.f pdgeequ.f pdgehd2.f pdgehrd.f
   pdgelq2.f pdgelqf.f pdgels.f  pdgeql2.f pdgeqlf.f pdgeqpf.f pdgeqr2.f 
   pdgeqrf.f pdgerfs.f pdgerq2.f pdgerqf.f                               
   pdgesv.f  pdgesvd.f pdgesvx.f pdgetf2.f pdgetrf.f pdgetri.f pdgetrs.f 
   pdggqrf.f pdggrqf.f pdlabrd.f 
   pdlacon.f pdlacp2.f pdlacpy.f pdlahrd.f pdlange.f pdlanhs.f pdlansy.f 
   pdlantr.f pdlapiv.f pdlapv2.f pdlaqge.f pdlaqsy.f pdlarf.f  pdlarfb.f 
   pdlarfg.f pdlarft.f pdlase2.f pdlaset.f pdlascl.f pdlassq.f pdlaswp.f 
   pdlatra.f pdlatrd.f pdlatrs.f pdlauu2.f pdlauum.f pdorg2l.f pdorg2r.f 
                       pdorgl2.f pdorglq.f pdorgql.f pdorgqr.f pdorgr2.f 
   pdorgrq.f           pdorm2l.f pdorm2r.f pdormbr.f pdormhr.f pdorml2.f 
   pdormlq.f pdormql.f pdormqr.f pdormr2.f pdormrq.f pdormtr.f pdpocon.f 
   pdpbsv.f  pdpbtrf.f pdpbtrs.f pdpbtrsv.f 
   pdptsv.f  pdpttrf.f pdpttrs.f pdpttrsv.f 
   pdpoequ.f pdporfs.f pdposv.f  pdposvx.f pdpotf2.f pdpotrf.f pdpotri.f 
   pdpotrs.f pdrscl.f  pdstein.f pdsyev.f  pdsyevd.f pdsyevx.f 
   pdsygs2.f pdsygst.f 
   pdsygvx.f pdsyngst.f pdsyntrd.f pdsyttrd.f pdsytd2.f pdsytrd.f pdtrti2.f 
   pdtrtri.f pdtrtrs.f 
   pdlaevswp.f 
   pdlarzb.f pdlarzt.f pdlarz.f pdlatrz.f pdtzrzf.f pdormr3.f pdormrz.f 
   pdlahqr.f pdlaconsb.f pdlacp3.f pdlawil.f 
   pdlasmsub.f pdlared2d.f pdlamr1d.f dlaref.f dlamsh.f dlasorte.f dsteqr2.f 
   ddbtf2.f  ddbtrf.f  ddttrf.f ddttrsv.f dpttrsv.f dtrmvt.f  pdsyevr.f
   bdlaapp.f bdlaexc.f bdtrexc.f dlaqr6.f pdtrord.f
   pdtrsen.f pdgebal.f pdhseqr.f pdlamve.f pdlaqr0.f pdlaqr1.f pdlaqr2.f
   pdlaqr3.f pdlaqr4.f pdlaqr5.f pdrot.f)

set (ZLASRC  
   pzdbsv.f  pzdbtrf.f pzdbtrs.f pzdbtrsv.f 
   pzdtsv.f  pzdttrf.f pzdttrs.f pzdttrsv.f 
   pzgbsv.f  pzgbtrf.f pzgbtrs.f 
   pzgebd2.f pzgebrd.f pzgecon.f pztrcon.f pzgeequ.f pzgehd2.f pzgehrd.f
   pzgelq2.f pzgelqf.f pzgels.f  pzgeql2.f pzgeqlf.f pzgeqpf.f pzgeqr2.f 
   pzgeqrf.f pzgerfs.f pzgerq2.f pzgerqf.f                               
   pzgesv.f  pzgesvd.f pzgesvx.f pzgetf2.f pzgetrf.f pzgetri.f pzgetrs.f 
   pzggqrf.f 
   pzggrqf.f pzheev.f  pzheevd.f pzheevx.f pzhegs2.f pzhegst.f pzhegvx.f 
   pzhengst.f pzhentrd.f pzhettrd.f pzhetd2.f 
   pzhetrd.f pzlabrd.f pzlacon.f pzlacgv.f pzlacp2.f pzlacpy.f pzlahrd.f 
   pzlahqr.f pzlaconsb.f pzlasmsub.f pzlacp3.f pzlawil.f 
   pzlange.f pzlanhe.f pzlanhs.f pzlansy.f pzlantr.f pzlapiv.f pzlapv2.f 
   pzlaqge.f pzlaqsy.f pzlarf.f  pzlarfb.f pzlarfc.f pzlarfg.f pzlarft.f 
   pzlascl.f pzlase2.f pzlaset.f pzlassq.f pzlaswp.f pzlatra.f pzlatrd.f 
   pzlattrs.f 
   pzlatrs.f pzlauu2.f pzlauum.f pzpocon.f pzpoequ.f pzporfs.f pzposv.f  
   pzpbsv.f  pzpbtrf.f pzpbtrs.f pzpbtrsv.f 
   pzptsv.f  pzpttrf.f pzpttrs.f pzpttrsv.f 
   pzposvx.f pzpotf2.f pzpotrf.f pzpotri.f pzpotrs.f pzdrscl.f pzstein.f 
   pztrevc.f pztrti2.f pztrtri.f pztrtrs.f pzung2l.f pzung2r.f 
             pzungl2.f pzunglq.f pzungql.f pzungqr.f pzungr2.f pzungrq.f 
             pzunm2l.f pzunm2r.f pzunmbr.f pzunmhr.f pzunml2.f pzunmlq.f 
   pzunmql.f pzunmqr.f pzunmr2.f pzunmrq.f pzunmtr.f 
   pzlaevswp.f 
   pzlarzb.f pzlarzt.f pzlarz.f  pzlarzc.f pzlatrz.f pztzrzf.f 
   pzunmr3.f pzunmrz.f pzmax1.f pdzsum1.f pzlamr1d.f 
   zdbtf2.f  zdbtrf.f  zdttrf.f zdttrsv.f  zpttrsv.f zsteqr2.f ztrmvt.f 
   zlamsh.f  zlaref.f  zlanv2.f zlahqr2.f pzheevr.f)
   
set(src 
      ${ALLAUX} ${SCLAUX} ${DZLAUX}
      ${SLASRC} ${CLASRC} ${DLASRC} ${ZLASRC} )
set(src-C 
      ${ALLAUX-C} )
