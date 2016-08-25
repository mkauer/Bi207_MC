CDECK  ID>, GEOMPAR.
*-- +KEEP,HITPAR.
c      CHARACTER*4 NAMESH(7),CHNMSV1(3),CHNMSV2(2)
c      DIMENSION ORIGH(7),FACT(7),NBITSH(7),NBITSV1(3),NBITSv2(2)
*--
c      DATA NAMESH/'X   ','Y   ','Z   ',
c     +       'TOF ','E   ','IPRT','ELOS'/
c      DATA NBITSH/7*32/ , ORIGH/3*1000.,4*0./
c      DATA FACT/3*10.,100.,1.0E8,1.,1.0E8/
c      DATA CHNMSV1 /'CALO','SBC ','SBR '/
c      DATA NBITSV1 /16,16,16/
c      DATA CHNMSV2 /'SCTB','SBP '/
c      DATA NBITSV2 /16,16/


*--
*--       SENSITIVE DETECTOR SETS
*--       -----------------------
*--
*--   SCINTILLATORS   (DETECTOR SET 'EMCL')

*--   -------------
*--   1 X       X COORDINATE OF THE TRACK (M.R.S.)
*--             AT THE ENTRY INTO THE VOLUME OR AT CREATION POINT ( cm )
*--   2 Y       Y COORDINATE OF THE TRACK (M.R.S.)
*--             AT THE ENTRY INTO THE VOLUME OR AT CREATION POINT ( cm )
*--   3 Z       Z COORDINATE OF THE TRACK (M.R.S.)
*--             AT THE ENTRY INTO THE VOLUME OR AT CREATION POINT ( cm )
*--   8 TOF     T.O.F. OF THE TRACK AT THE ENTRY INTO THE VOLUME OR
*--                                 AT CREATION POINT ( ns )
*--   9 E       KINETIC ENERGY OF THE TRACK AT THE ENTRY INTO THE VOLUME OR
*--                                         AT CREATION POINT ( GeV )
*--   10  ELOS  ENERGY LOSS OF THE TRACK INSIDE THE VOLUME ( GeV )
*--             WARNING:  IF MORE THAN ONE PARTICLE DEPOSIT NENERGY
*--                       IN THE BLOCK, SOME OF THE INFORMATIONS MAY NOT
*--                       BE VALID (X,Y,Z,TOF,E)
*--                       IN THIS SITUATION ELOSS IS THE TOTAL ENERGY DEPOSITED
*--                       IN THE BLOCK AND TOF MAY NOT BE THE T.O.F. OF THE
*--                       FIRST PARTICLE HITING THE BLOCK!
*--

CDECK  ID>, SUGEOM.
	subroutine UGEOM

       REAL motvol(3)
       REAL scintvol(3)
       REAL alvol(3)
c       REAL vacvol(3)
       REAL esrvol(3)

* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
* gconst.inc
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
*
* gconst.inc
*

       real AARY(2)
       real ZARY(2)
       real WMAT(2)
      
*-- Mother dimentions
       DATA motvol /50.,50.,50./
*-- Aluminum dimentions
       DATA alvol /20.,20.,0.1/  ! 0.2 cm = 2 mm
*-- Scintillator dimensions
c       DATA scintvol /2.5,2.5,1./  ! 5x5x2 cm
       DATA scintvol /10.,10.,7.5/  ! 20x20x15 cm (simulate large hex block)
*-- Vacuum dimensions
c       DATA vacvol /1.,1.,1./
*-- Wrap dimensions
c       DATA esrvol /20.,20.,0.00225/  ! 0.0045 = 45 um ESR
       DATA esrvol /20.,20.,0.0006/   ! 0.0012 = 12 um Al Mylar
       
*-- Material mixture definition PURE POLYSTYRENE :  C6H5CH=CH2
      
       DATA INPAR /-2/ ! -n means numbers by atomic proportions (not weights)
       DATA AARY /1.00794,12.011/
       DATA ZARY /1.0, 6.0/
       DATA WMAT /8., 8./
       DATA DENS /1.032/

      DATA MDSC/1/
      DATA MDFOIL/2/
      DATA MDAIR/3/
      DATA MDAL/4/
c      DATA MDVAC/5/
      DATA MDESR/5/


c        REAL PAR(10)
        INTEGER IVOL


*-- Standard material definition	
	CALL GMATE
*-- My material definition
        CALL GSMIXT(21,'POLY',AARY,ZARY,DENS,INPAR,WMAT)

*-- Media definition
        FIELDM = 0.
        IFIELD = 0.
        STEMAX = 1.E10
        DEEMAX = 0.01
        TFMAX = 5.
        EPSIL = 0.00001
        STMIN = 0.0001

	call GSTMED(MDSC,'PLAST',21,1,IFIELD,FIELDM,TMAXFD,
     +          STEMAX,DEEMAX,EPSIL,STMIN,0,0)

	call GSTMED(MDAIR,'AIRMED',15,0,IFIELD,FIELDM,TMAXFD,
     +          STEMAX,DEEMAX,EPSIL,STMIN,0,0)

c        call GSTMED(MDVAC,'VACUUM',16,0,IFIELD,FIELDM,TMAXFD,
c     +          STEMAX,DEEMAX,EPSIL,STMIN,0,0)

        call GSTMED(MDAL,'ALUM',9,0,IFIELD,FIELDM,TMAXFD,
     +          STEMAX,DEEMAX,EPSIL,STMIN,0,0)
       
        call GSTMED(MDESR,'ESR',9,0,IFIELD,FIELDM,TMAXFD,
     +          STEMAX,DEEMAX,EPSIL,STMIN,0,0)
       

c-- Book the volumes
        call gsvolu ('AIRV','BOX ',MDAIR,motvol,3,IVOL)
        call gsvolu ('SCNT','BOX ',MDSC,scintvol,3,IVOL)
        call gsvolu ('ESRV','BOX ',MDESR,esrvol,3,IVOL)
c        call gsvolu ('ALVO','BOX ',MDAL,alvol,3,IVOL)
c        call gsvolu ('VACV','BOX ',MDVAC,vacvol,3,IVOL)

c-- Place the volumes: foil+sc_walls inside mother volume
        call gspos ('SCNT',1,'AIRV',0.0,0.0,0.0,0,ONLY)
        call gspos ('ESRV',1,'AIRV',0.0,0.0,7.8,0,ONLY)
c        call gspos ('ALVO',1,'AIRV',0.0,0.0,5.2,0,ONLY)
c        call gspos ('VACV',1,'AIRV',0.0,0.0,2.2,0,ONLY)

c-- Print parameters for volume + medium
        call gpvolu(0)
        call gptmed(0)

c-- close the geometry description
	call GGCLOS
	print *,'UGINIT finished'

	END
c*****************************************************
CDECK  ID>, PROBA.
c------------------------------------------------
CDECK  ID>, MAIN.
	Program SNOBB
c        PARAMETER (NGBANK=1000000,NHBOOK=300000)
	PARAMETER (NGBANK=100000000,NHBOOK=30000000)

	COMMON /GCBANK/Q(NGBANK)
	COMMON /PAWC/H(NHBOOK)

	call GZEBRA(NGBANK)

	call HLIMIT(-NHBOOK)


      CALL TIMEST(1000000.)




c	IEVNTU=1
	call UGINIT
	print *,'--> Init OK <--'
	call GRUN
	
	call UGLAST
	STOP
      
c     ===============================================  
      CLOSE(11) ! close the output.dat file
      CLOSE(12) ! close the eloss.dat file
c     ===============================================  	
   
      END
      
*-- End of main procedure
c------------------------------------------------	


      subroutine UGINIT
*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
      COMMON /GCRLUX/ LUX,INTRLUX,K1RLUX,K2RLUX
      INTEGER         LUX,INTRLUX,K1RLUX,K2RLUX
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gccuts.inc
      COMMON/GCCUTS/CUTGAM,CUTELE,CUTNEU,CUTHAD,CUTMUO,BCUTE,BCUTM
     +             ,DCUTE ,DCUTM ,PPCUTM,TOFMAX,GCUTS(5)
C
*
* gccuts.inc
*

         PARAMETER (lrecl=8192)
	 COMMON/QUEST/IQUEST(100)

*-- Initialization of GEANT 
	call GINIT
	
*--   IQUEST(10) IS THE MAXIMUM NUMBER OF RECORDS
         IQUEST(10)=100000000

*-- Files opening
c     ===============================================  
         OPEN(unit=11,file="output.dat",status="old")
         OPEN(unit=12,file="eloss.dat",status="old")
c     ===============================================        
         call hropen(1,'myfile','results.hbook','NQ',lrecl,istat)
         
         IF (ISTAT.ne.0) THEN
            print *, 
     +     'HBOOK direct access file opening error: ', ISTAT
         ENDIF	
*-- Book ntuple

       call hbookn(10,'ntuple2',2,' ',5000,chtags)

c     ===============================================
       call hbook1(100,'vme_sim',4000,0,4000,0)
c     ===============================================       

       DO i = 1,2
          nreg(i) = 0
          nlost(i) = 0
       ENDDO

       cutele = 1.e-5
       cutgam = 1.e-5


c--	CALL HCDIR('//PAWC',' ')


	call GZINIT
	
	call GPART
*-- User geometry definition
	call UGEOM
	
	call GPHYSI


	END
*----------------------------------------	
	
	subroutine GUKINE
	
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gcflag.inc
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
*
* gcflag.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gckine.inc
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
*
* gckine.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
* gconst.inc
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
*
* gconst.inc
*

*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
      DATA MDSC/1/
      DATA MDFOIL/2/
      DATA MDAIR/3/
      DATA MDAL/4/
c      DATA MDVAC/5/
      DATA MDESR/5/

c+SEQ, GCONST.
c*********  GENBB commonblocks
	common/genevent/tevst,npfull,npgeant(100),pmoment(3,100),
     +                  ptime(100)
	common/genbbpar/nevents,ievstart,irndmst,iwrfile,chfile
        common/currentev/icurrent
        character chdspin*4
	common/enrange/ebb1,ebb2,toallevents,levelE,chdspin

      DIMENSION UBUF(83),IHEAD(10)
      CHARACTER*8  CHNUCL
      LOGICAL OK
      CHARACTER*20 CHNPAR
c*********

	COMMON/FOILPAR/dxmo,dmo,dsmo,fpos
	DIMENSION RNDM(3),PAR(3),ATT(10)
	INTEGER NATT,NPAR
	CHARACTER*4 CHNAME

	COMMON/GENBBDATA/I2BBS,CHNUCLIDE,ILEVEL,MODEBB,ISTART
	INTEGER I2BBS,ILEVEL,MODEBB,ISTART
	CHARACTER*16 CHNUCLIDE

        REAL ctet
        REAL stet
        REAL phi

c-- Number of events to generate = 4 million
c     ===============================================
        nevent = 4000000
c     ===============================================

c-- init hits data
	npart=0
	nsc=0
	IEVNTU=IEVENT
	IF (MOD(IEVENT,100000).eq.0) print *,'Proceeding event N ',IEVENT
        IF(IEVNTU.EQ.1) THEN
*-- Init genbb routines
        IER=1
        DO WHILE (IER.NE.0)
c     ===============================================
c     BYPASS THE USER PROMT FOR INT/EXT AND ISOTOPE
c     ===============================================
c           CALL GENBBDIA(I2BBS,CHNUCLIDE,ILEVEL,MODEBB,ISTART)
           I2BBS=2
           CHNUCLIDE='Bi207'
c     ===============================================
           CALL GENBBSUB(I2BBS,CHNUCLIDE,ILEVEL,MODEBB,ISTART,IER)
           print *,'\n Using isotope ==> ',CHNUCLIDE,'\n'
        ENDDO

*-- Or produce next decay
      ELSE
         CALL GENBBSUB(I2BBS,CHNUCLIDE,ILEVEL,MODEBB,ISTART,IER)
      ENDIF

c     ===============================================
c     SOURCE POSITION
c     ===============================================
      PAR(1)=0.0
      PAR(2)=0.0
      PAR(3)=13.5 !scriptme
c     ===============================================

	   call GSVERT(PAR,0,0,0,0,NVERT)
	DO i=1,npfull
	   IKINE=npgeant(I)
*-- Save particle time as user word 1
	   PMOMENT(1,I)=PMOMENT(1,I)/1000.
	   PMOMENT(2,I)=PMOMENT(2,I)/1000.
	   PMOMENT(3,I)=PMOMENT(3,I)/1000.
	   CALL GSKINE(PMOMENT(1,I),IKINE,NVERT,PTIME(I),1,NT)
	   IF(IKINE.eq.0) THEN
	      PMS=0
	   ELSEIF(IKINE.eq.2.or.IKINE.eq.3) THEN
	      PMS=EMASS
	   ELSE
	      CALL GFPART(IKINE,CHNPAR,ITRTYP,PMS,CHARGE,TLIFE,UB,NWB)
	   ENDIF
	      npart=npart+1
	      IGPART(npart)=IKINE
	      GV(1,npart)=PAR(1)
	      GV(2,npart)=PAR(2)
	      GV(3,npart)=PAR(3)
	      CALL UCOPY(PMOMENT(1,I),GP(1,npart),3)
	      GE(npart)=SQRT(PMOMENT(1,I)**2+PMOMENT(2,I)**2+PMOMENT(3,I)**2+PMS**2)-PMS

	ENDDO

	IF(IDEBUG.ne.0) THEN
		call GPRINT('VERT',0)
		call GPRINT('KINE',0)
	ENDIF

c	IEVNTU=IEVNTU+1
	END

*-- Subroutine GUTREV
      subroutine GUTREV
c--  Zero ntuple before looping over steps
*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
       do i=1,8
           n_tuple(i)=0.
c           n_tuple2(i)=0.
        enddo
        n_tuple2(1)=0;
        n_tuple2(2)=0;
        call gtreve
        END

*-- Subroutine GUSTEP
*-- automatic call after each step

	subroutine GUSTEP
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gckine.inc
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
*
* gckine.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
* gctrak.inc
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
      PARAMETER (MAXME1=30)
      COMMON/GCTPOL/POLAR(3), NAMEC1(MAXME1)
C
*
* gctrak.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
* gctmed.inc
      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCKOV
C
*
* gctmed.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gcflag.inc
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
*
* gcflag.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
* gcsets.inc
      COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20)
C
*
* gcsets.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
* gcking.inc
* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
* gckmax.inc
      INTEGER MXGKIN
      PARAMETER (MXGKIN=100)
*
* gckmax.inc
*

      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
      INTEGER       KCASE,NGKINE ,IFLGK,MXPHOT,NGPHOT
      REAL          GKIN,TOFD,XPHOT
C
      PARAMETER (MXPHOT=800)
      COMMON/GCKIN2/NGPHOT,XPHOT(11,MXPHOT)
C
      COMMON/GCKIN3/GPOS(3,MXGKIN)
      REAL          GPOS
C
*
* gcking.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gcjloc.inc
      COMMON/GCJLOC/NJLOC(2),JTM,JMA,JLOSS,JPROB,JMIXT,JPHOT,JANNI
     +                  ,JCOMP,JBREM,JPAIR,JDRAY,JPFIS,JMUNU,JRAYL
     +                  ,JMULOF,JCOEF,JRANG
C
      INTEGER       NJLOC   ,JTM,JMA,JLOSS,JPROB,JMIXT,JPHOT,JANNI
     +                  ,JCOMP,JBREM,JPAIR,JDRAY,JPFIS,JMUNU,JRAYL
     +                  ,JMULOF,JCOEF,JRANG
C
      COMMON/GCJLCK/NJLCK(2),JTCKOV,JABSCO,JEFFIC,JINDEX,JCURIN
     +                      ,JPOLAR,JTSTRA,JTSTCO,JTSTEN,JTASHO
C
      EQUIVALENCE (JLASTV,JTSTEN)
C
      INTEGER       NJLCK,JTCKOV,JABSCO,JEFFIC,JINDEX,JCURIN
     +                   ,JPOLAR,JLASTV,JTSTRA,JTSTCO,JTSTEN
     +                   ,JTASHO
C
*
* gcjloc.inc
*

      DATA MDSC/1/
      DATA MDFOIL/2/
      DATA MDAIR/3/
      DATA MDAL/4/
c      DATA MDVAC/5/
      DATA MDESR/5/

*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
*---

      if(numed.eq.MDVAC)then
         if(IPART.eq.3) ISTOP=1
      endif

c--   Books the energy deposited or lost in each material it has gone through
        if(numed.eq.MDSC) n_tuple(1)=n_tuple(1)+destep ! Scintillator 
c       if(numed.eq.MDFOIL) n_tuple(2)=n_tuple(2)+destep ! Selenium82
	if(numed.eq.MDAIR) n_tuple(2)=n_tuple(2)+destep ! gas (air for now)

c--   Considering the secondary particles  

        if(ngkine.ne.0) then
              do i=1,ngkine
                 call GSKING(i)
              enddo
        endif	
c-- Debug each particle
c	call GPCXYZ	

	END
	
*-- subroutine guout
*-- User routine called at the end of each event
*-- outputs all the parameters
	subroutine GUOUT
* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
* gclist.inc
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
*
* gclist.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gcflag.inc
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
*
* gcflag.inc
*

*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
      DATA MDSC/1/
      DATA MDFOIL/2/
      DATA MDAIR/3/
      DATA MDAL/4/
c      DATA MDVAC/5/
      DATA MDESR/5/

*--
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
* gctrak.inc
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
      PARAMETER (MAXME1=30)
      COMMON/GCTPOL/POLAR(3), NAMEC1(MAXME1)
C
*
* gctrak.inc
*

* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gckine.inc
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
*
* gckine.inc
*


      
      REAL adc_value
      REAL res
      REAL npeout
      REAL resoutsquared
      REAL mean_pe
      INTEGER pois_n
      REAL gaus_n
      REAL smear_ph
      REAL pois_dist
      REAL eloss
      
      if(n_tuple(1).gt.0.0001) then ! if greater then 100 keV (n_tuple(1) in GeV)
         
c     ===============================================
         res = 16.0 !scriptme
c     ===============================================
         
         res = res/100.
         
         npeout = (2.354*2.354)/(res*res) ! 6% = 1539 photons
         n_tuple2(1) = n_tuple(1)*1000.0 ! converts from GeV to MeV
         mean_pe = n_tuple2(1)*npeout ! converts from MeV to photo electrons
         CALL RNPSSN(mean_pe,pois_n,IERR)
         CALL RNORML(gaus_n,1)
         smear_ph=(1.0*pois_n)+gaus_n
         n_tuple2(2) = (smear_ph/npeout)*1000.0 ! converts from photons to MeV to keV
         adc_value = n_tuple2(2) ! multiply by some factor to put into realistic adc format
         
         WRITE(11,*)adc_value
         
         
c     try storing the energy lost in the air
c     ===============================================
         eloss = n_tuple(2) * 1000.0 ! converts from GeV to MeV
         eloss = eloss * 1000.0 ! converts MeV to keV
         
c         WRITE(12,*)eloss
         
         
         
c     FASTER IF I DON'T FILL THESE? --> not really
c         call HFN(10,n_tuple2)
c         call hfill(100,adc_value,0.,1.)
         
         
      endif
      END
      
*-- User last routine	
	subroutine UGLAST
*-- My Data 
	PARAMETER (MAXSC=100,MAXPART=100)
	COMMON /HDATA/istat,icycle,numeropart,momentpart,
     &                nlost,nreg,IEVNTU 
	COMMON /NTDATA/n_tuple
	COMMON /GDATA/npart,IGPART(MAXPART),GE(MAXPART),
     &    GP(3,MAXPART),GV(3,MAXPART)


	INTEGER IPRT,NSC,IW,ICLN,IROW,IEVNTU
	INTEGER npart,IGPAR
	REAL GE
        REAL GP
        REAL GV
	REAL E
        REAL T

	INTEGER istat,icycle,numeropart 
	DIMENSION nlost(2),nreg(2)
	REAL momentpart
        REAL n_tuple(10)
        REAL n_tuple2(2)
        CHARACTER*10 chtags(2)

        DATA chtags /'MeV','keV'/
c        DATA chtags /'Esc','Egas','Xc','Yc','Zc','Dist',
c     +       'Xvert','Yvert','Zvert','1000pe'/
c     DATA chtags /'Esc','Egas'/
c------------------------------------------------
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
* gcflag.inc
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
*
* gcflag.inc
*

	print *,' enter UGLAST'

        call hrout(0,icycle,' ')
        call hrend('myfile')
        close(1)                  

        call glast
        
c        DO i = 10,19
           write(*,*) 'lost events ----', nlost(2)
           write(*,*) 'registrated events ----', nreg(2)
           EF=real(nreg(2))/real(nevent)
           write(*,*) 'EFFICIENCY --- ', EF    
c        ENDDO

	END
*--

      SUBROUTINE ZGAUSS(SIGMA,XM,X)
      EXTERNAL GRNDM
      CALL RNORMX(R,1,GRNDM)
      X=XM+R*SIGMA
      RETURN
      END
*__________________________________________________________
CDECK  ID>, GRNDMNEMO.
      SUBROUTINE GRNDM(RVEC,LEN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       To generate a vector RVECV of LEN random numbers         *
C.    *         using the CERN Library routine RANLUX                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, many GEANT routines                  *
C.    *                                                                *
C.    ******************************************************************
C.
        CALL RANLUX(RVEC,LEN)
      RETURN
      END
CDECK  ID>, GRNDMQNEMO.
      SUBROUTINE GRNDMQ(IS1,IS2,ISEQ,CHOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       To set/retrieve the seed of the random number generator  *
C.    *                                                                *
C.    *  CHOPT  - A character specifying the action which              *
C.    *           GRNDMQ should take.                                  *
C.    *    'S'  - To initialyse random number generator                *
C.    *    'Q'  - To initialyse seed                                   *
C.    *    'G'  - To get status                                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, many GEANT routines                  *
C.    *       Author    V.I.Tretyak          *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
      COMMON /GCRLUX/ LUX,INTRLUX,K1RLUX,K2RLUX
      INTEGER         LUX,INTRLUX,K1RLUX,K2RLUX
      INTEGER IS1,IS2,ISEQ
      CHARACTER*(*) CHOPT
      CHARACTER*12  CCHOPT

      CCHOPT = CHOPT

      IF(INDEX(CHOPT,'S').NE.0) THEN
c        print *,' ISEQ=',ISEQ
c        IF(ISEQ.EQ.1) THEN
c	  LUX = IS1
c	  INTRLUX = IS2
c        END IF
	print *,'LUX=',LUX,' INTRLUX =',INTRLUX,
     &	' K1RLUX=',K1RLUX,' K2RLUX=',K2RLUX
        CALL RLUXGO(LUX,INTRLUX,K1RLUX,K2RLUX)
      END IF
	
      IF(INDEX(CHOPT,'Q').NE.0) THEN
        IF(ISEQ.LE.4 .AND. ISEQ.GT.0) THEN
	  LUX = ISEQ
	ELSE  
	  LUX = 3
	ENDIF
	
c	INTRLUX = seed_from_time()
C        INTRLUX = time()
	print *,' Q  LUX=',LUX,' INTRLUX =',INTRLUX,
     &	' K1RLUX=',K1RLUX,' K2RLUX=',K2RLUX
             
      END IF

      IF(INDEX(CHOPT,'G').NE.0) THEN
      CALL RLUXAT(LUX,INTRLUX,K1RLUX,K2RLUX)
        IS1 = K1RLUX
        IS2 = K2RLUX
      END IF
C
      END
*------------------------------------------------------------------------------
