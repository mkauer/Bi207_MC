*
* $Id: gxint.F,v 1.2 1997/01/07 10:25:42 cernlib Exp $
*
* $Log: gxint.F,v $
* Revision 1.2  1997/01/07 10:25:42  cernlib
* Remove #ifdef CERNLIB_MAIN; this shall be done via Imakefile.
*
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      PROGRAM GXINT
*
*     GEANT main program. To link with the MOTIF user interface
*     the routine GPAWPP(NWGEAN,NWPAW) should be called, whereas
*     the routine GPAW(NWGEAN,NWPAW) gives access to the basic
*     graphics version.
*
      PARAMETER (NWGEAN=3000000,NWPAW=1000000)
      COMMON/GCBANK/GEANT(NWGEAN)
      COMMON/PAWC/PAW(NWPAW)
*
      CALL GPAW(NWGEAN,NWPAW)
*
      END
      SUBROUTINE QNEXT
      END
      SUBROUTINE CZOPEN
      END
      SUBROUTINE CZTCP
      END
      SUBROUTINE CZCLOS
      END
      SUBROUTINE CZPUTA
      END
