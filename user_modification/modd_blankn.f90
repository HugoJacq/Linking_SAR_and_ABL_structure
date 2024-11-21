!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_BLANK_n
!     #################
!
!!****  *MODD_BLANK$n* -  Declarative module for MesoNH developpers namelist
!!
!!    PURPOSE
!!    -------
!!
!!      Offer dummy real, integer, logical and character variables for
!!    test and debugging purposes.
!!
!!**  METHOD
!!    ------
!!
!!      Eight dummy real, integer, logical and character*80 variables are
!!    defined and passed through the namelist read operations. None of the
!!    MesoNH routines uses any of those variables. When a developper choses
!!    to introduce temporarily a parameter to some subroutine, he has to
!!    introduce a USE MODD_BLANK statement into that subroutine. Then he
!!    can use any of the variables defined here and change them easily via
!!    the namelist input.
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!      K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    Original 25/04/96
!!      updated     17/11/00  (P Jabouille) Use dummy array
!!      updated     26/10/21  (Q.Rodier) Use for n model (grid-nesting)
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPDUMMY, JPMODELMAX
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE BLANK_t
!
  LOGICAL :: LDUMMY1
  LOGICAL :: LDUMMY2
  LOGICAL :: LDUMMY3
  LOGICAL :: LDUMMY4
  LOGICAL :: LDUMMY5
  LOGICAL :: LDUMMY6
  LOGICAL :: LDUMMY7
  LOGICAL :: LDUMMY8
!
  CHARACTER(len=80) :: CDUMMY1
  CHARACTER(len=80) :: CDUMMY2
  CHARACTER(len=80) :: CDUMMY3
  CHARACTER(len=80) :: CDUMMY4
  CHARACTER(len=80) :: CDUMMY5
  CHARACTER(len=80) :: CDUMMY6
  CHARACTER(len=80) :: CDUMMY7
  CHARACTER(len=80) :: CDUMMY8
!
  INTEGER :: NDUMMY1
  INTEGER :: NDUMMY2
  INTEGER :: NDUMMY3
  INTEGER :: NDUMMY4
  INTEGER :: NDUMMY5
  INTEGER :: NDUMMY6
  INTEGER :: NDUMMY7
  INTEGER :: NDUMMY8
!
  REAL :: XDUMMY1
  REAL :: XDUMMY2
  REAL :: XDUMMY3
  REAL :: XDUMMY4
  REAL :: XDUMMY5
  REAL :: XDUMMY6
  REAL :: XDUMMY7
  REAL :: XDUMMY8
!
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY1 ! used for THW_FLX in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY2 ! used for LM (mixing length) in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY3 ! used for UW_HFLX in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY4 ! used for RCONSW_FLX in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY5 ! used for UV_FLX in OUT files  
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY6 ! used for U_VAR in OUT files 
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY7 ! used for V_VAR in OUT files 
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY8 ! used for W_VAR in OUT files 
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY9 ! used for VW_HFLX in OUT files 
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY10 ! used for UW_VFLX in OUT files 
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY11 ! used for VW_VFLX in OUT files
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY12 ! used for SVT in OUT files. NOT USED
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY13 ! used for WSVT_VFLX in OUT files
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY14 ! used for USVT_VFLX in OUT files
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY15 ! used for VSVT_VFLX in OUT files
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY16 ! used for sv2_hvar in OUT files
  REAL, DIMENSION(:,:,:,:), POINTER :: TDUMMY17 ! used for sv2_vvar in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY18 ! used for Rt2_vvar in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY19 ! used for Rt2_hvar in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY20 ! used for THl2_vvar in OUT files
  REAL, DIMENSION(:,:,:), POINTER :: TDUMMY21 ! used for THl2_hvar in OUT files
!	Note : the sum of vertical and horizontal contribution is saved in OUT files : R_VAR = R_VVAR + R_HVAR
!	Note 2 : in mnh code, R_VVAR is called RTOT_VVAR. Mistake ?
!
END TYPE BLANK_t
!
TYPE(BLANK_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: BLANK_MODEL
!
LOGICAL, POINTER :: LDUMMY1=>NULL()
LOGICAL, POINTER :: LDUMMY2=>NULL()
LOGICAL, POINTER :: LDUMMY3=>NULL()
LOGICAL, POINTER :: LDUMMY4=>NULL()
LOGICAL, POINTER :: LDUMMY5=>NULL()
LOGICAL, POINTER :: LDUMMY6=>NULL()
LOGICAL, POINTER :: LDUMMY7=>NULL()
LOGICAL, POINTER :: LDUMMY8=>NULL()
!
CHARACTER(len=80), POINTER :: CDUMMY1=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY2=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY3=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY4=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY5=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY6=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY7=>NULL()
CHARACTER(len=80), POINTER :: CDUMMY8=>NULL()
!
INTEGER, POINTER :: NDUMMY1=>NULL()
INTEGER, POINTER :: NDUMMY2=>NULL()
INTEGER, POINTER :: NDUMMY3=>NULL()
INTEGER, POINTER :: NDUMMY4=>NULL()
INTEGER, POINTER :: NDUMMY5=>NULL()
INTEGER, POINTER :: NDUMMY6=>NULL()
INTEGER, POINTER :: NDUMMY7=>NULL()
INTEGER, POINTER :: NDUMMY8=>NULL()
!
REAL, POINTER :: XDUMMY1=>NULL()
REAL, POINTER :: XDUMMY2=>NULL()
REAL, POINTER :: XDUMMY3=>NULL()
REAL, POINTER :: XDUMMY4=>NULL()
REAL, POINTER :: XDUMMY5=>NULL()
REAL, POINTER :: XDUMMY6=>NULL()
REAL, POINTER :: XDUMMY7=>NULL()
REAL, POINTER :: XDUMMY8=>NULL()
!
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY1=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY2=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY3=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY4=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY5=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY6=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY7=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY8=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY9=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY10=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY11=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY12=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY13=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY14=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY15=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY16=>NULL()
REAL, DIMENSION(:,:,:,:),   POINTER :: TDUMMY17=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY18=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY19=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY20=>NULL()
REAL, DIMENSION(:,:,:),   POINTER :: TDUMMY21=>NULL()

!
CONTAINS
!
SUBROUTINE BLANK_GOTO_MODEL(KFROM,KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!    
LDUMMY1=>BLANK_MODEL(KTO)%LDUMMY1
LDUMMY2=>BLANK_MODEL(KTO)%LDUMMY2
LDUMMY3=>BLANK_MODEL(KTO)%LDUMMY3
LDUMMY4=>BLANK_MODEL(KTO)%LDUMMY4
LDUMMY5=>BLANK_MODEL(KTO)%LDUMMY5
LDUMMY6=>BLANK_MODEL(KTO)%LDUMMY6
LDUMMY7=>BLANK_MODEL(KTO)%LDUMMY7
LDUMMY8=>BLANK_MODEL(KTO)%LDUMMY8

CDUMMY1=>BLANK_MODEL(KTO)%CDUMMY1
CDUMMY2=>BLANK_MODEL(KTO)%CDUMMY2
CDUMMY3=>BLANK_MODEL(KTO)%CDUMMY3
CDUMMY4=>BLANK_MODEL(KTO)%CDUMMY4
CDUMMY5=>BLANK_MODEL(KTO)%CDUMMY5
CDUMMY6=>BLANK_MODEL(KTO)%CDUMMY6
CDUMMY7=>BLANK_MODEL(KTO)%CDUMMY7
CDUMMY8=>BLANK_MODEL(KTO)%CDUMMY8
!
NDUMMY1=>BLANK_MODEL(KTO)%NDUMMY1
NDUMMY2=>BLANK_MODEL(KTO)%NDUMMY2
NDUMMY3=>BLANK_MODEL(KTO)%NDUMMY3
NDUMMY4=>BLANK_MODEL(KTO)%NDUMMY4
NDUMMY5=>BLANK_MODEL(KTO)%NDUMMY5
NDUMMY6=>BLANK_MODEL(KTO)%NDUMMY6
NDUMMY7=>BLANK_MODEL(KTO)%NDUMMY7
NDUMMY8=>BLANK_MODEL(KTO)%NDUMMY8
!
XDUMMY1=>BLANK_MODEL(KTO)%XDUMMY1
XDUMMY2=>BLANK_MODEL(KTO)%XDUMMY2
XDUMMY3=>BLANK_MODEL(KTO)%XDUMMY3
XDUMMY4=>BLANK_MODEL(KTO)%XDUMMY4
XDUMMY5=>BLANK_MODEL(KTO)%XDUMMY5
XDUMMY6=>BLANK_MODEL(KTO)%XDUMMY6
XDUMMY7=>BLANK_MODEL(KTO)%XDUMMY7
XDUMMY8=>BLANK_MODEL(KTO)%XDUMMY8
!
BLANK_MODEL(KFROM)%TDUMMY1=>TDUMMY1
BLANK_MODEL(KFROM)%TDUMMY2=>TDUMMY2
BLANK_MODEL(KFROM)%TDUMMY3=>TDUMMY3
BLANK_MODEL(KFROM)%TDUMMY4=>TDUMMY4
BLANK_MODEL(KFROM)%TDUMMY5=>TDUMMY5
BLANK_MODEL(KFROM)%TDUMMY6=>TDUMMY6
BLANK_MODEL(KFROM)%TDUMMY7=>TDUMMY7
BLANK_MODEL(KFROM)%TDUMMY8=>TDUMMY8
BLANK_MODEL(KFROM)%TDUMMY9=>TDUMMY9
BLANK_MODEL(KFROM)%TDUMMY10=>TDUMMY10
BLANK_MODEL(KFROM)%TDUMMY11=>TDUMMY11
BLANK_MODEL(KFROM)%TDUMMY12=>TDUMMY12
BLANK_MODEL(KFROM)%TDUMMY13=>TDUMMY13
BLANK_MODEL(KFROM)%TDUMMY14=>TDUMMY14
BLANK_MODEL(KFROM)%TDUMMY15=>TDUMMY15
BLANK_MODEL(KFROM)%TDUMMY16=>TDUMMY16
BLANK_MODEL(KFROM)%TDUMMY17=>TDUMMY17
BLANK_MODEL(KFROM)%TDUMMY18=>TDUMMY18
BLANK_MODEL(KFROM)%TDUMMY19=>TDUMMY19
BLANK_MODEL(KFROM)%TDUMMY20=>TDUMMY20
BLANK_MODEL(KFROM)%TDUMMY21=>TDUMMY21
!
TDUMMY1=>BLANK_MODEL(KTO)%TDUMMY1 ! used for THW_FLX in OUT files
TDUMMY2=>BLANK_MODEL(KTO)%TDUMMY2 ! used for LM (mixing length) in OUT files
TDUMMY3=>BLANK_MODEL(KTO)%TDUMMY3 ! used for UW_HFLX in OUT files
TDUMMY4=>BLANK_MODEL(KTO)%TDUMMY4 ! used for RCONSW_FLX in OUT files
TDUMMY5=>BLANK_MODEL(KTO)%TDUMMY5 ! used for UV_FLX in OUT files  
TDUMMY6=>BLANK_MODEL(KTO)%TDUMMY6 ! used for U_VAR in OUT files 
TDUMMY7=>BLANK_MODEL(KTO)%TDUMMY7 ! used for V_VAR in OUT files 
TDUMMY8=>BLANK_MODEL(KTO)%TDUMMY8 ! used for W_VAR in OUT files 
TDUMMY9=>BLANK_MODEL(KTO)%TDUMMY9 ! used for VW_HFLX in OUT files 
TDUMMY10=>BLANK_MODEL(KTO)%TDUMMY10 ! used for UW_VFLX in OUT files 
TDUMMY11=>BLANK_MODEL(KTO)%TDUMMY11 ! used for VW_VFLX in OUT files
TDUMMY12=>BLANK_MODEL(KTO)%TDUMMY12
TDUMMY13=>BLANK_MODEL(KTO)%TDUMMY13 
TDUMMY14=>BLANK_MODEL(KTO)%TDUMMY14 
TDUMMY15=>BLANK_MODEL(KTO)%TDUMMY15 
TDUMMY16=>BLANK_MODEL(KTO)%TDUMMY16 
TDUMMY17=>BLANK_MODEL(KTO)%TDUMMY17 
TDUMMY18=>BLANK_MODEL(KTO)%TDUMMY18 
TDUMMY19=>BLANK_MODEL(KTO)%TDUMMY19
TDUMMY20=>BLANK_MODEL(KTO)%TDUMMY20 
TDUMMY21=>BLANK_MODEL(KTO)%TDUMMY21 
!
END SUBROUTINE BLANK_GOTO_MODEL
!
END MODULE MODD_BLANK_n
