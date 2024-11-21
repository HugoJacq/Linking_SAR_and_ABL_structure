!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_CONDSAMP
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_CONDSAMP  Declaration of conditional sampling tracers 
!
!!    AUTHOR
!!    ------
!	           : C.Lac                               
!	Creation   : 01/06/2011
!
!!    MODIFICATIONS
!!    -------------
!!      F.Brient                  * Tracer emission from the top
!!                                   of the boundary-layer * 05/2019
!-------------------------------------------------------------------------------
!
!
!*    0. DECLARATIONS
!        ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
LOGICAL            :: LCONDSAMP = .FALSE.  ! Switch to activate conditional sampling
!
INTEGER, PARAMETER :: JPCSMAX = 4     
!
INTEGER                            :: NCONDSAMP            ! Number of conditional
                                                           ! sampling tracers
REAL,         DIMENSION(JPCSMAX)   :: XRADIO               ! Radioactive decay period
REAL,         DIMENSION(JPCSMAX)   :: XSCAL                ! Scaling factor
REAL                               :: XHEIGHT_BASE         ! Distance below the
                              !         cloud base where the 2nd tracer is emitted
REAL                               :: XDEPTH_BASE          ! Depth in which the
                              !         2nd tracer is emitted
REAL                               :: XHEIGHT_TOP          ! Distance above the
                              !         cloud top  where the 3rd tracer is emitted
REAL                               :: XDEPTH_TOP           ! Depth in which the
                              !         3rd tracer is emitted
INTEGER                            :: NFINDTOP             ! Options for
                              !         the method for identifying the altitude 
                              !         where the 3rd tracer is emitted
REAL                               :: XTHVP                ! Threshold for 
                              !         identifying the PBL top based on virtual
                              !         potential temperature (IF NFINDTOP==2)
LOGICAL                            :: LTPLUS               ! Options for
                              !         allowing an emission of tracers one layer
                              !         below the cloud base and one level above
                              !         the PBL top (when the layers of emission 
                              !         are not detected)
LOGICAL				    :: LBISURF	! Options to use 2 tracer at surface        
CHARACTER(LEN=4)		    :: CBISURF	! Type of emission when LBISURF=TRUE 
REAL						:: XSST_CRIT ! Critical SST provided by user
!
END MODULE MODD_CONDSAMP
