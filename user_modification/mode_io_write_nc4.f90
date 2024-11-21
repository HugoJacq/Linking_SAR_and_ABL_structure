!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet may 2016  : use NetCDF Fortran module
!  J. Escobar  14/12/2017: correction for MNH_INT=8
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/12/2018: split of mode_netcdf into multiple modules/files
!  P. Wautelet 10/01/2019: replace handle_err by IO_Err_handle_nc4 for better netCDF error messages
!  P. Wautelet 11/01/2019: NVERB_INFO->NVERB_WARNING for zero size fields
!  P. Wautelet 01/02/2019: IO_Coordvar_write_nc4: bug: use of non-associated pointers (PIOCDF%DIM_Nx_y)
!  P. Wautelet 05/03/2019: rename IO subroutines and modules
!  P. Wautelet 12/07/2019: add support for 1D array of dates
!  P. Wautelet 10/09/2019: IO_Coordvar_write_nc4: split communication and file write operations
!                          + no more process coordination for Z-split files
!  P. Wautelet 18/09/2019: correct support of 64bit integers (MNH_INT=8)
!  P. Wautelet 19/09/2019: temporary workaround for netCDF bug if MNH_INT=8 (if netCDF fortran < 4.4.5)
!  P. Wautelet 11/02/2020: add 'dims' attribute in IO_Write_field_header_split_nc4
!  P. Wautelet 25/06/2020: remove workaround for netCDF bug (see 19/09/2019)
!  P. Wautelet 14/09/2020: IO_Coordvar_write_nc4: do not store 'time' coordinate in diachronic files
!  P. Wautelet 22/09/2020: add ldimreduced field to allow reduction in the number of dimensions of fields (used by 2D simulations)
!  P. Wautelet 10/11/2020: new data structures for netCDF dimensions
!  P. Wautelet 26/11/2020: add IO_Field_create_nc4 subroutine + use it for all IO_Field_write_nc4_*
!  P. Wautelet 04/12/2020: add IO_Field_partial_write_nc4 subroutines
!  P. Wautelet 11/01/2021: add coordinates for dimension variables in diachronic files
!  P. Wautelet 14/01/2021: add IO_Field_write_nc4_N4, IO_Field_partial_write_nc4_N2,
!                          IO_Field_partial_write_nc4_N3 and IO_Field_partial_write_nc4_N4 subroutines
!  P. Wautelet 30/03/2021: budgets: LES cartesian subdomain limits are defined in the physical domain
!  P. Wautelet 22/03/2022: correct time_les_avg and time_les_avg_bounds coordinates
!  P. Wautelet    06/2022: reorganize flyers
!  P. Wautelet 21/06/2022: bugfix: time_budget was not computed correctly (tdtexp -> tdtseg)
!  P. Wautelet 13/01/2023: IO_Coordvar_write_nc4: add optional dummy argument TPDTMODELN to force written model time
!  P. Wautelet 14/12/2023: add lossy compression for output files
!-----------------------------------------------------------------
#ifdef MNH_IOCDF4
module mode_io_write_nc4

use modd_field,        only: tfieldmetadata
use modd_io,           only: gsmonoproc, tfiledata
use modd_parameters,   only: NMNHNAMELGTMAX
use modd_precision,    only: CDFINT, MNHINT_NF90, MNHREAL32, MNHREAL_MPI, MNHREAL_NF90

use mode_io_tools_nc4, only: IO_Mnhname_clean, IO_Vdims_fill_nc4, IO_Dim_find_create_nc4, IO_Strdimid_get_nc4, IO_Err_handle_nc4
use mode_msg

use NETCDF,            only: NF90_CHAR, NF90_FLOAT, NF90_INT1,                          &
                             NF90_GLOBAL, NF90_NOERR,                                   &
                             NF90_DEF_VAR, NF90_DEF_VAR_DEFLATE, NF90_DEF_VAR_QUANTIZE, &
                             NF90_GET_ATT, NF90_INQ_VARID,                              &
                             NF90_INQUIRE_ATTRIBUTE, NF90_PUT_ATT, NF90_PUT_VAR

implicit none

private

public :: IO_Coordvar_write_nc4, IO_Header_write_nc4, IO_Field_header_split_write_nc4
public :: IO_Field_create_nc4, IO_Field_write_nc4, IO_Field_partial_write_nc4

INTERFACE IO_Field_write_nc4
   MODULE PROCEDURE IO_Field_write_nc4_X0,IO_Field_write_nc4_X1, &
                    IO_Field_write_nc4_X2,IO_Field_write_nc4_X3, &
                    IO_Field_write_nc4_X4,IO_Field_write_nc4_X5, &
                    IO_Field_write_nc4_X6,                       &
                    IO_Field_write_nc4_N0,IO_Field_write_nc4_N1, &
                    IO_Field_write_nc4_N2,IO_Field_write_nc4_N3, &
                    IO_Field_write_nc4_N4,                       &
                    IO_Field_write_nc4_L0,IO_Field_write_nc4_L1, &
                    IO_Field_write_nc4_C0,IO_Field_write_nc4_C1, &
                    IO_Field_write_nc4_T0,IO_Field_write_nc4_T1
END INTERFACE IO_Field_write_nc4

interface IO_Field_partial_write_nc4
   module procedure IO_Field_partial_write_nc4_X1, IO_Field_partial_write_nc4_X2, &
                    IO_Field_partial_write_nc4_X3, IO_Field_partial_write_nc4_X4, &
                    IO_Field_partial_write_nc4_N2, IO_Field_partial_write_nc4_N3, &
                    IO_Field_partial_write_nc4_N4
end interface IO_Field_partial_write_nc4

integer,parameter :: NSTRINGCHUNKSIZE = 16 !Dimension of the chunks of strings
                                           !(to limit the number of dimensions for strings)

integer(kind=CDFINT),parameter :: SHUFFLE = 1 !Set to 1 for (usually) better compression
integer(kind=CDFINT),parameter :: DEFLATE = 1

contains

subroutine IO_Field_header_split_write_nc4( tpfile, tpfield, knblocks )
use modd_field,      only: TYPEREAL
use modd_parameters, only: jphext

use mode_tools_ll,   only: Get_globaldims_ll

type(tfiledata),       intent(in) :: tpfile
class(tfieldmetadata), intent(in) :: tpfield
integer,               intent(in) :: knblocks

character(len=len(tpfield%cmnhname))  :: yvarname
integer                               :: iimax, ijmax
integer(kind=CDFINT)              :: istatus
integer(kind=CDFINT)              :: incid
integer(kind=CDFINT)              :: ivarid
integer(kind=CDFINT),dimension(3) :: ishape

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_header_split_write_nc4', 'called for field '//trim( tpfield%cmnhname ) )

if ( tpfield%ntype /= TYPEREAL ) then
  call Print_msg( NVERB_ERROR, 'IO', 'IO_Field_header_split_write_nc4', 'invalid ntype for field '//trim( tpfield%cmnhname ) )
  return
end if

! Get the Netcdf file ID
incid = tpfile%nncid

call IO_Mnhname_clean( tpfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( incid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then

  if ( tpfile%lncreduce_float_precision ) then
    istatus = NF90_DEF_VAR( incid, yvarname, NF90_FLOAT, ivarid )
  else
    istatus = NF90_DEF_VAR( incid, yvarname, MNHREAL_NF90, ivarid )
  end if

  if ( tpfield%ndims /= 3 ) call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_header_split_write_nc4', &
                  trim( tpfile%cname )//': '//trim( yvarname )//': NDIMS should be 3' )

  if ( tpfield%cdir /= 'XY' ) call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_header_split_write_nc4', &
                  trim( tpfile%cname )//': '//trim( yvarname )//': CDIR should be XY' )

  call Get_globaldims_ll( iimax, ijmax )
  ishape(1) = int( iimax + 2 * jphext, kind = CDFINT )
  ishape(2) = int( ijmax + 2 * jphext, kind = CDFINT )
  ishape(3) = knblocks
  call IO_Field_attr_write_nc4( tpfile, tpfield, ivarid, .false., kshape = ishape )

  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_DEF_VAR', trim(yvarname) )

  istatus = NF90_PUT_ATT( incid, ivarid,'split_variable', 'yes')
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'split_variable for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'split_mode', 'Z')
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'split_mode for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'split_nblocks', knblocks )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'split_nblocks for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'split_nfiles', tpfile%nsubfiles_ioz )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'split_nfiles for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'split_distribution', 'round-robin' )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'split_distribution for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'ndims', tpfield%ndims )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'ndims for '//trim( tpfield%cmnhname ) )

  istatus = NF90_PUT_ATT( incid, ivarid,'dims', ishape )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'dims for '//trim( tpfield%cmnhname ) )

  if ( tpfield%ltimedep ) then
    istatus = NF90_PUT_ATT( incid, ivarid,'time_dependent', 'yes' )
  else
    istatus = NF90_PUT_ATT( incid, ivarid,'time_dependent', 'no' )
  end if
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_header_split_write_nc4', 'NF90_PUT_ATT', &
                                                     'time_dependent for '//trim( tpfield%cmnhname ) )
else
  call Print_msg( NVERB_WARNING, 'IO', 'IO_Field_header_split_write_nc4', &
                  trim( tpfile%cname )//': '//trim( yvarname )//' already defined' )
end if

end subroutine IO_Field_header_split_write_nc4

SUBROUTINE IO_Field_attr_write_nc4(TPFILE,TPFIELD,KVARID,OEXISTED,KSHAPE,HCALENDAR,OISCOORD)
!
USE MODD_CONF,   ONLY: CPROGRAM, LCARTESIAN
USE MODD_CONF_n, ONLY: CSTORAGE_TYPE
use modd_field,  only: NMNHDIM_ARAKAWA, TYPEINT, TYPEREAL
!
TYPE(TFILEDATA),                              INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),                        INTENT(IN) :: TPFIELD
INTEGER(KIND=CDFINT),                         INTENT(IN) :: KVARID
LOGICAL,                                      INTENT(IN) :: OEXISTED !True if variable was already defined
INTEGER(KIND=CDFINT), DIMENSION(:), OPTIONAL, INTENT(IN) :: KSHAPE
CHARACTER(LEN=*),                   OPTIONAL, INTENT(IN) :: HCALENDAR
LOGICAL,                            OPTIONAL, INTENT(IN) :: OISCOORD   ! Is a coordinate variable (->do not write coordinates attribute)
!
CHARACTER(LEN=:),                   ALLOCATABLE :: YCOORDS
INTEGER(KIND=CDFINT)                            :: INCID
INTEGER(KIND=CDFINT)                            :: istatus
LOGICAL                                         :: GISCOORD
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','called for field '//TRIM(TPFIELD%CMNHNAME))
!
IF(LEN_TRIM(TPFIELD%CSTDNAME)==0 .AND. LEN_TRIM(TPFIELD%CLONGNAME)==0) THEN
  CALL PRINT_MSG(NVERB_WARNING,'IO','IO_Field_attr_write_nc4','at least long_name or standard_name must be provided &
  &to respect CF-convention for variable '//TRIM(TPFIELD%CMNHNAME))
ENDIF
!
IF (TPFIELD%NDIMS>1 .AND. .NOT.PRESENT(KSHAPE)) &
  CALL PRINT_MSG(NVERB_FATAL,'IO','IO_Field_attr_write_nc4','KSHAPE not provided for '//TRIM(TPFIELD%CMNHNAME))
!
IF (PRESENT(OISCOORD)) THEN
  GISCOORD = OISCOORD
ELSE
  GISCOORD = .FALSE.
END IF
!
INCID = TPFILE%NNCID
!
! Standard_name attribute definition (CF convention)
IF(LEN_TRIM(TPFIELD%CSTDNAME)==0) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','TPFIELD%CSTDNAME not set for variable '//TRIM(TPFIELD%CMNHNAME))
ELSE
  istatus = NF90_PUT_ATT(INCID, KVARID,'standard_name', TRIM(TPFIELD%CSTDNAME))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','standard_name for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! Long_name attribute definition (CF convention)
IF(LEN_TRIM(TPFIELD%CLONGNAME)==0) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','TPFIELD%CLONGNAME not set for variable '//TRIM(TPFIELD%CMNHNAME))
ELSE
  istatus = NF90_PUT_ATT(INCID, KVARID,'long_name', TRIM(TPFIELD%CLONGNAME))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','long_name for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! Canonical units attribute definition (CF convention)
IF(LEN_TRIM(TPFIELD%CUNITS)==0) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','TPFIELD%CUNITS not set for variable '//TRIM(TPFIELD%CMNHNAME))
ELSE
  istatus = NF90_PUT_ATT(INCID, KVARID,'units', TRIM(TPFIELD%CUNITS))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','units for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! GRID attribute definition
IF(TPFIELD%NGRID<0) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','TPFIELD%NGRID not set for variable '//TRIM(TPFIELD%CMNHNAME))
!Do not write GRID attribute if NGRID=0
ELSE IF (TPFIELD%NGRID>0) THEN
  istatus = NF90_PUT_ATT(INCID, KVARID, 'grid', TPFIELD%NGRID)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','grid for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! COMMENT attribute definition
IF(LEN_TRIM(TPFIELD%CCOMMENT)==0) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','TPFIELD%CCOMMENT not set for variable '//TRIM(TPFIELD%CMNHNAME))
ELSE
  istatus = NF90_PUT_ATT(INCID, KVARID,'comment', TRIM(TPFIELD%CCOMMENT))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','comment for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! Calendar (CF convention)
IF(PRESENT(HCALENDAR)) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','CALENDAR provided for variable '//TRIM(TPFIELD%CMNHNAME))
  istatus = NF90_PUT_ATT(INCID, KVARID,'calendar', TRIM(HCALENDAR))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','calendar for ' &
                                                   //trim(TPFIELD%CMNHNAME))
ENDIF
!
! Coordinates (CF convention)
IF (.NOT.GISCOORD) THEN
  !0D: nothing to do
  !1D: no direct correspondance with latitude(_x)/longitude(_x) 2D variables => nothing to do
  IF (.NOT.LCARTESIAN .AND. TPFIELD%NDIMS>1 .AND. TPFIELD%NGRID/=0) THEN
    IF (TPFIELD%CDIR=='XY') THEN
      if (       kshape(1) == tpfile%tncdims%tdims( NMNHDIM_ARAKAWA(tpfield%ngrid,1) )%nlen &
           .and. kshape(2) == tpfile%tncdims%tdims( NMNHDIM_ARAKAWA(tpfield%ngrid,2) )%nlen ) then
        SELECT CASE(TPFIELD%NGRID)
          CASE (0) !Not on Arakawa grid
            !Nothing to do
          CASE (1) !Mass point
            YCOORDS='latitude longitude'
          CASE (2) !u point
            YCOORDS='latitude_u longitude_u'
          CASE (3) !v point
            YCOORDS='latitude_v longitude_v'
          CASE (4) !w point
            YCOORDS='latitude longitude'
          CASE (5) !xi vorticity point (=f point =uv point)
            YCOORDS='latitude_f longitude_f'
          CASE (6) !eta vorticity point (=uw point)
            YCOORDS='latitude_u longitude_u'
          CASE (7) !zeta vorticity point (=vw point)
            YCOORDS='latitude_v longitude_v'
          CASE (8) !fw point (=uvw point)
            YCOORDS='latitude_f longitude_f'
          CASE DEFAULT
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_Field_attr_write_nc4','invalid NGRID for variable '//TRIM(TPFIELD%CMNHNAME))
        END SELECT
        !
        istatus = NF90_PUT_ATT(INCID, KVARID,'coordinates',YCOORDS)
        IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','coordinates')
        DEALLOCATE(YCOORDS)
      ELSE
        CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_attr_write_nc4','coordinates not implemented for variable ' &
                                                                    //TRIM(TPFIELD%CMNHNAME))
      END IF
    ELSE
      !No YCOORDS for CDIR/='XY'
    END IF
  END IF
ENDIF
!
IF(TPFIELD%NTYPE==TYPEINT .AND. TPFIELD%NDIMS>0) THEN
  IF (TPFIELD%NFILLVALUE>=TPFIELD%NVALIDMIN .AND. TPFIELD%NFILLVALUE<=TPFIELD%NVALIDMAX) &
    CALL PRINT_MSG(NVERB_WARNING,'IO','IO_Field_attr_write_nc4','_FillValue is not outside of valid_min - valid_max'// &
                                                                'interval for variable '//TRIM(TPFIELD%CMNHNAME))
  !
  ! Fillvalue (CF/COMODO convention)
  ! Remarks: * the attribute '_FillValue' is also recognized by the netCDF library
  !            and is used when pre-filling a variable
  !          * it cannot be modified if some data has already been written (->check OEXISTED)
  !BUG: NF90_PUT_ATT does not work for NF90_INT64 and _FillValue attribute if netCDF-fortran version < 4.4.5 (bug in netCDF)
  !     (see https://github.com/Unidata/netcdf-fortran/issues/62)
  IF(.NOT.OEXISTED) THEN
    istatus = NF90_PUT_ATT(INCID, KVARID,'_FillValue', TPFIELD%NFILLVALUE)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','_FillValue')
  END IF
  !
  ! Valid_min/max (CF/COMODO convention)
  istatus = NF90_PUT_ATT(INCID, KVARID,'valid_min', TPFIELD%NVALIDMIN)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','valid_min')
  !
  istatus = NF90_PUT_ATT(INCID, KVARID,'valid_max',TPFIELD%NVALIDMAX)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','valid_max')
ENDIF
!
IF(TPFIELD%NTYPE==TYPEREAL .AND. TPFIELD%NDIMS>0) THEN
  IF (TPFIELD%XFILLVALUE>=TPFIELD%XVALIDMIN .AND. TPFIELD%XFILLVALUE<=TPFIELD%XVALIDMAX) &
    CALL PRINT_MSG(NVERB_WARNING,'IO','IO_Field_attr_write_nc4','_FillValue is not outside of valid_min - valid_max'// &
                                                                'interval for variable '//TRIM(TPFIELD%CMNHNAME))
  !
  ! Fillvalue (CF/COMODO convention)
  ! Remarks: * the attribute '_FillValue' is also recognized by the netCDF library
  !            and is used when pre-filling a variable
  !          * it cannot be modified if some data has already been written (->check OEXISTED)
  IF(.NOT.OEXISTED) THEN
    IF (TPFILE%LNCREDUCE_FLOAT_PRECISION) THEN
      istatus = NF90_PUT_ATT(INCID, KVARID,'_FillValue', REAL(TPFIELD%XFILLVALUE,KIND=MNHREAL32))
    ELSE
      istatus = NF90_PUT_ATT(INCID, KVARID,'_FillValue', TPFIELD%XFILLVALUE)
    END IF
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','_FillValue')
  END IF
  !
  ! Valid_min/max (CF/COMODO convention)
  IF (TPFILE%LNCREDUCE_FLOAT_PRECISION) THEN
    istatus = NF90_PUT_ATT(INCID, KVARID,'valid_min', REAL(TPFIELD%XVALIDMIN,KIND=MNHREAL32))
  ELSE
    istatus = NF90_PUT_ATT(INCID, KVARID,'valid_min', TPFIELD%XVALIDMIN)
  END IF
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','valid_min')
  !
  IF (TPFILE%LNCREDUCE_FLOAT_PRECISION) THEN
    istatus = NF90_PUT_ATT(INCID, KVARID,'valid_max', REAL(TPFIELD%XVALIDMAX,KIND=MNHREAL32))
  ELSE
    istatus = NF90_PUT_ATT(INCID, KVARID,'valid_max',TPFIELD%XVALIDMAX)
  END IF
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_attr_write_nc4','NF90_PUT_ATT','valid_max')
ENDIF
!
END SUBROUTINE IO_Field_attr_write_nc4


subroutine IO_Field_create_nc4( tpfile, tpfield, kshape, hcalendar, oiscoord, kvarid, oisempty )

use modd_field,     only: NMNHDIM_TIME, TYPECHAR, TYPEDATE, TYPEINT, TYPELOG, TYPEREAL, TYPEUNDEF
use modd_precision, only: MNHINT_NF90, MNHREAL_NF90

type(tfiledata),       intent(in)            :: tpfile
class(tfieldmetadata), intent(in)            :: tpfield
integer, dimension(:), intent(in),  optional :: kshape
character(len=*),      intent(in),  optional :: hcalendar
logical,               intent(in),  optional :: oiscoord   ! Is a coordinate variable (->do not write coordinates attribute)
integer(kind=CDFINT),  intent(out), optional :: kvarid
logical,               intent(out), optional :: oisempty

character(len = Len( tpfield%cmnhname ))        :: yvarname
integer                                         :: idims
integer                                         :: idx
integer                                         :: ji
integer(kind=CDFINT)                            :: istatus
integer(kind=CDFINT)                            :: itype
integer(kind=CDFINT)                            :: ivarid
integer(kind=CDFINT), dimension(:), allocatable :: ivdims
integer(kind=CDFINT), dimension(:), allocatable :: ivdimstmp
integer(kind=CDFINT), dimension(:), allocatable :: ishape
logical                                         :: gexisted !True if variable was already defined

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_create_nc4', Trim( tpfile%cname ) // ': creating ' // Trim( tpfield%cmnhname ) )

gexisted = .false.

call IO_Mnhname_clean( tpfield%cmnhname, yvarname )

if ( .not. Present( kshape ) .and. tpfield%ndims > 0 ) then
  !kshape not provided => ndimlist has to be previously populated
  idims = tpfield%ndims
  if ( tpfield%ltimedep ) idims = idims + 1
  Allocate( ivdims(idims) )
  Allocate( ishape(idims) )

  do ji = 1, tpfield%ndims
    idx = tpfield%ndimlist(ji)
    if ( idx > tpfile%tncdims%nmaxdims .or. idx < 1 )                                        &
      call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_create_nc4', Trim( tpfile%cname ) // ': ' &
                      // Trim( tpfield%cmnhname ) // ': invalid ndimlist'  )
    ivdims(ji) = tpfile%tncdims%tdims(idx)%nid
    ishape(ji) = tpfile%tncdims%tdims(idx)%nlen
  end do

  !Set the last dimension if variable is time dependent
  if ( tpfield%ltimedep ) then
    ivdims(idims) = tpfile%tncdims%tdims(NMNHDIM_TIME)%nid
    ishape(idims) = 1
  end if
else
  !Guess the dimensions from the shape of the field
  if ( tpfield%ndims == 0 ) then
    Allocate( ishape(0) )
  else
    Allocate( ishape(Size( kshape )) )
    ishape(:) = kshape(:)
  end if

  !Get the netCDF dimensions
  if ( tpfield%ntype /= TYPECHAR ) then
    call IO_Vdims_fill_nc4( tpfile, tpfield, ishape, ivdims )
  else
    if ( tpfield%ndims == 0 ) then
      idims = 1
    else if ( tpfield%ndims == 1 ) then
      idims = 2
    else
      call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_create_nc4', Trim( tpfile%cname ) // ': ' &
                      // Trim( tpfield%cmnhname ) // ': ndims value not supported for character strings'  )
    end if

    Allocate( ivdims(idims) )
    ivdims(1) = IO_Strdimid_get_nc4( tpfile, Int( kshape(1), kind=CDFINT ) )

    if ( idims == 2 ) then
      call IO_Vdims_fill_nc4( tpfile, tpfield, [ Int( kshape(2), kind=CDFINT ) ], ivdimstmp )
      ivdims(2) = ivdimstmp(1)
      Deallocate( ivdimstmp )
    end if
  end if
end if

if ( Present( oisempty ) ) oisempty = .false.
do ji = 1, Size( ishape )
  if ( ishape(ji) == 0 ) then
    if ( Present( oisempty ) ) then
      call Print_msg( NVERB_WARNING, 'IO', 'IO_Field_create_nc4','ignoring variable with a zero size (' // Trim( yvarname ) // ')' )
      oisempty = .true.
      return
    else
      call Print_msg( NVERB_ERROR, 'IO', 'IO_Field_create_nc4','variable with a zero size (' // Trim( yvarname ) // ')' )
    end if
  end if
end do

! The variable should not already exist but who knows ?
istatus = NF90_INQ_VARID( tpfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  select case( tpfield%ntype )
    case ( TYPEINT )
      itype = MNHINT_NF90

    case ( TYPELOG )
      itype = NF90_INT1

    case ( TYPEREAL )
      if ( tpfile%lncreduce_float_precision .and. tpfield%ndims >= 1  ) then
        itype = NF90_FLOAT
      else
        itype = MNHREAL_NF90
      end if

    case ( TYPECHAR )
      itype = NF90_CHAR

    case ( TYPEDATE )
      itype = MNHREAL_NF90
      if ( .not. Present( hcalendar ) ) &
        call Print_msg( NVERB_ERROR, 'IO', 'IO_Field_create_nc4', Trim( tpfield%cmnhname ) // ': hcalendar not provided' )

    case ( TYPEUNDEF )
      call Print_msg( NVERB_ERROR, 'IO', 'IO_Field_create_nc4', Trim( tpfield%cmnhname ) // ': ntype is TYPEUNDEF' )
      return

    case default
      call Print_msg( NVERB_ERROR, 'IO', 'IO_Field_create_nc4', Trim( tpfield%cmnhname ) // ': invalid ntype' )
      return
  end select

  ! Define the variable
  istatus = NF90_DEF_VAR( tpfile%nncid, yvarname, itype, ivdims, ivarid )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_create_nc4', 'NF90_DEF_VAR', Trim( yvarname ) )

  ! Add compression if asked for
  if ( tpfile%lnccompress .and. tpfield%ntype == TYPEREAL .and. tpfield%ndims >= 1 ) then
    istatus = NF90_DEF_VAR_DEFLATE( tpfile%nncid, ivarid, SHUFFLE, DEFLATE, tpfile%nnccompress_level )
    if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_create_nc4', 'NF90_DEF_VAR_DEFLATE', Trim( yvarname ) )

    if ( tpfile%lnccompress_lossy ) then
      istatus = NF90_DEF_VAR_QUANTIZE( tpfile%nncid, ivarid, tpfile%nnccompress_lossy_algo, tpfile%nnccompress_lossy_nsd )
      if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'IO_Field_create_nc4', &
                                                           'NF90_DEF_VAR_QUANTIZE', Trim( yvarname ) )
    end if
  end if
else
  gexisted = .true.
  call Print_msg( NVERB_WARNING, 'IO', 'IO_Field_create_nc4', Trim( tpfile%cname ) // ': ' &
                  // Trim( yvarname ) // ' already defined' )
end if

! Write metadata
call IO_Field_attr_write_nc4( tpfile, tpfield, ivarid, gexisted, kshape = ishape, hcalendar = hcalendar, oiscoord = oiscoord )

if ( Present( kvarid ) ) kvarid = ivarid

end subroutine IO_Field_create_nc4


SUBROUTINE IO_Field_write_nc4_X0(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
REAL,                  INTENT(IN) :: PFIELD
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0

call IO_Field_create_nc4( tpfile, tpfield, kvarid = ivarid )

! Write the data
istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X0','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)

END SUBROUTINE IO_Field_write_nc4_X0


SUBROUTINE IO_Field_write_nc4_X1(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),TARGET,INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
REAL,DIMENSION(:),     INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( pfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X1','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_X1


SUBROUTINE IO_Field_write_nc4_X2(TPFILE,TPFIELD,PFIELD,KRESP,KVERTLEVEL,KZFILE,OISCOORD)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
REAL,DIMENSION(:,:),   INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,               INTENT(OUT):: KRESP
INTEGER,OPTIONAL,      INTENT(IN) :: KVERTLEVEL ! Number of the vertical level (needed for Z-level split files)
INTEGER,OPTIONAL,      INTENT(IN) :: KZFILE     ! Number of the Z-level split file
LOGICAL,OPTIONAL,      INTENT(IN) :: OISCOORD   ! Is a coordinate variable (->do not write coordinates attribute)
!
CHARACTER(LEN=4)               :: YSUFFIX
INTEGER(KIND=CDFINT)           :: istatus
INTEGER(KIND=CDFINT)           :: IVARID
logical                        :: gisempty
CLASS(TFIELDMETADATA), pointer :: TZFIELD
TYPE(TFILEDATA),       POINTER :: TZFILE
!
KRESP = 0
!
call IO_Select_split_file( tpfile, tpfield, tzfile, tzfield, kvertlevel, kzfile )
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X2',TRIM(TZFILE%CNAME)//': writing '//TRIM(TZFIELD%CMNHNAME))
!
call IO_Field_create_nc4( tzfile, tzfield, kshape = Shape( pfield ), oiscoord = oiscoord, kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TZFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X2','NF90_PUT_VAR',trim(TZFIELD%CMNHNAME),KRESP)
end if

if ( Present( kvertlevel ) ) deallocate( tzfield )

END SUBROUTINE IO_Field_write_nc4_X2


SUBROUTINE IO_Field_write_nc4_X3(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
REAL,DIMENSION(:,:,:), INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X3',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( pfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X3','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_X3


SUBROUTINE IO_Field_write_nc4_X4(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),           INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),     INTENT(IN) :: TPFIELD
REAL,DIMENSION(:,:,:,:),   INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                   INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X4',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( pfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X4','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_X4


SUBROUTINE IO_Field_write_nc4_X5(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),           INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),     INTENT(IN) :: TPFIELD
REAL,DIMENSION(:,:,:,:,:), INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                   INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X5',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( pfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X5','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_X5


SUBROUTINE IO_Field_write_nc4_X6(TPFILE,TPFIELD,PFIELD,KRESP)
!
TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),       INTENT(IN) :: TPFIELD
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                     INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_X6',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( pfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, PFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_X6','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_X6


SUBROUTINE IO_Field_write_nc4_N0(TPFILE,TPFIELD,KFIELD,KRESP)
!
#if 0
use modd_field,          only: NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_LEVEL
USE MODD_IO,             ONLY: LPACK, L2D
USE MODD_PARAMETERS_ll,  ONLY: JPHEXT, JPVEXT
#else
use modd_field,          only: NMNHDIM_LEVEL
USE MODD_PARAMETERS_ll,  ONLY: JPVEXT
#endif
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
INTEGER,               INTENT(IN) :: KFIELD
INTEGER,               INTENT(OUT):: KRESP
!
integer              :: iidx
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_N0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kvarid = ivarid )

! Write the data
istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, KFIELD)
IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_N0','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
!
! Use IMAX, JMAX, KMAX to define DIM_NI, DIM_NJ, DIM_LEVEL
! /!\ Can only work if IMAX, JMAX or KMAX are written before any array
!
#if 0
if ( tpfield%cmnhname == 'IMAX' .and. tpfile%tncdims%tdims(NMNHDIM_NI)%nid == -1 ) then
  call IO_Dim_find_create_nc4( tpfile, kfield + 2 * jphext, iidx, 'X' )
end if
if ( tpfield%cmnhname == 'JMAX' .and. tpfile%tncdims%tdims(NMNHDIM_NJ)%nid == -1 ) then
  if ( lpack .and. l2d ) then
    call IO_Dim_find_create_nc4( tpfile, 1,                   iidx, 'Y' )
  else
    call IO_Dim_find_create_nc4( tpfile, kfield + 2 * jphext, iidx, 'Z' )
  end if
end if
#endif
if ( tpfield%cmnhname == 'KMAX' .and. tpfile%tncdims%tdims(NMNHDIM_LEVEL)%nid == -1 ) then
  call IO_Dim_find_create_nc4( tpfile, Int( kfield + 2 * JPVEXT, kind = CDFINT ), iidx, 'Z' )
end if

END SUBROUTINE IO_Field_write_nc4_N0


SUBROUTINE IO_Field_write_nc4_N1(TPFILE,TPFIELD,KFIELD,KRESP)
!
#if 0
USE MODD_IO,             ONLY: LPACK, L2D
USE MODD_PARAMETERS_ll,  ONLY: JPHEXT, JPVEXT
#else
USE MODD_PARAMETERS_ll,  ONLY: JPVEXT
#endif
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
INTEGER, DIMENSION(:), INTENT(IN) :: KFIELD
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_N1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( kfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, KFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_N1','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_N1


SUBROUTINE IO_Field_write_nc4_N2(TPFILE,TPFIELD,KFIELD,KRESP)
!
TYPE(TFILEDATA),TARGET,INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
INTEGER,DIMENSION(:,:),INTENT(IN) :: KFIELD   ! array containing the data field
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_N2',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( kfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, KFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_N2','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_N2


SUBROUTINE IO_Field_write_nc4_N3(TPFILE,TPFIELD,KFIELD,KRESP)
!
TYPE(TFILEDATA),TARGET,  INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),   INTENT(IN) :: TPFIELD
INTEGER,DIMENSION(:,:,:),INTENT(IN) :: KFIELD   ! array containing the data field
INTEGER,                 INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_N3',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( kfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, KFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_N3','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_N3


SUBROUTINE IO_Field_write_nc4_N4(TPFILE,TPFIELD,KFIELD,KRESP)
!
TYPE(TFILEDATA),TARGET,    INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),     INTENT(IN) :: TPFIELD
INTEGER,DIMENSION(:,:,:,:),INTENT(IN) :: KFIELD   ! array containing the data field
INTEGER,                   INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_N4',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( kfield ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, KFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_N4','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_N4


SUBROUTINE IO_Field_write_nc4_L0(TPFILE,TPFIELD,OFIELD,KRESP)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
LOGICAL,               INTENT(IN) :: OFIELD
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER              :: IFIELD
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_L0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kvarid = ivarid )

!Convert LOGICAL to INTEGER (LOGICAL format not supported by netCDF files)
IF (OFIELD) THEN
  IFIELD = 1
ELSE
  IFIELD = 0
END IF

! Write the data
istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, IFIELD)
IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_L0','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)

END SUBROUTINE IO_Field_write_nc4_L0


SUBROUTINE IO_Field_write_nc4_L1(TPFILE,TPFIELD,OFIELD,KRESP)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
LOGICAL, DIMENSION(:), INTENT(IN) :: OFIELD
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER, DIMENSION(SIZE(OFIELD)) :: IFIELD
INTEGER(KIND=CDFINT)             :: istatus
INTEGER(KIND=CDFINT)             :: IVARID
logical                          :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_L1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
call IO_Field_create_nc4( tpfile, tpfield, kshape = Shape( ofield ), kvarid = ivarid, oisempty = gisempty )

!Convert LOGICAL to INTEGER (LOGICAL format not supported by netCDF files)
WHERE (OFIELD)
  IFIELD = 1
ELSEWHERE
  IFIELD = 0
END WHERE

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, IFIELD)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_L1','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_L1


SUBROUTINE IO_Field_write_nc4_C0(TPFILE,TPFIELD,HFIELD,KRESP)
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
CHARACTER(LEN=*),      INTENT(IN) :: HFIELD
INTEGER,               INTENT(OUT):: KRESP
!
INTEGER(KIND=CDFINT)          :: istatus
INTEGER(KIND=CDFINT)          :: IVARID
INTEGER                       :: ILEN
CHARACTER(LEN=:), ALLOCATABLE :: YFIELD
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_C0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0

!Store the character string in a string of a size multiple of NSTRINGCHUNKSIZE
!This is done to limit the number of dimensions in the netCDF file
ILEN = ((LEN_TRIM(HFIELD)+NSTRINGCHUNKSIZE-1)/NSTRINGCHUNKSIZE)*NSTRINGCHUNKSIZE
!If the string is empty, create it anyway with a non-zero size (to prevent problems later)
IF (ILEN==0) ILEN = NSTRINGCHUNKSIZE

!PW: si ndimlist populated, pas besoin de passer kshape...
call IO_Field_create_nc4( tpfile, tpfield, kshape = [ ilen ], kvarid = ivarid )

ALLOCATE(CHARACTER(LEN=ILEN)::YFIELD)
YFIELD(1:LEN_TRIM(HFIELD))=TRIM(HFIELD)
YFIELD(LEN_TRIM(HFIELD)+1:)=' '

! Write the data
istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, YFIELD)
IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_C0','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)

END SUBROUTINE IO_Field_write_nc4_C0


SUBROUTINE IO_Field_write_nc4_C1(TPFILE,TPFIELD,HFIELD,KRESP)
!  Modif
!    J.Escobar : 25/04/2018 : missing 'IF ALLOCATED(IVDIMSTMP)' DEALLOCATE
!----------------------------------------------------------------
TYPE(TFILEDATA),              INTENT(IN)  :: TPFILE
CLASS(TFIELDMETADATA),        INTENT(IN) :: TPFIELD
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN)  :: HFIELD
INTEGER,                      INTENT(OUT) :: KRESP
!
INTEGER(KIND=CDFINT),PARAMETER :: IONE = 1
!
INTEGER(KIND=CDFINT) :: istatus
INTEGER(KIND=CDFINT) :: IVARID
INTEGER(KIND=CDFINT) :: ILEN, ISIZE
logical              :: gisempty
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_C1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0

ILEN  = LEN(HFIELD)
ISIZE = SIZE(HFIELD)

call IO_Field_create_nc4( tpfile, tpfield, kshape = Int ([ ilen, isize ], kind = Kind( 1 ) ), kvarid = ivarid, oisempty = gisempty )

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, HFIELD(1:ISIZE)(1:ILEN), START=(/IONE,IONE/), COUNT=(/ILEN,ISIZE/))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_C1','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_C1


SUBROUTINE IO_Field_write_nc4_T0(TPFILE,TPFIELD,TPDATA,KRESP)
!
USE MODD_TIME_n,     ONLY: TDTMOD
USE MODD_TYPE_DATE
!
USE MODE_DATETIME
!
TYPE(TFILEDATA),       INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA), INTENT(IN) :: TPFIELD
TYPE (DATE_TIME),      INTENT(IN) :: TPDATA
INTEGER,               INTENT(OUT):: KRESP
!
CHARACTER(LEN=40)                  :: YUNITS
INTEGER(KIND=CDFINT)               :: istatus
INTEGER(KIND=CDFINT)               :: IVARID
REAL                               :: ZDELTATIME !Distance in seconds since reference date and time
CLASS(TFIELDMETADATA), ALLOCATABLE :: TZFIELD
TYPE(DATE_TIME)                    :: TZREF
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_T0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
Allocate( TZFIELD, source = TPFIELD )
!
! Model beginning date (TDTMOD%TDATE) is used as the reference date
! Reference time is set to 0.
IF (.NOT.ASSOCIATED(TDTMOD)) THEN
  CALL PRINT_MSG(NVERB_WARNING,'IO','IO_Field_write_nc4_T0',TRIM(TPFILE%CNAME)// &
                 ': '//TRIM(TZFIELD%CMNHNAME)//': DTMOD is not associated and not known. Reference date set to 2000/01/01')
  TZREF%nyear  = 2000
  TZREF%nmonth = 1
  TZREF%nday   = 1
  TZREF%xtime  = 0.
ELSE
  TZREF = TDTMOD
  TZREF%xtime = 0.
END IF
WRITE(YUNITS,'( "seconds since ",I4.4,"-",I2.2,"-",I2.2," 00:00:00 +0:00" )') &
      TZREF%nyear, TZREF%nmonth, TZREF%nday
TZFIELD%CUNITS = TRIM(YUNITS)
!
call IO_Field_create_nc4( tpfile, tzfield, kvarid = ivarid, hcalendar = 'standard' )
!
! Compute the temporal distance from reference
CALL DATETIME_DISTANCE(TZREF,TPDATA,ZDELTATIME)

! Write the data
istatus = NF90_PUT_VAR(TPFILE%NNCID, IVARID, ZDELTATIME)
IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_T0','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)

END SUBROUTINE IO_Field_write_nc4_T0


SUBROUTINE IO_Field_write_nc4_T1(TPFILE,TPFIELD,TPDATA,KRESP)
!
USE MODD_TIME_n,     ONLY: TDTMOD
USE MODD_TYPE_DATE
!
USE MODE_DATETIME
!
TYPE(TFILEDATA),                INTENT(IN) :: TPFILE
CLASS(TFIELDMETADATA),          INTENT(IN) :: TPFIELD
TYPE (DATE_TIME), DIMENSION(:), INTENT(IN) :: TPDATA
INTEGER,                        INTENT(OUT):: KRESP
!
CHARACTER(LEN=40)                  :: YUNITS
INTEGER                            :: JI
INTEGER(KIND=CDFINT)               :: istatus
INTEGER(KIND=CDFINT)               :: IVARID
logical                            :: gisempty
REAL, DIMENSION(:),    ALLOCATABLE :: ZDELTATIME !Distance in seconds since reference date and time
CLASS(TFIELDMETADATA), ALLOCATABLE :: TZFIELD
TYPE(DATE_TIME)                    :: TZREF
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Field_write_nc4_T1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
!
KRESP = 0
!
Allocate( TZFIELD, source = TPFIELD )
!
! Model beginning date (TDTMOD%TDATE) is used as the reference date
! Reference time is set to 0.
IF (.NOT.ASSOCIATED(TDTMOD)) THEN
  CALL PRINT_MSG(NVERB_WARNING,'IO','IO_Field_write_nc4_T1',TRIM(TPFILE%CNAME)// &
                 ': '//TRIM(TZFIELD%CMNHNAME)//': DTMOD is not associated and not known. Reference date set to 2000/01/01')
  TZREF%nyear  = 2000
  TZREF%nmonth = 1
  TZREF%nday   = 1
  TZREF%xtime  = 0.
ELSE
  TZREF = TDTMOD
  TZREF%xtime = 0.
END IF
WRITE(YUNITS,'( "seconds since ",I4.4,"-",I2.2,"-",I2.2," 00:00:00 +0:00" )') &
      TZREF%nyear, TZREF%nmonth, TZREF%nday
TZFIELD%CUNITS = TRIM(YUNITS)
!
call IO_Field_create_nc4( tpfile, tzfield, kshape = Shape( tpdata), kvarid = ivarid, hcalendar = 'standard', oisempty = gisempty )
!
! Compute the temporal distances from reference
ALLOCATE( ZDELTATIME( SIZE( TPDATA ) ) )

DO JI = 1, SIZE( TPDATA )
  CALL DATETIME_DISTANCE( TZREF, TPDATA(JI ), ZDELTATIME(JI) )
END DO

! Write the data
if ( .not. gisempty ) then
  istatus = NF90_PUT_VAR( TPFILE%NNCID, IVARID, ZDELTATIME(:) )
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_Field_write_nc4_T1','NF90_PUT_VAR',trim(TPFIELD%CMNHNAME),KRESP)
end if

END SUBROUTINE IO_Field_write_nc4_T1


subroutine IO_Field_partial_write_nc4_X1( tpfile, tpfield, pfield, koffset, kresp )

type(tfiledata),                intent(in)  :: tpfile
class(tfieldmetadata),          intent(in)  :: tpfield
real,             dimension(:), intent(in)  :: pfield   ! array containing the data field
integer,          dimension(1), intent(in)  :: koffset
integer,                        intent(out) :: kresp

character(len=NMNHNAMELGTMAX)      :: yvarname
integer(kind=CDFINT)               :: istatus
integer(kind=CDFINT)               :: ivarid
integer(kind=CDFINT), dimension(1) :: istarts

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_X1',&
                Trim( tpfile%cname ) // ': writing ' // Trim( tpfield%cmnhname ) )

kresp = 0

call IO_Mnhname_clean( tpfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tpfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_X1', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( pfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tpfile%nncid, ivarid, pfield(:), start = istarts(:), count = Int( Shape( pfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_X1', 'NF90_PUT_VAR', Trim( tpfield%cmnhname ), kresp )
end if

end subroutine IO_Field_partial_write_nc4_X1


subroutine IO_Field_partial_write_nc4_X2( tpfile, tpfield, pfield, koffset, kresp, kvertlevel, kzfile )

type(tfiledata),                  intent(in)  :: tpfile
class(tfieldmetadata),            intent(in)  :: tpfield
real,             dimension(:,:), intent(in)  :: pfield   ! array containing the data field
integer,          dimension(2),   intent(in)  :: koffset
integer,                          intent(out) :: kresp
integer,                optional, intent(in)  :: kvertlevel ! Number of the vertical level (needed for Z-level split files)
integer,                optional, intent(in)  :: kzfile     ! Number of the Z-level split file

character(len=4)                    :: ysuffix
character(len=NMNHNAMELGTMAX)       :: yvarname
integer(kind=CDFINT)                :: istatus
integer(kind=CDFINT)                :: ivarid
integer(kind=CDFINT),  dimension(2) :: istarts
class(tfieldmetadata), pointer      :: tzfield
type(tfiledata),       pointer      :: tzfile

kresp = 0

call IO_Select_split_file( tpfile, tpfield, tzfile, tzfield, kvertlevel, kzfile )

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_X2',&
                Trim( tzfile%cname ) // ': writing ' // Trim( tzfield%cmnhname ) )

call IO_Mnhname_clean( tzfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tzfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_X2', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( pfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tzfile%nncid, ivarid, pfield(:,:), start = istarts(:), count = Int( Shape( pfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_X2', 'NF90_PUT_VAR', Trim( tzfield%cmnhname ), kresp )
end if

if ( Present( kvertlevel ) ) deallocate( tzfield )

end subroutine IO_Field_partial_write_nc4_X2


subroutine IO_Field_partial_write_nc4_X3( tpfile, tpfield, pfield, koffset, kresp )

type(tfiledata),                    intent(in)  :: tpfile
class(tfieldmetadata),              intent(in)  :: tpfield
real,             dimension(:,:,:), intent(in)  :: pfield   ! array containing the data field
integer,          dimension(3),     intent(in)  :: koffset
integer,                            intent(out) :: kresp

character(len=NMNHNAMELGTMAX)      :: yvarname
integer(kind=CDFINT)               :: istatus
integer(kind=CDFINT)               :: ivarid
integer(kind=CDFINT), dimension(3) :: istarts

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_X3',&
                Trim( tpfile%cname ) // ': writing ' // Trim( tpfield%cmnhname ) )

kresp = 0

call IO_Mnhname_clean( tpfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tpfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_X3', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( pfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tpfile%nncid, ivarid, pfield(:,:,:), start = istarts(:), count = Int( Shape( pfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_X3', 'NF90_PUT_VAR', Trim( tpfield%cmnhname ), kresp )
end if

end subroutine IO_Field_partial_write_nc4_X3


subroutine IO_Field_partial_write_nc4_X4( tpfile, tpfield, pfield, koffset, kresp )

type(tfiledata),                      intent(in)  :: tpfile
class(tfieldmetadata),                intent(in)  :: tpfield
real,             dimension(:,:,:,:), intent(in)  :: pfield   ! array containing the data field
integer,          dimension(4),       intent(in)  :: koffset
integer,                              intent(out) :: kresp

character(len=NMNHNAMELGTMAX)      :: yvarname
integer(kind=CDFINT)               :: istatus
integer(kind=CDFINT)               :: ivarid
integer(kind=CDFINT), dimension(4) :: istarts

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_X4',&
                Trim( tpfile%cname ) // ': writing ' // Trim( tpfield%cmnhname ) )

kresp = 0

call IO_Mnhname_clean( tpfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tpfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_X4', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( pfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tpfile%nncid, ivarid, pfield(:,:,:,:), start = istarts(:), count = Int( Shape( pfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_X4', 'NF90_PUT_VAR', Trim( tpfield%cmnhname ), kresp )
end if

end subroutine IO_Field_partial_write_nc4_X4


subroutine IO_Field_partial_write_nc4_N2( tpfile, tpfield, kfield, koffset, kresp, kvertlevel, kzfile )

type(tfiledata),                  intent(in)  :: tpfile
class(tfieldmetadata),            intent(in)  :: tpfield
integer,          dimension(:,:), intent(in)  :: kfield   ! array containing the data field
integer,          dimension(2),   intent(in)  :: koffset
integer,                          intent(out) :: kresp
integer,                optional, intent(in)  :: kvertlevel ! Number of the vertical level (needed for Z-level split files)
integer,                optional, intent(in)  :: kzfile     ! Number of the Z-level split file

character(len=4)                    :: ysuffix
character(len=NMNHNAMELGTMAX)       :: yvarname
integer(kind=CDFINT)                :: istatus
integer(kind=CDFINT)                :: ivarid
integer(kind=CDFINT),  dimension(2) :: istarts
class(tfieldmetadata), pointer      :: tzfield
type(tfiledata),       pointer      :: tzfile

kresp = 0

call IO_Select_split_file( tpfile, tpfield, tzfile, tzfield, kvertlevel, kzfile )

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_N2',&
                Trim( tzfile%cname ) // ': writing ' // Trim( tzfield%cmnhname ) )

call IO_Mnhname_clean( tzfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tzfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_N2', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( kfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tzfile%nncid, ivarid, kfield(:,:), start = istarts(:), count = Int( Shape( kfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_N2', 'NF90_PUT_VAR', Trim( tzfield%cmnhname ), kresp )
end if

if ( Present( kvertlevel ) ) deallocate( tzfield )

end subroutine IO_Field_partial_write_nc4_N2


subroutine IO_Field_partial_write_nc4_N3( tpfile, tpfield, kfield, koffset, kresp, kvertlevel, kzfile )

type(tfiledata),                    intent(in)  :: tpfile
class(tfieldmetadata),              intent(in)  :: tpfield
integer,          dimension(:,:,:), intent(in)  :: kfield   ! array containing the data field
integer,          dimension(3),     intent(in)  :: koffset
integer,                            intent(out) :: kresp
integer,                  optional, intent(in)  :: kvertlevel ! Number of the vertical level (needed for Z-level split files)
integer,                  optional, intent(in)  :: kzfile     ! Number of the Z-level split file

character(len=4)                    :: ysuffix
character(len=NMNHNAMELGTMAX)       :: yvarname
integer(kind=CDFINT)                :: istatus
integer(kind=CDFINT)                :: ivarid
integer(kind=CDFINT),  dimension(3) :: istarts
class(tfieldmetadata), pointer      :: tzfield
type(tfiledata),       pointer      :: tzfile

kresp = 0

call IO_Select_split_file( tpfile, tpfield, tzfile, tzfield, kvertlevel, kzfile )

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_N3',&
                Trim( tzfile%cname ) // ': writing ' // Trim( tzfield%cmnhname ) )

call IO_Mnhname_clean( tzfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tzfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_N3', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( kfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tzfile%nncid, ivarid, kfield(:,:,:), start = istarts(:), count = Int( Shape( kfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_N3', 'NF90_PUT_VAR', Trim( tzfield%cmnhname ), kresp )
end if

if ( Present( kvertlevel ) ) deallocate( tzfield )

end subroutine IO_Field_partial_write_nc4_N3


subroutine IO_Field_partial_write_nc4_N4( tpfile, tpfield, kfield, koffset, kresp, kvertlevel, kzfile )

type(tfiledata),                      intent(in)  :: tpfile
class(tfieldmetadata),                intent(in)  :: tpfield
integer,          dimension(:,:,:,:), intent(in)  :: kfield   ! array containing the data field
integer,          dimension(4),       intent(in)  :: koffset
integer,                              intent(out) :: kresp
integer,                    optional, intent(in)  :: kvertlevel ! Number of the vertical level (needed for Z-level split files)
integer,                    optional, intent(in)  :: kzfile     ! Number of the Z-level split file

character(len=4)                    :: ysuffix
character(len=NMNHNAMELGTMAX)       :: yvarname
integer(kind=CDFINT)                :: istatus
integer(kind=CDFINT)                :: ivarid
integer(kind=CDFINT),  dimension(4) :: istarts
class(tfieldmetadata), pointer      :: tzfield
type(tfiledata),       pointer      :: tzfile

kresp = 0

call IO_Select_split_file( tpfile, tpfield, tzfile, tzfield, kvertlevel, kzfile )

call Print_msg( NVERB_DEBUG, 'IO', 'IO_Field_partial_write_nc4_N4',&
                Trim( tzfile%cname ) // ': writing ' // Trim( tzfield%cmnhname ) )

call IO_Mnhname_clean( tzfield%cmnhname, yvarname )

istatus = NF90_INQ_VARID( tzfile%nncid, yvarname, ivarid )
if ( istatus /= NF90_NOERR ) then
  call Print_msg( NVERB_FATAL, 'IO', 'IO_Field_partial_write_nc4_N4', 'variable ' // Trim( yvarname ) &
                  // ' not yet created (IO_Field_create not yet called?)' )
end if

! Write the data
if ( Size( kfield ) > 0 ) then
  istarts(:) = koffset(:) + 1
  istatus = NF90_PUT_VAR( tzfile%nncid, ivarid, kfield(:,:,:,:), start = istarts(:), count = Int( Shape( kfield ), kind = CDFINT ) )
  if (istatus /= NF90_NOERR) &
    call IO_Err_handle_nc4( istatus, 'IO_Field_partial_write_nc4_N4', 'NF90_PUT_VAR', Trim( tzfield%cmnhname ), kresp )
end if

if ( Present( kvertlevel ) ) deallocate( tzfield )

end subroutine IO_Field_partial_write_nc4_N4


subroutine IO_Coordvar_write_nc4( tpfile, hprogram_orig, tpdtmodeln )
use modd_aircraft_balloon
use modd_budget,     only: cbutype, lbu_icp, lbu_jcp, lbu_kcp, nbuih, nbuil, nbujh, nbujl, nbukh, nbukl, nbukmax, &
                           nbustep, nbutotwrite
use modd_conf,       only: cprogram, l2d, lcartesian
use modd_conf_n,     only: cstorage_type
use modd_dim_n,      only: nkmax
use modd_dyn_n,      only: xtstep
use modd_field,      only: NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NI_U, NMNHDIM_NJ_U, NMNHDIM_NI_V, NMNHDIM_NJ_V, &
                           NMNHDIM_LEVEL, NMNHDIM_LEVEL_W, NMNHDIM_TIME,                                   &
                           NMNHDIM_BUDGET_CART_NI,   NMNHDIM_BUDGET_CART_NJ,   NMNHDIM_BUDGET_CART_NI_U,   &
                           NMNHDIM_BUDGET_CART_NJ_U, NMNHDIM_BUDGET_CART_NI_V, NMNHDIM_BUDGET_CART_NJ_V,   &
                           NMNHDIM_BUDGET_CART_LEVEL, NMNHDIM_BUDGET_CART_LEVEL_W,                         &
                           NMNHDIM_BUDGET_MASK_LEVEL, NMNHDIM_BUDGET_MASK_LEVEL_W,                         &
                           NMNHDIM_BUDGET_TIME, NMNHDIM_BUDGET_LES_TIME, NMNHDIM_BUDGET_LES_AVG_TIME,      &
                           NMNHDIM_BUDGET_LES_LEVEL,                                                       &
                           NMNHDIM_SPECTRA_2PTS_NI, NMNHDIM_SPECTRA_2PTS_NJ, NMNHDIM_SPECTRA_LEVEL,        &
                           NMNHDIM_SERIES_LEVEL, NMNHDIM_SERIES_LEVEL_W, NMNHDIM_SERIES_TIME,              &
                           NMNHDIM_PROFILER_TIME, NMNHDIM_STATION_TIME,                                    &
                           tfieldlist
use modd_grid,       only: xlatori, xlonori
use modd_grid_n,     only: lsleve, xxhat, xxhatm, xyhat, xyhatm, xzhat, xzhatm, xxhat_ll, xyhat_ll, xxhatm_ll, xyhatm_ll
use modd_les,        only: cles_level_type, cspectra_level_type, nlesn_iinf, nlesn_isup, nlesn_jinf, nlesn_jsup, &
                           nles_k, nles_levels, nspectra_k, nspectra_levels,                                     &
                           xles_altitudes, xspectra_altitudes
use modd_les_n,      only: nles_dtcount, nles_mean_end, nles_mean_start, nles_mean_step, nles_mean_times, &
                           nles_times, nspectra_ni, nspectra_nj, tles_dates, xles_times
use modd_netcdf,     only: tdimnc
use modd_parameters, only: jphext, JPVEXT
use modd_profiler_n, only: lprofiler, tprofilers_time
use modd_series,     only: lseries
use modd_series_n,   only: nsnbstept, tpsdates
use modd_station_n,  only: lstation, tstations_time
use modd_time,       only: tdtseg
use modd_time_n,     only: tdtcur
use modd_type_date,  only: date_time

use mode_field,      only: Find_field_id_from_mnhname
use mode_gridproj,   only: Sm_latlon
use mode_nest_ll,    only: Get_model_number_ll, Go_tomodel_ll

type(tfiledata),            intent(in) :: tpfile
character(len=*), optional, intent(in) :: hprogram_orig !To emulate a file coming from this program
type(date_time),  optional, intent(in) :: tpdtmodeln    !Time of model (to force model date written in file)

character(len=:),                         allocatable :: ystdnameprefix
character(len=:),                         allocatable :: yprogram
integer                                               :: iiu, iju
integer                                               :: id, iid, iresp
integer                                               :: imi
integer                                               :: ji
integer                                               :: jt
integer                                               :: jtb, jte
integer(kind=cdfint)                                  :: incid
logical                                               :: gchangemodel
logical                                               :: gdealloc
logical,                         pointer              :: gsleve
real,            dimension(:),   pointer              :: zxhat, zyhat, zzhat
real,            dimension(:),   pointer              :: zxhatm, zyhatm, zzhatm !Coordinates at mass points in the transformed space
real,            dimension(:),            allocatable :: zles_levels
real,            dimension(:),            allocatable :: zspectra_levels
real,            dimension(:,:), pointer              :: zlat, zlon
type(tdimnc),                    pointer              :: tzdim_ni, tzdim_nj, tzdim_ni_u, tzdim_nj_u, tzdim_ni_v, tzdim_nj_v
type(date_time), dimension(:),            allocatable :: tzdates
type(date_time), dimension(:,:),          allocatable :: tzdates_bound

real, dimension(:),   pointer :: zxhat_glob,  zyhat_glob
real, dimension(:),   pointer :: zxhatm_glob, zyhatm_glob
!These variables are save: they are populated once for the master Z-split file and freed after the last file has been written
real, dimension(:,:), pointer, save :: zlatm_glob  => null(), zlonm_glob  => null()
real, dimension(:,:), pointer, save :: zlatu_glob  => null(), zlonu_glob  => null()
real, dimension(:,:), pointer, save :: zlatv_glob  => null(), zlonv_glob  => null()
real, dimension(:,:), pointer, save :: zlatf_glob  => null(), zlonf_glob  => null()
type(tfiledata)                                       :: tzfile


call Print_msg( NVERB_DEBUG, 'IO', 'IO_Coordvar_write_nc4', 'called for ' // Trim( tpfile%cname ) )

tzfile = tpfile
! Disable lossy compression to not loose precision for coordinates variables
tzfile%lnccompress_lossy = .false.

zxhat  => null()
zyhat  => null()
zzhat  => null()
zxhatm => null()
zyhatm => null()
zzhatm => null()
zxhat_glob  => null()
zyhat_glob  => null()
zxhatm_glob => null()
zyhatm_glob => null()

gchangemodel = .false.

if ( Present( hprogram_orig ) ) then
  yprogram = hprogram_orig
else
  yprogram = cprogram
endif

! Get the Netcdf file ID
incid = tzfile%nncid

call Get_model_number_ll( imi )

if ( tzfile%nmodel > 0 ) then
  call Find_field_id_from_mnhname( 'XHAT', iid, iresp )
  zxhat => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'YHAT', iid, iresp )
  zyhat => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'XHATM', iid, iresp )
  zxhatm => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'YHATM', iid, iresp )
  zyhatm => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'ZHAT', iid, iresp )
  zzhat => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'ZHATM', iid, iresp )
  zzhatm => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'XHAT_ll', iid, iresp )
  zxhat_glob => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'YHAT_ll', iid, iresp )
  zyhat_glob => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'XHATM_ll', iid, iresp )
  zxhatm_glob => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'YHATM_ll', iid, iresp )
  zyhatm_glob => tfieldlist(iid)%tfield_x1d(tzfile%nmodel)%data
  call Find_field_id_from_mnhname( 'SLEVE', iid, iresp )
  gsleve => tfieldlist(iid)%tfield_l0d(tzfile%nmodel)%data

  if ( imi /= tzfile%nmodel ) then
    !This is necessary to have correct domain sizes (used by Gather_xxfield)
    call Go_tomodel_ll( tzfile%nmodel, iresp )
    gchangemodel = .true.
  end if
else
  zxhat  => xxhat
  zyhat  => xyhat
  zzhat  => xzhat
  zxhatm => xxhatm
  zyhatm => xyhatm
  zzhatm => xzhatm
  zxhat_glob  => xxhat_ll
  zyhat_glob  => xyhat_ll
  zxhatm_glob => xxhatm_ll
  zyhatm_glob => xyhatm_ll
  gsleve => lsleve
end if

iiu = Size( zxhat )
iju = Size( zyhat )

if ( lcartesian ) then
  ystdnameprefix = 'plane'
else
  ystdnameprefix = 'projection'
endif

if ( Associated( tzfile%tncdims ) ) then
  tzdim_ni   => tzfile%tncdims%tdims(NMNHDIM_NI)
  tzdim_nj   => tzfile%tncdims%tdims(NMNHDIM_NJ)
  tzdim_ni_u => tzfile%tncdims%tdims(NMNHDIM_NI_U)
  tzdim_nj_u => tzfile%tncdims%tdims(NMNHDIM_NJ_U)
  tzdim_ni_v => tzfile%tncdims%tdims(NMNHDIM_NI_V)
  tzdim_nj_v => tzfile%tncdims%tdims(NMNHDIM_NJ_V)
else
  tzdim_ni   => Null()
  tzdim_nj   => Null()
  tzdim_ni_u => Null()
  tzdim_nj_u => Null()
  tzdim_ni_v => Null()
  tzdim_nj_v => Null()
end if

call Write_hor_coord1d( tzdim_ni,   'x-dimension of the grid', &
                        trim(ystdnameprefix)//'_x_coordinate',               'X', 0.,   jphext, jphext, zxhatm_glob )
call Write_hor_coord1d( tzdim_nj,   'y-dimension of the grid', &
                        trim(ystdnameprefix)//'_y_coordinate',               'Y', 0.,   jphext, jphext, zyhatm_glob )
call Write_hor_coord1d( tzdim_ni_u, 'x-dimension of the grid at u location', &
                        trim(ystdnameprefix)//'_x_coordinate_at_u_location', 'X', -0.5, jphext, 0,      zxhat_glob  )
call Write_hor_coord1d( tzdim_nj_u, 'y-dimension of the grid at u location', &
                        trim(ystdnameprefix)//'_y_coordinate_at_u_location', 'Y', 0.,   jphext, jphext, zyhatm_glob )
call Write_hor_coord1d( tzdim_ni_v, 'x-dimension of the grid at v location', &
                        trim(ystdnameprefix)//'_x_coordinate_at_v_location', 'X', 0.,   jphext, jphext, zxhatm_glob )
call Write_hor_coord1d( tzdim_nj_v, 'y-dimension of the grid at v location', &
                        trim(ystdnameprefix)//'_y_coordinate_at_v_location', 'Y', -0.5, jphext, 0,      zyhat_glob  )

!Warning: the following block has to be reenabled if IO_Coordvar_write_nc4 is again called for Z-split files
#if 0
!Deallocate only if it is a non Z-split file or the last Z-split subfile
gdealloc = .false.
if ( Associated( tzfile%tmainfile ) ) then
  if ( tzfile%cname == tzfile%tmainfile%tfiles_ioz(tzfile%tmainfile%nsubfiles_ioz)%tfile%cname ) gdealloc = .true.
else if ( tzfile%nsubfiles_ioz == 0 .and. .not. Associated( tzfile%tmainfile ) ) then
  gdealloc = .true.
end if
#else
gdealloc = .true.
#endif

if ( .not. lcartesian ) then
  !Compute latitude/longitude for the Arakawa points
  Allocate( zlat(iiu, iju), zlon(iiu, iju) )

  !If the file is a Z-split subfile, coordinates are already collected
  if ( .not. associated( tzfile%tmainfile ) ) then
    call Gather_hor_coord2d( zxhatm, zyhatm, zlatm_glob, zlonm_glob )
    call Gather_hor_coord2d( zxhat,  zyhatm, zlatu_glob, zlonu_glob )
    call Gather_hor_coord2d( zxhatm, zyhat,  zlatv_glob, zlonv_glob )
    call Gather_hor_coord2d( zxhat,  zyhat,  zlatf_glob, zlonf_glob )
  end if

  ! Mass point
  call Write_hor_coord2d( zlatm_glob, zlonm_glob, 'latitude',  'longitude')
  ! u point
  call Write_hor_coord2d( zlatu_glob, zlonu_glob, 'latitude_u','longitude_u')
  ! v point
  call Write_hor_coord2d( zlatv_glob, zlonv_glob, 'latitude_v','longitude_v')
  ! xi vorticity point (=f point =uv point)
  call Write_hor_coord2d( zlatf_glob, zlonf_glob, 'latitude_f','longitude_f')

  Deallocate( zlat, zlon )

  !The zlat/lon._glob were allocated in Gather_hor_coord2d calls
  !Deallocate only if it is non Z-split file or the last Z-split subfile
  if ( gdealloc ) Deallocate( zlatm_glob, zlonm_glob, zlatu_glob, zlonu_glob, zlatv_glob, zlonv_glob, zlatf_glob, zlonf_glob )
end if

if ( tzfile%lmaster ) then !vertical coordinates in the transformed space are the same on all processes
  if ( Trim( yprogram ) /= 'PGD' .and. Trim( yprogram ) /= 'NESPGD' .and. Trim( yprogram ) /= 'ZOOMPG' &
      .and. .not. ( Trim( yprogram ) == 'REAL' .and. cstorage_type == 'SU') ) then !condition to detect prep_surfex

    call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_LEVEL),  'position z in the transformed space',              '', &
                          'altitude',                0.,  JPVEXT, JPVEXT, ZZHATM )
    call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_LEVEL_W),'position z in the transformed space at w location','', &
                          'altitude_at_w_location', -0.5, JPVEXT, 0,      ZZHAT )
  END IF
END IF

!Write time scale
if ( tzfile%lmaster ) then !time scale is the same on all processes
  if ( Trim( yprogram ) /= 'PGD' .and. Trim( yprogram ) /= 'NESPGD' .and. Trim( yprogram ) /= 'ZOOMPG' &
      .and. .not. ( Trim( yprogram ) == 'REAL' .and. cstorage_type == 'SU' ) ) then !condition to detect prep_surfex
    if ( tzfile%ctype /= 'MNHDIACHRONIC' ) then
      if ( Present( tpdtmodeln ) ) then
        call Write_time_coord( tzfile%tncdims%tdims(nmnhdim_time), 'time axis', [ tpdtmodeln ] )
      else if ( Associated( tdtcur ) ) then
        call Write_time_coord( tzfile%tncdims%tdims(nmnhdim_time), 'time axis', [ tdtcur ] )
      end if
    end if
  end if
end if

if ( tzfile%lmaster ) then
  !Write coordinates used in diachronic files
  if ( tzfile%ctype == 'MNHDIACHRONIC' ) then
    if ( cbutype == 'CART' .or. cbutype == 'SKIP' ) then
      !Coordinates for the budgets in cartesian boxes
      if ( .not. lbu_icp )                                                                                                  &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NI), 'x-dimension of the budget cartesian box',    &
                         trim(ystdnameprefix)//'_x_coordinate', 'X', 0., 0, 0, zxhatm_glob(nbuil + jphext : nbuih + jphext) )
      if ( .not. lbu_jcp )                                                                                                  &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NJ), 'y-dimension of the budget cartesian box',    &
                         trim(ystdnameprefix)//'_y_coordinate', 'Y', 0., 0, 0, zyhatm_glob(nbujl + jphext : nbujh + jphext) )
      if ( .not. lbu_icp )                                                                    &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NI_U),               &
                                'x-dimension of the budget cartesian box at u location',      &
                                trim(ystdnameprefix)//'_x_coordinate_at_u_location',          &
                                'X', -0.5, 0, 0, zxhat_glob (nbuil + jphext : nbuih + jphext) )
      if ( .not. lbu_jcp )                                                                    &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NJ_U),               &
                                'y-dimension of the budget cartesian box at u location',      &
                                trim(ystdnameprefix)//'_y_coordinate_at_u_location',          &
                                'Y', 0.,   0, 0, zyhatm_glob(nbujl + jphext : nbujh + jphext) )
      if ( .not. lbu_icp )                                                                    &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NI_V),               &
                                'x-dimension of the budget cartesian box at v location',      &
                                trim(ystdnameprefix)//'_x_coordinate_at_v_location',          &
                                'X', 0.,   0, 0, zxhatm_glob(nbuil + jphext : nbuih + jphext) )
      if ( .not. lbu_jcp )                                                                    &
        call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_NJ_V),               &
                                'y-dimension of the budget cartesian box at v location',      &
                                trim(ystdnameprefix)//'_y_coordinate_at_v_location',          &
                                'Y', -0.5, 0, 0, zyhat_glob (nbujl + jphext : nbujh + jphext) )
      if ( .not. lbu_kcp )                                                                                      &
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_LEVEL),                                  &
                              'position z in the transformed space of the budget cartesian box',                &
                              '', 'altitude',               0.,   0, 0, zzhatm(nbukl + JPVEXT : nbukh + JPVEXT) )
      if ( .not. lbu_kcp )                                                                                     &
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_CART_LEVEL_W),                               &
                              'position z in the transformed space at w location of the budget cartesian box', &
                              '', 'altitude_at_w_location', -0.5, 0, 0, zzhat (nbukl + JPVEXT : nbukh + JPVEXT) )
    else if ( cbutype == 'MASK' ) then
      !Coordinates for the budgets masks
      if ( nbukmax > 0 )                                                                                        &
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_MASK_LEVEL),                                  &
                              'position z in the transformed space of the budget mask',                         &
                              '', 'altitude',               0.,   0, 0, zzhatm(nbukl + JPVEXT : nbukh + JPVEXT) )
      if ( nbukmax > 0 )                                                                                        &
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_MASK_LEVEL_W),                                &
                              'position z in the transformed space at w location of the budget mask',           &
                              '', 'altitude',               -0.5, 0, 0, zzhat (nbukl + JPVEXT : nbukh + JPVEXT) )

      !NMNHDIM_BUDGET_MASK_NBUMASK: not a true dimension
    end if !cbutype

    !Write time_budget coordinate + its boundaries
    if ( nbutotwrite > 0 ) then
      Allocate( tzdates(nbutotwrite) )
      Allocate( tzdates_bound(2, nbutotwrite) )

      do jt = 1, nbutotwrite
        tzdates(jt)%nyear  = tdtseg%nyear
        tzdates(jt)%nmonth = tdtseg%nmonth
        tzdates(jt)%nday   = tdtseg%nday
        tzdates(jt)%xtime  = tdtseg%xtime + nbustep * ( ( jt - 1 )  + 0.5  ) * xtstep

        tzdates_bound(1, jt)%nyear  = tdtseg%nyear
        tzdates_bound(1, jt)%nmonth = tdtseg%nmonth
        tzdates_bound(1, jt)%nday   = tdtseg%nday
        tzdates_bound(1, jt)%xtime  = tdtseg%xtime + nbustep * ( jt - 1 ) * xtstep

        tzdates_bound(2, jt)%nyear  = tdtseg%nyear
        tzdates_bound(2, jt)%nmonth = tdtseg%nmonth
        tzdates_bound(2, jt)%nday   = tdtseg%nday
        tzdates_bound(2, jt)%xtime  = tdtseg%xtime + nbustep * jt * xtstep
      end do

      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_TIME), 'time axis for budgets', tzdates, tzdates_bound )

      Deallocate( tzdates_bound )
      Deallocate( tzdates )
    end if

    !Coordinates for the number of LES budget time samplings
    if ( nles_times > 0 ) &
      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_LES_TIME), 'time axis for LES budgets', tles_dates )

    !Coordinates for the number of LES budget time averages
    !Condition also on nles_times to not create this coordinate when not used (no time average if nles_times=0)
    if ( nles_times > 0 .and. nles_mean_times > 0 ) then
      Allocate( tzdates(nles_mean_times) )
      Allocate( tzdates_bound(2, nles_mean_times) )

      do jt = 1, nles_mean_times
        jtb = ( nles_mean_start + ( jt - 1 ) * nles_mean_step ) / nles_dtcount
        jte = MIN( jtb + nles_mean_step / nles_dtcount, nles_mean_end / nles_dtcount, nles_times )
        ! jtb could be 0 if nles_mean_start is smaller than the first LES measurement
        ! For example, it occurs if xles_temp_mean_start is smaller than xles_temp_sampling (if xles_temp_mean_start=0.)
        ! Do this correction only after computation of jte
        if ( jtb < 1 ) jtb = 1

        tzdates(jt)%nyear  = tdtseg%nyear
        tzdates(jt)%nmonth = tdtseg%nmonth
        tzdates(jt)%nday   = tdtseg%nday
        tzdates(jt)%xtime        = tdtseg%xtime + ( xles_times(jtb) + xles_times(jte) ) / 2.
        !Not necessary:  call Datetime_correctdate( tzdates(jt ) )

        tzdates_bound(1, jt)%nyear  = tdtseg%nyear
        tzdates_bound(1, jt)%nmonth = tdtseg%nmonth
        tzdates_bound(1, jt)%nday   = tdtseg%nday
        tzdates_bound(1, jt)%xtime  = tdtseg%xtime + xles_times(jtb)

        tzdates_bound(2, jt)%nyear  = tdtseg%nyear
        tzdates_bound(2, jt)%nmonth = tdtseg%nmonth
        tzdates_bound(2, jt)%nday   = tdtseg%nday
        tzdates_bound(2, jt)%xtime  = tdtseg%xtime + xles_times(jte)
      end do
      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_LES_AVG_TIME), 'time axis for LES budget time averages', &
                             tzdates, tzdates_bound )

      Deallocate( tzdates_bound )
      Deallocate( tzdates )
    end if

    !Coordinates for the number of vertical levels for local LES budgets
    if ( nles_k > 0 ) then
      if ( cles_level_type == 'K' ) then
        Allocate( zles_levels(nles_k) )
        do ji = 1, nles_k
          zles_levels(ji) = zzhatm(nles_levels(ji) + JPVEXT)
        end do
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_LES_LEVEL),           &
                              'position z in the transformed space of the LES budgets', &
                              '', 'altitude', 0., 0, 0, zles_levels(:)                  )
        Deallocate( zles_levels )
      else if ( cles_level_type == 'Z' ) then
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_BUDGET_LES_LEVEL), 'altitude levels for the LES budgets', &
                              'altitude', '', 0., 0, 0, xles_altitudes(1:nles_k) )
      else
        call Print_msg( NVERB_ERROR, 'IO', 'IO_Coordvar_write_nc4','invalid cles_level_type' )
      end if
    end if

    !NMNHDIM_BUDGET_LES_SV: not a true dimension

    !Coordinates for the number of horizontal wavelengths for non-local LES budgets (2 points correlations)
    if ( nspectra_ni > 0 ) &
      call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_SPECTRA_2PTS_NI), 'x-dimension of the LES budget cartesian box',    &
                            trim(ystdnameprefix)//'_x_coordinate', 'X', 0., 0, 0,                                              &
                            zxhatm_glob(nlesn_iinf(imi) + jphext : nlesn_isup(imi) + jphext) )
    if ( nspectra_nj > 0 .and. .not. l2d ) &
      call Write_hor_coord1d( tzfile%tncdims%tdims(NMNHDIM_SPECTRA_2PTS_NJ), 'y-dimension of the LES budget cartesian box',    &
                              trim(ystdnameprefix)//'_y_coordinate', 'Y', 0., 0, 0,                                            &
                              zyhatm_glob(nlesn_jinf(imi) + jphext : nlesn_jsup(imi) + jphext) )


    !NMNHDIM_SPECTRA_SPEC_NI, NMNHDIM_SPECTRA_SPEC_NJ: not true dimensions: spectra wavelengths

    !Coordinates for the number of vertical levels for non-local LES budgets
    if ( nspectra_k > 0 ) then
      if ( cspectra_level_type == 'K' ) then
        Allocate( zspectra_levels(nspectra_k) )
        do ji = 1, nspectra_k
          zspectra_levels(ji) = zzhatm(nspectra_levels(ji) + JPVEXT)
        end do
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_SPECTRA_LEVEL),                        &
                              'position z in the transformed space of the non-local LES budgets', &
                              '', 'altitude', 0., 0, 0, zspectra_levels(:)                        )
        Deallocate( zspectra_levels )
      else if ( cspectra_level_type == 'Z' ) then
        call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_SPECTRA_LEVEL), 'altitude levels for the non-local LES budgets', &
                              'altitude', '', 0., 0, 0, xspectra_altitudes(1 : nspectra_k) )
      else
        call Print_msg( NVERB_ERROR, 'IO', 'IO_Coordvar_write_nc4','invalid cspectra_level_type' )
      end if
    end if

    !Coordinates for the number of profiler times
    if ( lprofiler ) &
      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_PROFILER_TIME), 'time axis for profilers', tprofilers_time%tpdates )

    !Coordinates for the number of station times
    if ( lstation ) &
      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_STATION_TIME), 'time axis for stations', tstations_time%tpdates )

    !Dimension for the number of series times
    if ( lseries .and. nsnbstept > 0 ) then
      call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_SERIES_LEVEL),                   &
                            'position z in the transformed space of the temporal series', &
                            '', 'altitude', 0., 0, 0, zzhatm(1 + JPVEXT : nkmax + JPVEXT) )
      call Write_ver_coord( tzfile%tncdims%tdims(NMNHDIM_SERIES_LEVEL_W),                                 &
                            'position z in the transformed space at w location of the temporal series',   &
                            '', 'altitude_at_w_location', -0.5, 0, 0, zzhat (1 + JPVEXT : nkmax + JPVEXT) )
      call Write_time_coord( tzfile%tncdims%tdims(NMNHDIM_SERIES_TIME), 'time axis for temporal series', tpsdates )
    end if

    if ( lflyer ) then
      ! Remark: to work flyer data must be on the file master rank
      ! This is currently ensured in WRITE_AIRCRAFT_BALLOON subroutine
      do ji = 1, nballoons
        if ( associated( tballoons(ji)%tballoon ) ) then
          call Write_flyer_time_coord( tballoons(ji)%tballoon )
        else
          call Print_msg( NVERB_ERROR, 'IO', 'IO_Coordvar_write_nc4','tballoon not associated' )
        end if
      end do

      do ji = 1, naircrafts
        if ( associated( taircrafts(ji)%taircraft ) ) then
          call Write_flyer_time_coord( taircrafts(ji)%taircraft )
        else
          call Print_msg( NVERB_ERROR, 'IO', 'IO_Coordvar_write_nc4','taircraft not associated' )
        end if
      end do
    end if

  end if !MNHDIACHRONIC

end if


if ( gchangemodel ) call Go_tomodel_ll( imi, iresp )


contains

#if 0
!Not used anymore
subroutine Gather_hor_coord1d( haxis, pcoords_loc, pcoords_glob )
  use mode_allocbuffer_ll, only: Allocbuffer_ll
  use mode_gather_ll,      only: Gather_xxfield

  character(len=*),            intent(in)  :: haxis
  real, dimension(:),          intent(in)  :: pcoords_loc
  real, dimension(:), pointer, intent(out) :: pcoords_glob

  character(len=2) :: ydir
  integer          :: ierr
  logical          :: galloc

  if ( haxis == 'X' ) then
    ydir = 'XX'
  else if ( haxis == 'Y' ) then
    ydir = 'YY'
  else
    call Print_msg( NVERB_FATAL, 'IO', 'Gather_hor_coord1d', 'invalid haxis ('//trim(haxis)//')' )
  end if

  ! Allocate pcoords_glob
  if ( gsmonoproc ) then ! sequential execution
    allocate( pcoords_glob( size( pcoords_loc) ) )
  else if ( tzfile%nsubfiles_ioz > 0 ) then
    !If there are Z-split subfiles, all subfile writers need the coordinates
    call Allocbuffer_ll( pcoords_glob, pcoords_loc, ydir, galloc )
  else if ( .not. tzfile%lmaster ) then
    allocate( pcoords_glob(0 ) ) !to prevent false positive with valgrind
  else !Master process
    call Allocbuffer_ll( pcoords_glob, pcoords_loc, ydir, galloc )
  end if

  !Gather coordinates
  if ( gsmonoproc ) then ! sequential execution
      pcoords_glob(: ) = pcoords_loc(: )
  else ! multiprocesses execution
      call Gather_xxfield( ydir, pcoords_loc, pcoords_glob, tzfile%nmaster_rank, tzfile%nmpicomm )
  endif

  !If the file has Z-split subfiles, broadcast the coordinates to all processes
  if ( tzfile%nsubfiles_ioz > 0 ) &
    call MPI_BCAST( pcoords_glob, size( pcoords_glob ), MNHREAL_MPI, tzfile%nmaster_rank - 1,  tzfile%nmpicomm, ierr )

end subroutine Gather_hor_coord1d
#endif


subroutine Gather_hor_coord2d( px, py, plat_glob, plon_glob )
  use mode_allocbuffer_ll, only: Allocbuffer_ll
  use mode_gather_ll,      only: Gather_xyfield

  real,dimension(:), intent(in) :: px
  real,dimension(:), intent(in) :: py
  real, dimension(:,:), pointer, intent(out) :: plat_glob
  real, dimension(:,:), pointer, intent(out) :: plon_glob

  integer          :: ierr
  logical :: galloc1, galloc2

  call Sm_latlon( xlatori, xlonori,                             &
                  spread( source = px, dim = 2, ncopies = iju), &
                  spread( source = py, dim = 1, ncopies = iiu), &
                  zlat, zlon )

  ! Allocate coordinate arrays
  if ( gsmonoproc ) then ! sequential execution
    allocate( plat_glob( size( zlat, 1 ), size( zlat, 2 ) ) )
    allocate( plon_glob( size( zlon, 1 ), size( zlon, 2 ) ) )
  else if ( tzfile%nsubfiles_ioz > 0 ) then
    !If there are Z-split subfiles, all subfile writers need the coordinates
    call Allocbuffer_ll( plat_glob, zlat, 'XY', galloc1 )
    call Allocbuffer_ll( plon_glob, zlon, 'XY', galloc2 )
  else if ( .not. tzfile%lmaster ) then
    allocate( plat_glob( 0, 0 ), plon_glob( 0, 0 ) ) !to prevent false positive with valgrind
  else !Master process
    call Allocbuffer_ll( plat_glob, zlat, 'XY', galloc1 )
    call Allocbuffer_ll( plon_glob, zlon, 'XY', galloc2 )
  end if

  !Gather coordinates
  if ( gsmonoproc ) then ! sequential execution
      plat_glob(:, : ) = zlat(:, : )
      plon_glob(:, : ) = zlon(:, : )
  else ! multiprocesses execution
      call Gather_xyfield( zlat, plat_glob, tzfile%nmaster_rank, tzfile%nmpicomm )
      call Gather_xyfield( zlon, plon_glob, tzfile%nmaster_rank, tzfile%nmpicomm )
  endif

  !If the file has Z-split subfiles, broadcast the coordinates to all processes
  if ( tzfile%nsubfiles_ioz > 0 ) then
    call MPI_BCAST( plat_glob, size( plat_glob ), MNHREAL_MPI, tzfile%nmaster_rank - 1,  tzfile%nmpicomm, ierr )
    call MPI_BCAST( plon_glob, size( plon_glob ), MNHREAL_MPI, tzfile%nmaster_rank - 1,  tzfile%nmpicomm, ierr )
  end if
end subroutine Gather_hor_coord2d


subroutine Write_hor_coord1d(TDIM,HLONGNAME,HSTDNAME,HAXIS,PSHIFT,KBOUNDLOW,KBOUNDHIGH,PCOORDS)
  USE MODE_ALLOCBUFFER_ll, ONLY: ALLOCBUFFER_ll
  USE MODE_GATHER_ll,      ONLY: GATHER_XXFIELD

  TYPE(tdimnc),               INTENT(IN) :: TDIM
  CHARACTER(LEN=*),           INTENT(IN) :: HLONGNAME
  CHARACTER(LEN=*),           INTENT(IN) :: HSTDNAME
  CHARACTER(LEN=*),           INTENT(IN) :: HAXIS
  REAL,                       INTENT(IN) :: PSHIFT
  INTEGER,                    INTENT(IN) :: KBOUNDLOW
  INTEGER,                    INTENT(IN) :: KBOUNDHIGH
  REAL, DIMENSION(:), TARGET, INTENT(IN) :: PCOORDS

  CHARACTER(LEN=64)             :: YRANGE
  CHARACTER(LEN=:),ALLOCATABLE  :: YVARNAME
  INTEGER                       :: IRESP
  INTEGER                       :: ISIZE
  INTEGER(KIND=CDFINT)          :: IVARID
  INTEGER(KIND=CDFINT)          :: IVDIM
  INTEGER(KIND=CDFINT)          :: ISTATUS

  IF (tzfile%LMASTER) THEN
    isize    = tdim%nlen
    yvarname = Trim( tdim%cname )
    ivdim    = tdim%nid

    ISTATUS = NF90_INQ_VARID(INCID, YVARNAME, IVARID)
    IF (ISTATUS /= NF90_NOERR) THEN
      ! Define the coordinate variable
      ISTATUS = NF90_DEF_VAR(INCID, YVARNAME, MNHREAL_NF90, IVDIM, IVARID)
      IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_DEF_VAR',trim(YVARNAME))
    ELSE
      CALL PRINT_MSG(NVERB_ERROR,'IO','Write_hor_coord1d',TRIM(YVARNAME)//' already defined')
    END IF

    ! Write metadata
    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'long_name',HLONGNAME)
    IF (ISTATUS /= NF90_NOERR) &
      CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','long_name for '//trim(YVARNAME))

    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'standard_name',HSTDNAME)
    IF (ISTATUS /= NF90_NOERR) &
      CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','standard_name for '//trim(YVARNAME))

    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'units','m')
    IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','units for '//trim(YVARNAME))

    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'axis',HAXIS)
    IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','axis for '//trim(YVARNAME))

    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'c_grid_axis_shift',PSHIFT)
    IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','c_grid_axis_shift for ' &
                                                     //trim(YVARNAME))

    WRITE(YRANGE,'( I0,":",I0 )') 1+KBOUNDLOW,ISIZE-KBOUNDHIGH
    ISTATUS = NF90_PUT_ATT(INCID, IVARID, 'c_grid_dynamic_range',TRIM(YRANGE))
    IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_ATT','c_grid_dynamic_range for ' &
                                                     //trim(YVARNAME))

    ! Write the data
    ISTATUS = NF90_PUT_VAR(INCID, IVARID, PCOORDS)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'Write_hor_coord1d','NF90_PUT_VAR',trim(YVARNAME),IRESP)
  END IF
end subroutine Write_hor_coord1d


subroutine Write_hor_coord2d( plat, plon, hlat, hlon )
  real,dimension(:,:), intent(in) :: plat
  real,dimension(:,:), intent(in) :: plon
  character(len=*),    intent(in) :: hlat
  character(len=*),    intent(in) :: hlon

  if ( tzfile%lmaster ) then
    call Find_field_id_from_mnhname( hlat, id, iresp )
    call IO_Field_write_nc4_x2( tzfile, tfieldlist(id ), plat, iresp, oiscoord = .true. )
    call Find_field_id_from_mnhname( hlon, id, iresp )
    call IO_Field_write_nc4_x2( tzfile, tfieldlist(id ), plon, iresp, oiscoord = .true. )
  end if
end subroutine Write_hor_coord2d


SUBROUTINE WRITE_VER_COORD(TDIM,HLONGNAME,HSTDNAME,HCOMPNAME,PSHIFT,KBOUNDLOW,KBOUNDHIGH,PCOORDS)
  TYPE(tdimnc),          INTENT(IN) :: TDIM
  CHARACTER(LEN=*),      INTENT(IN) :: HLONGNAME
  CHARACTER(LEN=*),      INTENT(IN) :: HSTDNAME
  CHARACTER(LEN=*),      INTENT(IN) :: HCOMPNAME
  REAL,                  INTENT(IN) :: PSHIFT
  INTEGER,               INTENT(IN) :: KBOUNDLOW
  INTEGER,               INTENT(IN) :: KBOUNDHIGH
  REAL,DIMENSION(:),     INTENT(IN) :: PCOORDS

  CHARACTER(LEN=64)             :: YRANGE
  CHARACTER(LEN=:),ALLOCATABLE  :: YVARNAME
  INTEGER                       :: IRESP
  INTEGER                       :: ISIZE
  INTEGER                       :: JI
  INTEGER(KIND=CDFINT)          :: IVARID
  INTEGER(KIND=CDFINT)          :: IVDIM
  INTEGER(KIND=CDFINT)          :: istatus

  isize    = tdim%nlen
  yvarname = Trim( tdim%cname )
  ivdim    = tdim%nid

  istatus = NF90_INQ_VARID(INCID, YVARNAME, IVARID)
  IF (istatus /= NF90_NOERR) THEN
    ! Define the coordinate variable
    istatus = NF90_DEF_VAR(INCID, YVARNAME, MNHREAL_NF90, IVDIM, IVARID)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_DEF_VAR',trim(YVARNAME))
  ELSE
    CALL PRINT_MSG(NVERB_ERROR,'IO','WRITE_VER_COORD',TRIM(YVARNAME)//' already defined')
  END IF

  ! Write metadata
  istatus = NF90_PUT_ATT(INCID, IVARID, 'long_name',HLONGNAME)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','long_name for '//trim(YVARNAME))
  if ( Len_trim( hstdname ) > 0 ) then
    istatus = NF90_PUT_ATT(INCID, IVARID, 'standard_name',HSTDNAME)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','standard_name for '//trim(YVARNAME))
  end if
  istatus = NF90_PUT_ATT(INCID, IVARID, 'units','m')
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','units for '//trim(YVARNAME))
  istatus = NF90_PUT_ATT(INCID, IVARID, 'axis','Z')
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','axis for '//trim(YVARNAME))
  istatus = NF90_PUT_ATT(INCID, IVARID, 'positive','up')
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','positive for '//trim(YVARNAME))
  istatus = NF90_PUT_ATT(INCID, IVARID, 'c_grid_axis_shift',PSHIFT)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','c_grid_axis_shift for ' &
                                                   //trim(YVARNAME))
  WRITE(YRANGE,'( I0,":",I0 )') 1+KBOUNDLOW,ISIZE-KBOUNDHIGH
  istatus = NF90_PUT_ATT(INCID, IVARID, 'c_grid_dynamic_range',TRIM(YRANGE))
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','c_grid_dynamic_range for ' &
                                                   //trim(YVARNAME))
  !
  if ( Len_trim( hcompname ) > 0 ) then
    IF (GSLEVE) THEN
      !Remark: ZS, ZSMT and ZTOP in the formula are the same for mass point or flux point
      istatus = NF90_PUT_ATT(INCID, IVARID,'formula_terms','s: '//TRIM(YVARNAME)//                   &
                                          ' height: ZTOP oro_ls: ZSMT oro: ZS len1: LEN1 len2: LEN2')
      IF (istatus /= NF90_NOERR) &
        CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','formula_terms for '//trim(YVARNAME))
      istatus = NF90_PUT_ATT(INCID, IVARID, 'formula_definition','z(n,k,j,i)=s(k)'//                                      &
                            '+ oro_ls(j,i)*sinh((height/len1)**1.35-(s(k)/len1)**1.35)/sinh((s(k)/len1)**1.35)'//        &
                            '+(oro(j,i)-oro_ls(j,i))*sinh((height/len2)**1.35-(s(k)/len2)**1.35)/sinh((s(k)/len2)**1.35)')
      IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','formula_definition for ' &
                                                      //trim(YVARNAME))
    ELSE
      !Remark: ZS and ZTOP in the formula are the same for mass point or flux point
      istatus = NF90_PUT_ATT(INCID, IVARID, 'formula_terms','s: '//TRIM(YVARNAME)//' height: ZTOP orog: ZS')
      IF (istatus /= NF90_NOERR) &
        CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','formula_terms for '//trim(YVARNAME))
      istatus = NF90_PUT_ATT(INCID, IVARID, 'formula_definition','z(n,k,j,i)=s(k)*(height-orog(j,i))/height+orog(j,i)')
      IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','formula_definition for ' &
                                                      //trim(YVARNAME))
    ENDIF
    !
    istatus = NF90_PUT_ATT(INCID, IVARID, 'computed_standard_name',HCOMPNAME)
    IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_ATT','computed_standard_name for ' &
                                                     //trim(YVARNAME))
  end if

  ! Write the data
  istatus = NF90_PUT_VAR(INCID, IVARID, PCOORDS)
  IF (istatus /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'WRITE_VER_COORD','NF90_PUT_VAR',trim(YVARNAME))

END SUBROUTINE WRITE_VER_COORD

subroutine Write_time_coord( tdim, hlongname, tpdates, tpdates_bound )
  use modd_field,      only: NMNHDIM_PAIR, tfieldlist
  use modd_time_n,     only: tdtmod
  use modd_type_date,  only: date_time

  use mode_datetime,   only: Datetime_distance
  use mode_field,      only: Find_field_id_from_mnhname

  type(tdimnc),                                       intent(in) :: tdim
  character(len=*),                                   intent(in) :: hlongname
  type(date_time), dimension(:),                      intent(in) :: tpdates
  type(date_time), dimension(:,:),          optional, intent(in) :: tpdates_bound !Boundaries of the date intervals

  character(len=40)                                 :: yunits
  character(len=:),                     allocatable :: yvarname
  integer                                           :: jt
  integer(kind=CDFINT)                              :: ivarid
  integer(kind=CDFINT)                              :: ivdim
  integer(kind=CDFINT)                              :: istatus
  integer(kind=CDFINT), dimension(2)                :: ivdims
  real,                 dimension(:),   allocatable :: zdeltatimes_1d
  real,                 dimension(:,:), allocatable :: zdeltatimes_2d
  type(date_time)                                   :: tzref

  if ( Associated( tdtmod ) ) then
    ! Model beginning date (tdtmod%tdate) is used as the reference date
    ! Reference time is set to 0.
    tzref = tdtmod
    tzref%xtime = 0.
  else
    tzref%nyear  = 2000
    tzref%nmonth = 1
    tzref%nday   = 1
    tzref%xtime  = 0.
  end if

  yvarname = Trim( tdim%cname )
  ivdim    = tdim%nid

  istatus = NF90_INQ_VARID( incid, yvarname, ivarid )
  if ( istatus /= NF90_NOERR ) then
    ! Define the coordinate variable
    istatus = NF90_DEF_VAR( incid, yvarname, mnhreal_nf90, ivdim, ivarid )
    if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_DEF_VAR', Trim( yvarname ) )
  else
    call Print_msg( NVERB_ERROR, 'IO', 'Write_time_coord', Trim( yvarname ) // ' already defined' )
  end if

  ! Write metadata
  istatus = NF90_PUT_ATT( incid, ivarid, 'long_name', Trim( hlongname ) )
  if ( istatus /= NF90_NOERR ) &
    call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', 'long_name for ' // Trim( yvarname ) )
  istatus = NF90_PUT_ATT( incid, ivarid, 'standard_name','time' )
  IF ( istatus /= NF90_NOERR ) &
    call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', 'standard_name for ' // Trim( yvarname ) )
  Write( yunits, '( "seconds since ", i4.4, "-", i2.2, "-", i2.2, " 00:00:00 +0:00" )' ) &
         tzref%nyear, tzref%nmonth, tzref%nday
  istatus = NF90_PUT_ATT( incid, ivarid, 'units', yunits )
  if ( istatus /= NF90_NOERR ) &
    call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', 'units for ' // Trim( yvarname ) )
  istatus = NF90_PUT_ATT( incid, ivarid, 'axis', 'T' )
  if ( istatus /= NF90_NOERR ) &
    call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', 'axis for ' // Trim( yvarname ) )
  istatus = NF90_PUT_ATT( incid, ivarid, 'calendar', 'standard' )
  if ( istatus /= NF90_NOERR ) &
    call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', 'calendar for ' // Trim( yvarname ) )
  if ( Present( tpdates_bound ) ) then
    istatus = NF90_PUT_ATT( incid, ivarid, 'bounds', Trim( yvarname ) // '_bounds' )
    if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_ATT', &
                                                         'bounds for ' // Trim( yvarname ) )
  end if

  ! Compute the temporal distance from reference
  Allocate( zdeltatimes_1d(Size( tpdates )) )
  do jt = 1, Size( tpdates )
    call Datetime_distance( tzref, tpdates(jt), zdeltatimes_1d(jt) )
  end do

  ! Write the data
  istatus = NF90_PUT_VAR( incid, ivarid, zdeltatimes_1d )
  if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_VAR', Trim( yvarname ) )

  Deallocate( zdeltatimes_1d )

  !Write the date interval boundaries (if provided)
  if ( Present( tpdates_bound ) ) then
    yvarname = Trim( tdim%cname ) // '_bounds'
    istatus = NF90_INQ_VARID( incid, yvarname, ivarid )
    if ( istatus /= NF90_NOERR ) then
      ! Define the coordinate variable
      ivdims(1) = tzfile%tncdims%tdims(NMNHDIM_PAIR)%nid
      ivdims(2) = tdim%nid
      istatus = NF90_DEF_VAR( incid, yvarname, MNHREAL_NF90, ivdims, ivarid )
      if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_DEF_VAR', Trim( yvarname ) )
    else
      call Print_msg( NVERB_ERROR, 'IO', 'Write_time_coord', Trim( yvarname ) // ' already defined' )
    end if

    ! Compute the temporal distance from reference
    Allocate( zdeltatimes_2d(2, Size( tpdates_bound, 2 )) )
    do jt = 1, Size( tpdates_bound, 2 )
      call Datetime_distance( tzref, tpdates_bound(1, jt), zdeltatimes_2d(1, jt) )
      call Datetime_distance( tzref, tpdates_bound(2, jt), zdeltatimes_2d(2, jt) )
    end do

    ! Write the data
    istatus = NF90_PUT_VAR( incid, ivarid, zdeltatimes_2d(:,:) )
    if ( istatus /= NF90_NOERR ) call IO_Err_handle_nc4( istatus, 'Write_time_coord', 'NF90_PUT_VAR', Trim( yvarname ) )
  end if

end subroutine Write_time_coord


subroutine Write_flyer_time_coord( tpflyer )
  use NETCDF

  use modd_aircraft_balloon
  use modd_parameters,       only: NBUNAMELGTMAX, XUNDEF

  use mode_aircraft_balloon, only: Aircraft_balloon_longtype_get
  use mode_io_tools_nc4,     only: IO_Mnhname_clean


  class(tflyerdata), intent(in) :: tpflyer

  character(len=NBUNAMELGTMAX) :: ytype
  character(len=NBUNAMELGTMAX) :: ytype_clean
  integer                      :: istatus
  integer(kind=CDFINT)         :: icatid
  integer(kind=CDFINT)         :: isubcatid
  integer(kind=CDFINT)         :: idimid
  type(tdimnc),        pointer :: tzdim

  !Do it only if correct model level and has really flown
  if ( tpflyer%nmodel == imi ) then
  ! Do the second if only if the first one is OK (if not, xx may be not allocated)
  if ( Count( tpflyer%xx /= XUNDEF) > 1 ) then
    Allocate( tzdim )

    istatus = NF90_INQ_NCID( tzfile%nncid, 'Flyers', icatid )
    if ( istatus /= NF90_NOERR ) then
      call Print_msg( NVERB_ERROR, 'IO', 'Write_flyer_time_coord', &
                      Trim( tzfile%cname ) // ': group Flyers not found' )
    end if

    call Aircraft_balloon_longtype_get( tpflyer, ytype )
    call IO_Mnhname_clean( ytype, ytype_clean )
    istatus = NF90_INQ_NCID( icatid, Trim( ytype_clean ), isubcatid )
    if ( istatus /= NF90_NOERR ) then
      call Print_msg( NVERB_ERROR, 'IO', 'Write_flyer_time_coord', &
                      Trim( tzfile%cname ) // ': group ' // Trim( ytype_clean ) // ' not found' )
    end if

    istatus = NF90_INQ_NCID( isubcatid, Trim( tpflyer%cname ), incid )
    if ( istatus /= NF90_NOERR ) then
      call Print_msg( NVERB_ERROR, 'IO', 'Write_flyer_time_coord', &
                      Trim( tzfile%cname ) // ': group '// Trim( tpflyer%cname ) // ' not found' )
    end if

    istatus = NF90_INQ_DIMID( incid, 'time_flyer', idimid )
    if ( istatus /= NF90_NOERR ) then
      call Print_msg( NVERB_ERROR, 'IO', 'Write_flyer_time_coord', &
                      Trim( tzfile%cname ) // ': group ' // Trim( tpflyer%cname ) // ' time_flyer dimension not found' )
    end if

    tzdim%cname = 'time_flyer'
    istatus = NF90_INQUIRE_DIMENSION( incid, idimid, len = tzdim%nlen )
    tzdim%nid = idimid

    !Remark: incid is used in Write_time_coord
    call Write_time_coord( tzdim, 'time axis for flyer', tpflyer%tflyer_time%tpdates )

    Deallocate( tzdim )

    !Restore file identifier to root group
    incid = tzfile%nncid
  end if
  end if

end subroutine Write_flyer_time_coord

END SUBROUTINE IO_Coordvar_write_nc4


SUBROUTINE IO_Header_write_nc4(TPFILE)
!
TYPE(TFILEDATA), INTENT(IN)  :: TPFILE ! File structure
!
INTEGER(KIND=CDFINT) :: ISTATUS
!
IF (TRIM(TPFILE%CFORMAT)/='NETCDF4' .AND. TRIM(TPFILE%CFORMAT)/='LFICDF4') RETURN
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_Header_write_nc4','called for file '//TRIM(TPFILE%CNAME))
!
IF (TPFILE%LMASTER)  THEN
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'Conventions', 'CF-1.10 COMODO-1.4')
  IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_FILE_WRITE_HEADER','NF90_PUT_ATT','Conventions')

#if (MNH_REAL == 8)
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'MNH_REAL', '8')
#else
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'MNH_REAL', '4')
#endif
  IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_FILE_WRITE_HEADER','NF90_PUT_ATT','MNH_REAL')

#if (MNH_INT == 4)
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'MNH_INT', '4')
#else
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'MNH_INT', '8')
#endif
  IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_FILE_WRITE_HEADER','NF90_PUT_ATT','MNH_INT')

  if ( tpfile%ldimreduced ) then
    istatus = NF90_PUT_ATT( tpfile%nncid, NF90_GLOBAL, 'MNH_REDUCE_DIMENSIONS_IN_FILES', '1')
  else
    istatus = NF90_PUT_ATT( tpfile%nncid, NF90_GLOBAL, 'MNH_REDUCE_DIMENSIONS_IN_FILES', '0')
  endif

!title

  !history
  CALL IO_History_append_nc4(TPFILE)

!institution

!source

!comment

!references
END IF
!
END SUBROUTINE IO_Header_write_nc4


SUBROUTINE IO_History_append_nc4(TPFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN)  :: TPFILE ! File structure
!
INTEGER,PARAMETER :: YEAR=1, MONTH=2, DAY=3, HH=5, MM=6, SS=7
!
CHARACTER(len=5)             :: YZONE
CHARACTER(LEN=:),ALLOCATABLE :: YCMD, YHISTORY, YHISTORY_NEW, YHISTORY_PREV
INTEGER                      :: ILEN_CMD, ILEN_PREV
INTEGER(KIND=CDFINT)         :: ILEN_NC
INTEGER(KIND=CDFINT)         :: ISTATUS
INTEGER,DIMENSION(8)         :: IDATETIME
!
IF (TRIM(TPFILE%CFORMAT)/='NETCDF4' .AND. TRIM(TPFILE%CFORMAT)/='LFICDF4') RETURN
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_History_append_nc4','called for file '//TRIM(TPFILE%CNAME))
!
IF (TPFILE%LMASTER)  THEN
  !Check if history attribute already exists in file and read it
  ISTATUS = NF90_INQUIRE_ATTRIBUTE(TPFILE%NNCID, NF90_GLOBAL, 'history', LEN=ILEN_NC)
  ILEN_PREV = int( ILEN_NC, kind=kind(ILEN_PREV) )
  IF (ISTATUS == NF90_NOERR) THEN
    ALLOCATE(CHARACTER(LEN=ILEN_PREV) :: YHISTORY_PREV)
    ISTATUS = NF90_GET_ATT(TPFILE%NNCID, NF90_GLOBAL, 'history', YHISTORY_PREV)
    IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(ISTATUS,'IO_History_append_nc4','NF90_GET_ATT','history')
    YHISTORY_PREV = YHISTORY_PREV
  ELSE
    ILEN_PREV = 0
  END IF

  !Get date and time
  call DATE_AND_TIME( VALUES=IDATETIME, ZONE=YZONE )
  call GET_COMMAND(lENGTH=ILEN_CMD)
  IF (ILEN_CMD>0) THEN
    ALLOCATE(CHARACTER(ILEN_CMD) :: YCMD)
    CALL GET_COMMAND(COMMAND=YCMD)
  END IF

  !Prepare new history entry
  ALLOCATE(CHARACTER(26+ILEN_CMD) :: YHISTORY_NEW)
  WRITE(YHISTORY_NEW, '( I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,A5,": ", A )') &
        IDATETIME(YEAR),IDATETIME(MONTH),IDATETIME(DAY),IDATETIME(HH),IDATETIME(MM),IDATETIME(SS),YZONE, YCMD

  !Write full history
  IF (ILEN_PREV == 0) THEN
    YHISTORY = YHISTORY_NEW
  ELSE
    YHISTORY = YHISTORY_NEW//NEW_LINE('A')//YHISTORY_PREV
  END IF
  ISTATUS = NF90_PUT_ATT(TPFILE%NNCID, NF90_GLOBAL, 'history', YHISTORY  )
  IF (ISTATUS /= NF90_NOERR) CALL IO_Err_handle_nc4(istatus,'IO_History_append_nc4','NF90_PUT_ATT','history')
END IF

END SUBROUTINE IO_History_append_nc4


subroutine IO_Select_split_file( tpfile, tpfield, tpfileout, tpfieldout, kvertlevel, kzfile )
use modd_field,          only: NMNHDIM_LEVEL, NMNHDIM_LEVEL_W, NMNHDIM_NLWB, NMNHDIM_NSWB, NMNHDIM_NOTLISTED, NMNHDIM_UNKNOWN

type(tfiledata),  target,         intent(in)  :: tpfile
class(tfieldmetadata), target,    intent(in)  :: tpfield
class(tfieldmetadata), pointer,   intent(out) :: tpfieldout
type(tfiledata),  pointer,        intent(out) :: tpfileout
integer,                optional, intent(in)  :: kvertlevel ! Number of the vertical level (needed for Z-level split files)
integer,                optional, intent(in)  :: kzfile     ! Number of the Z-level split file

character(len=4) :: ysuffix
integer          :: ivertdim
integer          :: jdim

if ( Present( kvertlevel ) ) then
  if ( kvertlevel > 9999 ) call Print_msg( NVERB_FATAL, 'IO', 'IO_Select_split_file','too many vertical levels' )
  if ( .not. Present( kzfile ) ) call Print_msg( NVERB_FATAL, 'IO', 'IO_Select_split_file', 'kzfile argument not provided' )
  if ( kzfile > tpfile%nsubfiles_ioz ) call Print_msg( NVERB_FATAL, 'IO', 'IO_Select_split_file', 'kzfile value too high' )

  Write( ysuffix, '( i4.4 )' ) kvertlevel
  tpfileout => tpfile%tfiles_ioz(kzfile)%tfile
  !Copy the values of tpfield to the pointer tpfieldout
  Allocate( tpfieldout, source = tpfield )
  tpfieldout%cmnhname  = Trim( tpfieldout%cmnhname ) // ysuffix
  if ( Len_trim( tpfieldout%cstdname  ) > 0 )  tpfieldout%cstdname  = Trim( tpfieldout%cstdname  ) // '_at_level_' // ysuffix
  if ( Len_trim( tpfieldout%clongname ) > 0 )  tpfieldout%clongname = Trim( tpfieldout%clongname ) // ' at level ' // ysuffix
  tpfieldout%ndims = 2

  !Remove the vertical dimension
  ivertdim = 0
  do jdim = 1, size(tpfieldout%ndimlist)
    if (      tpfieldout%ndimlist(jdim) == NMNHDIM_LEVEL     &
         .or. tpfieldout%ndimlist(jdim) == NMNHDIM_LEVEL_W   &
         .or. tpfieldout%ndimlist(jdim) == NMNHDIM_NLWB      &
         .or. tpfieldout%ndimlist(jdim) == NMNHDIM_NSWB      &
         .or. tpfieldout%ndimlist(jdim) == NMNHDIM_NOTLISTED ) then
      if ( jdim < size(tpfieldout%ndimlist) ) &
        tpfieldout%ndimlist(jdim:size(tpfieldout%ndimlist)-1) = tpfieldout%ndimlist(jdim+1:)
      tpfieldout%ndimlist(size(tpfieldout%ndimlist)) = NMNHDIM_UNKNOWN
      ivertdim = jdim
      exit !max 1 vertical dimension
    end if
  end do

  if ( ivertdim == 0 ) then
    call Print_msg( NVERB_FATAL, 'IO', 'IO_Select_split_file', 'vertical dimension not found for field ' // trim(tpfield%cmnhname) )
  else if ( ivertdim /= 3 ) then
    !Note: this situation is unexpected for the moment. If it is OK, please remove this warning
    call Print_msg( NVERB_WARNING, 'IO', 'IO_Select_split_file', 'removed vertical dimension was not at 3rd position for field ' &
                    // trim(tpfield%cmnhname) )
  end if
else
  tpfileout  => tpfile
  tpfieldout => tpfield
endif

end subroutine IO_Select_split_file


end module mode_io_write_nc4
#else
!
! External dummy subroutines
!
subroutine IO_Coordvar_write_nc4(a, b)
use mode_msg
integer :: a, b
CALL PRINT_MSG(NVERB_ERROR,'IO','IO_Coordvar_write_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine IO_Coordvar_write_nc4
!
subroutine IO_Field_write_nc4(a, b, c, d, e, f, g)
use mode_msg
integer :: a, b, c, d, e, f, g
CALL PRINT_MSG(NVERB_ERROR,'IO','IO_Field_write_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine IO_Field_write_nc4
!
subroutine IO_Header_write_nc4(a)
use mode_msg
integer :: a
CALL PRINT_MSG(NVERB_ERROR,'IO','IO_Header_write_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine IO_Header_write_nc4
!
#endif
