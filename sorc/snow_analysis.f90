MODULE M_Snow_Analysis

USE NETCDF
USE M_DA
!USE MPI
Use, Intrinsic :: IEEE_ARITHMETIC	

CONTAINS

subroutine calculate_IMS_fSCA(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC, IVEGSRC, &  
                                num_assim_steps, dT_Asssim,  & 
                                max_num_nearIMS, num_subgrd_ims_cels, &
                                SFC_FORECAST_PREFIX, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH)
                                                        
        !----------------------------------------------------------------------
        ! Input arguments: 
        ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
        ! IY, IM, ID, IH = year, month, day, hour of current model step   
        ! MYRANK: rank/id of the MPI process
        ! ...
        !
        ! Inputs, read from file:
        ! RLA, RLO: lat lon information for the tile
        !                   
        !----------------------------------------------------------------------
        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &
                                  IY, IM, ID, IH, LENSFC, IVEGSRC, num_assim_steps 
        CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
!  SFC_FORECAST_PREFIX = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Forecast/C768/snow."        
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
! If (IDIM == 96) then          
!         num_subgrd_ims_cels = 627   ! number of ims sub-grids            
!         bkgst_srch_rad = 240.	      ! Km radius of background search for obs location
! elseif (IDIM == 128) then
!         num_sub = 627               
!         max_distance = 240.		!Km 
! elseif (IDIM == 768) then
!         num_sub = 30
!         max_distance = 27.          !Km  
        REAL, intent(In)    :: dT_Asssim    
        INTEGER, intent(in) :: max_num_nearIMS, num_subgrd_ims_cels 
        CHARACTER(LEN=5)    :: TILE_NUM
        Character(LEN=3)    :: rank_str
        INTEGER             :: IERR     
        INTEGER             :: LANDMASK(LENSFC)
        CHARACTER(len=250)   :: ims_inp_file, ims_inp_file_indices
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile              
        REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:)
        REAL                 :: lat_min, lat_max, lon_min, lon_max      
        Real                 :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
        
        REAL                :: SNDFCS(LENSFC), SWEFCS(LENSFC), VETFCS(LENSFC)
        REAL                :: SNODENS_Grid(LENSFC), snodens
        REAL                :: SNO_IMS_at_Grid(LENSFC), SNUP_Array(LENSFC)

        REAL                :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC) !, OROG_UF(LENSFC)
        INTEGER             :: Num_Ims
        INTEGER             :: jndx, zndx, ncol, nrow
        Integer                 :: ims_assm_hour

        CHARACTER(len=250)       :: forc_inp_file, da_out_file
        CHARACTER(LEN=3)         :: RANKCH 

        Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice
    ! for mpi par
        INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
        INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiReal_size, rsize
        REAL               :: tmp
        INTEGER            :: istep, IY_loc, IM_loc, ID_loc, IH_loc
        REAL               :: IH_real
        INTEGER, PARAMETER :: PRINTRANK = 4 
! CSD-todo, should be same as in sfcsub. Share properly
        real, parameter :: nodata_val = -9999.9

!=============================================================================================
! 1. initialise vars,set-up processors, and read lat/lon from orog files.
!=============================================================================================
        !obs_srch_rad = 250. ! radius of observation search
        ims_assm_hour = 18 
        ! noah models specific? Needed to ID glaciers.
        if (IVEGSRC == 2) then   ! sib
                veg_type_landice=13
        else
                veg_type_landice=15
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
        IF (myrank ==PRINTRANK) PRINT*,"IMS fSCA: total num proc ", NPROCS, " Num tiles : ", NUM_TILES

        Np_ext = MOD(NPROCS, NUM_TILES)  ! extra/inactive procs
        if (MYRANK >  NPROCS - Np_ext - 1) goto 999
        Np_til = NPROCS / NUM_TILES  ! num proc. per tile 
        p_tN = MOD(MYRANK, NUM_TILES)  ! tile for proc.
        p_tRank = MYRANK / NUM_TILES  ! proc. rank within tile
        N_sA = LENSFC / Np_til  ! sub array length per proc
        N_sA_Ext = LENSFC - N_sA * Np_til ! extra grid cells
        if(p_tRank == 0) then 
                mp_start = 1
        else
                mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        endif
        mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc        
        If(myrank == PRINTRANK )PRINT*,"snowDA: sub array length ", N_sA, " extra sub array: ", N_sA_Ext
! if (p_tN /= 4 ) goto 999

! must have IMS obs. If no obs files: go to end 
        if ((IMS_SNOWCOVER_PATH(1:8).eq.'        ') .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
                print*, "Observation paths don't exist!, skipping"
                goto 999
        end if

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
        CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,LENSFC) !OROG_UF,
        IF (p_tRank==PRINTRANK) PRINT*,"IMS fSCA on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
        
        RLO_Tile = RLO ! copy so that RLO used later is not modified
        Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
        lat_min = MAX(MINVAL(RLA) - 1., -90.) ! CSD - why is the limit 1?
        lat_max = MIN(MAXVAL(RLA) + 1., 90.)
        lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
        lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)      
        if (p_tN==3)  then  
                lon_min = 125.       ! lon_min is left, lon_max is right, not necessary min/max value
                lon_max = -145.
        endif   
        if ((p_tRank==PRINTRANK) ) then !.and. print_deb) then
                print*, TILE_NUM, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
        endif

! If multiple time steps are simulated in a given time period (window) given by num_assim_steps * dT_Assim
! Note: the DA outputs from the last time step are returned        
        IH_real = IH; IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH; ! these variables change inside loop below
    Do istep = 1, num_assim_steps        
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc
        write(fvs_tile, "(I3)") IDIM

!=============================================================================================
! 2. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
!=============================================================================================

       ! READ THE INPUT SURFACE DATA ON THE CUBED-SPHERE TILE p_tN. A
       ! Also get vegtype (VETFCS) to identify glacier.    
        if (SFC_FORECAST_PREFIX(1:8).eq.'        ') then
                ! FNBGSI = "./fnbgsi." // RANKCH
                WRITE(RANKCH, '(I3.3)') (p_tN+1)
                forc_inp_file = "./fnbgsi." // RANKCH
        else
                forc_inp_file = TRIM(SFC_FORECAST_PREFIX)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"."//TRIM(h_str)// &
                                     "0000.sfc_data."//TILE_NUM//".nc"
        end if

        if (MYRANK==PRINTRANK) PRINT *, 'reading model backgroundfile', trim(forc_inp_file) 
                                     
        Call READ_Forecast_Data_atPath(forc_inp_file, veg_type_landice, LENSFC, SWEFCS, SNDFCS, &
                                       VETFCS, LANDMASK)

        ! grid snow density
        SNODENS_Grid = SWEFCS/SNDFCS
        WHERE (SNODENS_Grid < 0.0001) SNODENS_Grid = 0.

        if (COUNT (LANDMASK==1 .and. SNDFCS> 0.01) > 0) then 
                ! mean density over snow-covered land
                snodens = SUM(SNODENS_Grid, Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01) 
                if (MYRANK==PRINTRANK) & 
                        PRINT *, 'snowDA: mean density ', snodens 
        else 
                snodens = 0.1  ! default value if have no snow in tile
                if (MYRANK==PRINTRANK) & 
                        PRINT *, 'snowDA: no snow in current tiles, using default density ', snodens 
        endif
        ! PRINT *, 'snowDA:density ', MYRANK, snodens 
        ! for grid cells with no valid density, fill in the average snodens
        Where(.not.(LANDMASK==1 .and. SNDFCS>0.01 )) SNODENS_Grid = snodens
        !If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        !If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        !If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp

        If (MYRANK==PRINTRANK) PRINT *, 'reading IMS Snow cover'
        ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
        ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                                    TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
        ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
        ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                                    ".IMS.Indices."//TRIM(TILE_NUM)//".nc"                       
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading IMS file', trim(ims_inp_file) 
        Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
                                                    MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
        if((p_tRank==0) .and. (p_tN==PRINTRANK) .and. print_deb) then
                PRINT*, "SNCOV from rank: ", MYRANK
                PRINT*, SNCOV_IMS
        endif 
        if (myrank==PRINTRANK) PRINT*,'Finished reading SNCOV, converting to snow depth' 
 
        ! SNUP array will be used later, to test whether SCF > 100%
        call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC, SNO_IMS_at_Grid, SNUP_Array)

        SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_Grid ! convert SWE to SND
        
        if ((p_tN==PRINTRANK) .and. (p_tRank==0) .and. print_deb) then
                PRINT*, "SNCOV obs at each grid cell from rank: ", MYRANK
                PRINT*, SNO_IMS_at_Grid
        endif
        if (myrank==PRINTRANK) PRINT*,'Finished converting SNCOV observations'

        ! write outputs 
        Write(rank_str, '(I3.3)') (MYRANK+1)
        ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/"
        da_out_file = "./IMSfSCA."// &  !
                                  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !  
        call Write_fSCA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
                              LANDMASK, SNCOV_IMS, SNO_IMS_at_Grid)  !, anl_fSCA) !updated snocov

998 CONTINUE
        Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==PRINTRANK) PRINT*,'Finished OI DA at datetime: ', y_str, m_str, d_str, h_str
        
    End do

999 CONTINUE
        !PRINT*,'Finished OI DA ON RANK: ', MYRANK
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        !STOP

        RETURN

 END subroutine calculate_IMS_fSCA

  ! SWE threshold for 100% snow cover
 subroutine get_SWE_Threshold(VETFCS_in, SNUP)
        
        IMPLICIT NONE
        !
        Real, Intent(In)        :: VETFCS_in
        Real, Intent(Out)       :: SNUP
        
        INTEGER            :: VETFCS
        REAL               :: snupx(30)

        !This is for the IGBP veg classification scheme.
        ! SWE at which snow cover reaches 100%
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

        VETFCS = INT(VETFCS_in)
        If (VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
        
        SNUP = snupx(VETFCS) * 1000.0 ! mm
           
        RETURN
            
 END Subroutine get_SWE_Threshold
 

 Subroutine update_snow_cover_fraction(LENSFC, SNOANL, VETFCS_in, anl_fSCA)

	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: LENSFC
	REAL, intent(In)   :: SNOANL(LENSFC), VETFCS_in(LENSFC)
	REAL, intent(Out)   ::  anl_fSCA(LENSFC)
	INTEGER                :: VETFCS(LENSFC)

	REAL               :: snupx(30), SNEQV(LENSFC), SNUP, SALP, RSNOW
	Integer 		   :: indx, vtype_int

	!This is for the IGBP veg classification scheme.
	snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020, 	&
			0.020, 0.060, 0.040, 0.020, 0.010, 0.020,			&
			0.020, 0.020, 0.013, 0.013, 0.010, 0.020,			&
			0.020, 0.020, 0.000, 0.000, 0.000, 0.000,			&
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

	SNEQV = 0.001 * SNOANL   ! units mm->m
	SALP = -4.0
	
	VETFCS = INT(VETFCS_in)
	Where(VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
	
	Do indx=1, LENSFC
		SNUP = snupx(VETFCS(indx))
		if (SNUP == 0.) then
			print*, " 0.0 snup value, check vegclasses", vtype_int
			Stop
		endif

		IF (SNEQV(indx) .LT. SNUP) THEN
			RSNOW = SNEQV(indx)/SNUP
			anl_fSCA(indx) = 1. - (EXP(SALP*RSNOW) - RSNOW*EXP(SALP))
		ELSE
			anl_fSCA(indx) = 1.0
		ENDIF

		if (SNEQV(indx) < 0.00001)  anl_fSCA(indx) = 0.0	

	End do
	
	RETURN

 End Subroutine update_snow_cover_fraction


  ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_fSCA_Outputs(output_file, idim, jdim, lensfc, myrank,   &
                                 landmask, SNCOV_IMS, SNO_IMS_at_Grid)  !, anl_fSCA) !updated snocov
        !------------------------------------------------------------------
        !------------------------------------------------------------------
        implicit none

        CHARACTER(LEN=*), Intent(In)      :: output_file
        integer, intent(in)         :: idim, jdim, lensfc, myrank
        integer, intent(in)         :: landmask(lensfc)
        Real, intent(in)            ::  SNCOV_IMS(lensfc), SNO_IMS_at_Grid(lensfc)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        integer                     :: error, i, ncid
        integer                     :: dim_x, dim_y, dim_time
        integer                     :: id_x, id_y, id_time
        integer       :: id_imscov, id_imssnd, id_landmask
 
        real(kind=4)                :: times
        real(kind=4), allocatable   :: x_data(:), y_data(:)
        real(kind=8), allocatable   :: dum2d(:,:)

        include "mpif.h"

        ! print*
        ! print*,"Process ", myrank, "writing output data to: ",trim(output_file)

        !--- create the file
        error = NF90_CREATE(output_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
        call netcdf_err(error, 'CREATING FILE='//trim(output_file) )

        !--- define dimensions
        error = nf90_def_dim(ncid, 'xaxis_1', idim, dim_x)
        call netcdf_err(error, 'DEFINING XAXIS DIMENSION' )
        error = nf90_def_dim(ncid, 'yaxis_1', jdim, dim_y)
        call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
        error = nf90_def_dim(ncid, 'Time', 1, dim_time)
        call netcdf_err(error, 'DEFINING TIME DIMENSION' )

        !--- define fields
        error = nf90_def_var(ncid, 'xaxis_1', NF90_FLOAT, dim_x, id_x)
        call netcdf_err(error, 'DEFINING XAXIS_1 FIELD' )
        error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
        call netcdf_err(error, 'DEFINING XAXIS_1 LONG NAME' )
        error = nf90_put_att(ncid, id_x, "units", "none")
        call netcdf_err(error, 'DEFINING XAXIS_1 UNITS' )
        error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
        call netcdf_err(error, 'WRITING XAXIS_1 FIELD' )

        error = nf90_def_var(ncid, 'yaxis_1', NF90_FLOAT, dim_y, id_y)
        call netcdf_err(error, 'DEFINING YAXIS_1 FIELD' )
        error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
        call netcdf_err(error, 'DEFINING YAXIS_1 LONG NAME' )
        error = nf90_put_att(ncid, id_y, "units", "none")
        call netcdf_err(error, 'DEFINING YAXIS_1 UNITS' )
        error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
        call netcdf_err(error, 'WRITING YAXIS_1 FIELD' )

        error = nf90_def_var(ncid, 'Time', NF90_FLOAT, dim_time, id_time)
        call netcdf_err(error, 'DEFINING TIME FIELD' )
        error = nf90_put_att(ncid, id_time, "long_name", "Time")
        call netcdf_err(error, 'DEFINING TIME LONG NAME' )
        error = nf90_put_att(ncid, id_time, "units", "time level")
        call netcdf_err(error, 'DEFINING TIME UNITS' )
        error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
        call netcdf_err(error, 'WRITING TIME FIELD' )

        dims_3d(1) = dim_x
        dims_3d(2) = dim_y
        dims_3d(3) = dim_time

        error = nf90_def_var(ncid, 'LandMask', NF90_INT, dims_3d, id_landmask)
        call netcdf_err(error, 'DEFINING LandMask' )
        error = nf90_put_att(ncid, id_landmask, "long_name", "Masl: 1=non glacier land")
        call netcdf_err(error, 'DEFINING LandMask LONG NAME' )
        error = nf90_put_att(ncid, id_landmask, "units", "binary")
        call netcdf_err(error, 'DEFINING LandMask UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        error = nf90_def_var(ncid, 'imsSND', NF90_DOUBLE, dims_3d, id_imssnd)
        call netcdf_err(error, 'DEFINING imsSND' )
        error = nf90_put_att(ncid, id_imssnd, "long_name", "IMS Snow Depth")
        call netcdf_err(error, 'DEFINING imsSND LONG NAME' )
        error = nf90_put_att(ncid, id_imssnd, "units", "mm")
        call netcdf_err(error, 'DEFINING imsSND UNITS' )

        error = nf90_enddef(ncid, header_buffer_val,4,0,4)
        call netcdf_err(error, 'DEFINING HEADER' )

        allocate(x_data(idim))
        do i = 1, idim
        x_data(i) = float(i)
        enddo
        allocate(y_data(jdim))
        do i = 1, jdim
        y_data(i) = float(i)
        enddo
        times = 1.0

        error = nf90_put_var( ncid, id_x, x_data)
        call netcdf_err(error, 'WRITING XAXIS RECORD' )
        error = nf90_put_var( ncid, id_y, y_data)
        call netcdf_err(error, 'WRITING YAXIS RECORD' )
        error = nf90_put_var( ncid, id_time, times)
        call netcdf_err(error, 'WRITING TIME RECORD' )

        allocate(dum2d(idim,jdim))
        dims_strt(1:3) = 1
        dims_end(1) = idim
        dims_end(2) = jdim
        dims_end(3) = 1
        
        dum2d = reshape(landmask, (/idim,jdim/))
        error = nf90_put_var( ncid, id_landmask, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING LandMask RECORD' )

        dum2d = reshape(SNCOV_IMS, (/idim, jdim/))
        error = nf90_put_var(ncid, id_imscov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD')

        dum2d = reshape(SNO_IMS_at_Grid, (/idim, jdim/))        
        error = nf90_put_var(ncid, id_imssnd, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsSND RECORD')


        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 End subroutine Write_fSCA_Outputs

 SUBROUTINE Find_Nearest_GridIndices_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, &
                                            num_src, num_tar, RLA_cg, RLO_cg, RLA, RLO, index_ens_atGrid)
                            
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"
        
        Integer, Intent(In)          :: Myrank, NUM_TILES, p_tN, p_tRank, Np_til, num_src, num_tar
        Real, Intent(In)                :: RLA_cg(num_src), RLO_cg(num_src), RLA(num_tar), RLO(num_tar)
        Integer, Intent(Out)         :: index_ens_atGrid(num_tar)
        Real                        :: RLA_cg_rad(num_src), RLO_cg_rad(num_src)
        Real                        :: RLA_rad(num_tar), RLO_rad(num_tar)
        
        INTEGER                     :: indx, min_indx
        Real                        :: distArr(num_src), haversinArr(num_src)
        Real                        :: d_latArr(num_src), d_lonArr(num_src)
        Real(16), Parameter         :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter         :: pi_div_180 = PI_16/180.0
        Real, Parameter                 :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiInt_size, isize, IERR

        !Np_til ! num proc. per tile p_tRank ! proc. rank within tile !p_tN  ! tile for proc.
        N_sA = num_tar / Np_til  ! sub array length per proc
        N_sA_Ext = num_tar - N_sA * Np_til ! extra grid cells
        if(p_tRank == 0) then 
                mp_start = 1
        else
                mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        endif
        mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc
                
        index_ens_atGrid = -1  
        ! at each target point compute its distance from source RLA/RLO pairs and find the position of the minimum      
        RLA_rad =  pi_div_180 * RLA
        RLO_rad =  pi_div_180 * RLO
        RLA_cg_rad =  pi_div_180 * RLA_cg
        RLO_cg_rad =  pi_div_180 * RLO_cg       
        ! https://en.wikipedia.org/wiki/Haversine_formula
        Do indx = mp_start, mp_end    ! num_tar 
                d_latArr = (RLA_rad(indx) - RLA_cg_rad) / 2.
                d_lonArr = (RLO_rad(indx) - RLO_cg_rad) / 2.
                haversinArr = sin(d_latArr)**2 + cos(RLA_rad(indx)) * cos(RLA_cg_rad) * sin(d_lonArr)**2
                WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
                distArr = 2 * earth_rad * asin(sqrt(haversinArr))               
                min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
                index_ens_atGrid(indx) = min_indx
        end do

        isize = SIZEOF(N_sA) 
        Call MPI_TYPE_SIZE(MPI_INTEGER, mpiInt_size, IERR) 
        If (isize == 2 ) then 
                mpiInt_size = MPI_INTEGER2
        elseif (isize == 4 ) then 
                mpiInt_size = MPI_INTEGER4
        elseif (isize == 8 ) then 
                mpiInt_size = MPI_INTEGER8
        else
                PRINT*," Possible mismatch between Fortran Int ", isize," and Mpi Int ", mpiInt_size
                Stop
        endif
    
        if (MYRANK > (NUM_TILES - 1) ) then
                call MPI_SEND(index_ens_atGrid(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
                                                MYRANK*1000, MPI_COMM_WORLD, IERR)
        else !if (MYRANK == p_tN ) then  
                Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                        dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                        send_proc = MYRANK +  pindex * NUM_TILES
                        call MPI_RECV(index_ens_atGrid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
                                                send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                enddo
        endif
    !ToDO: better way to do this?
        ! now share the whole array
        if (MYRANK < NUM_TILES ) then   !if (MYRANK == p_tN ) then      
                Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
                        rec_proc = MYRANK +  pindex * NUM_TILES
                        call MPI_SEND(index_ens_atGrid, num_tar, mpiInt_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
                enddo
        else 
                call MPI_RECV(index_ens_atGrid, num_tar, mpiInt_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        endif
             
    RETURN
        
 END SUBROUTINE Find_Nearest_GridIndices_Parallel

 SUBROUTINE READ_Forecast_Data(MYRANK, LENSFC, veg_type_landice, SWEFCS, SNDFCS, VETFCS, LANDMASK)
    
        IMPLICIT NONE

        include "mpif.h"

        INTEGER, INTENT(IN)       :: MYRANK, LENSFC, veg_type_landice
        REAL, INTENT(OUT)         :: SWEFCS(LENSFC),SNDFCS(LENSFC),VETFCS(LENSFC)
        INTEGER, INTENT(OUT)      :: LANDMASK(LENSFC) 

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR, i

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)

        WRITE(RANKCH, '(I3.3)') (MYRANK+1)

        FNBGSI = "./fnbgsi." // RANKCH
        if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        ERROR=NF90_OPEN(TRIM(FNBGSI),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNBGSI) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM,JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/)) 

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING vtype' )
        VETFCS = RESHAPE(DUMMY, (/LENSFC/))    

        ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING slmsk' )
        SLMASK = RESHAPE(DUMMY, (/LENSFC/))    

        do i = 1, LENSFC 
           ! if land, but not land ice, set mask to 1.
           if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
                ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
                LANDMASK(i) = 1 
           else 
                LANDMASK(i) = 0
           endif 
        enddo

        ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data

 SUBROUTINE READ_Forecast_Data_atPath(forc_inp_path, veg_type_landice, LENSFC, SWEFCS, SNDFCS, VETFCS, LANDMASK) !VEGFCS, !SRFLAG)
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: forc_inp_path
        INTEGER, INTENT(IN)               :: LENSFC, veg_type_landice
        REAL, INTENT(OUT)                 :: SWEFCS(LENSFC), SNDFCS(LENSFC), VETFCS(LENSFC)  !VEGFCS(LENSFC), 
        INTEGER, INTENT(OUT)              :: LANDMASK(LENSFC) 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID, i
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! WRITE(RANKCH, '(I3.3)') (MYRANK+1)
        ! FNBGSI = "./fnbgsi." // RANKCH
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(forc_inp_path)
        
        INQUIRE(FILE=trim(forc_inp_path), EXIST=file_exists)

        if (.not. file_exists) then 
                print *, 'READ_Forecast_Data_atPath error,file does not exist', &   
                        trim(forc_inp_path) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
                

        ERROR=NF90_OPEN(TRIM(forc_inp_path), NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(forc_inp_path) )
        ! ERROR=NF90_OPEN(TRIM(FNBGSI),NF90_NOWRITE,NCID)
        ! CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNBGSI) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM,JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING vtype' )
        VETFCS = RESHAPE(DUMMY, (/LENSFC/))    

        ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING slmsk' )
        SLMASK = RESHAPE(DUMMY, (/LENSFC/))    

        do i = 1, LENSFC 
           ! if land, but not land ice, set mask to 1.
           if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
                ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
                LANDMASK(i) = 1 
           else 
                LANDMASK(i) = 0
           endif 
        enddo

        ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Forecast_Data_atPath

 SUBROUTINE READ_Analysis_Data(anl_inp_path, LENSFC, SWEFCS, SNDFCS) 
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: anl_inp_path
        INTEGER, INTENT(IN)       :: LENSFC     !MYRANK, 
        REAL, INTENT(OUT)         :: SWEFCS(LENSFC), SNDFCS(LENSFC)

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)

        ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(anl_inp_path)

        ERROR=NF90_OPEN(TRIM(anl_inp_path),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(anl_inp_path) )

        ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
        CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

        ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
        CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

        IF ((IDIM*JDIM) /= LENSFC) THEN
        PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM, JDIM))

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        SWEFCS = RESHAPE(DUMMY, (/LENSFC/))

        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        SNDFCS = RESHAPE(DUMMY, (/LENSFC/)) 
 
        
        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Analysis_Data
    
 SUBROUTINE READ_LAT_LON_CoarseRes(inp_path,RLA,RLO,IDIM,JDIM,IJDIM)
    
     IMPLICIT NONE
    
     include "mpif.h"
    
     INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM
         
         CHARACTER(LEN=*), Intent(In)      :: inp_path
    
     REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
    
     INTEGER                :: ERROR, NCID
     INTEGER                :: I, II, J, JJ
     INTEGER                :: ID_DIM, ID_VAR, NX, NY
    
     REAL, ALLOCATABLE      :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
    
     !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
    
     if (print_deb) then
        PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(inp_path)
     endif
    
     ERROR=NF90_OPEN(TRIM(inp_path),NF90_NOWRITE,NCID)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_path) )
    
     ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )
    
     ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )
    
     IF ((NX/2) /= IDIM .OR. (NY/2) /= JDIM) THEN
       PRINT*,'FATAL ERROR: DIMENSIONS IN FILE: ',(NX/2),(NY/2)
       PRINT*,'DO NOT MATCH GRID DIMENSIONS: ',IDIM,JDIM
       CALL MPI_ABORT(MPI_COMM_WORLD, 130)
     ENDIF
    
     ALLOCATE(GEOLON(NX+1,NY+1))
     ALLOCATE(GEOLAT(NX+1,NY+1))
    
     ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )
    
     ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )
    
     ALLOCATE(DUMMY(IDIM,JDIM))
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLON(II,JJ)
       ENDDO
     ENDDO
    
     RLO = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLON)
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLAT(II,JJ)
       ENDDO
     ENDDO
    
     RLA = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLAT, DUMMY)
    
     ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_LAT_LON_CoarseRes

 SUBROUTINE READ_LAT_LON_OROG_atRank(MYRANK, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,IJDIM) !OROG_UF,
    
    !--------------------------------------------------------------
    ! READ LATITUDE, LONGITUDE, FILTERED OROGRAPHY, AND
    ! UNFILTERED OROGRAPHY FOR THE CUBED-SPHERE TILE FROM
    ! THE "GRID" FILE.
    !--------------------------------------------------------------
    
     IMPLICIT NONE
    
     include "mpif.h"
    
     INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM
    
     CHARACTER(LEN=5), INTENT(OUT) :: TILE_NUM
    
     REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
     REAL, INTENT(OUT)      :: OROG(IJDIM)   !,OROG_UF(IJDIM)
    
     CHARACTER(LEN=50)      :: FNOROG, FNGRID
     CHARACTER(LEN=3)       :: RANKCH
    
     INTEGER                :: ERROR, NCID, NCID_OROG
     INTEGER                :: I, II, J, JJ, MYRANK
     INTEGER                :: ID_DIM, ID_VAR, NX, NY
    
     REAL, ALLOCATABLE         :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
     REAL(KIND=4), ALLOCATABLE :: DUMMY4(:,:)
    
     !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
    
     WRITE(RANKCH, '(I3.3)') (MYRANK+1)
    
     FNGRID = "./fngrid." // RANKCH
    
     if (print_deb) then
        PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(FNGRID)
     endif
    
     ERROR=NF90_OPEN(TRIM(FNGRID),NF90_NOWRITE,NCID)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNGRID) )
    
     ERROR=NF90_INQ_DIMID(NCID, 'nx', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NX)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NX' )
    
     ERROR=NF90_INQ_DIMID(NCID, 'ny', ID_DIM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY ID' )
    
     ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NY)
     CALL NETCDF_ERR(ERROR, 'ERROR READING NY' )
    
     IF ((NX/2) /= IDIM .OR. (NY/2) /= JDIM) THEN
       PRINT*,'FATAL ERROR: DIMENSIONS IN FILE: ',(NX/2),(NY/2)
       PRINT*,'DO NOT MATCH GRID DIMENSIONS: ',IDIM,JDIM
       CALL MPI_ABORT(MPI_COMM_WORLD, 130)
     ENDIF
    
     ALLOCATE(GEOLON(NX+1,NY+1))
     ALLOCATE(GEOLAT(NX+1,NY+1))
    
     ERROR=NF90_INQ_VARID(NCID, 'x', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLON)
     CALL NETCDF_ERR(ERROR, 'ERROR READING X RECORD' )
    
     ERROR=NF90_INQ_VARID(NCID, 'y', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, GEOLAT)
     CALL NETCDF_ERR(ERROR, 'ERROR READING Y RECORD' )
    
     ALLOCATE(DUMMY(IDIM,JDIM))
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLON(II,JJ)
       ENDDO
     ENDDO
    
     RLO = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLON)
    
     DO J = 1, JDIM
       DO I = 1, IDIM
         II = 2*I
         JJ = 2*J
         DUMMY(I,J) = GEOLAT(II,JJ)
       ENDDO
     ENDDO
    
     RLA = RESHAPE(DUMMY, (/IJDIM/))
    
     DEALLOCATE(GEOLAT, DUMMY)
    
     ERROR=NF90_INQ_VARID(NCID, 'tile', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING TILE ID' )
     ERROR=NF90_GET_VAR(NCID, ID_VAR, TILE_NUM)
     CALL NETCDF_ERR(ERROR, 'ERROR READING TILE RECORD' )
    
     ERROR = NF90_CLOSE(NCID)
    
     FNOROG = "./fnorog." // RANKCH
    
     if (print_deb) PRINT*, "READ FV3 OROG INFO FROM: "//TRIM(FNOROG)
    
     ERROR=NF90_OPEN(TRIM(FNOROG),NF90_NOWRITE,NCID_OROG)
     CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(FNOROG) )
    
     ALLOCATE(DUMMY4(IDIM,JDIM))
    
!      ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_raw', ID_VAR)
!      CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw ID' )
!      ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
!      CALL NETCDF_ERR(ERROR, 'ERROR READING orog_raw RECORD' )
!      OROG_UF = RESHAPE(DUMMY4, (/IJDIM/))
    
     ERROR=NF90_INQ_VARID(NCID_OROG, 'orog_filt', ID_VAR)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt ID' )
     ERROR=NF90_GET_VAR(NCID_OROG, ID_VAR, DUMMY4)
     CALL NETCDF_ERR(ERROR, 'ERROR READING orog_filt RECORD' )
     OROG = RESHAPE(DUMMY4, (/IJDIM/))
    
     DEALLOCATE(DUMMY4)
    
     ERROR = NF90_CLOSE(NCID_OROG)
    
 END SUBROUTINE READ_LAT_LON_OROG_atRank

 ! 10.14.20: subroutine copied from UEBFortran at
 ! https://github.com/dtarb/UEBFortran/blob/master/Model/snowxv.f90
 !                 Update time for each time step
 SUBROUTINE UPDATEtime(YEAR,MONTH,DAY,HOUR,DT)
        
        IMPLICIT NONE

        INTEGER   :: YEAR, MONTH, DAY, DMON(12), DM, I       ! 30/03/2004 ITB 
        !INTEGER   :: LYEAR  ! 30/03/2004 ITB  
        Real      :: hour, dt  ! DGT Dec 10, 2004.  Fixing ITB errors 
 
        DATA (DMON(I),I=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
        HOUR=HOUR+DT
        DM=DMON(MONTH)
!  check for leap years 
        if(month .eq. 2)dm=lyear(year)
10   continue
        IF(HOUR.GE.24.0) THEN
          HOUR=HOUR-24.0
          DAY=DAY+1
          go to 10
        ENDIF
20   continue
        IF(DAY.GT.DM) THEN
          DAY=day - dm
          MONTH=MONTH+1
          IF(MONTH.GT.12) THEN
                MONTH=1
                YEAR=YEAR+1
                  DM=DMON(MONTH)
                if(month .eq. 2)dm=lyear(year)
                endif
                go to 20
        ENDIF
        RETURN

 END SUBROUTINE UPDATEtime

 ! 10.14.20: subroutine copied from UEBFortran at
 ! https://github.com/dtarb/UEBFortran/blob/master/Model/snowxv.f90
!    function to return number of days in February checking for leap years
 function lyear(year) result(Alyear)  
        
        IMPLICIT NONE

        Integer  :: year, Alyear
        IF(MOD(YEAR,4).GT.0 .or. &
                (mod(year,100) .eq.0 .and. mod(year,400) .ne. 0)) THEN
        ! Leap years are every 4 years 
        ! - except for years that are multiples of centuries (e.g. 1800, 1900)
        ! - except again that when the century is divisible by 4 (e.g. 1600, 2000)
        !   then it is a leap year 
                Alyear=28
        ELSE
                Alyear=29
        ENDIF
          
 end function lyear

 ! copied from UEB Fortran at 
 ! https://github.com/dtarb/UEBFortran/blob/master/Model/functions.f90
 !================================================================
    SUBROUTINE JULDAT (I, M, K, H, TJD)
        !THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
        !time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT time VALUE
        !CAN BE IN ANY UT-LIKE time SCALE (UTC, UT1, TT, ETC.) - OUTPUT
        !JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
        !VAN FLANDERN.
        !SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
        !I = YEAR (IN)
        !M = MONTH NUMBER (IN)
        !K = DAY OF MONTH (IN)
        !H = UT HOURS (IN)
        !TJD = JULIAN DATE (OUT)
        Implicit None
        DOUBLE PRECISION H,TJD,JD
        Integer:: I,M,K
        !JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
        JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12-3*((I+4900+(M-14)/12)/100)/4
        TJD = JD - 0.5D0 + H/24.D0
        RETURN
    END SUBROUTINE JULDAT
    !================================================================
    !================================================================
    SUBROUTINE CALDAT (TJD, I, M, K, H)
        !THIS SUBROUTINE COMPUTES CALENDAR DATE AND time, GIVEN JULIAN
        !DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
        !(UTC, UT1, TT, ETC.) - OUTPUT time VALUE WILL HAVE SAME BASIS.
        !OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
        !VAN FLANDERN.
        !SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
        !
        !TJD = JULIAN DATE (IN)
        !I = YEAR (OUT)
        !M = MONTH NUMBER (OUT)
        !K = DAY OF MONTH (OUT)
        !H = UT HOURS (OUT)
        Implicit None
        DOUBLE PRECISION TJD,H,DJD,DMOD,JD
        Integer L,N,I,M,K
        DJD = TJD + 0.5D0
        JD = DJD
        H = DMOD (DJD,1.D0)*24 ! 24.D0
        !JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
        L = JD + 68569
        N = 4*L/146097
        L = L - (146097*N+3)/4
        !I=YEAR, M=MONTH, K=DAY
        I = 4000*(L+1)/1461001
        L = L - 1461*I/4 + 31
        M = 80*L/2447
        K = L - 2447*M/80
        L = M / 11
        M = M + 2 - 12*L
        I = 100*(N-49) + I + L
        RETURN

    END SUBROUTINE CALDAT    

 END MODULE M_Snow_Analysis
 
