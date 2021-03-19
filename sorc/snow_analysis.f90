MODULE M_Snow_Analysis

USE NETCDF
USE M_DA
Use, Intrinsic :: IEEE_ARITHMETIC       

CONTAINS

 subroutine Snow_Analysis_OI(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  
                                LENSFC, IVEGSRC,  &
                                L_horz , h_ver, obs_tolerance, &
                                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, &
                                ims_max_ele, num_subgrd_ims_cels, &
                                assim_snowdepth_stn, assim_snowdepth_ims, &
                                GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, &
                                IMS_INDEXES_PATH, SFC_FORECAST_PREFIX, &
                                SNDANL)
                                                        
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
        ! Outputs:
        ! SNDANL: snow depth analysis (units?)
        ! 
        ! Draper - changes to snow variables to clarify names, removed unnecesary conversions
        !          SWE - snow water equivalent 
        !          SND - snow depth 
        !                   
        !----------------------------------------------------------------------
        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &
                                  IY, IM, ID, IH, LENSFC, IVEGSRC
        LOGICAL, intent(in)    :: assim_snowdepth_stn, assim_snowdepth_ims 
        CHARACTER(LEN=*), Intent(In)   :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
        CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX
!  GHCND_SNOWDEPTH_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/GHCND/",
!  IMS_SNOWCOVER_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/",
!  IMS_INDEXES_PATH="/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/",
!  SFC_FORECAST_PREFIX = "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/fnbgsio/snow."
        Real, intent(In)    :: L_horz , h_ver, obs_tolerance
        Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
        INTEGER, intent(in) :: max_num_nearStn, num_subgrd_ims_cels 
        REAL, intent(Out)   :: SNDANL(LENSFC) 
        CHARACTER(LEN=5)    :: TILE_NUM
        Character(LEN=1)    :: rank_str
        INTEGER                     :: IERR     
        REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC)  
        REAL                :: SNDFCS(LENSFC), SWEFCS(LENSFC), SWEANL(LENSFC), SCFANL(LENSFC) 
        REAL                :: VETFCS(LENSFC), SNUP_Array(LENSFC)
        INTEGER             :: LANDMASK(LENSFC)
        CHARACTER(len=250)   :: dim_name, ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
        REAL, ALLOCATABLE    :: SNDOBS_stn(:), SNDFCS_at_stn(:), SNDANL_at_stn(:)
        REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROGFCS_at_stn(:)  
        REAL                 :: lat_min, lat_max, lon_min, lon_max      
        Real                 :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
        Real                 :: SND_IMS_at_Grid(LENSFC), SWE_IMS_at_Grid(LENSFC)

        INTEGER :: num_stn, Num_Ims
        INTEGER :: jndx, zndx, ncol, nrow
        Integer, Allocatable   :: index_back_at_nearStn(:), index_back_at_nearIMS(:) !loc_near_Obs(:), 
        Integer                :: num_loc, num_loc_1, num_loc_2
        Real, Parameter         :: stdev_back_depth = 30., stdev_obsv_stn = 40., stdev_obsv_ims = 80. ! mm 

        Real(dp), Allocatable    :: B_cov_mat(:,:), b_cov_vect(:)
        Real(dp), Allocatable    :: O_cov_mat(:,:), W_wght_vect(:)
        Real, Allocatable        :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:)
        Real, Allocatable        :: obs_Innov(:), OmB_innov_at_stn(:)

        CHARACTER(len=250)       :: forc_inp_file, da_out_file  !anl_inp_path,
        CHARACTER(LEN=3)         :: RANKCH 

        Integer, ALLOCATABLE    :: index_back_atObs(:)   ! the location of background corresponding obs
        Integer, ALLOCATABLE    :: array_index_back_atObs(:,:)   ! the location of background corresponding obs

        Real               :: snodens, SNODENS_Grid(LENSFC)
        LOGICAL            :: assim_sncov_thisGridCell    !assimilate sncov, 

        Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice
    ! for mpi par
        INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
        INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiReal_size, rsize
        REAL               :: tmp, gross_thold 
        INTEGER, PARAMETER :: PRINTRANK = 4 
! CSD-todo, should be same as in sfcsub. Share properly
        real, parameter :: nodata_val = -9999.9
        real, parameter :: nodata_tol = 0.0001

!=============================================================================================
! 1. initialise vars,set-up processors, and read lat/lon from orog files.
!=============================================================================================

        if( (.not. assim_snowdepth_stn) .and. (.not. assim_snowdepth_ims) ) then
                print*, "snowDA: No obs types selecteed, existing snow anal"  
                return
        end if
        !initialse output with nodata
        SNDANL= nodata_val

        !obs_srch_rad = 250. ! radius of observation search
        ! noah models specific? Needed to ID glaciers.
        if (IVEGSRC == 2) then   ! sib
                veg_type_landice=13
        else
                veg_type_landice=15
        endif

!total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
        IF (myrank ==PRINTRANK) PRINT*,"snowDA: OI total num proc ", NPROCS, " Num tiles : ", NUM_TILES 

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

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
        CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,LENSFC) !OROG_UF,
        if (print_deb) PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
        
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
        if ((p_tRank==0) .and. print_deb) then
                print*, TILE_NUM, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
        endif

        write(y_str, "(I4)") IY
        write(m_str, "(I0.2)") IM
        write(d_str, "(I0.2)") ID
        write(h_str, "(I0.2)") IH
        write(fvs_tile, "(I3)") IDIM

! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 

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

        if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading model backgroundfile', trim(forc_inp_file) 
                                     
        Call READ_Forecast_Data_atPath(forc_inp_file, veg_type_landice, LENSFC, SWEFCS, SNDFCS, &
                                       VETFCS, LANDMASK)

        ! get the snow density  = SWE/SND
        ! where snow present, use forecsts
        ! where snow not present, use average from snow forecasts over land

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
        PRINT *, 'snowDA:density ', MYRANK, snodens 

        ! for grid cells with no valid density, fill in the average snodens
        Where( SNDFCS<0.01 ) SNODENS_Grid = snodens

        If (p_tRank==0)  print*, "snowDA: tile ", p_tN, ' mean snow density', snodens
        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "snowDA: tile ", p_tN,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "snowDA: tile ", p_tN,  ' mean SND', tmp
        
!=============================================================================================
! 3. Read observations
!=============================================================================================

! 3a. Read station obs (of snow depth or SWE)

        if (assim_snowdepth_stn) then 
             ghcnd_inp_file = TRIM(GHCND_SNOWDEPTH_PATH)//"/GHCND.SNWD."// &
                                         TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"
             if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
             dim_name = "Site_Id"    
             ! here: all valid stn obs, within lat/lon box (can be outside of tile though) 
             call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name, &
                        lat_min, lat_max, lon_min, lon_max, & 
                        num_stn, SNDOBS_stn,              &
                        Lat_stn, Lon_stn, MYRANK) 

            if ((p_tRank==0) .and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if ((p_tRank==0) .and. (p_tN==2) .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNDOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
            endif
          PRINT*,'snowDA: read ', num_stn, ' obs on RANK', MYRANK

       endif

! CSD beyond this point, there should be no specific mention of the station data source

! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 

       if (assim_snowdepth_ims) then
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
            ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"/IMS.SNCOV."// &
                                       TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
            ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"/C"//TRIM(ADJUSTL(fvs_tile))//&
                                                       ".IMS.Indices."//TRIM(TILE_NUM)//".nc"                       
             if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading IMS file', trim(ims_inp_file) 
            Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
                                                       MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
            if((p_tRank==0) .and. (p_tN==1) .and. print_deb) then
                    PRINT*, "SNCOV from rank: ", MYRANK
                    PRINT*, SNCOV_IMS
            endif 

           ! SNUP array will be used later, to test whether SCF > 100%
            call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC, SWE_IMS_at_Grid, SNUP_Array)

            SND_IMS_at_Grid = SWE_IMS_at_Grid/SNODENS_Grid

            if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then
                    PRINT*, "SNCOV obs at each grid cell from rank: ", MYRANK
                    PRINT*, SND_IMS_at_Grid
            endif

            ! QC: ims -remove if fail elevation check
            !         -remove if model and IMS have snow cover

            where (OROG >  ims_max_ele) SND_IMS_at_Grid=nodata_val
            where ( (SWEFCS >= SNUP_Array)  .and. ( SNCOV_IMS > 0.99) ) SND_IMS_at_Grid=nodata_val
            if (myrank==PRINTRANK) PRINT*,'snowDA: Finished  reading and converting SNCOV observations'
        
        endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!      snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNDFCS_at_stn(num_stn))
            ALLOCATE(OROGFCS_at_stn(num_stn)) 
            ALLOCATE(index_back_atObs(num_stn)) 
            ALLOCATE(SNDANL_at_stn(num_stn))
            ALLOCATE(array_index_back_atObs(2,num_stn)) 
            ALLOCATE(OmB_innov_at_stn(num_stn)) 
 
            ! threshold for gross error check 
            gross_thold =  obs_tolerance * sqrt(stdev_back_depth**2 + stdev_obsv_stn**2)
            ! QC: remove non-land obs, and do gross-error check
            Call Observation_Operator_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, & 
                                RLA, RLO, OROG, Lat_stn, Lon_stn, SNDOBS_stn,   &
                                LENSFC, num_stn, bkgst_srch_rad, SNDFCS, LANDMASK, gross_thold,  &
                                SNDFCS_at_stn, OROGFCS_at_stn, index_back_atObs )  
            ! from here, invalid obs have no sndfcs_at_stn = NaN, and index_back_atObs = -1 
            ! Could speed up below by removing invalid obs
            ! At least, should add explicit assim_flag.

            OmB_innov_at_stn = SNDOBS_stn - SNDFCS_at_stn

            if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then
                    PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                    PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
                    PRINT*, "Model elevation at station locations from rank: ", MYRANK
                    PRINT*, OROGFCS_at_stn
                    PRINT*, "Background Indices at obs points"
                    PRINT*, index_back_atObs
                    PRINT*, "Background snow depth at station locations from rank: ", MYRANK
                    PRINT*, SNDFCS_at_stn  
                    PRINT*, "O - B (innovation at obs points)"
                    PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==PRINTRANK) PRINT*,'snowDA: Finished observation operator for station data'         
        endif ! num_stn > 0

!=============================================================================================
! 5.  obs QC goes here
!=============================================================================================

!CSDCSD - todo. Add QC here.

! QC steps:
! if station elevation >1500, discard stn_obs
! is model eleveation > threshold, discard IMS obs * 
! **GHCN has incomplete station info, and currently we're not reading it in. 
!   and not doing the station elevation check.
! Should we be reading it in, and discarding stations with no elevation? 

! min/max limits on station obs
! temperature check on all obs  (later if no temperature data with GHCN?)

! if abs (model - obs ) elevation > ??,  discard obs  (threshold, try 200 - 400 m-ish)

! QC obs for land cover (below QC's update cell, but not obs) 
! gross error check * 
! screen stn obs for IMS snow cover ???
! screen snow cover-derived snow depth (IMS) if model has snow  * 

!=============================================================================================
! 6. Perform the DA update, by looping over each grid cell 
!=============================================================================================
 
        !obs_srch_rad = 250. ! radius of observation search
        if (myrank==PRINTRANK) PRINT*,'snowDA: starting DA loop'
        Do jndx = mp_start, mp_end   
           ! QC: only update this grid cell if it is land.
           if( LANDMASK(jndx) == 1 ) then 
                num_loc_1 = 0
                num_loc_2 = 0
                assim_sncov_thisGridCell = .FALSE.
                if (print_deb) print*, "proc ", myrank, " grid: ", jndx
                if(num_stn>0) then 
                ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                ! QC: gross error check is done in this call.
                        call nearest_Observations_Locations(RLA(jndx), RLO(jndx),                    &
                                        Lat_stn, Lon_stn,  num_stn, obs_srch_rad, max_num_nearStn,   &
                                        stdev_back_depth, stdev_obsv_stn,                            &
                                        SNDFCS_at_stn, SNDOBS_stn,                                   &
                                        index_back_at_nearStn,  num_loc_1) !,      &LENSFC,
                        if (print_deb) print*, "number of stn sndpth obs ", num_loc_1
                endif         
                if( assim_snowdepth_ims  .and.   &
                    (.NOT. IEEE_IS_NAN(SNDFCS(jndx))) .AND. &
                    (.NOT. IEEE_IS_NAN(SND_IMS_at_Grid(jndx))) .AND. &
                    (  SND_IMS_at_Grid(jndx) > (nodata_val + nodata_tol) ) )  THEN  
                        num_loc_2 = 1 
                        assim_sncov_thisGridCell = .TRUE.                
                endif
                num_loc = num_loc_1 + num_loc_2                

                if((num_loc > 0) ) then  ! skip DA if no obs, use background
                        ! get background states
                        Allocate(back_at_Obs(num_loc))
                        Allocate(obs_Array(num_loc))
                        Allocate(Lat_Obs(num_loc))
                        Allocate(Lon_Obs(num_loc))
                        Allocate(orogfcs_Obs(num_loc))                        
                        if(num_loc_1 > 0) then
                                Do zndx = 1, num_loc_1     
                                        back_at_Obs(zndx) = SNDFCS_at_stn(index_back_at_nearStn(zndx))
                                        obs_Array(zndx) = SNDOBS_stn(index_back_at_nearStn(zndx))
                                        Lat_Obs(zndx) = Lat_stn(index_back_at_nearStn(zndx))
                                        Lon_Obs(zndx) = Lon_stn(index_back_at_nearStn(zndx))
                                        orogfcs_Obs(zndx) = OROGFCS_at_stn(index_back_at_nearStn(zndx)) 
                                End Do
                        End if

                        if(assim_sncov_thisGridCell) then                                
                                back_at_Obs(num_loc) = SNDFCS(jndx)
                                obs_Array(num_loc) = SND_IMS_at_Grid(jndx)
                                Lat_Obs(num_loc) = RLA(jndx)   
                                Lon_Obs(num_loc) = RLO(jndx) 
                                orogfcs_Obs(num_loc) = OROG(jndx)
                        endif
                        ! compute covariances
                        Allocate(B_cov_mat(num_loc, num_loc))
                        Allocate(b_cov_vect(num_loc))
                        Allocate(O_cov_mat(num_loc, num_loc))
                        Allocate(W_wght_vect(num_loc))   

                        call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), SNDFCS(jndx),    &
                                Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc,   &
                                stdev_back_depth, stdev_obsv_stn, stdev_obsv_ims,      &
                                L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                                assim_sncov_thisGridCell,                          &
                                B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)

                        Allocate(obs_Innov(num_loc))

                        call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
                                SNDFCS(jndx), SNDANL(jndx), obs_Innov)
                        if ((p_tN==4) .and. (p_tRank==0) .and. print_deb) then  !
                                print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
                                PRINT*, " background at obs pts: "
                                PRINT*, back_at_Obs     
                                PRINT*, "Observed"
                                PRINT*,  obs_Array
                                PRINT*, "Obs innovation: "
                                PRINT*, obs_Innov
                                PRINT*, "Weight vector: "
                                PRINT*, W_wght_vect     
                                print*, "forec: ", SNDFCS(jndx), " anl: ", SNDANL(jndx)
                        endif           
                        !free mem
                        DEALLOCATE(back_at_Obs, obs_Array)
                        DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs, obs_Innov)
                        DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
                else ! no obs were available,  keep background
                        SNDANL(jndx) = SNDFCS(jndx) 
                endif
                if (allocated(index_back_at_nearStn))  Deallocate(index_back_at_nearStn) 
                if (allocated(index_back_at_nearIMS))  Deallocate(index_back_at_nearIMS)                 
           endif ! not a land cell
        End do
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==PRINTRANK) PRINT*, 'snowDA: Finished DA loops'
        
!=============================================================================================
! 7. Clean up, write outputs 
!=============================================================================================


! collect results onto main tasks, if necessary

! ToDO: Better way to handle this? ! CSD - I'll fix this later.
! Real data type size corresponding to mpi
        rsize = SIZEOF(snodens)
        Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
        If (rsize == 4 ) then 
                mpiReal_size = MPI_REAL4
        elseif (rsize == 8 ) then 
                mpiReal_size = MPI_REAL8
        elseif (rsize == 16 ) then 
                mpiReal_size = MPI_REAL16
        else
                PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real ", mpiReal_size
                Stop
        endif
        ! send analyses arrays to 'tile-level root' proc.               
        if (MYRANK > (NUM_TILES - 1) ) then
                call MPI_SEND(SNDANL(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                                          MYRANK, MPI_COMM_WORLD, IERR) 
        else    !if(p_tRank == 0) then  
                Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                        dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                        send_proc = MYRANK +  pindex * NUM_TILES
                        call MPI_RECV(SNDANL(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
                                          send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                enddo
        endif
        if (myrank==PRINTRANK) PRINT*,'snowDA: Finished Data copy'

        if (MYRANK > NUM_TILES - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

         
        ! avoid -ve anl
        Where( (SNDANL < 0.) .and. (SNDANL>(nodata_val+nodata_tol)) ) SNDANL = 0.
        if (print_deb) then
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNDANL
        endif

! Calculate variables for diagnistic output file
! (only SNDANL is passed back to calling program) 

! Grid variable, analysis SWE and SCF
! note that in sfcsub, density in locations without current snow  is new snow density
! SWE results will be slightly different
        SWEANL = nodata_val
        Where( LANDMASK==1 ) SWEANL  = SNDANL*SNODENS_Grid
        
        !call calculate_snow_cover_fraction(LENSFC, SWEANL, VETFCS, nodata_val, nodata_tol, SCFANL)

! get index for station obs in original array 
 
        array_index_back_atObs = -9  
        do jndx = 1, num_stn 
           if   ( index_back_atObs(jndx) > 0. ) then 
              array_index_back_atObs(1,jndx) =  index_back_atObs(jndx) / JDIM + 1
              array_index_back_atObs(2,jndx) =  mod(index_back_atObs(jndx), JDIM) 
              if (array_index_back_atObs(2,jndx) ==0 ) array_index_back_atObs(2,jndx) = JDIM
           endif
        enddo


! station depth from analysis
! get analysis at obs locations, to O-A
        gross_thold = 999999999. ! set to large, so don't filter out any analyses
        if (num_stn > 0 ) then
                Call Observation_Operator_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, &
                            RLA, RLO, OROG, Lat_stn, Lon_stn, SNDOBS_stn,   &
                            LENSFC, num_stn, bkgst_srch_rad, SNDANL, LANDMASK, gross_thold,  & 
                            SNDANL_at_stn, OROGFCS_at_stn, index_back_atObs )
        endif

        ! write outputs 
        Write(rank_str, '(I1.1)') (MYRANK+1)
        da_out_file = "./SNOANLOI.tile"//rank_str//".nc"  !

        ! use array_index to detect stn obs not assimilated (later, add an explicit flag)
        call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
                              RLO, RLA, & 
                              SWEFCS, SWEANL, SNDFCS, SNDANL, LANDMASK,&  
                              SNCOV_IMS, SND_IMS_at_Grid, & 
                              num_stn, Lat_stn, Lon_stn, array_index_back_atObs, & 
                              SNDOBS_stn, SNDFCS_at_stn, SNDANL_at_stn) 

998 CONTINUE
    ! clean up
        if (allocated(SNDOBS_stn))      DEALLOCATE(SNDOBS_stn)
        if (allocated(SNDFCS_at_stn))   DEALLOCATE(SNDFCS_at_stn)
        if (allocated(SNDANL_at_stn))   DEALLOCATE(SNDANL_at_stn)
        if (allocated(OmB_innov_at_stn))   DEALLOCATE(OmB_innov_at_stn)
        if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
        if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
        if (allocated(OROGFCS_at_stn))  DEALLOCATE(OROGFCS_at_stn) 
        if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)
        if (allocated(array_index_back_atObs)) DEALLOCATE(array_index_back_atObs)
999 CONTINUE
        PRINT*,'snowDA: Finished OI DA ON RANK: ', MYRANK
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        !STOP

        RETURN

 END subroutine Snow_Analysis_OI

 Subroutine calculate_snow_cover_fraction(LENSFC, SWE, VETFCS_in, nodata_val, nodata_tol, fSCA)
! Routine to calculate snow cover fraction from SWE. Taken from Noah model code 
! Assuming -ve SWE = nodata (sloppy) 

        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in) :: LENSFC
        REAL, intent(In)   :: SWE(LENSFC), VETFCS_in(LENSFC)
        REAL, intent(in)   :: nodata_val, nodata_tol
        REAL, intent(Out)   ::  fSCA(LENSFC)
        INTEGER                :: VETFCS(LENSFC)

        REAL               :: snupx(30), SNEQV(LENSFC), SNUP, SALP, RSNOW
        Integer                    :: indx, vtype_int

        !This is for the IGBP veg classification scheme.
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                        0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                        0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                        0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                        0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

        SNEQV = 0.001 * SWE   ! units mm->m
        SALP = -4.0
        
        VETFCS = INT(VETFCS_in)
        Where(VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
        
        Do indx=1, LENSFC
                SNUP = snupx(VETFCS(indx))
                if (SNUP == 0.) then
                        print*, " 0.0 snup value, check vegclasses", vtype_int
                        Stop
                endif

                IF (SNEQV(indx) <  (nodata_val+nodata_tol) ) THEN 
                          fSCA(indx) = nodata_val
                ELSEIF (SNEQV(indx) < 0.00001) THEN 
                          fSCA(indx) = 0.0
                ELSEIF (SNEQV(indx) .LT. SNUP) THEN
                        RSNOW = SNEQV(indx)/SNUP
                        fSCA(indx) = 1. - (EXP(SALP*RSNOW) - RSNOW*EXP(SALP))
                ELSE
                        fSCA(indx) = 1.0
                ENDIF

        End do
        
        RETURN

 End Subroutine calculate_snow_cover_fraction

 ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_DA_Outputs(output_file, idim, jdim, lensfc, myrank,   &
                                 lons_grid, lats_grid, & 
                                 snoforc, snoanl, snwdforc, snwdanal, landmask,  &
                                 SNCOV_IMS, SND_IMS_at_Grid, &
                                 num_stn, Lat_atObs, Lon_atObs, index_stn,OBS_stn, FCS_at_stn, ANL_at_stn )  
!**
        !------------------------------------------------------------------
        ! Write DA ouputs: 
        ! forecast SWE
        ! analysis SWE
        ! analysis Snow Depth
        ! innovation at grid
        !------------------------------------------------------------------
        implicit none

        CHARACTER(LEN=*), Intent(In)      :: output_file
        integer, intent(in)         :: idim, jdim, lensfc, num_stn
        real, intent(in)            :: snoforc(lensfc), snoanl(lensfc) , lons_grid(lensfc), lats_grid(lensfc)
        real, intent(in)            ::  snwdforc(lensfc), snwdanal(lensfc)
        integer, intent(in)         :: landmask(lensfc)
        integer, intent(in)         :: index_stn(2,num_stn)
        Real, intent(in)            :: OBS_stn(num_stn), FCS_at_stn(num_stn), ANL_at_stn(num_stn)
        Real, intent(in)            :: Lat_atObs(num_stn), Lon_atObs(num_stn)
        Real, intent(in)            :: SND_IMS_at_Grid(lensfc), SNCOV_IMS(lensfc)  !, anl_fSCA(lensfc)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        integer                     :: error, i, ncid
        integer                     :: dim_x, dim_y, dim_time, dim_stn
        integer                     :: id_x, id_y, id_time
        integer       :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_imssno, id_imscov  
        integer       :: id_latstn, id_lonstn, id_obsstn, id_forcstn, id_analstn, id_landmask
        integer       :: id_iindex, id_jindex, id_lon, id_lat
        
        integer                     :: myrank

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
        ! obs points
        error = nf90_def_dim(ncid, 'obs_points', num_stn, dim_stn)
        call netcdf_err(error, 'DEFINING obs_points' )

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

        error = nf90_def_var(ncid, 'Lon', NF90_DOUBLE, dims_3d, id_lon)
        call netcdf_err(error, 'DEFINING Longitude' )
        error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
        call netcdf_err(error, 'DEFINING Lon LONG NAME' )
        error = nf90_put_att(ncid, id_lon, "units", "degrees East")
        call netcdf_err(error, 'DEFINING Lon UNITS' )

        error = nf90_def_var(ncid, 'Lat', NF90_DOUBLE, dims_3d, id_lat)
        call netcdf_err(error, 'DEFINING Lattitude' )
        error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
        call netcdf_err(error, 'DEFINING Lat LONG NAME' )
        error = nf90_put_att(ncid, id_lat, "units", "degrees North")
        call netcdf_err(error, 'DEFINING Lat UNITS' )

        error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_3d, id_swe_forc)
        call netcdf_err(error, 'DEFINING SWE Forecast' )
        error = nf90_put_att(ncid, id_swe_forc, "long_name", "Forecast Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_swe_forc, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE Forecast UNITS' )

        error = nf90_def_var(ncid, 'SWE_Analysis', NF90_DOUBLE, dims_3d, id_swe)
        call netcdf_err(error, 'DEFINING SWE_Analysis' )
        error = nf90_put_att(ncid, id_swe, "long_name", "Analysis Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE LONG NAME' )
        error = nf90_put_att(ncid, id_swe, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE UNITS' )

        error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_3d, id_snwdf)
        call netcdf_err(error, 'DEFINING SNF Forecast' )
        error = nf90_put_att(ncid, id_snwdf, "long_name", "Forecast Snow Depth")
        call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_snwdf, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

        error = nf90_def_var(ncid, 'SND_Analysis', NF90_DOUBLE, dims_3d, id_snwd)
        call netcdf_err(error, 'DEFINING SND Analyis' )
        error = nf90_put_att(ncid, id_snwd, "long_name", "Analysis Snow Depth")
        call netcdf_err(error, 'DEFINING SND Analysis LONG NAME' )
        error = nf90_put_att(ncid, id_snwd, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Analysis UNITS' )

        error = nf90_def_var(ncid, 'imsSND', NF90_DOUBLE, dims_3d, id_imssno)
        call netcdf_err(error, 'DEFINING imsSND' )
        error = nf90_put_att(ncid, id_imssno, "long_name", "IMS derived snow depth") 
        call netcdf_err(error, 'DEFINING imsSND LONG NAME' )
        error = nf90_put_att(ncid, id_imssno, "units", "mm")
        call netcdf_err(error, 'DEFINING imsSND UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        ! obs points
        error = nf90_def_var(ncid, 'latitude@MetaData', NF90_DOUBLE, dim_stn, id_latstn)
        call netcdf_err(error, 'DEFINING  latitude@MetaData' )
        error = nf90_put_att(ncid, id_latstn, "long_name", "Latitude at Observation Points")
        call netcdf_err(error, 'DEFINING  latitude@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_latstn, "units", "deg")
        call netcdf_err(error, 'DEFINING  latitude@MetaData UNITS' )

        error = nf90_def_var(ncid, 'longitude@MetaData', NF90_DOUBLE, dim_stn, id_lonstn)
        call netcdf_err(error, 'DEFINING longitude@MetaData' )
        error = nf90_put_att(ncid, id_lonstn, "long_name", "Longitude at Observation Points")
        call netcdf_err(error, 'DEFINING longitude@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_lonstn, "units", "deg")
        call netcdf_err(error, 'DEFINING longitude@MetaData UNITS' )
        
        error = nf90_def_var(ncid, 'IDIM_index@MetaData', NF90_INT, dim_stn, id_iindex)
        call netcdf_err(error, 'DEFINING IDIM_index@MetaData' )
        error = nf90_put_att(ncid, id_iindex, "long_name", "IDIM index at Observation Points")
        call netcdf_err(error, 'DEFINING IDIM_indexe@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_iindex, "units", "-")
        call netcdf_err(error, 'DEFINING IDIM_index@MetaData UNITS' )

        error = nf90_def_var(ncid, 'JDIM_index@MetaData', NF90_INT, dim_stn, id_jindex)
        call netcdf_err(error, 'DEFINING JDIM_index@MetaData' )
        error = nf90_put_att(ncid, id_jindex, "long_name", "JDIM index at Observation Points")
        call netcdf_err(error, 'DEFINING JDIM_indexe@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_jindex, "units", "-")
        call netcdf_err(error, 'DEFINING JDIM_index@MetaData UNITS' )

        error = nf90_def_var(ncid, 'snwdph@ObsValue', NF90_DOUBLE, dim_stn, id_obsstn)
        call netcdf_err(error, 'DEFINING snwdph@ObsValue' )
        error = nf90_put_att(ncid, id_obsstn, "long_name", "Observed at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue LONG NAME' )
        error = nf90_put_att(ncid, id_obsstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue UNITS' )
        
        error = nf90_def_var(ncid, 'snwdph@hofx', NF90_DOUBLE, dim_stn, id_forcstn) 
        call netcdf_err(error, 'DEFINING snwdph@hofx' )
        error = nf90_put_att(ncid, id_forcstn, "long_name", "Forecast at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@hofx LONG NAME' )
        error = nf90_put_att(ncid, id_forcstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@hofx UNITS' )

        error = nf90_def_var(ncid, 'snwdph@anal', NF90_DOUBLE, dim_stn, id_analstn)
        call netcdf_err(error, 'DEFINING snwdph@anal' )
        error = nf90_put_att(ncid, id_analstn, "long_name", "Analysis at Observation Points") 
        call netcdf_err(error, 'DEFINING snwdph@anal LONG NAME' )
        error = nf90_put_att(ncid, id_analstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@anal UNITS' )

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

        dum2d = reshape(lons_grid, (/idim,jdim/))
        error = nf90_put_var( ncid, id_lon, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING Longitude RECORD' )

        dum2d = reshape(lats_grid, (/idim,jdim/))
        error = nf90_put_var( ncid, id_lat, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING Lattitude RECORD' )

        dum2d = reshape(snoforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe_forc, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Forecast RECORD' )

        dum2d = reshape(snoanl, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Analysis RECORD' ) 

        dum2d = reshape(snwdforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwdf, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND Forecast RECORD' )

        dum2d = reshape(snwdanal, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwd, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND analysis RECORD' )

        dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )

         dum2d = reshape(SND_IMS_at_Grid, (/idim,jdim/))
         error = nf90_put_var( ncid, id_imssno, dum2d, dims_strt, dims_end)
         call netcdf_err(error, 'WRITING imsSND RECORD' )
        
        ! obs points (obs, hofx, omb) 
        error = nf90_put_var( ncid, id_latstn, Lat_atObs)
        call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
        call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_iindex, index_stn(1,:))
        call netcdf_err(error, 'WRITING idim index RECORD' )

        error = nf90_put_var( ncid, id_jindex, index_stn(2,:))
        call netcdf_err(error, 'WRITING jdim index RECORD' )

        error = nf90_put_var( ncid, id_obsstn, OBS_stn)
        call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_forcstn, FCS_at_stn)
        call netcdf_err(error, 'WRITING SNDFCS_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_analstn, ANL_at_stn)
        call netcdf_err(error, 'WRITING anal_atObsPts RECORD' )

        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 End subroutine Write_DA_Outputs

 SUBROUTINE READ_Forecast_Data_atPath(forc_inp_path, veg_type_landice, LENSFC, SWEFCS, SNDFCS, VETFCS, LANDMASK)  
!**
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: forc_inp_path
        INTEGER, INTENT(IN)               :: LENSFC, veg_type_landice
        REAL, INTENT(OUT)                 :: SWEFCS(LENSFC), SNDFCS(LENSFC), VETFCS(LENSFC) 
        INTEGER, INTENT(OUT)              :: LANDMASK(LENSFC) 

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID, i
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists

        INQUIRE(FILE=trim(forc_inp_path), EXIST=file_exists)

        if (.not. file_exists) then 
                print *, 'READ_Forecast_Data_atPath error,file does not exist', &   
                        trim(forc_inp_path) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif

        ERROR=NF90_OPEN(TRIM(forc_inp_path), NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(forc_inp_path) )

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

 SUBROUTINE READ_LAT_LON_OROG_atRank(MYRANK, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,IJDIM) 
!** 
    
    !--------------------------------------------------------------
    ! READ LATITUDE, LONGITUDE, FILTERED OROGRAPHY,
    ! ON THE CUBED-SPHERE TILE FROM
    ! THE "GRID" FILE.
    !--------------------------------------------------------------
    
     IMPLICIT NONE
    
     include "mpif.h"
    
     INTEGER, INTENT(IN)    :: IDIM, JDIM, IJDIM
    
     CHARACTER(LEN=5), INTENT(OUT) :: TILE_NUM
    
     REAL, INTENT(OUT)      :: RLA(IJDIM),RLO(IJDIM)
     REAL, INTENT(OUT)      :: OROG(IJDIM)  
    
     CHARACTER(LEN=50)      :: FNOROG, FNGRID
     CHARACTER(LEN=3)       :: RANKCH
    
     INTEGER                :: ERROR, NCID, NCID_OROG
     INTEGER                :: I, II, J, JJ, MYRANK
     INTEGER                :: ID_DIM, ID_VAR, NX, NY
    
     REAL, ALLOCATABLE         :: DUMMY(:,:), GEOLAT(:,:), GEOLON(:,:)
     REAL(KIND=4), ALLOCATABLE :: DUMMY4(:,:)
     LOGICAL                   :: file_exists
    
     WRITE(RANKCH, '(I3.3)') (MYRANK+1)
    
     FNGRID = "./fngrid." // RANKCH
    
     if (print_deb) then
        PRINT*, "READ FV3 GRID INFO FROM: "//TRIM(FNGRID)
     endif

    INQUIRE(FILE=trim(FNGRID), EXIST=file_exists)

    if (.not. file_exists) then 
            print *, 'READ_LAT_LON_OROG_atRank error, grid file does not exist', &   
                    trim(FNGRID) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10)
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

     INQUIRE(FILE=trim(FNOROG), EXIST=file_exists)

     if (.not. file_exists) then 
            print *, 'READ_LAT_LON_OROG_atRank error orogfile does not exist', &   
                    trim(FNGRID) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10)
     endif
    
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

 END MODULE M_Snow_Analysis
 
