MODULE M_Snow_Analysis
    
USE NETCDF
USE M_DA
USE M_UFSLAND_SNOW_UPDATE, ufsu_NETCDF_ERR=>NETCDF_ERR
USE M_SCF_Den

!USE MPI
Use, Intrinsic :: IEEE_ARITHMETIC	

CONTAINS

 subroutine OI_Snow_Analysis_NOAH(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  ! num_assim_steps, dT_Asssim,  & 
                                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
                                L_horz , h_ver, obs_tolerance, &
                                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                                ims_max_ele, num_subgrd_ims_cels, &
                            STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, & 
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, SFC_FORECAST_PREFIX, output_prefix, &   ! CURRENT_ANALYSIS_PREFIX,  &
                                SNDANL, exclude_obs_at_grid) !incr_at_Grid_out, SNOFCS_out, , SWEANL_out, LANDMASK)
                                                        
        !----------------------------------------------------------------------
        ! Input arguments: 
        ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
        ! IY, IM, ID, IH = year, month, day, hour of current model step   
        ! MYRANK: rank/id of the MPI process
        ! ...
        !
        ! Inputs, read from file:
        ! RLA, RLO: lat lon information for the tile
        ! SNOFCS(LENSFC): forecast  snowdepth, or snow water equivalent forecat, as requested
        !
        ! Outputs:
        ! SNOANL: snow analysis (of SWE or SND, see below)
        ! 
        ! Draper - changes to snow variables to clarify names, removed unnecesary conversions
        !          SWE - snow water equivalent 
        !          SND - snow depth 
        !          SNO - generic variable to refer to either (determined by value of 
        !                snow_OI_type
        !                   
        !----------------------------------------------------------------------
        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &  !SNOW_OI_TYPE, 
                                  IY, IM, ID, IH, LENSFC, IVEGSRC !, num_assim_steps 
        !LOGICAL, intent(in)    :: assim_SnowPack_obs, assim_SnowCov_obs, ims_correlated
        !CHARACTER(len=4), intent(in)   :: stn_var ! should probably be called control_var
        CHARACTER(LEN=*), Intent(In)   :: STN_OBS_PREFIX, &
                                          STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, & 
                                          IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
        CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX, output_prefix   !, CURRENT_ANALYSIS_PREFIX
 
        REAL, intent(In)    :: PERCENT_OBS_WITHHELD  !, dT_Asssim 
        Real, intent(In)    :: L_horz , h_ver, obs_tolerance
        Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
        INTEGER, intent(in) :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels 

        REAL, intent(Out)   :: SNDANL(LENSFC)  !, incr_at_Grid_out(LENSFC), 
                               !, SNOFCS_out(LENSFC), SWEANL_out(LENSFC)
        LOGICAL             :: exclude_obs_at_grid

        INTEGER             :: LANDMASK(LENSFC)
        CHARACTER(LEN=5)    :: TILE_NUM
        Character(LEN=3)    :: rank_str
        INTEGER             :: IERR     
        REAL        :: RLA(LENSFC), RLO(LENSFC), RLO_Tile(LENSFC), OROG(LENSFC)  !, OROG_UF(LENSFC)
        REAL        :: SNDFCS(LENSFC), SWEFCS(LENSFC), SWEANL(LENSFC)
        ! REAL        :: SNDANL_Cur(LENSFC), SWEANL_Cur(LENSFC)       !, SNOANL_Cur(LENSFC)
        REAL        :: VETFCS(LENSFC), SNUP_Array(LENSFC)           ! SNOFCS(LENSFC), 
        
        CHARACTER(len=250)   :: ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
        REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:)               
        REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROGFCS_at_stn(:)  
        REAL                 :: lat_min, lat_max, lon_min, lon_max      
        Real                 :: SNCOV_IMS(LENSFC)  ! ims resampled at each grid
        Real                 :: SNO_IMS_at_Grid(LENSFC)
        Real                 :: IMS_Foot_Print(LENSFC)  ! grid cells where IMS was assimilated

        INTEGER :: num_stn, Num_Ims, num_Eval !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon
        !Real    :: bkgst_srch_rad   ! radius_of_influence for selecting state at observation point
        INTEGER :: jndx, zndx, ncol, nrow
        Integer, Allocatable   :: index_back_at_nearStn(:), index_back_at_nearIMS(:) !loc_near_Obs(:), 
        Integer                :: num_loc, num_loc_1, num_loc_2
        Real, Parameter         :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 80. ! mm 
        real                    :: stdev_obsv, stdev_back
        ! Integer                 :: ims_assm_hour
        !Real                           :: obs_tolerance, ims_max_ele

        Real(dp), Allocatable    :: B_cov_mat(:,:), b_cov_vect(:)
        Real(dp), Allocatable    :: O_cov_mat(:,:), W_wght_vect(:)
        Real, Allocatable        :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:)
        REAL                     :: incr_at_Grid(LENSFC)    ! increment at grid
        Real, Allocatable        :: obs_Innov(:), OmB_innov_at_stn(:)

        CHARACTER(len=250)       :: forc_inp_file, da_out_file, anl_inp_path
        CHARACTER(LEN=3)         :: RANKCH 

        REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), incr_atEvalPts(:), SNOANL_atEvalPts(:)   !, SNOANL_Cur_atEvalPts(:)  !evalution points 
        REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
        Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
        Integer, ALLOCATABLE    :: index_back_atObs(:)   ! the location of background corresponding obs

        Real               :: snodens, SNODENS_Grid(LENSFC)
        !LOGICAL            :: assim_snpack_stn, assim_SWE    !use swe as model state var? (instead of snow depth)
        LOGICAL            :: assim_sncov_thisGridCell    !    assim_sncov, assimilate sncov, 

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

        !initialse output with nodata
        ! SNOANL= nodata_val
        IMS_Foot_Print = IEEE_VALUE(IMS_Foot_Print, IEEE_QUIET_NAN)

        stdev_obsv = stdev_obsv_depth
        stdev_back = stdev_back_depth  
        !obs_srch_rad = 250. ! radius of observation search
        ! ims_assm_hour = 18 
        ! noah models specific? Needed to ID glaciers.
        if (IVEGSRC == 2) then   ! sib
                veg_type_landice=13
        else
                veg_type_landice=15
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
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

! must have all obs
        if ((STN_OBS_PREFIX(1:8).eq.'        ') .OR. (IMS_SNOWCOVER_PATH(1:8).eq.'        ') &
             .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
                print*, "One or more observation paths don't exist!, skipping the OI DA"
                goto 999
        end if

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
        CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG,TILE_NUM,IDIM,JDIM,LENSFC) !OROG_UF,
        IF (p_tRank==PRINTRANK) PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
        
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
    ! Do istep = 1, num_assim_steps        
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc
        write(fvs_tile, "(I3)") IDIM

! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 
        num_Eval = 0

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

! 10.22.20 Cur_Analysis: existing SNODEP based analysis
    !     anl_inp_path = TRIM(CURRENT_ANALYSIS_PREFIX)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." //TRIM(h_str)// &
    !                        "0000.sfcanl_data."//TILE_NUM//".nc"
    !     Call READ_Analysis_Data(anl_inp_path, LENSFC, SWEANL_Cur, SNDANL_Cur)   ! p_tN,

        ! initialise analysis to forecast (mostly to retain non-land snow states, which are not updated)
        SWEANL = SWEFCS
        SNDANL = SNDFCS
        incr_at_Grid = 0 

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
        If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp
        
!=============================================================================================
! 3. Read observations
!=============================================================================================

! 3a. Read station obs (of snow depth or SWE)

        !if (assim_snpack_stn) then 
             ghcnd_inp_file = TRIM(STN_OBS_PREFIX)//"GHCND.SNWD."// &
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"
             if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
            !  dim_name = "Site_Id"    
             ! here: all valid stn obs, within lat/lon box (can be outside of tile though) 
             call Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, TRIM(STN_DIM_NAME), &
                        lat_min, lat_max, lon_min, lon_max, & 
                        num_stn, SNOOBS_stn,              &
                        Lat_stn, Lon_stn, MYRANK) 
            ! ghcnd_inp_file = TRIM(STN_OBS_PREFIX)// &   !GHCND.SNWD TRIM(h_str)//
            !                 TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//".nc"
            ! if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
            ! Call Observation_Read_GHCND_IODA(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
            !                 TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
            !                 num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn)

            if (p_tRank==0) then   !.and. (p_tN==PRINTRANK).and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if ((p_tRank==0) .and. (p_tN==PRINTRANK) .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNOOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
            endif
            if (myrank==0) PRINT*,'Finished reading station data'

        !endif

! CSD beyond this point, there should be no specific mention of the station data source

! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 

        ! if (assim_sncov) then
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
        
        ! endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!      snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        ! if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNOFCS_at_stn(num_stn))
            ALLOCATE(OROGFCS_at_stn(num_stn)) 
            ALLOCATE(index_back_atObs(num_stn)) 
            ALLOCATE(OmB_innov_at_stn(num_stn)) 
             ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
            num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
            if (num_Eval > 0) then 
                ALLOCATE(index_back_atEval(num_Eval)) 
                ALLOCATE(Obs_atEvalPts(num_Eval)) 
                ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
                ALLOCATE(Lat_atEvalPts(num_Eval))
                ALLOCATE(Lon_atEvalPts(num_Eval)) 
                ALLOCATE(incr_atEvalPts(num_Eval))
                ALLOCATE(SNOANL_atEvalPts(num_Eval))                  
                if(p_tRank == 0) then 
                        PRINT*, "Tile ", p_tN+1, " ", num_Eval, ' points for evaluation excluded from DA'       
                endif  
            endif 

    ! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)
    !           for JEDI, probably also want to add snow depth derived from snow cover here.
            Call Observation_Operator_tiles_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, & 
                                RLA, RLO, OROG, Lat_stn, Lon_stn,   &
        LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNDFCS, SNOOBS_stn, LANDMASK,  &
        SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval,  &
                Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            if ((p_tN==PRINTRANK) .and. (p_tRank==0) .and. print_deb) then 
                    PRINT*, "Background Indices at eval points"
                    PRINT*, index_back_atEval       
                    PRINT*, "Obs at Eval Points" 
                    PRINT*, Obs_atEvalPts   
                    PRINT*, "Forecast at Eval Points"
                    PRINT*, SNOFCS_atEvalPts                             
                    PRINT*, "Lat at Eval Points"
                    PRINT*, Lat_atEvalPts
                    PRINT*, "Lon at Eval Points"
                    PRINT*, Lon_atEvalPts
            endif

            OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn

            if ((p_tN==PRINTRANK) .and. (p_tRank==0) .and. print_deb) then
                    PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                    PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
                !     PRINT*, "Lat at Obs Points"
                !     PRINT*, Lat_stn
                !     PRINT*, "Lon at Obs Points"
                !     PRINT*, Lon_stn
                    PRINT*, "Model elevation at station locations from rank: ", MYRANK
                    PRINT*, OROGFCS_at_stn
                    PRINT*, "Background Indices at obs points"
                    PRINT*, index_back_atObs
                    PRINT*, "Background at station locations from rank: ", MYRANK
                    PRINT*, SNOFCS_at_stn  
                !     PRINT*, "Obs at obs stns" 
                !     PRINT*, SNOOBS_stn   
                    PRINT*, "O - B (innovation at obs points)"
                    PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished observation operator for station data'         
        ! endif ! num_stn > 0

!=============================================================================================
! 5.  obs QC goes here
!=============================================================================================
! As of 2.4.22: Most qcs are being done at obs read time and inside the DA loop 
! for efficiency; 
! In the future consider consolidating QC's here

!CSDCSD - todo. Add QC here.

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
        if (myrank==PRINTRANK) PRINT*,'Starting DA loop'
        Do jndx = mp_start, mp_end     !LENSFC/2, LENSFC ! 442369, 442370   !
           ! QC: only update this grid cell if it is land.
                if( LANDMASK(jndx) == 1 ) then 
                    num_loc_1 = 0
                    num_loc_2 = 0
                    assim_sncov_thisGridCell = .FALSE.
                    if (print_deb) print*, "proc ", myrank, " grid: ", jndx
                    if (num_stn>0) then 
                    ! CSD - speed up by skipping over all grid cells where model and IMS agree on no snow? 
                    ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                    ! QC: gross error check is done in this call.
                            call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
                                            Lat_stn, Lon_stn,  num_stn, obs_srch_rad, max_num_nearStn,   &
                                            stdev_back, stdev_obsv, obs_tolerance,                 &
                                            SNOFCS_at_stn, SNOOBS_stn,                                              &
                                            index_back_at_nearStn,  num_loc_1) !,      &LENSFC,
                            if (print_deb) print*, "number of stn sndpth obs ", num_loc_1
                    endif         
                    !     
                    if( (.NOT. IEEE_IS_NAN(SNDFCS(jndx))) .AND. &     !(IH_loc == ims_assm_hour) .AND. &
                        (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(jndx))) .AND. &
                        ! (OROG(jndx) <= ims_max_ele) .AND. &
                        ( .not.(SWEFCS(jndx) >= SNUP_Array(jndx) .AND. & 
                                SNCOV_IMS(jndx) >= 0.99) ) ) then
                            num_loc_2 = 1 
                            assim_sncov_thisGridCell = .TRUE.                
                    endif
                    num_loc = num_loc_1 + num_loc_2                
                    ! if assim_sncov=false >> num_loc_1=num_loc
                    ! QC: only update this grid cell if it is land.
                    if(num_loc > 0) then            !((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then 
                        ! get background states
                        Allocate(back_at_Obs(num_loc))
                        Allocate(obs_Array(num_loc))
                        Allocate(Lat_Obs(num_loc))
                        Allocate(Lon_Obs(num_loc))
                        Allocate(orogfcs_Obs(num_loc))                        
                        if(num_loc_1 > 0) then
                                Do zndx = 1, num_loc_1     
                                    back_at_Obs(zndx) = SNOFCS_at_stn(index_back_at_nearStn(zndx))
                                    obs_Array(zndx) = SNOOBS_stn(index_back_at_nearStn(zndx))
                                    Lat_Obs(zndx) = Lat_stn(index_back_at_nearStn(zndx))
                                    Lon_Obs(zndx) = Lon_stn(index_back_at_nearStn(zndx))
                                    orogfcs_Obs(zndx) = OROGFCS_at_stn(index_back_at_nearStn(zndx)) 
                                End Do
                        End if
                        ! Append IMS-derived snow depth to the obs array 
                        if(assim_sncov_thisGridCell) then   
                                IMS_Foot_Print(jndx) = 1.0                          
                                back_at_Obs(num_loc) = SNDFCS(jndx)
                                obs_Array(num_loc) = SNO_IMS_at_Grid(jndx)
                                Lat_Obs(num_loc) = RLA(jndx)   
                                Lon_Obs(num_loc) = RLO(jndx) 
                                orogfcs_Obs(num_loc) = OROG(jndx)
                        endif
                        
                        if(exclude_obs_at_grid .and. (num_loc_1 > 1)) then 
                            Allocate(B_cov_mat(num_loc-1, num_loc-1))
                            Allocate(b_cov_vect(num_loc-1))
                            Allocate(O_cov_mat(num_loc-1, num_loc-1))
                            Allocate(W_wght_vect(num_loc-1))  
                            Allocate(obs_Innov(num_loc-1)) 
                            ! compute covariances
                            call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx),  &  !SNDFCS(jndx),    &
                                Lat_Obs(2:num_loc), Lon_Obs(2:num_loc), orogfcs_Obs(2:num_loc),  &
                                num_loc-1, Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                                L_horz, h_ver, assim_sncov_thisGridCell,             &
                                B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                        
                            call Snow_DA_OI(back_at_Obs(2:num_loc), obs_Array(2:num_loc), &
                                num_loc-1, W_wght_vect,            &
                                SNDFCS(jndx), incr_at_Grid(jndx), SNDANL(jndx), obs_Innov)
                        else 
                            Allocate(B_cov_mat(num_loc, num_loc))
                            Allocate(b_cov_vect(num_loc))
                            Allocate(O_cov_mat(num_loc, num_loc))
                            Allocate(W_wght_vect(num_loc))  
                            Allocate(obs_Innov(num_loc)) 
                            call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), &  !SNDFCS(jndx),    &
                                    Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc,   &
                                    Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                                    L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                                    assim_sncov_thisGridCell,                          &
                                    B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                        
                            call Snow_DA_OI(back_at_Obs, obs_Array, num_loc, W_wght_vect,            &
                            SNDFCS(jndx), incr_at_Grid(jndx), SNDANL(jndx), obs_Innov)
                        endif

                        if ((p_tN==PRINTRANK) .and. (p_tRank==0) .and. print_deb) then  !
                                print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1, "total obs", num_loc
                                PRINT*, " background at obs pts: "
                                PRINT*, back_at_Obs     
                                PRINT*, "Observed"
                                PRINT*,  obs_Array
                                PRINT*, "Obs innovation: "
                                PRINT*, obs_Innov
                                PRINT*, "Weight vector: "
                                PRINT*, W_wght_vect  
                        endif                        
                        if ((p_tN==PRINTRANK) .and. (p_tRank==0) .and. print_deb) then  !
                                print*, "forec: ", SNDFCS(jndx)   
                                print*, "increment SND and SNCOV: ", incr_at_Grid(jndx), " anl: ", SNDANL(jndx)
                        endif           
                        !free mem
                        DEALLOCATE(back_at_Obs, obs_Array)
                        DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs, obs_Innov)
                        DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                            
                    ! else 
                        ! no obs available: keep background (copied earlier)
                        ! SNOANL(jndx) = SNOFCS(jndx) 
                    endif
                    if (allocated(index_back_at_nearStn))  Deallocate(index_back_at_nearStn) 
                    if (allocated(index_back_at_nearIMS))  Deallocate(index_back_at_nearIMS)                 
                endif ! not a land cell
        End do
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==PRINTRANK) PRINT*, ' Finished DA loops'
        
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
            call MPI_SEND(incr_at_Grid(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                                        MYRANK*100, MPI_COMM_WORLD, IERR)
        else    !if(p_tRank == 0) then  
            Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                send_proc = MYRANK +  pindex * NUM_TILES
                call MPI_RECV(SNDANL(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
                                    send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(incr_at_Grid(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,      &
                                    send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            enddo
        endif
        if (myrank==PRINTRANK) PRINT*,'Finished Data copy'

        if (MYRANK > NUM_TILES - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

! fill in SWE and SND arrays       
        ! avoid -ve anl
        Where(SNDANL < 0.) SNDANL = 0.
        if (print_deb) then
                PRINT*, "Weighted increment SWE/snwd from rank: ", MYRANK
            PRINT*, incr_at_Grid       
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNDANL
        endif
        ! d = swe/sndi
        ! this is for writing outputs at observation and evaluation points
        WHERE (LANDMASK==1) SWEANL = SNDANL * SNODENS_Grid

!TODO Compute updated snocov 
        ! Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

        ! copy values at eval points
        incr_atEvalPts = IEEE_VALUE(incr_atEvalPts, IEEE_QUIET_NAN)
        SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
        Do jndx = 1, num_Eval
                if (index_back_atEval(jndx) > 0)  then   ! 
                    incr_atEvalPts(jndx) = incr_at_Grid(index_back_atEval(jndx))
                    SNOANL_atEvalPts(jndx) = SNDANL(index_back_atEval(jndx))      
                endif
        End do
        ! write outputs 
        Write(rank_str, '(I3.3)') (MYRANK+1)
        ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis/"
        da_out_file = TRIM(output_prefix)// &
                                  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"  !  
        call Write_DA_Outputs(da_out_file, IDIM, JDIM, LENSFC, MYRANK, &
                              SWEFCS, SWEANL, SNDFCS, SNDANL, LANDMASK, &  !
                              num_stn, Lat_stn, Lon_stn, SNOOBS_stn, SNOFCS_at_stn, OmB_innov_at_stn,  &                                
                              incr_at_Grid, SNCOV_IMS, IMS_Foot_Print, &
                              num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
                              SNOFCS_atEvalPts, incr_atEvalPts, SNOANL_atEvalPts)  !, anl_fSCA) !updated snocov

998 CONTINUE
        ! clean up
             if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
             if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)
             if (allocated(OmB_innov_at_stn))   DEALLOCATE(OmB_innov_at_stn)
             if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
             if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
             if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)
             if (allocated(OROGFCS_at_stn))  DEALLOCATE(OROGFCS_at_stn) 
             if (allocated(index_back_atEval)) DEALLOCATE(index_back_atEval)
             if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
             if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts)
             if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
             if (allocated(SNOFCS_atEvalPts)) DEALLOCATE(SNOFCS_atEvalPts)
             if (allocated(incr_atEvalPts))  DEALLOCATE(incr_atEvalPts)
             if (allocated(SNOANL_atEvalPts)) DEALLOCATE(SNOANL_atEvalPts)    

        ! Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        ! IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==PRINTRANK) PRINT*,'Finished OI DA at datetime: ', y_str, m_str, d_str, h_str
        
    ! End do

999 CONTINUE
        !PRINT*,'Finished OI DA ON RANK: ', MYRANK
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        !STOP

        RETURN

 END subroutine OI_Snow_Analysis_NOAH

  subroutine OI_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  !num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
                L_horz , h_ver, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, &   !SFC_FORECAST_PREFIX, CURRENT_ANALYSIS_PREFIX,  &
            vector_restart_prefix, vector_noda_prefix, static_prefix, output_prefix, &
            snowUpdateOpt, PRINTRANK, print_debg_info, &            !fv3_index, vector_inputs, &
            SNDANL_out, &   !SNDFCS_out, SWEANL_out, & incr_at_Grid_out, 
            Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)
                                                        
        !----------------------------------------------------------------------
        ! Input arguments: 
        ! IDIM * JDIM = LENSFC: number of grid cells in tile = xdim * ydim   
        ! IY, IM, ID, IH = year, month, day, hour of current model step   
        ! MYRANK: rank/id of the MPI process
        ! ...
        !
        ! Inputs, read from file:
        ! RLA, RLO: lat lon information for the tile
        ! SNOFCS(LENSFC): forecast  snowdepth, or snow water equivalent forecat, as requested
        !
        ! Outputs:
        ! SNOANL: snow analysis (of SWE or SND, see below)
        ! 
        ! Draper - changes to snow variables to clarify names, removed unnecesary conversions
        !          SWE - snow water equivalent 
        !          SND - snow depth 
        !          SNO - generic variable to refer to either (determined by value of 
        !                snow_OI_type
        !                   
        !----------------------------------------------------------------------
        IMPLICIT NONE
        !
        include 'mpif.h'
        
        integer, parameter :: dp = kind(1.d0)

        INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &  !SNOW_OI_TYPE, 
                                  IY, IM, ID, IH, LENSFC, IVEGSRC   !, num_assim_steps 
        LOGICAL, intent(in)    :: assim_SnowPack_obs, assim_SnowCov_obs    !, ims_correlated
        !CHARACTER(len=4), intent(in)   :: stn_var ! should probably be called control_var
        CHARACTER(LEN=*), Intent(In)   :: STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME
        CHARACTER(LEN=*), Intent(In)   :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
        ! CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX, CURRENT_ANALYSIS_PREFIX
        CHARACTER(LEN=*), Intent(In)   :: vector_restart_prefix, vector_noda_prefix, &
                                          static_prefix, output_prefix

        REAL, intent(In)    :: PERCENT_OBS_WITHHELD   !, dT_Asssim 
        Real, intent(In)    :: L_horz , h_ver, obs_tolerance
        Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
        INTEGER, intent(in) :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels 
        INTEGER, intent(in) :: snowUpdateOpt, PRINTRANK
        REAL, intent(out)   :: SNDANL_out(LENSFC)
        LOGICAL             :: print_debg_info   !, fv3_indexvector_inputs, 

         ! for mpi par
        INTEGER   :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
        Integer   :: LENSFC_proc, begloc, endloc
        LOGICAL             :: exclude_obs_at_grid
        
        CHARACTER(LEN=5)    :: TILE_NUM
        ! Character(LEN=3)    :: rank_str
        INTEGER             :: IERR     
        REAL                :: SNDFCS(LENSFC_proc), SWEFCS(LENSFC_proc), &
                               SNDANL(LENSFC_proc), SWEANL(LENSFC_proc)
        REAL                :: SNDFCS_out(LENSFC), SWEFCS_out(LENSFC), &
                               incr_at_Grid_out(LENSFC), SWEANL_out(LENSFC), &
                               SCF_Grid(LENSFC_proc), SCF_Grid_out(LENSFC), &
                               SNO_IMS_at_Grid_out(LENSFC)
        REAL                :: SNDnoDA(LENSFC), SWEnoDA(LENSFC)       !, SNOANL_Cur(LENSFC)
        REAL, allocatable   :: SNDnoDA_atEval(:), SWEnoDA_atEval(:) !, SND_ml_noDA_atEval(:)   
        
        REAL   :: RLA(LENSFC_proc), RLO(LENSFC_proc), RLO_Tile(LENSFC_proc), OROG(LENSFC_proc)  !, OROG_UF(LENSFC)
        REAL            :: VETFCS(LENSFC_proc), SNUP_Array(LENSFC_proc)           ! SNOFCS(LENSFC), 
        Integer         :: tile_xy(LENSFC_proc), Idim_xy(LENSFC_proc), Jdim_xy(LENSFC_proc)
        INTEGER         :: LANDMASK(LENSFC_proc)  !, DAMASK(LENSFC)

        CHARACTER(len=250)   :: ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
        CHARACTER(len=4)     :: y_str, m_str, d_Str, h_str, fvs_tile
        REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:), SNOANL_atStn(:)  !, SNOFCS_at_stn_ML(:)                
        REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROG_at_stn(:)  
        REAL                 :: lat_min, lat_max, lon_min, lon_max      
        Real                 :: SNCOV_IMS(LENSFC_proc)  ! ims resampled at each grid
        Real                 :: SNO_IMS_at_Grid(LENSFC_proc)
        Real                 :: IMS_Foot_Print(LENSFC_proc)  ! grid cells where IMS was assimilated

        INTEGER :: num_stn, Num_Ims, num_Eval !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon
        !Real    :: bkgst_srch_rad   ! radius_of_influence for selecting state at observation point
        INTEGER :: jndx, zndx, ncol, nrow
        Integer, Allocatable   :: index_at_nearStn(:)  !, index_at_nearIMS(:) !loc_near_Obs(:), 
        Integer                :: num_loc, num_loc_1, num_loc_2
        Real, Parameter   :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 80. ! mm 
        real                    :: stdev_obsv, stdev_back
        ! Integer                 :: ims_assm_hour

        !Real                           :: obs_tolerance, ims_max_ele

        Real(dp), Allocatable    :: B_cov_mat(:,:), b_cov_vect(:)
        Real(dp), Allocatable    :: O_cov_mat(:,:), W_wght_vect(:)
        Real, Allocatable        :: back_at_Obs(:), obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:)
        REAL                     :: incr_at_Grid(LENSFC_proc)    ! increment at grid
        Real, Allocatable        :: obs_Innov(:), OmB_innov_at_stn(:)

        !ML
        real                     :: incr_at_Grid_ML(LENSFC_proc), SNDANL_ML(LENSFC_proc)

        CHARACTER(len=250)       :: forc_inp_file, da_out_file, noda_inp_path
        CHARACTER(LEN=4)         :: RANKCH 
        CHARACTER(LEN=500)       :: static_filename

        REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), incr_atEvalPts(:), SNOANL_atEvalPts(:) !, SNOANL_Cur_atEvalPts(:)  !evalution points 
        REAL, ALLOCATABLE  :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
        REAL, ALLOCATABLE  :: Orog_obs_atEvalPts(:), Orog_fcs_atEvalPts(:)
        Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
        Integer                 :: NEXC 
        Integer, ALLOCATABLE    :: index_back_atObs(:), Index_Obs_Excluded(:)   ! the location of background corresponding obs
!1.4 track obserations assimilated
        Integer            :: index_obs_assmilated(max_num_nearStn+1, LENSFC_proc)

        Real               :: snodens, SNODENS_Grid(LENSFC_proc)
        !LOGICAL            :: assim_snpack_stn, assim_SWE    !use swe as model state var? (instead of snow depth)
        LOGICAL            :: assim_sncov_thisGridCell    !    assim_sncov, assimilate sncov, 

        Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice
        INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, &
                              dest_Aoffset_end, arLen, pindex
        INTEGER            :: mpiReal_size, rsize, isize, mpiInt_size, ixy
        REAL               :: tmp
        INTEGER            :: istep, IY_loc, IM_loc, ID_loc, IH_loc
        REAL               :: IH_real
        ! INTEGER, PARAMETER :: PRINTRANK = 4 
! CSD-todo, should be same as in sfcsub. Share properly
        real, parameter        :: nodata_val = -9999.9

!   Snow partion options for Noah-MP
        integer, parameter     :: DO_NOTHING = 0
        integer, parameter     :: TOP_LAYER = 1
        integer, parameter     :: BOTTOM_LAYER = 2
        integer, parameter     :: ALL_LAYERS = 3

        integer, parameter     :: lsm_type = 2

        ! noah mp
        type(noahmp_type)      :: noahmp
        type(observation_type) :: obs
        ! type noahmp_type
        allocate(noahmp%swe                (LENSFC_proc))
        allocate(noahmp%snow_depth         (LENSFC_proc))
        allocate(noahmp%active_snow_layers (LENSFC_proc))
        allocate(noahmp%swe_previous       (LENSFC_proc))
        allocate(noahmp%snow_soil_interface(LENSFC_proc,7))
        allocate(noahmp%temperature_snow   (LENSFC_proc,3))
        allocate(noahmp%snow_ice_layer     (LENSFC_proc,3))
        allocate(noahmp%snow_liq_layer     (LENSFC_proc,3))
        allocate(noahmp%temperature_soil   (LENSFC_proc))  
        ! end type noahmp_type 

!=============================================================================================
! 1. initialise vars,set-up processors, and read lat/lon from orog files.
!=============================================================================================

        !initialse output with nodata
        ! SNOANL= nodata_val
        IMS_Foot_Print = IEEE_VALUE(IMS_Foot_Print, IEEE_QUIET_NAN)
        ! DAMASK = 0

        stdev_obsv = stdev_obsv_depth
        stdev_back = stdev_back_depth  
        !obs_srch_rad = 250. ! radius of observation search
        ! ims_assm_hour = 18 
        ! noah models specific? Needed to ID glaciers.
        if (IVEGSRC == 2) then   ! sib
                veg_type_landice=13
        else
                veg_type_landice=15
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
        IF (myrank ==PRINTRANK) PRINT*,"snowDA: OI total num proc ", NPROCS, " Num tiles : ", NUM_TILES
  
! if (p_tN /= 4 ) goto 999

! must have all obs
        if ((STN_OBS_PREFIX(1:8).eq.'        ') .OR. (IMS_SNOWCOVER_PATH(1:8).eq.'        ') &
             .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
                print*, "One or more observation paths don't exist!, skipping the OI DA"
                goto 999
        end if

        write(fvs_tile, "(I0.2)") IDIM

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
        static_filename = trim(static_prefix)//"/ufs-land_C"//trim(fvs_tile)//"_static_fields.nc"

        Call ReadTileInfo(static_filename, LENSFC_proc, veg_type_landice, mp_start, mp_end, &
                tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 
    
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
                print*, "group ", p_tN, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
        endif

! If multiple time steps are simulated in a given time period (window) given by num_assim_steps * dT_Assim
! Note: the DA outputs from the last time step are returned        
        IH_real = IH; IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH; ! these variables change inside loop below
    ! Do istep = 1, num_assim_steps        
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc

! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 
        num_Eval = 0

!=============================================================================================
! 2. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
!=============================================================================================

       ! READ THE INPUT SURFACE DATA from vector
        ! "ufs_land_restart.", yyyy, "-", mm, "-", dd, "_", hh, "-", nn, "-", ss, ".nc"
        forc_inp_file=TRIM(vector_restart_prefix)//"/ufs_land_restart."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
        ! Read no da outputs
        noda_inp_path=TRIM(vector_noda_prefix)//"/ufs_land_output."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
        Call ReadRestartNoahMP(myrank, forc_inp_file, noda_inp_path, &
                               mp_start, mp_end, LENSFC_proc, &
                               LENSFC, noahmp, SNDnoDA, SWEnoDA, SCF_Grid)  !, SNDFCS, SWEFCS)

        ! initialise analysis to forecast (mostly to retain non-land snow states, which are not updated)
        SWEFCS = noahmp%swe
        SNDFCS = noahmp%snow_depth
        SWEANL = SWEFCS
        SNDANL = SNDFCS
        incr_at_Grid = 0 

        ! grid snow density
        ! SNODENS_Grid = SWEFCS/SNDFCS
        call calc_density(LENSFC_proc, lsm_type, LANDMASK, SWEFCS, SNDFCS, &
                        noahmp%temperature_soil, SNODENS_Grid)
        
        WHERE (SNODENS_Grid < 0.0001) SNODENS_Grid = 0.

        if (COUNT (LANDMASK==1 .and. SNDFCS> 0.001) > 0) then   ! 
            ! mean density over snow-covered land
            snodens = SUM(SNODENS_Grid, Mask = (LANDMASK==1 .and. SNDFCS>0.001 )) &  !
                        / COUNT (LANDMASK==1 .and. SNDFCS> 0.001)                    !
            ! if (MYRANK==PRINTRANK) & 
            PRINT *, 'proc', myrank, ' snowDA: mean density ', snodens 
        else 
            snodens = 0.1  ! default value if have no snow in array
                ! if (MYRANK==PRINTRANK) & 
            PRINT *, 'proc ', myrank, ' snowDA: no snow in current array, using default density ', snodens 
        endif
        ! for grid cells with no valid density, fill in the average snodens
        Where(.not.(LANDMASK==1 .and. SNDFCS>0.01 )) SNODENS_Grid = snodens                     !

        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.001 )) &          !
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.001)                 !
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.001 )) &              !
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.001)                  !
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' mean SND', tmp
        
!=============================================================================================
! 3. Read observations
!=============================================================================================
! 3a. Read station obs (of snow depth or SWE)
        if (assim_SnowPack_obs) then 
            ghcnd_inp_file = TRIM(STN_OBS_PREFIX)// &   !GHCND.SNWD TRIM(h_str)//
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//".nc"
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
             
            Call Observation_Read_GHCND_IODA(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
                            TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
                            num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn, &
                            NEXC, Index_Obs_Excluded)
            ! Call Observation_Read_GHCND_All_excNaN(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
            !                 TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
            !                 num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn)

            if (myrank==PRINTRANK) then   !.and. (p_tN==PRINTRANK).and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNOOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
                    PRINT*, "Elevation at station locations from rank: ", MYRANK
                    PRINT*, OROG_at_stn
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished reading station data'
        endif

! CSD beyond this point, there should be no specific mention of the station data source
! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 

        if (assim_SnowCov_obs) then
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
            ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
            WRITE(RANKCH, '(I0.1)') myrank                            
            ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                            ".IMS.Indices."//trim(ADJUSTL(RANKCH))//".nc"                       
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading IMS file', trim(ims_inp_file) 
            ! Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
            !                         MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
            Call Read_IMS_and_Resample_toTarget(ims_inp_file, IMS_INDEXES_PATH, & !fv3_index, 
                    LENSFC_proc, num_subgrd_ims_cels, IDIM, JDIM, &
                    tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
            if(myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV from rank: ", MYRANK
                PRINT*, SNCOV_IMS
            endif 
            if (myrank==PRINTRANK) PRINT*,'Finished reading SNCOV, converting to snow depth' 
           ! SNUP array will be used later, to test whether SCF > 100%
            ! call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC_proc, SNO_IMS_at_Grid, SNUP_Array)
            ! SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_Grid ! convert SWE to SND

            call calcSD_noahmp(LENSFC_proc, SNCOV_IMS, VETFCS, &
                                   SNODENS_Grid, SNO_IMS_at_Grid)

            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV derived snodepth at each grid cell from rank: ", MYRANK
                PRINT*, SNO_IMS_at_Grid
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished converting SNCOV observations'
        
        endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!      snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNOFCS_at_stn(num_stn))
            ALLOCATE(SNOANL_atStn(num_stn))
            ALLOCATE(index_back_atObs(num_stn)) 
            ALLOCATE(OmB_innov_at_stn(num_stn)) 
             ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
            num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
            if (num_Eval > 0) then 
                ALLOCATE(index_back_atEval(num_Eval)) 
                ALLOCATE(Obs_atEvalPts(num_Eval)) 
                ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
                ALLOCATE(Lat_atEvalPts(num_Eval))
                ALLOCATE(Lon_atEvalPts(num_Eval)) 

                ALLOCATE(Orog_obs_atEvalPts(num_Eval)) 
                ALLOCATE(Orog_fcs_atEvalPts(num_Eval)) 

                ALLOCATE(incr_atEvalPts(num_Eval))
                ALLOCATE(SNOANL_atEvalPts(num_Eval)) 
                ALLOCATE(SNDnoDA_atEval(num_Eval)) 
                ALLOCATE(SWEnoDA_atEval(num_Eval))  
                ! ALLOCATE(SND_ml_noDA_atEval(num_Eval))   
                ! ALLOCATE(SNOANL_Cur_atEvalPts(num_Eval))                  
                ! if(p_tRank == PRINTRANK) then 
                        PRINT*,"Proc ", myrank," ",num_Eval,' points for evaluation excluded from DA'       
                ! endif   
            endif 

    ! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)
    ! for JEDI, probably also want to add snow depth derived from snow cover here.
            ! Call Observation_Operator_tiles_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, & 
            !                     RLA, RLO, OROG, Lat_stn, Lon_stn,   &
            !                     LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNDFCS, SNOOBS_stn, LANDMASK,  &
            !                     SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval,  &
            !                     Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            Call Observation_Operator_Parallel(Myrank, NPROCS, LENSFC, LENSFC_proc, &
                 num_stn, num_Eval, RLA, RLO, Lat_stn, Lon_stn, OROG_at_stn, bkgst_srch_rad, &
                 SNDFCS, SNOOBS_stn, LANDMASK, SNOFCS_at_stn, index_back_atObs, index_back_atEval, &
                 Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts, Orog_obs_atEvalPts)
            if (myrank==PRINTRANK) PRINT*,'Finished stn observation operator'
 
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then 
                    PRINT*, "Background Indices at eval points"
                    PRINT*, index_back_atEval       
                    PRINT*, "Obs at Eval Points" 
                    PRINT*, Obs_atEvalPts   
                    PRINT*, "Forecast at Eval Points"
                    PRINT*, SNOFCS_atEvalPts                             
                    PRINT*, "Lat at Eval Points"
                    PRINT*, Lat_atEvalPts
                    PRINT*, "Lon at Eval Points"
                    PRINT*, Lon_atEvalPts
            endif

            OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn

            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
            !     PRINT*, "Lat at Obs Points"
            !     PRINT*, Lat_stn
            !     PRINT*, "Lon at Obs Points"
            !     PRINT*, Lon_stn
                ! PRINT*, "Model elevation at station locations from rank: ", MYRANK
                ! PRINT*, OROGFCS_at_stn
                PRINT*, "Background Indices at obs points"
                PRINT*, index_back_atObs
                PRINT*, "Background at station locations from rank: ", MYRANK
                PRINT*, SNOFCS_at_stn  
            !     PRINT*, "Obs at obs stns" 
            !     PRINT*, SNOOBS_stn   
                PRINT*, "O - B (innovation at obs points)"
                PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished observation operator for station data'         
        endif ! num_stn > 0

!=============================================================================================
! 5.  obs QC goes here
!=============================================================================================
! As of 2.4.22: Most qcs are being done at obs read time and inside the DA loop 
! for efficiency; 
! In the future consider consolidating QC's here

!CSDCSD - todo. Add QC here.

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
        if (myrank==PRINTRANK) PRINT*,'Starting DA loop'
        Do jndx = 1, LENSFC_proc    !mp_start, mp_end     
            ! if (myrank == 4 .and. (jndx == 3943 .or. jndx == 4033) ) then  
            ! QC: only update this grid cell if it is land.
            if( LANDMASK(jndx) == 1 ) then 
                num_loc_1 = 0
                num_loc_2 = 0
                assim_sncov_thisGridCell = .FALSE.
                if (print_debg_info) print*, "proc ", myrank, " grid: ", jndx
                if (num_stn>0) then 
                ! CSD - speed up by skipping over all grid cells where model and IMS agree on no snow? 
                ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                ! QC: gross error check is done in this call.
                    call nearest_Observations_Locations(RLA(jndx), RLO(jndx),    &
                            Lat_stn, Lon_stn, num_stn, obs_srch_rad, max_num_nearStn,   &
                            stdev_back, stdev_obsv, obs_tolerance,                 &
                            SNOFCS_at_stn, SNOOBS_stn,                                              &
                            index_at_nearStn,  num_loc_1) !,      &LENSFC,
                    if (print_debg_info) print*, "number of stn sndpth obs ", num_loc_1
                    
                endif         
                !     
                if( assim_SnowCov_obs  .AND. &  !.and. (IH_loc == ims_assm_hour)
                    (.NOT. IEEE_IS_NAN(SNDFCS(jndx))) .AND. &
                    (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(jndx))) .AND. &
                    ! (OROG(jndx) <= ims_max_ele) .AND. &
                    ( .not.(SCF_Grid(jndx) >= 0.99 .AND. & 
                            SNCOV_IMS(jndx) >= 0.99) ) ) then
                        num_loc_2 = 1 
                        assim_sncov_thisGridCell = .TRUE.                
                endif
                num_loc = num_loc_1 + num_loc_2                
                ! if assim_sncov=false >> num_loc_1=num_loc
                ! QC: only update this grid cell if it is land.
                if(num_loc > 0) then            !((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then 
                    ! get background states
                    Allocate(back_at_Obs(num_loc))
                    Allocate(obs_Array(num_loc))
                    Allocate(Lat_Obs(num_loc))
                    Allocate(Lon_Obs(num_loc))
                    Allocate(orogfcs_Obs(num_loc))                        
                    if(num_loc_1 > 0) then
                        index_obs_assmilated(1, jndx) = num_loc_1
                        Do zndx = 1, num_loc_1     
                            index_obs_assmilated(zndx+1, jndx) = index_at_nearStn(zndx)
                            back_at_Obs(zndx) = SNOFCS_at_stn(index_at_nearStn(zndx))
                            obs_Array(zndx) = SNOOBS_stn(index_at_nearStn(zndx))
                            Lat_Obs(zndx) = Lat_stn(index_at_nearStn(zndx))
                            Lon_Obs(zndx) = Lon_stn(index_at_nearStn(zndx))
                            orogfcs_Obs(zndx) = OROG_at_stn(index_at_nearStn(zndx)) 
                        End Do
                    End if
                    ! Append IMS-derived snow depth to the obs array 
                    if(assim_sncov_thisGridCell) then   
                        IMS_Foot_Print(jndx) = 1.0                          
                        back_at_Obs(num_loc) = SNDFCS(jndx)
                        obs_Array(num_loc) = SNO_IMS_at_Grid(jndx)
                        Lat_Obs(num_loc) = RLA(jndx)   
                        Lon_Obs(num_loc) = RLO(jndx) 
                        orogfcs_Obs(num_loc) = OROG(jndx)  
                    endif
                    ! compute covariances
                    if(exclude_obs_at_grid .and. (num_loc_1 > 1)) then
                        Allocate(B_cov_mat(num_loc-1, num_loc-1))
                        Allocate(b_cov_vect(num_loc-1))
                        Allocate(O_cov_mat(num_loc-1, num_loc-1))
                        Allocate(W_wght_vect(num_loc-1)) 
                        Allocate(obs_Innov(num_loc-1))  
                        call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), &    !SNDFCS(jndx),    &
                                Lat_Obs(2:num_loc), Lon_Obs(2:num_loc), orogfcs_Obs(2:num_loc), &
                                num_loc-1,   &
                                Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                                L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                                assim_sncov_thisGridCell,                          &
                                B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                   
                        call Snow_DA_OI(back_at_Obs(2:num_loc), obs_Array(2:num_loc), num_loc-1, &
                            W_wght_vect, SNDFCS(jndx), incr_at_Grid(jndx), SNDANL(jndx), obs_Innov)
                    else
                        Allocate(B_cov_mat(num_loc, num_loc))
                        Allocate(b_cov_vect(num_loc))
                        Allocate(O_cov_mat(num_loc, num_loc))
                        Allocate(W_wght_vect(num_loc)) 
                        Allocate(obs_Innov(num_loc))  
                        call compute_covariances(RLA(jndx), RLO(jndx), OROG(jndx), &    !SNDFCS(jndx),    &
                                Lat_Obs, Lon_Obs, orogfcs_Obs, num_loc,   &
                                Stdev_back, Stdev_Obsv, Stdev_Obsv_ims,      &
                                L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                                assim_sncov_thisGridCell,                          &
                                B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                   
                        call Snow_DA_OI(back_at_Obs, obs_Array, num_loc,   &
                            W_wght_vect, SNDFCS(jndx), incr_at_Grid(jndx), SNDANL(jndx), obs_Innov)
                    endif
                    if (myrank==PRINTRANK .and. print_debg_info) then  !
                            print*, "proc ", myrank, "loop ", jndx, &
                            "num depth obs ", num_loc_1, "total obs", num_loc
                            PRINT*, " background at obs pts: "
                            PRINT*, back_at_Obs     
                            PRINT*, "Observed"
                            PRINT*,  obs_Array
                            PRINT*, "Obs innovation: "
                            PRINT*, obs_Innov
                            PRINT*, "Weight vector: "
                            PRINT*, W_wght_vect  
                    endif  
                    if (myrank==PRINTRANK .and. print_debg_info) then  !
                        print*, "forec: ", SNDFCS(jndx)   
                        print*, "increment SND and SNCOV: ", incr_at_Grid(jndx), " anl: ", SNDANL(jndx)
                    endif
                 
                    !free mem
                    DEALLOCATE(back_at_Obs, obs_Array)
                    DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs, obs_Innov)
                    DEALLOCATE(B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)                            
                ! else 
                    ! no obs available: keep background (copied earlier)
                    ! SNOANL(jndx) = SNOFCS(jndx) 
                endif
                if (allocated(index_at_nearStn))  Deallocate(index_at_nearStn) 
                ! if (allocated(index_at_nearIMS))  Deallocate(index_at_nearIMS)                 
            endif ! non-glacier land
        ! endif
        End do
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==PRINTRANK) PRINT*, ' Finished DA loops'
        

!=============================================================================================
! 7. Clean up, write outputs 
!=============================================================================================
! collect results onto main tasks, if necessary
! fill in SWE and SND arrays       
        ! avoid -ve anl
        Where(SNDANL < 0.) SNDANL = 0.
        if (print_debg_info) then
            PRINT*, "Weighted increment SWE/snwd from rank: ", MYRANK
            PRINT*, incr_at_Grid       
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNDANL
        endif
        ! d = swe/sndi
        ! this is for writing outputs at observation and evaluation points
        WHERE (LANDMASK==1) SWEANL = SNDANL * SNODENS_Grid
! !Compute updated snocov 
        !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)

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
! outputs
        ! Enddo 
        if (MYRANK > 0) then  
            call MPI_SEND(SNDFCS, LENSFC_proc, mpiReal_size, 0, &
                            MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(SWEFCS, LENSFC_proc, mpiReal_size, 0, &
                            40000*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(incr_at_Grid, LENSFC_proc, mpiReal_size, 0, &
                            100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SNDANL, LENSFC_proc, mpiReal_size, 0, &
                            200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SWEANL, LENSFC_proc, mpiReal_size, 0, &
                            300*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SCF_Grid, LENSFC_proc, mpiReal_size, 0, &
                            400*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SNO_IMS_at_Grid, LENSFC_proc, mpiReal_size, 0, &
                            500*(MYRANK+1), MPI_COMM_WORLD, IERR)
        endif
        if (MYRANK == 0) then   
            SNDFCS_out(1:LENSFC_proc) = SNDFCS
            SWEFCS_out(1:LENSFC_proc) = SWEFCS
            incr_at_Grid_out(1:LENSFC_proc) = incr_at_Grid            
            SNDANL_out(1:LENSFC_proc) = SNDANL   
            SWEANL_out(1:LENSFC_proc) = SWEANL        
            SCF_Grid_out(1:LENSFC_proc) = SCF_Grid     
            SNO_IMS_at_Grid_out(1:LENSFC_proc) = SNO_IMS_at_Grid       
            Do ixy =  1, NPROCS - 1  ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1 
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(SNDFCS_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SWEFCS_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 40000*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(incr_at_Grid_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)                
                call MPI_RECV(SNDANL_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SWEANL_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SCF_Grid_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 400*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SNO_IMS_at_Grid_out(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 500*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
    ! copy values at eval points
            incr_atEvalPts = IEEE_VALUE(incr_atEvalPts, IEEE_QUIET_NAN)
            SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
            SNDnoDA_atEval = IEEE_VALUE(SNDnoDA_atEval, IEEE_QUIET_NAN)
            SWEnoDA_atEval = IEEE_VALUE(SWEnoDA_atEval, IEEE_QUIET_NAN)
            Orog_fcs_atEvalPts = IEEE_VALUE(Orog_fcs_atEvalPts, IEEE_QUIET_NAN)
            ! SND_ml_noDA_atEval = IEEE_VALUE(SND_ml_noDA_atEval, IEEE_QUIET_NAN)
            ! SNOANL_Cur_atEvalPts = IEEE_VALUE(SNOANL_Cur_atEvalPts, IEEE_QUIET_NAN)
            Do jndx = 1, num_Eval
                if (index_back_atEval(jndx) > 0)  then   !                     
                    ! SND_ml_noDA_atEval(jndx) = noahmp%snow_depth(index_back_atEval(jndx))
                    SNDnoDA_atEval(jndx) = SNDnoDA(index_back_atEval(jndx))
                    SWEnoDA_atEval(jndx) = SWEnoDA(index_back_atEval(jndx))
                    ! SNOANL_Cur_atEvalPts(jndx) = SNDANL_Cur(index_back_atEval(jndx))
                    incr_atEvalPts(jndx) = incr_at_Grid_out(index_back_atEval(jndx))
                    SNOANL_atEvalPts(jndx) = SNDANL_out(index_back_atEval(jndx))                 
                    Orog_fcs_atEvalPts(jndx) = OROG(index_back_atEval(jndx))
                endif
            End do  
            SNOANL_atStn = IEEE_VALUE(SNOANL_atStn, IEEE_QUIET_NAN)
            Do jndx = 1, num_stn
                if (index_back_atObs(jndx) > 0)  then   !       index_back_atObs()              
                    SNOANL_atStn(jndx) = SNDANL_out(index_back_atObs(jndx))                        
                endif
            End do         
        endif
        call MPI_BCAST(SNDANL_out, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        if (myrank==PRINTRANK) PRINT*,'Finished Data copy'

        select case (snowUpdateOpt)
            case (DO_NOTHING)
                write(*,*) " proc ", myrank, " Option 0: do nothing"
                noahmp%swe = SWEANL
                noahmp%snow_depth = SNDANL
            case (TOP_LAYER)
                write(*,*) " proc ",myrank, " Option 1: UpdateTopLayer"
                call UpdateTopLayer(myrank, LENSFC_proc, noahmp, incr_at_Grid, &
                        SNDANL, LANDMASK)
            case (BOTTOM_LAYER)
                write(*,*) " proc ",myrank, " Option 2: UpdateBottomLayer"
                call UpdateBottomLayer(myrank, LENSFC_proc, noahmp, incr_at_Grid, &
                        SNDANL, LANDMASK)
            case (ALL_LAYERS)
                write(*,*) " proc ",myrank, " Option 3: UpdateAllLayers"
                call UpdateAllLayers(myrank, LENSFC_proc, noahmp, incr_at_Grid, &
                        SNDANL, SNODENS_Grid, LANDMASK)
            case default
                write(*,*) "choose a valid partition option"
        end select 
        ! avoid -ve anl
        Where(noahmp%swe < 0) noahmp%swe = 0.
        Where(noahmp%snow_depth < 0) noahmp%snow_depth = 0.

        if (print_debg_info .and. myrank==PRINTRANK) then
            print*, "noah SWE ", noahmp%swe 
            print*, "noah snowdepth", noahmp%snow_depth
        endif

! write outputs 
        if (myrank==PRINTRANK) then 
            PRINT*,"proc ", myrank, "starting writing output"
        endif
        ! Write(rank_str, '(I3.3)') (MYRANK+1)
        da_out_file = trim(output_prefix)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"_18.nc"  !  
        call Write_DA_Outputs_restart_vector(MYRANK, NPROCS, da_out_file, forc_inp_file, &
            LENSFC, LENSFC_proc, num_stn, num_Eval, max_num_nearStn, &
            N_sA, N_sA_Ext, mpiReal_size, mpiInt_size, noahmp, &      !IDIM, JDIM, 
            SNDFCS_out, SWEFCS_out, incr_at_Grid_out, SNDANL_out, SWEANL_out, &       !SWEANL_Cur, SNDANL_Cur, LANDMASK, &  !
        Lat_stn, Lon_stn, OROG_at_stn, SNOOBS_stn, SNOFCS_at_stn, OmB_innov_at_stn,  &                                
                index_back_atObs, index_obs_assmilated,   &
                NEXC, Index_Obs_Excluded,   &
            SNOANL_atStn, SNCOV_IMS, IMS_Foot_Print, SCF_Grid_out, SNO_IMS_at_Grid_out, & ! , &
                Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
                SNOFCS_atEvalPts, incr_atEvalPts, SNOANL_atEvalPts, &
                SNDnoDA, SWEnoDA, SNDnoDA_atEval, SWEnoDA_atEval, index_back_atEval, &
                Orog_obs_atEvalPts, Orog_fcs_atEvalPts)   !, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov
        ! if (myrank==PRINTRANK) 
        PRINT*,"proc ", myrank, "finished writing output"
        ! endif
998 CONTINUE
        ! clean up
        if (allocated(SNDnoDA_atEval))   deallocate(SNDnoDA_atEval)
        if (allocated(SWEnoDA_atEval))  deallocate(SWEnoDA_atEval)
        ! if (allocated(SND_ml_noDA_atEval))  deallocate(SND_ml_noDA_atEval)        

        if (allocated(noahmp%swe))   deallocate(noahmp%swe)
        if (allocated(noahmp%snow_depth))  deallocate(noahmp%snow_depth)
        if (allocated(noahmp%active_snow_layers)) deallocate(noahmp%active_snow_layers)
        if (allocated(noahmp%swe_previous))  deallocate(noahmp%swe_previous)
        if (allocated(noahmp%snow_soil_interface))  deallocate(noahmp%snow_soil_interface)
        if (allocated(noahmp%temperature_snow))  deallocate(noahmp%temperature_snow)
        if (allocated(noahmp%snow_ice_layer))  deallocate(noahmp%snow_ice_layer)
        if (allocated(noahmp%snow_liq_layer))  deallocate(noahmp%snow_liq_layer)
        if (allocated(noahmp%temperature_soil))  deallocate(noahmp%temperature_soil)
        ! if (allocated(SNOFCS_at_stn_ML))   DEALLOCATE(SNOFCS_at_stn_ML)

        if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
        if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)
        if (allocated(SNOANL_atStn))   DEALLOCATE(SNOANL_atStn)
        
        if (allocated(OmB_innov_at_stn))   DEALLOCATE(OmB_innov_at_stn)
        if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
        if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
        if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs) 

        if (allocated(Index_Obs_Excluded)) DEALLOCATE(Index_Obs_Excluded)

        if (allocated(OROG_at_stn))  DEALLOCATE(OROG_at_stn) 
        if (allocated(index_back_atEval)) DEALLOCATE(index_back_atEval)
        if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
        if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts)
        if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
        if (allocated(Orog_obs_atEvalPts))   DEALLOCATE(Orog_obs_atEvalPts)
        if (allocated(Orog_fcs_atEvalPts))   DEALLOCATE(Orog_fcs_atEvalPts)

        if (allocated(SNOFCS_atEvalPts)) DEALLOCATE(SNOFCS_atEvalPts)
        if (allocated(incr_atEvalPts))  DEALLOCATE(incr_atEvalPts)
        if (allocated(SNOANL_atEvalPts)) DEALLOCATE(SNOANL_atEvalPts) 

        ! Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        ! IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==PRINTRANK) PRINT*,'Finished OI DA at datetime: ', y_str, m_str, d_str, h_str
        
    ! End do

999 CONTINUE
        !PRINT*,'Finished OI DA ON RANK: ', MYRANK
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        !STOP

        RETURN

 END subroutine OI_Snow_Analysis_NOAHMP

 subroutine EnKF_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, & !num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, & 
                L_horz , h_ver, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                ens_size, rcov_localize, ens_inflate,  &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
                vector_restart_prefix, vector_noda_prefix, static_prefix, output_prefix, &
                snowUpdateOpt, PRINTRANK, print_debg_info, &   !fv3_index, vector_inputs, &
                SNOANL_out, &   !SNDFCS_out, SWEANL_out, & incr_at_Grid_out, 
                Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)
   !----------------------------------------------------------------------
   ! Input arguments: 
         
    !----------------------------------------------------------------------
    IMPLICIT NONE
    !
    include 'mpif.h'

    integer, parameter :: dp = kind(1.d0)
   ! , enkf_ensrf
    INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &  
                              IY, IM, ID, IH, LENSFC, IVEGSRC   !, num_assim_steps 
    Integer, intent(in)    :: ens_size
    LOGICAL, intent(in)    :: assim_SnowPack_obs, assim_SnowCov_obs
    Logical, intent(in)    :: rcov_localize, ens_inflate               ! localize R Covariance
    !CHARACTER(len=4), intent(in)   :: stn_var ! should probably be called control_var
    CHARACTER(LEN=*), Intent(In)   :: STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME
    CHARACTER(LEN=*), Intent(In)   :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
    ! CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX, CURRENT_ANALYSIS_PREFIX
    CHARACTER(LEN=*), Intent(In)   :: vector_restart_prefix, vector_noda_prefix, &
                                      static_prefix, output_prefix

    REAL, intent(in)    :: PERCENT_OBS_WITHHELD   !, dT_Asssim
    Real, intent(In)    :: L_horz , h_ver, obs_tolerance
    Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
    INTEGER, intent(in) :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels
    INTEGER, intent(in) :: snowUpdateOpt, PRINTRANK
    REAL, intent(out)   :: SNOANL_out(LENSFC)
    LOGICAL             :: print_debg_info  !, fv3_index vector_inputs, 
        ! for mpi par
    INTEGER   :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
    Integer   :: LENSFC_proc, begloc, endloc
    LOGICAL             :: exclude_obs_at_grid

    CHARACTER(LEN=5)    :: TILE_NUM
    ! Character(LEN=3)    :: rank_str
    INTEGER	            :: IERR	
    REAL             :: SNDFCS(ens_size+1, LENSFC_proc), SWEFCS(ens_size+1, LENSFC_proc), &
                        SNDANL(ens_size+1, LENSFC_proc), SWEANL(ens_size+1, LENSFC_proc)
    REAL             :: SNDnoDA(LENSFC_proc), SWEnoDA(LENSFC_proc)       !, SNOANL_Cur(LENSFC)
    
    REAL   :: RLA(LENSFC_proc), RLO(LENSFC_proc), RLO_Tile(LENSFC_proc), OROG(LENSFC_proc)  !, OROG_UF(LENSFC)
    REAL            :: VETFCS(LENSFC_proc), SNUP_Array(LENSFC_proc)           ! SNOFCS(LENSFC), 
    Integer         :: tile_xy(LENSFC_proc), Idim_xy(LENSFC_proc), Jdim_xy(LENSFC_proc)
    INTEGER             :: LANDMASK(LENSFC_proc)  !, DAMASK(LENSFC)

    CHARACTER(len=250)   :: ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
    CHARACTER(len=4)     :: y_str, m_str, d_Str, h_str, fvs_tile
    REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:), SNOANL_atStn(:)  !, SNOFCS_at_stn_ML(:)                
    REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROG_at_stn(:)  
    REAL                 :: lat_min, lat_max, lon_min, lon_max      
    Real                 :: SNCOV_IMS(LENSFC_proc)  ! ims resampled at each grid
    Real                 :: SNO_IMS_at_Grid(ens_size+1, LENSFC_proc)
    Real                 :: IMS_Foot_Print(LENSFC_proc)  ! grid cells where IMS was assimilated

    INTEGER                :: num_stn, Num_Ims, num_Eval !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon

    Integer, Allocatable   :: index_at_nearStn(:)  !, index_at_nearIMS(:) !loc_near_Obs(:), 
    Integer                :: num_loc, num_loc_1, num_loc_2
    Real, Parameter        :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 80. ! mm 
    real                   :: stdev_obsv, stdev_back
    ! Integer                :: ims_assm_hour
	INTEGER                :: jndx, zndx, ncol, nrow

    !ens
    Real, ALLOCATABLE    :: SNOFCS_atObs_ens(:,:)  !SWEFCS_Ens(:,:), SNDFCS_Ens(:,:), , back_at_Obs_ens(:,:)

!1.4 track obserations assimilated
    Integer                :: index_obs_assmilated(max_num_nearStn+1, LENSFC_proc)

    Real, Allocatable        :: obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:) !back_at_Obs(:), 
    REAL                     :: incr_at_Grid(ens_size+1, LENSFC_proc) !, incr_at_Grid_ensM(LENSFC)    ! increment at grid
    Real, Allocatable        :: obs_Innov(:)  !, OmB_innov_at_stn(:)

    CHARACTER(len=250)       :: forc_inp_file, da_out_file, noda_inp_path
    CHARACTER(LEN=4)         :: RANKCH 
    CHARACTER(LEN=500)       :: static_filename

    ! REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), incr_atEvalPts(:), SNOANL_atEvalPts(:) !, SNOANL_Cur_atEvalPts(:)  !evalution points 
    REAL, ALLOCATABLE    :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
    ! REAL, ALLOCATABLE  :: Orog_obs_atEvalPts(:), Orog_fcs_atEvalPts(:)
    Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
    Integer                 :: NEXC 
    Integer, ALLOCATABLE    :: index_back_atObs(:), Index_Obs_Excluded(:), &
                               index_obs_atEval(:)   ! the location of background corresponding obs
   
    Real               :: snodens, SNODENS_Grid(ens_size+1, LENSFC_proc)
    !LOGICAL            :: assim_snpack_stn, assim_SWE    !use swe as model state var? (instead of snow depth)
    LOGICAL            :: assim_sncov_thisGridCell    !    assim_sncov, assimilate sncov, 

    Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice

    INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, &
                            dest_Aoffset_end, arLen, pindex
    INTEGER            :: mpiReal_size, rsize, isize, mpiInt_size, ixy, ie
    REAL               :: tmp
    INTEGER            :: istep, IY_loc, IM_loc, ID_loc, IH_loc
    REAL               :: IH_real

        ! INTEGER, PARAMETER :: PRINTRANK = 4 
    ! CSD-todo, should be same as in sfcsub. Share properly
    Real, parameter    :: nodata_val = -9999.9

    integer, parameter :: lsm_type = 2

!   Snow partion options for Noah-MP
    integer, parameter     :: DO_NOTHING = 0
    integer, parameter     :: TOP_LAYER = 1
    integer, parameter     :: BOTTOM_LAYER = 2
    integer, parameter     :: ALL_LAYERS = 3

    ! noah mp
    Real  :: SCF_Grid(ens_size+1, LENSFC_proc)
    !Real  :: SCF_Grid_Land(ens_size+1, LENSFC_proc)
    Real     :: noahmp_ensm_swe(LENSFC_proc), noahmp_ensm_snow_depth(LENSFC_proc)

    type(noahmp_type)      :: noahmp(ens_size)
    ! type(observation_type) :: obs    
!   noahmp%model%weasd  = this%snow_water_equivalent
!   noahmp%model%snwdph = this%snow_depth * 1000.0  ! driver wants mm
!   noahmp%model%canopy = this%canopy_water
!   noahmp%model%tskin  = this%skin_temperature
!   noahmp%model%stc    = this%soil_temperature
!   noahmp%model%smc    = this%soil_moisture
!   noahmp%model%slc    = this%soil_liquid
!   noahmp%model%tsurf  = noahmp%model%tskin

!   ! type noahmp_type
!     snwdph  ! depth 
!     weasd     ! SWE 
!     snowc   !fractional snow cover
!     sncovr1 ! snow cover over land
!     tsurf ! surface skin temperature
!     stc(time, soil_levels, location) !soil temperature
!     tsnoxy(time, snow_levels, location) ! snow temperature
    Do ie = 1, ens_size
        allocate(noahmp(ie)%swe                (LENSFC_proc))
        allocate(noahmp(ie)%snow_depth         (LENSFC_proc))
        allocate(noahmp(ie)%active_snow_layers (LENSFC_proc))
        allocate(noahmp(ie)%swe_previous       (LENSFC_proc))
        allocate(noahmp(ie)%snow_soil_interface(LENSFC_proc,7))
        allocate(noahmp(ie)%temperature_snow   (LENSFC_proc,3))
        allocate(noahmp(ie)%snow_ice_layer     (LENSFC_proc,3))
        allocate(noahmp(ie)%snow_liq_layer     (LENSFC_proc,3))
        allocate(noahmp(ie)%temperature_soil   (LENSFC_proc)) 
    End do       
    ! end type noahmp_type 

    !=============================================================================================
    ! 1. initialise vars,set-up processors, and read lat/lon from orog files.
    !=============================================================================================

    !initialse output with nodata
    IMS_Foot_Print = IEEE_VALUE(IMS_Foot_Print, IEEE_QUIET_NAN)
    stdev_obsv = stdev_obsv_depth
    stdev_back = stdev_back_depth  
    !obs_srch_rad = 250. ! radius of observation search
    ! ims_assm_hour = 18 
    ! noah models specific? Needed to ID glaciers.
    if (IVEGSRC == 2) then   ! sib
            veg_type_landice=13
    else
            veg_type_landice=15
    endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    !total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
    IF (myrank ==PRINTRANK) PRINT*,"snowDA: EnSRF total num proc ", NPROCS, " Num tiles : ", NUM_TILES
  
! if (p_tN /= 4 ) goto 999

! must have all obs
    if ((STN_OBS_PREFIX(1:8).eq.'        ') .OR. (IMS_SNOWCOVER_PATH(1:8).eq.'        ') &
            .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
            print*, "One or more observation paths don't exist!, skipping the DA"
            goto 999
    end if

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
    write(fvs_tile, "(I0.2)") IDIM
    static_filename = trim(static_prefix)//"/ufs-land_C"//trim(fvs_tile)// &
                    "_static_fields.nc"
    Call ReadTileInfo(static_filename, LENSFC_proc, veg_type_landice, &
            mp_start, mp_end, &
            tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 
    ! make RLO copy before so that RLO (used later) is not modified
    RLO_Tile = RLO 
    Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
    ! 1 degree buffer around tile boundaries (but not to go out of valid coords)
    lat_min = MAX(MINVAL(RLA) - 1., -90.)
    lat_max = MIN(MAXVAL(RLA) + 1., 90.)
    lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
    lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)      
  
    if ((p_tRank==PRINTRANK) ) then !.and. print_deb) then
            print*, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
    endif

! If multiple time steps are simulated in a given time period (window) given by num_assim_steps * dT_Assim
! Note: the DA outputs from the last time step are returned        
        IH_real = IH; IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH; ! these variables change inside loop below
    ! Do istep = 1, num_assim_steps
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc

        ! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 
        num_Eval = 0

!=============================================================================================
! 2. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
!=============================================================================================

       ! READ THE INPUT SURFACE DATA from vector and no da outputs
        noda_inp_path=TRIM(vector_noda_prefix)//"/ufs_land_output."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
        Call ReadRestartNoahMP_Ens(myrank, LENSFC_proc, vector_restart_prefix, noda_inp_path, &
                        y_str, m_str, d_str, h_str, ens_size, mp_start, mp_end, &
                        noahmp, SNDnoDA, SWEnoDA, SCF_Grid)  !, SNDFCS, SWEFCS)
        ! initialise analysis to forecast (mostly to retain non-land snow states, which are not updated)
        SWEFCS(ens_size+1,:) = 0.0
        SNDFCS(ens_size+1,:) = 0.0
        ! SNODENS_Grid(ens_size+1,:) = 0.0
        SCF_Grid(ens_size+1,:) = 0.0
        !SCF_Grid_Land(ens_size+1,:) = 0.0
        Do ie = 1, ens_size
            SWEFCS(ie,:) = noahmp(ie)%swe(:)
            SNDFCS(ie,:) = noahmp(ie)%snow_depth(:)
            ! grid snow density
            call calc_density(LENSFC_proc, lsm_type, LANDMASK, SWEFCS(ie,:), SNDFCS(ie,:), &
                               noahmp(ie)%temperature_soil, SNODENS_Grid(ie,:))
            SWEFCS(ens_size+1,:) = SWEFCS(ens_size+1,:) + SWEFCS(ie,:)
            SNDFCS(ens_size+1,:) = SNDFCS(ens_size+1,:) + SNDFCS(ie,:)
            SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) + SCF_Grid(ie,:)
            !SCF_Grid_Land(ens_size+1,:) = SCF_Grid_Land(ens_size+1,:) + SCF_Grid_Land(ie,:)
        End do         
        ! SUM(SNDFCS(1:ens_size, jndx)) / ens_size
        SWEFCS(ens_size+1,:) = SWEFCS(ens_size+1,:) / ens_size
        SNDFCS(ens_size+1,:) = SNDFCS(ens_size+1,:) / ens_size
        SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) / ens_size
        !SCF_Grid_Land(ens_size+1,:) = SCF_Grid_Land(ens_size+1,:) / ens_size

        SWEANL(:,:) = SWEFCS(:,:)
        SNDANL(:,:) = SNDFCS(:,:)
        incr_at_Grid = 0.0

        tmp = SUM(SWEFCS(ens_size+1,:),  Mask = (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>=0 )) &          !
                         / COUNT (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>= 0)                 !
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' ensemble mean SWE', tmp
        tmp = SUM(SNDFCS(ens_size+1,:),  Mask = (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>=0 )) &              !
                         / COUNT (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>= 0)     
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' ensemble mean SND', tmp  
        SNODENS_Grid(ens_size+1,:) =  SWEFCS(ens_size+1,:)/SNDFCS(ens_size+1,:)

!=============================================================================================
! 3. Read observations
!=============================================================================================

! 3a. Read station obs (of snow depth or SWE)
        if (assim_SnowPack_obs) then 
            ghcnd_inp_file = TRIM(STN_OBS_PREFIX)// &   !GHCND.SNWD TRIM(h_str)//
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//".nc"
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
            
            Call Observation_Read_GHCND_IODA(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
                            TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
                            num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn, &
                            NEXC, Index_Obs_Excluded)
            ! Call Observation_Read_GHCND_All_excNaN(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
            !                 TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
            !                 num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn)

            if (myrank==PRINTRANK) then   !.and. (p_tN==PRINTRANK).and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNOOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
                    PRINT*, "Elevation at station locations from rank: ", MYRANK
                    PRINT*, OROG_at_stn
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished reading station data'
        endif
! CSD beyond this point, there should be no specific mention of the station data source

! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 
        if (assim_SnowCov_obs) then
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
            ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
            WRITE(RANKCH, '(I0.1)') myrank                            
            ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                        ".IMS.Indices."//trim(ADJUSTL(RANKCH))//".nc"                       
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading IMS file', trim(ims_inp_file) 
            ! Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
            !                         MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
            Call Read_IMS_and_Resample_toTarget(ims_inp_file, IMS_INDEXES_PATH, & !fv3_index, 
                    LENSFC_proc, num_subgrd_ims_cels, IDIM, JDIM, &
                    tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
            if(myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV from rank: ", MYRANK
                PRINT*, SNCOV_IMS
            endif 
            if (myrank==PRINTRANK) PRINT*,'Finished reading SNCOV, converting to snow depth' 
           ! SNUP array will be used later, to test whether SCF > 100%
            ! call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC_proc, &
            !                            SNO_IMS_at_Grid, SNUP_Array)
            ! SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_Grid(ens_size+1,:) ! convert SWE to SND
            SNO_IMS_at_Grid(ens_size+1,:) = 0.0
            Do ie = 1, ens_size
                call calcSD_noahmp(LENSFC_proc, SNCOV_IMS, VETFCS, &
                                   SNODENS_Grid(ie,:), SNO_IMS_at_Grid(ie,:))
                SNO_IMS_at_Grid(ens_size+1,:) = SNO_IMS_at_Grid(ens_size+1,:) + &
                                                SNO_IMS_at_Grid(ie,:)
            End do
            SNO_IMS_at_Grid(ens_size+1,:) = SNO_IMS_at_Grid(ens_size+1,:) / ens_size
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV derived snodepth at each grid cell from rank: ", MYRANK
                PRINT*, SNO_IMS_at_Grid
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished converting SNCOV observations'
        
        endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!     snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNOFCS_at_stn(num_stn))
            ALLOCATE(SNOANL_atStn(num_stn))
            ! Allocate(SNOFCS_at_stn_ML(num_stn))
            ! ALLOCATE(OROGFCS_at_stn(num_stn)) 
            ALLOCATE(index_back_atObs(num_stn)) 
        
	! ALLOCATE(OmB_innov_at_stn_ens(num_stn, ens_size)) 
        ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
            num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
            if (num_Eval > 0) then 
                ALLOCATE(index_back_atEval(num_Eval)) 
                ALLOCATE(Obs_atEvalPts(num_Eval)) 
                ! ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
                ALLOCATE(Lat_atEvalPts(num_Eval))
                ALLOCATE(Lon_atEvalPts(num_Eval))            
                ALLOCATE(index_obs_atEval(num_Eval)) 
                
                ! if(p_tRank == PRINTRANK) then 
                        PRINT*,"Proc ", myrank," ",num_Eval,' points for evaluation excluded from DA'       
                ! endif   
            endif 
! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)
!           for JEDI, probably also want to add snow depth derived from snow cover here.
            ! Call Observation_Operator_tiles_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, & 
            !                     RLA, RLO, OROG, Lat_stn, Lon_stn,   &
            !                     LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNDFCS, SNOOBS_stn, LANDMASK,  &
            !                     SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval,  &
            !                     Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            Call Observation_Operator_Parallel_vect(Myrank, NPROCS, LENSFC, LENSFC_proc, &
                 num_stn, num_Eval, RLA, RLO, Lat_stn, Lon_stn, OROG_at_stn, bkgst_srch_rad, &
                 SNDFCS(ens_size+1,:), SNOOBS_stn, LANDMASK, SNOFCS_at_stn, &
                 index_back_atObs, index_back_atEval, index_obs_atEval, &
                 Obs_atEvalPts, Lat_atEvalPts, Lon_atEvalPts) !SNOFCS_atEvalPts, , Orog_obs_atEvalPts)
            if (myrank==PRINTRANK) PRINT*,'Finished stn observation operator'
            ! if(snowUpdateOpt > 0) then 
            !     SNOFCS_at_stn_ML = IEEE_VALUE(SNOFCS_at_stn_ML, IEEE_QUIET_NAN) 
            !     do jndx = 1, num_stn
            !         if(index_back_atObs(jndx) > 0) then
            !             SNOFCS_at_stn_ML(jndx)= noahmp%snow_depth(index_back_atObs(jndx))
            !         endif
            !     enddo 
            ! endif
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then 
                PRINT*, "Background Indices at eval points"
                PRINT*, index_back_atEval       
                PRINT*, "Obs at Eval Points" 
                PRINT*, Obs_atEvalPts   
                ! PRINT*, "Forecast at Eval Points"
                ! PRINT*, SNOFCS_atEvalPts                             
                PRINT*, "Lat at Eval Points"
                PRINT*, Lat_atEvalPts
                PRINT*, "Lon at Eval Points"
                PRINT*, Lon_atEvalPts
            endif
            ! OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn
            ALLOCATE(SNOFCS_atObs_ens(ens_size, num_stn))
            SNOFCS_atObs_ens = IEEE_VALUE(SNOFCS_atObs_ens, IEEE_QUIET_NAN)
            Do jndx = 1, num_stn
                if (index_back_atObs(jndx) > 0) then
                    SNOFCS_atObs_ens(:, jndx) = SNDFCS(1:ens_size, index_back_atObs(jndx))
                endif
            end do
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
            !     PRINT*, "Lat at Obs Points"
            !     PRINT*, Lat_stn
            !     PRINT*, "Lon at Obs Points"
            !     PRINT*, Lon_stn
                    ! PRINT*, "Model elevation at station locations from rank: ", MYRANK
                    ! PRINT*, OROGFCS_at_stn
                PRINT*, "Background Indices at obs points"
                PRINT*, index_back_atObs
                 PRINT*, "Background at station locations from rank: ", MYRANK
                PRINT*, SNOFCS_at_stn  
            !     PRINT*, "Obs at obs stns" 
            !     PRINT*, SNOOBS_stn   
                ! PRINT*, "O - B (innovation at obs points)"
                ! PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished observation operator for station data'         
        endif ! num_stn > 0

!=============================================================================================
! 5.  obs QC goes here
!=============================================================================================
! As of 2.4.22: Most qcs are being done at obs read time and inside the DA loop 
! for efficiency; 
! In the future consider consolidating QC's here

!CSDCSD - todo. Add QC here.

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

        !Do nrow = 1, Num_Snotel !JDIM/2
        !obs_srch_rad = 250. ! radius of observation search
        if (myrank==PRINTRANK) PRINT*,'Starting DA loop'
        Do jndx = 1, LENSFC_proc    !mp_start, mp_end     
            ! if (myrank == 4 .and. (jndx == 3943 .or. jndx == 4033) ) then  
            ! QC: only update this grid cell if it is land.
            if( LANDMASK(jndx) == 1 ) then 
                num_loc_1 = 0
                num_loc_2 = 0
                assim_sncov_thisGridCell = .FALSE.
                if (print_debg_info) print*, "proc ", myrank, " grid: ", jndx
                if (num_stn>0) then 
                ! CSD - speed up by skipping over all grid cells where model and IMS agree on no snow? 
                ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                ! QC: gross error check is done in this call.
                ! 8.19.20: we are assuming here if deterministic forecast exists/not null
	            ! at a point, then ensemble members also have valid value

! 10.25.21 Note the QC based on forecast ensemble mean
! may need to use a different (perhaps larger) obs_tolerance
                    call nearest_Observations_Locations(RLA(jndx), RLO(jndx),          &
                            Lat_stn, Lon_stn, num_stn, obs_srch_rad, max_num_nearStn,  &
                            stdev_back, stdev_obsv, obs_tolerance,             &
                            SNOFCS_at_stn, SNOOBS_stn,                         &
                            index_at_nearStn,  num_loc_1) !,      &LENSFC,
                    if (print_debg_info) print*, "number of stn sndpth obs ", num_loc_1
                endif
! 10.25.21 Note the QC based on forecast ensemble mean  
                if( assim_SnowCov_obs .AND. &
                    (.NOT. IEEE_IS_NAN(SNDFCS(ens_size+1,jndx))) .AND. &   !.and. (IH_loc == ims_assm_hour) 
                    (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(ens_size+1, jndx))) .AND. &
                !   ( .not.(SWEFCS(ens_size+1,jndx) >= SNUP_Array(jndx) .AND. & 
                    ( .not.(SCF_Grid(ens_size+1,jndx) >= 0.99 .AND. & 
                            SNCOV_IMS(jndx) >= 0.99) )  &
! 10.25.21 CSD decided not to use this criteria anymore
                    !(OROG(jndx) <= ims_max_ele) .AND. &
                    ) then
                        num_loc_2 = 1 
                        assim_sncov_thisGridCell = .TRUE.                
                endif
                num_loc = num_loc_1 + num_loc_2                
                ! if assim_sncov=false >> num_loc_1=num_loc
                ! QC: only update this grid cell if it is land.	
                if( num_loc > 0 ) then            !((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then 
                    ! Allocate(back_at_Obs(num_loc))
                    Allocate(obs_Array(num_loc))
                    Allocate(Lat_Obs(num_loc))
                    Allocate(Lon_Obs(num_loc))
                    Allocate(orogfcs_Obs(num_loc))  
                    ! ghcnd
                    if(num_loc_1 > 0) then
                        index_obs_assmilated(1, jndx) = num_loc_1
                        Do zndx = 1, num_loc_1    
                            index_obs_assmilated(zndx+1, jndx) = index_at_nearStn(zndx) 
                            ! back_at_Obs(zndx) = SNOFCS_at_stn(index_at_nearStn(zndx))
                            obs_Array(zndx) = SNOOBS_stn(index_at_nearStn(zndx))
                            Lat_Obs(zndx) = Lat_stn(index_at_nearStn(zndx))
                            Lon_Obs(zndx) = Lon_stn(index_at_nearStn(zndx))
                            orogfcs_Obs(zndx) = OROG_at_stn(index_at_nearStn(zndx)) 
                        End Do
                    End if
                    ! Append IMS-derived snow depth to the obs array 
                    if(assim_sncov_thisGridCell) then   
                        IMS_Foot_Print(jndx) = 1.0                           
                        ! back_at_Obs(num_loc) = SNDFCS(jndx)
                        obs_Array(num_loc) = SNO_IMS_at_Grid(ens_size+1, jndx)
                        Lat_Obs(num_loc) = RLA(jndx)   
                        Lon_Obs(num_loc) = RLO(jndx) 
                        orogfcs_Obs(num_loc) = OROG(jndx)
                    endif	
                    ! EnKF 
                    if(exclude_obs_at_grid .and. (num_loc_1 > 1)) then 
                        Allocate(obs_Innov(num_loc-1))
                        Call snow_DA_EnKF_LocCov(RLA(jndx), RLO(jndx), OROG(jndx),   &
                            num_loc_1-1, num_loc-1, num_stn, jndx, ens_size, &
                            Lat_Obs(2:num_loc), Lon_Obs(2:num_loc), orogfcs_Obs(2:num_loc),      &
                            L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
                            assim_sncov_thisGridCell, rcov_localize,       & !ens_inflate, 
                            SNDFCS(1:ens_size, jndx),   & !LENSFC, 
                            index_at_nearStn(2:num_loc), SNOFCS_atObs_ens,	 &
                            stdev_obsv, Stdev_Obsv_ims,         &    !stdev_back,    &
                            obs_Array(2:num_loc),                          &
                            obs_Innov, incr_at_Grid(:, jndx), SNDANL(:, jndx)) 	
                    else
                        Allocate(obs_Innov(num_loc))
                        Call snow_DA_EnKF_LocCov(RLA(jndx), RLO(jndx), OROG(jndx),   &
                            num_loc_1, num_loc, num_stn, jndx, ens_size, &
                            Lat_Obs, Lon_Obs, orogfcs_Obs,                              &
                            L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
                            assim_sncov_thisGridCell, rcov_localize,       & !ens_inflate, 
                            SNDFCS(1:ens_size, jndx),   & !LENSFC, 
                            index_at_nearStn, SNOFCS_atObs_ens,	 &
                            stdev_obsv, Stdev_Obsv_ims,         &    !stdev_back,    &
                            obs_Array,                          &
                            obs_Innov, incr_at_Grid(:, jndx), SNDANL(:, jndx)) 	
                    endif
                    if (myrank==PRINTRANK .and. print_debg_info) then  !
                        print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1
                        print*, "total obs", num_loc, " ens size ", ens_size
                        PRINT*, "Ens background: "
                        PRINT*, SNDFCS(:, jndx)	
                        PRINT*, "Observed"
                        PRINT*,  obs_Array
                        PRINT*, "Obs innovation: "
                        PRINT*, obs_Innov
                        PRINT*, "Ens increment at grid: "
                        PRINT*, incr_at_Grid(:, jndx)
                        PRINT*, "Ens analysis at grid: "
                        PRINT*, SNDANL(:, jndx)
                    endif					
                    DEALLOCATE(obs_Array, obs_Innov)
                    DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs)
                ! else
                !     SNDANL(:, jndx) = SNDFCS(:,jndx) !anl_at_Grid_ens(1:ens_size, jndx) = SNDFCS_Ens(:, jndx)
                !     SNDANL(ens_size+1, jndx) = SUM(SNDFCS(1:ens_size, jndx)) / ens_size
                endif
                if (allocated(index_at_nearStn))  Deallocate(index_at_nearStn)  
            endif ! not a land cell  
	    End do
        !End do
        if (myrank==PRINTRANK) PRINT*, 'Finished DA loops'
	! PRINT*, 'Finished DA loops', ' proc:', myrank
!=============================================================================================
! 7. Clean up, write outputs 
!=============================================================================================
! collect results onto main tasks, if necessary
! fill in SWE and SND arrays       
        ! avoid -ve anl
        Where(SNDANL < 0.) SNDANL = 0.
        if (print_debg_info) then
            PRINT*, "Weighted increment SWE/snwd from rank: ", MYRANK
            PRINT*, incr_at_Grid       
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNDANL
        endif
        ! d = swe/sndi
        ! this is for writing outputs at observation and evaluation points
        ! SWEANL = SNDANL * SNODENS_Grid
        Do ie = 1, ens_size+1
            WHERE (LANDMASK==1) SWEANL(ie,:) = SNDANL(ie,:) * SNODENS_Grid(ie,:)
        Enddo
        
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
! outputs
        if (MYRANK > 0) then  
            call MPI_SEND(SNDANL(ens_size+1, :), LENSFC_proc, mpiReal_size, 0, &
                            500*(MYRANK+1), MPI_COMM_WORLD, IERR) 
        endif
        if (MYRANK == 0) then
            SNOANL_out(1:LENSFC_proc) = SNDANL(ens_size+1, :)         
            Do ixy =  1, NPROCS - 1  ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1 
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(SNOANL_out(dest_Aoffset:dest_Aoffset_end),arLen, mpiReal_size, &
                            ixy, 500*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        endif
        call MPI_BCAST(SNOANL_out, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        if (myrank==PRINTRANK) PRINT*,'Finished Data copy'
        ! PRINT*, 'Finished Data copy', ' proc:', myrank
        ! if (MYRANK > NUM_TILES - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998
        noahmp_ensm_swe = 0.0
        noahmp_ensm_snow_depth = 0.0
        Do ie = 1, ens_size             
            select case (snowUpdateOpt)
            case (DO_NOTHING)
                ! write(*,*) " proc ", myrank, " Option 0: do nothing"
                noahmp(ie)%swe = SWEANL(ie,:)
                noahmp(ie)%snow_depth = SNDANL(ie,:)
            case (TOP_LAYER)
                ! write(*,*) " proc ",myrank, " Option 1: UpdateTopLayer"
                call UpdateTopLayer(myrank, LENSFC_proc, noahmp(ie), &
                incr_at_Grid(ie,:), SNDANL(ie,:), LANDMASK)
            case (BOTTOM_LAYER)
                ! write(*,*) " proc ",myrank, " Option 2: UpdateBottomLayer"
                call UpdateBottomLayer(myrank, LENSFC_proc, noahmp(ie), &
                incr_at_Grid(ie,:), SNDANL(ie,:), LANDMASK)
            case (ALL_LAYERS)
                ! write(*,*) " proc ",myrank, " Option 3: UpdateAllLayers"
                ! call UpdateAllLayers_Ens(ie, myrank, LENSFC_proc, noahmp(ie), &
                ! incr_at_Grid(ie,:), SNDANL(ie,:), SNODENS_Grid(ie,:), LANDMASK)
                call UpdateAllLayers(myrank, LENSFC_proc, noahmp(ie), incr_at_Grid(ie,:), &
                                SNDANL(ie,:), SNODENS_Grid(ie,:), LANDMASK)
            case default
                write(*,*) "choose a valid partition option"
                stop
            end select 
            ! avoid -ve anl
            ! print*, " proc ", myrank, "done partitioning ens ", ie
            Where(noahmp(ie)%swe < 0) noahmp(ie)%swe = 0.
            Where(noahmp(ie)%snow_depth < 0) noahmp(ie)%snow_depth = 0.
            ! ens mean
            noahmp_ensm_swe = noahmp_ensm_swe + noahmp(ie)%swe
            noahmp_ensm_snow_depth = noahmp_ensm_snow_depth + noahmp(ie)%snow_depth
            ! print*, " proc ", myrank, "done loop ", ie
        enddo
        ! SUM(SNDFCS(1:ens_size, jndx)) / ens_size
        noahmp_ensm_swe = noahmp_ensm_swe / ens_size
        noahmp_ensm_snow_depth = noahmp_ensm_snow_depth / ens_size
        if (print_debg_info .and. myrank==PRINTRANK) then
            print*, "noah SWE ", noahmp_ensm_swe
            print*, "noah snowdepth", noahmp_ensm_snow_depth
        endif
! !Compute updated snocov 
        !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)
        ! update ensemble scf
        SNODENS_Grid(ens_size+1,:) =  noahmp_ensm_swe/noahmp_ensm_snow_depth
        SCF_Grid(ens_size+1,:) = 0.0
        Do ie = 1, ens_size 
            ! grid snow density
            ! call calc_density(LENSFC_proc, lsm_type, LANDMASK, noahmp(ie)%swe, &
            ! noahmp(ie)%snow_depth, noahmp(ie)%temperature_soil, SNODENS_Grid(ie,:))             
            call calcSCF_noahmp(LENSFC_proc, VETFCS, SNODENS_Grid(ie,:), &
                noahmp(ie)%snow_depth, SCF_Grid(ie, :)) 
            SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) + SCF_Grid(ie,:)
! update  sncovr1 "snow cover over land" 
        Enddo 
        SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) / ens_size
        
        ! write outputs	
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==PRINTRANK) then 
            PRINT*,"proc ", myrank, "starting writing output"
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        da_out_file = trim(output_prefix)// &
                      TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"_18.nc"  !     
        call Write_DA_Outputs_ens_vect(da_out_file, vector_restart_prefix, &
                y_str, m_str, d_str, h_str, &        
                MYRANK, NPROCS, LENSFC, LENSFC_proc, ens_size, num_stn, num_Eval, &
                max_num_nearStn,  &
                N_sA, N_sA_Ext, mpiReal_size, mpiInt_size,     &
                SNDnoDA, SWEnoDA, noahmp,   &
                SNDFCS,  SWEFCS, SNDANL, SWEANL, incr_at_Grid,   &  ! 
                index_obs_assmilated,   &
                SNCOV_IMS, IMS_Foot_Print, SNO_IMS_at_Grid(ens_size+1,:),  &
                SNODENS_Grid, SCF_Grid, &  ! SAVe this as snowc
                noahmp_ensm_swe, noahmp_ensm_snow_depth, &       
                Lat_stn, Lon_stn, OROG_at_stn, SNOOBS_stn, index_back_atObs, &
                NEXC, Index_Obs_Excluded, &
                Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, &
                index_back_atEval, index_obs_atEval)            
	! if (myrank==PRINTRANK) 
        PRINT*,"proc ", myrank, "finished writing output"
    ! endif
    998 CONTINUE
        ! clean up      
        Do ie = 1, ens_size
            if (allocated(noahmp(ie)%swe))   deallocate(noahmp(ie)%swe)
            if (allocated(noahmp(ie)%snow_depth))  deallocate(noahmp(ie)%snow_depth)
            if (allocated(noahmp(ie)%active_snow_layers)) deallocate(noahmp(ie)%active_snow_layers)
            if (allocated(noahmp(ie)%swe_previous))  deallocate(noahmp(ie)%swe_previous)
            if (allocated(noahmp(ie)%snow_soil_interface))  deallocate(noahmp(ie)%snow_soil_interface)
            if (allocated(noahmp(ie)%temperature_snow))  deallocate(noahmp(ie)%temperature_snow)
            if (allocated(noahmp(ie)%snow_ice_layer))  deallocate(noahmp(ie)%snow_ice_layer)
            if (allocated(noahmp(ie)%snow_liq_layer))  deallocate(noahmp(ie)%snow_liq_layer)
            if (allocated(noahmp(ie)%temperature_soil))  deallocate(noahmp(ie)%temperature_soil)
        End do 
        if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)	
        if (allocated(SNOANL_atStn))   DEALLOCATE(SNOANL_atStn)
        if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
        if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
        if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
        if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)
        if (allocated(Index_Obs_Excluded)) DEALLOCATE(Index_Obs_Excluded)
        if (allocated(OROG_at_stn))  DEALLOCATE(OROG_at_stn)        
        if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
        if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts)
        if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
        if (allocated(index_back_atEval)) DEALLOCATE(index_back_atEval)
        if (allocated(index_obs_atEval))   DEALLOCATE(index_obs_atEval)

        ! Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        ! IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==PRINTRANK) PRINT*,'Finished EnSRF DA at datetime: ', y_str, m_str, d_str, h_str
    ! End do

999 CONTINUE
    !PRINT*,'Finished OI DA ON RANK: ', MYRANK
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    !STOP

    RETURN

END subroutine EnKF_Snow_Analysis_NOAHMP

subroutine EnSRF_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, & !num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, & 
                L_horz , h_ver, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                ens_size, rcov_localize, ens_inflate,  &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
                vector_restart_prefix, vector_noda_prefix, static_prefix, output_prefix, &
                snowUpdateOpt, PRINTRANK, print_debg_info, &   !fv3_index, vector_inputs, &
                SNOANL_out, &   !SNDFCS_out, SWEANL_out, & incr_at_Grid_out, 
                Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)
   !----------------------------------------------------------------------
   ! Input arguments: 
         
    !----------------------------------------------------------------------
    IMPLICIT NONE
    !
    include 'mpif.h'

    integer, parameter :: dp = kind(1.d0)
   ! , enkf_ensrf
    INTEGER, intent(in)    :: NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &  
                              IY, IM, ID, IH, LENSFC, IVEGSRC   !, num_assim_steps 
    Integer, intent(in)    :: ens_size
    LOGICAL, intent(in)    :: assim_SnowPack_obs, assim_SnowCov_obs
    Logical, intent(in)    :: rcov_localize, ens_inflate               ! localize R Covariance
    !CHARACTER(len=4), intent(in)   :: stn_var ! should probably be called control_var
    CHARACTER(LEN=*), Intent(In)   :: STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME
    CHARACTER(LEN=*), Intent(In)   :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH
    ! CHARACTER(LEN=*), Intent(In)   :: SFC_FORECAST_PREFIX, CURRENT_ANALYSIS_PREFIX
    CHARACTER(LEN=*), Intent(In)   :: vector_restart_prefix, vector_noda_prefix, &
                                      static_prefix, output_prefix

    REAL, intent(in)    :: PERCENT_OBS_WITHHELD   !, dT_Asssim
    Real, intent(In)    :: L_horz , h_ver, obs_tolerance
    Real, intent(In)    :: obs_srch_rad, bkgst_srch_rad, ims_max_ele        
    INTEGER, intent(in) :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels
    INTEGER, intent(in) :: snowUpdateOpt, PRINTRANK
    REAL, intent(out)   :: SNOANL_out(LENSFC)
    LOGICAL             :: vector_inputs, print_debg_info, fv3_index
        ! for mpi par
    INTEGER   :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
    Integer   :: LENSFC_proc, begloc, endloc
    LOGICAL             :: exclude_obs_at_grid

    CHARACTER(LEN=5)    :: TILE_NUM
    ! Character(LEN=3)    :: rank_str
    INTEGER	            :: IERR	
    REAL             :: SNDFCS(ens_size+1, LENSFC_proc), SWEFCS(ens_size+1, LENSFC_proc), &
                        SNDANL(ens_size+1, LENSFC_proc), SWEANL(ens_size+1, LENSFC_proc)
    REAL             :: SNDnoDA(LENSFC_proc), SWEnoDA(LENSFC_proc)       !, SNOANL_Cur(LENSFC)
    
    REAL   :: RLA(LENSFC_proc), RLO(LENSFC_proc), RLO_Tile(LENSFC_proc), OROG(LENSFC_proc)  !, OROG_UF(LENSFC)
    REAL            :: VETFCS(LENSFC_proc), SNUP_Array(LENSFC_proc)           ! SNOFCS(LENSFC), 
    Integer         :: tile_xy(LENSFC_proc), Idim_xy(LENSFC_proc), Jdim_xy(LENSFC_proc)
    INTEGER             :: LANDMASK(LENSFC_proc)  !, DAMASK(LENSFC)

    CHARACTER(len=250)   :: ghcnd_inp_file, ims_inp_file, ims_inp_file_indices
    CHARACTER(len=4)     :: y_str, m_str, d_Str, h_str, fvs_tile
    REAL, ALLOCATABLE    :: SNOOBS_stn(:), SNOFCS_at_stn(:), SNOANL_atStn(:)  !, SNOFCS_at_stn_ML(:)                
    REAL, ALLOCATABLE    :: Lat_stn(:), Lon_stn(:), OROG_at_stn(:)  
    REAL                 :: lat_min, lat_max, lon_min, lon_max      
    Real                 :: SNCOV_IMS(LENSFC_proc)  ! ims resampled at each grid
    Real                 :: SNO_IMS_at_Grid(ens_size+1, LENSFC_proc)
    Real                 :: IMS_Foot_Print(LENSFC_proc)  ! grid cells where IMS was assimilated

    INTEGER                :: num_stn, Num_Ims, num_Eval !num_subgrd_ims_cels !Num_Ims_Lat, Num_Ims_Lon

    Integer, Allocatable   :: index_at_nearStn(:)  !, index_at_nearIMS(:) !loc_near_Obs(:), 
    Integer                :: num_loc, num_loc_1, num_loc_2
    Real, Parameter        :: Stdev_back_depth = 30., Stdev_Obsv_depth = 40., Stdev_Obsv_ims = 80. ! mm 
    real                   :: stdev_obsv, stdev_back
    ! Integer                :: ims_assm_hour
	INTEGER                :: jndx, zndx, ncol, nrow

    !ens
    Real, ALLOCATABLE    :: SNOFCS_atObs_ens(:,:)  !SWEFCS_Ens(:,:), SNDFCS_Ens(:,:), , back_at_Obs_ens(:,:)

!1.4 track obserations assimilated
    Integer                :: index_obs_assmilated(max_num_nearStn+1, LENSFC_proc)

    Real, Allocatable        :: obs_Array(:), Lat_Obs(:), Lon_Obs(:), orogfcs_Obs(:) !back_at_Obs(:), 
    REAL                     :: incr_at_Grid(ens_size+1, LENSFC_proc) !, incr_at_Grid_ensM(LENSFC)    ! increment at grid
    Real, Allocatable        :: obs_Innov(:)  !, OmB_innov_at_stn(:)

    CHARACTER(len=250)       :: forc_inp_file, da_out_file, noda_inp_path
    CHARACTER(LEN=4)         :: RANKCH 
    CHARACTER(LEN=500)       :: static_filename

    ! REAL, ALLOCATABLE  :: SNOFCS_atEvalPts(:), incr_atEvalPts(:), SNOANL_atEvalPts(:) !, SNOANL_Cur_atEvalPts(:)  !evalution points 
    REAL, ALLOCATABLE    :: Lat_atEvalPts(:), Lon_atEvalPts(:), Obs_atEvalPts(:)     !evalution points
    ! REAL, ALLOCATABLE  :: Orog_obs_atEvalPts(:), Orog_fcs_atEvalPts(:)
    Integer, ALLOCATABLE    :: index_back_atEval(:)     ! background locations at eval points 
    Integer                 :: NEXC 
    Integer, ALLOCATABLE    :: index_back_atObs(:), Index_Obs_Excluded(:), &
                               index_obs_atEval(:)   ! the location of background corresponding obs
   
    Real               :: snodens, SNODENS_Grid(ens_size+1, LENSFC_proc)
    !LOGICAL            :: assim_snpack_stn, assim_SWE    !use swe as model state var? (instead of snow depth)
    LOGICAL            :: assim_sncov_thisGridCell    !    assim_sncov, assimilate sncov, 

    Integer            :: veg_type_landice  ! 10.21.20: no assmn over land ice

    INTEGER            :: send_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, &
                            dest_Aoffset_end, arLen, pindex
    INTEGER            :: mpiReal_size, rsize, isize, mpiInt_size, ixy, ie
    REAL               :: tmp
    INTEGER            :: istep, IY_loc, IM_loc, ID_loc, IH_loc
    REAL               :: IH_real

    integer            :: print_i
        ! INTEGER, PARAMETER :: PRINTRANK = 4 
    ! CSD-todo, should be same as in sfcsub. Share properly
    Real, parameter    :: nodata_val = -9999.9

    integer, parameter :: lsm_type = 2

!   Snow partion options for Noah-MP
    integer, parameter     :: DO_NOTHING = 0
    integer, parameter     :: TOP_LAYER = 1
    integer, parameter     :: BOTTOM_LAYER = 2
    integer, parameter     :: ALL_LAYERS = 3

    ! noah mp
    Real  :: SCF_Grid(ens_size+1, LENSFC_proc)
    !Real  :: SCF_Grid_Land(ens_size+1, LENSFC_proc)
    Real     :: noahmp_ensm_swe(LENSFC_proc), noahmp_ensm_snow_depth(LENSFC_proc)

    type(noahmp_type)      :: noahmp(ens_size)
    ! type(observation_type) :: obs    
!   noahmp%model%weasd  = this%snow_water_equivalent
!   noahmp%model%snwdph = this%snow_depth * 1000.0  ! driver wants mm
!   noahmp%model%canopy = this%canopy_water
!   noahmp%model%tskin  = this%skin_temperature
!   noahmp%model%stc    = this%soil_temperature
!   noahmp%model%smc    = this%soil_moisture
!   noahmp%model%slc    = this%soil_liquid
!   noahmp%model%tsurf  = noahmp%model%tskin

!   ! type noahmp_type
!     snwdph  ! depth 
!     weasd     ! SWE 
!     snowc   !fractional snow cover
!     sncovr1 ! snow cover over land
!     tsurf ! surface skin temperature
!     stc(time, soil_levels, location) !soil temperature
!     tsnoxy(time, snow_levels, location) ! snow temperature
    Do ie = 1, ens_size
        allocate(noahmp(ie)%swe                (LENSFC_proc))
        allocate(noahmp(ie)%snow_depth         (LENSFC_proc))
        allocate(noahmp(ie)%active_snow_layers (LENSFC_proc))
        allocate(noahmp(ie)%swe_previous       (LENSFC_proc))
        allocate(noahmp(ie)%snow_soil_interface(LENSFC_proc,7))
        allocate(noahmp(ie)%temperature_snow   (LENSFC_proc,3))
        allocate(noahmp(ie)%snow_ice_layer     (LENSFC_proc,3))
        allocate(noahmp(ie)%snow_liq_layer     (LENSFC_proc,3))
        allocate(noahmp(ie)%temperature_soil   (LENSFC_proc)) 
    End do       
    ! end type noahmp_type 

    !=============================================================================================
    ! 1. initialise vars,set-up processors, and read lat/lon from orog files.
    !=============================================================================================

    !initialse output with nodata
    IMS_Foot_Print = IEEE_VALUE(IMS_Foot_Print, IEEE_QUIET_NAN)
    stdev_obsv = stdev_obsv_depth
    stdev_back = stdev_back_depth  
    !obs_srch_rad = 250. ! radius of observation search
    ! ims_assm_hour = 18 
    ! noah models specific? Needed to ID glaciers.
    if (IVEGSRC == 2) then   ! sib
            veg_type_landice=13
    else
            veg_type_landice=15
    endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    !total number of processors used = Multiple of 6: any extra processors sent to end of subroutine
    IF (myrank ==PRINTRANK) PRINT*,"snowDA: EnSRF total num proc ", NPROCS, " Num tiles : ", NUM_TILES
  
! if (p_tN /= 4 ) goto 999

! must have all obs
    if ((STN_OBS_PREFIX(1:8).eq.'        ') .OR. (IMS_SNOWCOVER_PATH(1:8).eq.'        ') &
            .OR. (IMS_INDEXES_PATH(1:8).eq.'        ')) then
            print*, "One or more observation paths don't exist!, skipping the DA"
            goto 999
    end if

! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE p_tN
    write(fvs_tile, "(I0.2)") IDIM
    static_filename = trim(static_prefix)//"/ufs-land_C"//trim(fvs_tile)// &
                    "_static_fields.nc"
    Call ReadTileInfo(static_filename, LENSFC_proc, veg_type_landice, &
            mp_start, mp_end, &
            tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 
    ! make RLO copy before so that RLO (used later) is not modified
    RLO_Tile = RLO 
    Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
    ! 1 degree buffer around tile boundaries (but not to go out of valid coords)
    lat_min = MAX(MINVAL(RLA) - 1., -90.)
    lat_max = MIN(MAXVAL(RLA) + 1., 90.)
    lon_min = MAX(MINVAL(RLO_Tile) - 1., -180.)
    lon_max = MIN(MAXVAL(RLO_Tile) + 1., 180.)      
  
    if ((p_tRank==PRINTRANK) ) then !.and. print_deb) then
            print*, " min/max lat/lon ", lat_min, lat_max, lon_min, lon_max
    endif

! If multiple time steps are simulated in a given time period (window) given by num_assim_steps * dT_Assim
! Note: the DA outputs from the last time step are returned        
        IH_real = IH; IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH; ! these variables change inside loop below
    ! Do istep = 1, num_assim_steps
                
        write(y_str, "(I4)") IY_loc
        write(m_str, "(I0.2)") IM_loc
        write(d_str, "(I0.2)") ID_loc
        write(h_str, "(I0.2)") IH_loc

        ! controls calling of obs operator for stn data. If remains 0 will not be called.
        num_stn = 0 
        num_Eval = 0

!=============================================================================================
! 2. Read model forecast here, as need VETFCS and snow density for IMS snow depth calc. (later, separate read routines) 
!=============================================================================================

       ! READ THE INPUT SURFACE DATA from vector and no da outputs
        noda_inp_path=TRIM(vector_noda_prefix)//"/ufs_land_output."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
        Call ReadRestartNoahMP_Ens(myrank, LENSFC_proc, vector_restart_prefix, noda_inp_path, &
                        y_str, m_str, d_str, h_str, ens_size, mp_start, mp_end, &
                        noahmp, SNDnoDA, SWEnoDA, SCF_Grid)  !, SNDFCS, SWEFCS)
        ! initialise analysis to forecast (mostly to retain non-land snow states, which are not updated)
        SWEFCS(ens_size+1,:) = 0.0
        SNDFCS(ens_size+1,:) = 0.0
        ! SNODENS_Grid(ens_size+1,:) = 0.0
        SCF_Grid(ens_size+1,:) = 0.0
        !SCF_Grid_Land(ens_size+1,:) = 0.0
        Do ie = 1, ens_size
            SWEFCS(ie,:) = noahmp(ie)%swe(:)
            SNDFCS(ie,:) = noahmp(ie)%snow_depth(:)
            ! grid snow density
            call calc_density(LENSFC_proc, lsm_type, LANDMASK, SWEFCS(ie,:), SNDFCS(ie,:), &
                               noahmp(ie)%temperature_soil, SNODENS_Grid(ie,:))
            SWEFCS(ens_size+1,:) = SWEFCS(ens_size+1,:) + SWEFCS(ie,:)
            SNDFCS(ens_size+1,:) = SNDFCS(ens_size+1,:) + SNDFCS(ie,:)
            SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) + SCF_Grid(ie,:)
            !SCF_Grid_Land(ens_size+1,:) = SCF_Grid_Land(ens_size+1,:) + SCF_Grid_Land(ie,:)
        End do         
        ! SUM(SNDFCS(1:ens_size, jndx)) / ens_size
        SWEFCS(ens_size+1,:) = SWEFCS(ens_size+1,:) / ens_size
        SNDFCS(ens_size+1,:) = SNDFCS(ens_size+1,:) / ens_size
        SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) / ens_size
        !SCF_Grid_Land(ens_size+1,:) = SCF_Grid_Land(ens_size+1,:) / ens_size

        SWEANL(:,:) = SWEFCS(:,:)
        SNDANL(:,:) = SNDFCS(:,:)
        incr_at_Grid = 0.0

        tmp = SUM(SWEFCS(ens_size+1,:),  Mask = (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>=0 )) &          !
                         / COUNT (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>= 0)                 !
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' ensemble mean SWE', tmp
        tmp = SUM(SNDFCS(ens_size+1,:),  Mask = (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>=0 )) &              !
                         / COUNT (LANDMASK==1 .and. SNDFCS(ens_size+1,:)>= 0)     
        ! If (p_tRank==0)  
        print*, "proc ", myrank,  ' ensemble mean SND', tmp  
        SNODENS_Grid(ens_size+1,:) =  SWEFCS(ens_size+1,:)/SNDFCS(ens_size+1,:)

!=============================================================================================
! 3. Read observations
!=============================================================================================

! 3a. Read station obs (of snow depth or SWE)
        if (assim_SnowPack_obs) then 
            ghcnd_inp_file = TRIM(STN_OBS_PREFIX)// &   !GHCND.SNWD TRIM(h_str)//
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//".nc"
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading GHCN file', trim(ghcnd_inp_file) 
            Call Observation_Read_GHCND_IODA(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
                            TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
                            num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn, &
                            NEXC, Index_Obs_Excluded)
            ! Call Observation_Read_GHCND_All_excNaN(ghcnd_inp_file, TRIM(STN_DIM_NAME),  &
            !                 TRIM(STN_VAR_NAME), TRIM(STN_ELE_NAME), &
            !                 num_stn, SNOOBS_stn, OROG_at_stn, Lat_stn, Lon_stn)

            if (myrank==PRINTRANK) then   !.and. (p_tN==PRINTRANK).and. print_deb) then
                    print*, "Tile ", p_tN, " num. Stn obs ", num_stn
            endif
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                    PRINT*, "Stn SND from rank: ", MYRANK
                    PRINT*, SNOOBS_stn
                    PRINT*, "Lat at Stn from rank: ", MYRANK
                    PRINT*, Lat_stn
                    PRINT*, "Lon at Stn from rank: ", MYRANK
                    PRINT*, Lon_stn
                    PRINT*, "Elevation at station locations from rank: ", MYRANK
                    PRINT*, OROG_at_stn
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished reading station data'
        endif
! CSD beyond this point, there should be no specific mention of the station data source

! 3b. Read remotely sensed snow cover, and convert to  snow depth or SWE. 
        if (assim_SnowCov_obs) then
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS/"
            ims_inp_file = TRIM(IMS_SNOWCOVER_PATH)//"IMS.SNCOV."// &
                            TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//TRIM(h_str)//".nc"                      !
            ! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/IMS_INDEXES/"
            WRITE(RANKCH, '(I0.1)') myrank                            
            ims_inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
                                        ".IMS.Indices."//trim(ADJUSTL(RANKCH))//".nc"                       
            if (MYRANK==PRINTRANK) PRINT *, 'snowDA: reading IMS file', trim(ims_inp_file) 
            ! Call Observation_Read_IMS_Full(ims_inp_file, ims_inp_file_indices, &
            !                         MYRANK, JDIM, IDIM, num_subgrd_ims_cels, SNCOV_IMS)
            Call Read_IMS_and_Resample_toTarget(ims_inp_file, IMS_INDEXES_PATH, & !fv3_index, 
                    LENSFC_proc, num_subgrd_ims_cels, IDIM, JDIM, &
                    tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
            if(myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV from rank: ", MYRANK
                PRINT*, SNCOV_IMS
            endif 
            if (myrank==PRINTRANK) PRINT*,'Finished reading SNCOV, converting to snow depth' 
           ! SNUP array will be used later, to test whether SCF > 100%
            ! call CalcSWEFromSnowCover(SNCOV_IMS, VETFCS, LENSFC_proc, &
            !                            SNO_IMS_at_Grid, SNUP_Array)
            ! SNO_IMS_at_Grid = SNO_IMS_at_Grid/SNODENS_Grid(ens_size+1,:) ! convert SWE to SND
            SNO_IMS_at_Grid(ens_size+1,:) = 0.0
            Do ie = 1, ens_size
                call calcSD_noahmp(LENSFC_proc, SNCOV_IMS, VETFCS, &
                                   SNODENS_Grid(ie,:), SNO_IMS_at_Grid(ie,:))
                SNO_IMS_at_Grid(ens_size+1,:) = SNO_IMS_at_Grid(ens_size+1,:) + &
                                                SNO_IMS_at_Grid(ie,:)
            End do
            SNO_IMS_at_Grid(ens_size+1,:) = SNO_IMS_at_Grid(ens_size+1,:) / ens_size
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "SNCOV derived snodepth at each grid cell from rank: ", MYRANK
                PRINT*, SNO_IMS_at_Grid
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished converting SNCOV observations'
        
        endif ! read_IMS 

!=============================================================================================
! 4. Get H(x): Read forecast snow fields from restart files, then interpolate to obs location.
!=============================================================================================

! 4a. read the forecast file on model grid : this was done earlier, as need veg type and 
!     snow density for IMS snow depth conversion

! 4b. get H(x) for station obs
        ! Get model states at obs points
        if (num_stn > 0) then ! skip if not reading in station data / no obs were available
            ALLOCATE(SNOFCS_at_stn(num_stn))
            ALLOCATE(SNOANL_atStn(num_stn))
            ! Allocate(SNOFCS_at_stn_ML(num_stn))
            ! ALLOCATE(OROGFCS_at_stn(num_stn)) 
            ALLOCATE(index_back_atObs(num_stn)) 
        
	! ALLOCATE(OmB_innov_at_stn_ens(num_stn, ens_size)) 
        ! using PERCENT_OBS_WITHHELD % of stn locations for evaluation
            num_Eval = floor(0.01 * PERCENT_OBS_WITHHELD * num_stn)  
            if (num_Eval > 0) then 
                ALLOCATE(index_back_atEval(num_Eval)) 
                ALLOCATE(Obs_atEvalPts(num_Eval)) 
                ! ALLOCATE(SNOFCS_atEvalPts(num_Eval)) 
                ALLOCATE(Lat_atEvalPts(num_Eval))
                ALLOCATE(Lon_atEvalPts(num_Eval))            
                ALLOCATE(index_obs_atEval(num_Eval)) 
                
                ! if(p_tRank == PRINTRANK) then 
                        PRINT*,"Proc ", myrank," ",num_Eval,' points for evaluation excluded from DA'       
                ! endif   
            endif 
! CSD todo: separate out eval call and obs operator call. model info (SNOOBS_stn shouldn't be in the obs operator)
!           for JEDI, probably also want to add snow depth derived from snow cover here.
            ! Call Observation_Operator_tiles_Parallel(Myrank, NUM_TILES, p_tN, p_tRank, Np_til, & 
            !                     RLA, RLO, OROG, Lat_stn, Lon_stn,   &
            !                     LENSFC, num_stn, num_Eval, bkgst_srch_rad, SNDFCS, SNOOBS_stn, LANDMASK,  &
            !                     SNOFCS_at_stn, OROGFCS_at_stn, index_back_atObs, index_back_atEval,  &
            !                     Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
            Call Observation_Operator_Parallel_vect(Myrank, NPROCS, LENSFC, LENSFC_proc, &
                 num_stn, num_Eval, RLA, RLO, Lat_stn, Lon_stn, OROG_at_stn, bkgst_srch_rad, &
                 SNDFCS(ens_size+1,:), SNOOBS_stn, LANDMASK, SNOFCS_at_stn, &
                 index_back_atObs, index_back_atEval, index_obs_atEval, &
                 Obs_atEvalPts, Lat_atEvalPts, Lon_atEvalPts) !SNOFCS_atEvalPts, , Orog_obs_atEvalPts)
            if (myrank==PRINTRANK) PRINT*,'Finished stn observation operator'

            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then 
                PRINT*, "Background Indices at eval points"
                PRINT*, index_back_atEval       
                PRINT*, "Obs at Eval Points" 
                PRINT*, Obs_atEvalPts   
                ! PRINT*, "Forecast at Eval Points"
                ! PRINT*, SNOFCS_atEvalPts                             
                PRINT*, "Lat at Eval Points"
                PRINT*, Lat_atEvalPts
                PRINT*, "Lon at Eval Points"
                PRINT*, Lon_atEvalPts
            endif
            ! OmB_innov_at_stn = SNOOBS_stn - SNOFCS_at_stn
            ALLOCATE(SNOFCS_atObs_ens(ens_size, num_stn))
            SNOFCS_atObs_ens = IEEE_VALUE(SNOFCS_atObs_ens, IEEE_QUIET_NAN)
            Do jndx = 1, num_stn
                if (index_back_atObs(jndx) > 0) then
                    SNOFCS_atObs_ens(:, jndx) = SNDFCS(1:ens_size, index_back_atObs(jndx))
                endif
            end do
            if (myrank==PRINTRANK .and. print_debg_info .and. print_deb) then
                PRINT*, "station Lat range from rank: ", MYRANK, MINVAL(Lat_stn), " ", MAXVAL(Lat_stn)
                PRINT*, "station Lon range from rank: ", MYRANK, MINVAL(Lon_stn), " ", MAXVAL(Lon_stn)
            !     PRINT*, "Lat at Obs Points"
            !     PRINT*, Lat_stn
            !     PRINT*, "Lon at Obs Points"
            !     PRINT*, Lon_stn
                    ! PRINT*, "Model elevation at station locations from rank: ", MYRANK
                    ! PRINT*, OROGFCS_at_stn
                PRINT*, "Background Indices at obs points"
                PRINT*, index_back_atObs
                 PRINT*, "Background at station locations from rank: ", MYRANK
                PRINT*, SNOFCS_at_stn  
            !     PRINT*, "Obs at obs stns" 
            !     PRINT*, SNOOBS_stn   
                ! PRINT*, "O - B (innovation at obs points)"
                ! PRINT*, OmB_innov_at_stn 
            endif
            if (myrank==PRINTRANK) PRINT*,'Finished observation operator for station data'         
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

        !Do nrow = 1, Num_Snotel !JDIM/2
        !obs_srch_rad = 250. ! radius of observation search
        if (myrank==PRINTRANK) PRINT*,'Starting DA loop'
        print_i = 99
        Do jndx = 1, LENSFC_proc    !mp_start, mp_end     
            ! if (myrank == 4 .and. (jndx == 3943 .or. jndx == 4033) ) then  
            ! QC: only update this grid cell if it is land.
            if( LANDMASK(jndx) == 1 ) then 
                num_loc_1 = 0
                num_loc_2 = 0
                assim_sncov_thisGridCell = .FALSE.
                if (print_debg_info) print*, "proc ", myrank, " grid: ", jndx
                if (num_stn>0) then 
                ! CSD - speed up by skipping over all grid cells where model and IMS agree on no snow? 
                ! currently: find station obs in radius, do gross error check, and limit to 50 obs
                ! QC: gross error check is done in this call.
                ! 8.19.20: we are assuming here if deterministic forecast exists/not null
	            ! at a point, then ensemble members also have valid value

! 10.25.21 Note the QC based on forecast ensemble mean
! may need to use a different (perhaps larger) obs_tolerance
                    call nearest_Observations_Locations(RLA(jndx), RLO(jndx),          &
                            Lat_stn, Lon_stn, num_stn, obs_srch_rad, max_num_nearStn,  &
                            stdev_back, stdev_obsv, obs_tolerance,             &
                            SNOFCS_at_stn, SNOOBS_stn,                         &
                            index_at_nearStn,  num_loc_1) !,      &LENSFC,
                    if (print_debg_info) print*, "number of stn sndpth obs ", num_loc_1
                endif
! 10.25.21 Note the QC based on forecast ensemble mean  
                if( assim_SnowCov_obs .AND. &    !.and. (IH_loc == ims_assm_hour) 
                    (.NOT. IEEE_IS_NAN(SNDFCS(ens_size+1,jndx))) .AND. &
                    (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(ens_size+1, jndx))) .AND. &
                !   ( .not.(SWEFCS(ens_size+1,jndx) >= SNUP_Array(jndx) .AND. & 
                    ( .not.(SCF_Grid(ens_size+1,jndx) >= 0.99 .AND. & 
                            SNCOV_IMS(jndx) >= 0.99) )  &
! 10.25.21 CSD decided not to use this criteria anymore
                    !(OROG(jndx) <= ims_max_ele) .AND. &
                    ) then
                        num_loc_2 = 1 
                        assim_sncov_thisGridCell = .TRUE.                
                endif
                num_loc = num_loc_1 + num_loc_2                
                ! if assim_sncov=false >> num_loc_1=num_loc
                ! QC: only update this grid cell if it is land.	
                if( num_loc > 0 ) then            !((num_loc > 0) .and. ( LANDMASK(jndx) == 1 )) then 
                    ! Allocate(back_at_Obs(num_loc))
                    Allocate(obs_Array(num_loc))
                    Allocate(Lat_Obs(num_loc))
                    Allocate(Lon_Obs(num_loc))
                    Allocate(orogfcs_Obs(num_loc))  
                    ! ghcnd
                    if(num_loc_1 > 0) then
                        index_obs_assmilated(1, jndx) = num_loc_1
                        Do zndx = 1, num_loc_1    
                            index_obs_assmilated(zndx+1, jndx) = index_at_nearStn(zndx) 
                            ! back_at_Obs(zndx) = SNOFCS_at_stn(index_at_nearStn(zndx))
                            obs_Array(zndx) = SNOOBS_stn(index_at_nearStn(zndx))
                            Lat_Obs(zndx) = Lat_stn(index_at_nearStn(zndx))
                            Lon_Obs(zndx) = Lon_stn(index_at_nearStn(zndx))
                            orogfcs_Obs(zndx) = OROG_at_stn(index_at_nearStn(zndx)) 
                        End Do
                    End if
                    ! Append IMS-derived snow depth to the obs array 
                    if(assim_sncov_thisGridCell) then   
                        IMS_Foot_Print(jndx) = 1.0                           
                        ! back_at_Obs(num_loc) = SNDFCS(jndx)
                        obs_Array(num_loc) = SNO_IMS_at_Grid(ens_size+1, jndx)
                        Lat_Obs(num_loc) = RLA(jndx)   
                        Lon_Obs(num_loc) = RLO(jndx) 
                        orogfcs_Obs(num_loc) = OROG(jndx)
                    endif	
                    ! EnSRF
                    if(exclude_obs_at_grid .and. (num_loc_1 > 1)) then
                        Allocate(obs_Innov(num_loc-1))
                        Call snow_DA_EnSRF( RLA(jndx), RLO(jndx), OROG(jndx),     &
                            num_loc_1-1, num_loc-1, num_stn, jndx, ens_size,   &
                            L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                            Lat_Obs(2:num_loc), Lon_Obs(2:num_loc), orogfcs_Obs(2:num_loc),     		    &
                            assim_sncov_thisGridCell, rcov_localize, ens_inflate,           &
                            SNDFCS(1:ens_size, jndx),   & !LENSFC, 
                            index_at_nearStn(2:num_loc), SNOFCS_atObs_ens,	 &
                            stdev_obsv, Stdev_Obsv_ims, stdev_back,         &
                            obs_Array(2:num_loc),                          &
                            obs_Innov, incr_at_Grid(:, jndx), SNDANL(:, jndx))         
                        ! print_i = print_i + 1
                    else
                        Allocate(obs_Innov(num_loc))
                        Call snow_DA_EnSRF( RLA(jndx), RLO(jndx), OROG(jndx),     &
                            num_loc_1, num_loc, num_stn, jndx, ens_size,   &
                            L_horz, h_ver,                                   &   !L_horz in Km, h_ver in m
                            Lat_Obs, Lon_Obs, orogfcs_Obs,     		    &
                            assim_sncov_thisGridCell, rcov_localize, ens_inflate,           &
                            SNDFCS(1:ens_size, jndx),   & !LENSFC, 
                            index_at_nearStn, SNOFCS_atObs_ens,	 &
                            stdev_obsv, Stdev_Obsv_ims, stdev_back,         &
                            obs_Array,                          &
                            obs_Innov, incr_at_Grid(:, jndx), SNDANL(:, jndx))         
                        ! print_i = print_i + 1
                    endif
                    if (myrank==PRINTRANK .and. print_debg_info) then  !
                        print*, "proc ", myrank, "loop ", jndx, "num depth obs ", num_loc_1
                        print*, "total obs", num_loc, " ens size ", ens_size
                        PRINT*, "Ens background: "
                        PRINT*, SNDFCS(:, jndx)	
                        PRINT*, "Observed"
                        PRINT*,  obs_Array
                        PRINT*, "Obs innovation: "
                        PRINT*, obs_Innov
                        PRINT*, "Ens increment at grid: "
                        PRINT*, incr_at_Grid(:, jndx)
                        PRINT*, "Ens analysis at grid: "
                        PRINT*, SNDANL(:, jndx)
                    endif					
                    DEALLOCATE(obs_Array, obs_Innov)
                    DEALLOCATE(Lat_Obs, Lon_Obs, orogfcs_Obs)
                ! else
                !     SNDANL(:, jndx) = SNDFCS(:,jndx) !anl_at_Grid_ens(1:ens_size, jndx) = SNDFCS_Ens(:, jndx)
                !     SNDANL(ens_size+1, jndx) = SUM(SNDFCS(1:ens_size, jndx)) / ens_size
                endif
                if (allocated(index_at_nearStn))  Deallocate(index_at_nearStn)  
            endif ! not a land cell  
	    End do
        !End do
        if (myrank==PRINTRANK) PRINT*, 'Finished DA loops'
	! PRINT*, 'Finished DA loops', ' proc:', myrank
!=============================================================================================
! 7. Clean up, write outputs 
!=============================================================================================
! collect results onto main tasks, if necessary
! fill in SWE and SND arrays       
        ! avoid -ve anl
        Where(SNDANL < 0.) SNDANL = 0.
        if (print_debg_info) then
            PRINT*, "Weighted increment SWE/snwd from rank: ", MYRANK
            PRINT*, incr_at_Grid       
            PRINT*, "Analysis SWE/ snwd  from rank: ", MYRANK
            PRINT*, SNDANL
        endif
        ! d = swe/sndi
        ! this is for writing outputs at observation and evaluation points
        ! SWEANL = SNDANL * SNODENS_Grid
        Do ie = 1, ens_size+1
            WHERE (LANDMASK==1) SWEANL(ie,:) = SNDANL(ie,:) * SNODENS_Grid(ie,:)
        Enddo
        
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
! outputs
        if (MYRANK > 0) then  
            call MPI_SEND(SNDANL(ens_size+1, :), LENSFC_proc, mpiReal_size, 0, &
                            500*(MYRANK+1), MPI_COMM_WORLD, IERR) 
        endif
        if (MYRANK == 0) then
            SNOANL_out(1:LENSFC_proc) = SNDANL(ens_size+1, :)         
            Do ixy =  1, NPROCS - 1  ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1 
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(SNOANL_out(dest_Aoffset:dest_Aoffset_end),arLen, mpiReal_size, &
                            ixy, 500*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        endif
        call MPI_BCAST(SNOANL_out, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        if (myrank==PRINTRANK) PRINT*,'Finished Data copy'
        ! PRINT*, 'Finished Data copy', ' proc:', myrank
        ! if (MYRANK > NUM_TILES - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998
        noahmp_ensm_swe = 0.0
        noahmp_ensm_snow_depth = 0.0
        Do ie = 1, ens_size             
            select case (snowUpdateOpt)
            case (DO_NOTHING)
                ! write(*,*) " proc ", myrank, " Option 0: do nothing"
                noahmp(ie)%swe = SWEANL(ie,:)
                noahmp(ie)%snow_depth = SNDANL(ie,:)
            case (TOP_LAYER)
                ! write(*,*) " proc ",myrank, " Option 1: UpdateTopLayer"
                call UpdateTopLayer(myrank, LENSFC_proc, noahmp(ie), &
                incr_at_Grid(ie,:), SNDANL(ie,:), LANDMASK)
            case (BOTTOM_LAYER)
                ! write(*,*) " proc ",myrank, " Option 2: UpdateBottomLayer"
                call UpdateBottomLayer(myrank, LENSFC_proc, noahmp(ie), &
                incr_at_Grid(ie,:), SNDANL(ie,:), LANDMASK)
            case (ALL_LAYERS)
                ! write(*,*) " proc ",myrank, " Option 3: UpdateAllLayers"
                ! call UpdateAllLayers_Ens(ie, myrank, LENSFC_proc, noahmp(ie), &
                ! incr_at_Grid(ie,:), SNDANL(ie,:), SNODENS_Grid(ie,:), LANDMASK)
                call UpdateAllLayers(myrank, LENSFC_proc, noahmp(ie), incr_at_Grid(ie,:), &
                                SNDANL(ie,:), SNODENS_Grid(ie,:), LANDMASK)
            case default
                write(*,*) "choose a valid partition option"
                stop
            end select 
            ! avoid -ve anl
            ! print*, " proc ", myrank, "done partitioning ens ", ie
            Where(noahmp(ie)%swe < 0) noahmp(ie)%swe = 0.
            Where(noahmp(ie)%snow_depth < 0) noahmp(ie)%snow_depth = 0.
            ! ens mean
            noahmp_ensm_swe = noahmp_ensm_swe + noahmp(ie)%swe
            noahmp_ensm_snow_depth = noahmp_ensm_snow_depth + noahmp(ie)%snow_depth
            ! print*, " proc ", myrank, "done loop ", ie
        enddo
        ! SUM(SNDFCS(1:ens_size, jndx)) / ens_size
        noahmp_ensm_swe = noahmp_ensm_swe / ens_size
        noahmp_ensm_snow_depth = noahmp_ensm_snow_depth / ens_size
        if (print_debg_info .and. myrank==PRINTRANK) then
            print*, "noah SWE ", noahmp_ensm_swe
            print*, "noah snowdepth", noahmp_ensm_snow_depth
        endif
! !Compute updated snocov 
        !Call update_snow_cover_fraction(LENSFC, SNOANL, VETFCS, anl_fSCA)
        ! update ensemble scf
        SNODENS_Grid(ens_size+1,:) =  noahmp_ensm_swe/noahmp_ensm_snow_depth
        SCF_Grid(ens_size+1,:) = 0.0
        Do ie = 1, ens_size 
            ! grid snow density
            ! call calc_density(LENSFC_proc, lsm_type, LANDMASK, noahmp(ie)%swe, &
            ! noahmp(ie)%snow_depth, noahmp(ie)%temperature_soil, SNODENS_Grid(ie,:))             
            call calcSCF_noahmp(LENSFC_proc, VETFCS, SNODENS_Grid(ie,:), &
                noahmp(ie)%snow_depth, SCF_Grid(ie, :)) 
            SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) + SCF_Grid(ie,:)
! update  sncovr1 "snow cover over land" 
        Enddo 
        SCF_Grid(ens_size+1,:) = SCF_Grid(ens_size+1,:) / ens_size
        
        ! write outputs	
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
        if (myrank==PRINTRANK) then 
            PRINT*,"proc ", myrank, "starting writing output"
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        da_out_file = trim(output_prefix)// &
                      TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"_18.nc"  !     
        call Write_DA_Outputs_ens_vect(da_out_file, vector_restart_prefix, &
                y_str, m_str, d_str, h_str, &        
                MYRANK, NPROCS, LENSFC, LENSFC_proc, ens_size, num_stn, num_Eval, &
                max_num_nearStn,  &
                N_sA, N_sA_Ext, mpiReal_size, mpiInt_size,     &
                SNDnoDA, SWEnoDA, noahmp,   &
                SNDFCS,  SWEFCS, SNDANL, SWEANL, incr_at_Grid,   &  ! 
                index_obs_assmilated,   &
                SNCOV_IMS, IMS_Foot_Print, SNO_IMS_at_Grid(ens_size+1, :),    &
                SNODENS_Grid, SCF_Grid, &  ! SAVe this as snowc
                noahmp_ensm_swe, noahmp_ensm_snow_depth, &       
                Lat_stn, Lon_stn, OROG_at_stn, SNOOBS_stn, index_back_atObs, &
                NEXC, Index_Obs_Excluded,  &
                Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, &
                index_back_atEval, index_obs_atEval)            
	! if (myrank==PRINTRANK) 
        PRINT*,"proc ", myrank, "finished writing output"
    ! endif
    998 CONTINUE
        ! clean up      
        Do ie = 1, ens_size
            if (allocated(noahmp(ie)%swe))   deallocate(noahmp(ie)%swe)
            if (allocated(noahmp(ie)%snow_depth))  deallocate(noahmp(ie)%snow_depth)
            if (allocated(noahmp(ie)%active_snow_layers)) deallocate(noahmp(ie)%active_snow_layers)
            if (allocated(noahmp(ie)%swe_previous))  deallocate(noahmp(ie)%swe_previous)
            if (allocated(noahmp(ie)%snow_soil_interface))  deallocate(noahmp(ie)%snow_soil_interface)
            if (allocated(noahmp(ie)%temperature_snow))  deallocate(noahmp(ie)%temperature_snow)
            if (allocated(noahmp(ie)%snow_ice_layer))  deallocate(noahmp(ie)%snow_ice_layer)
            if (allocated(noahmp(ie)%snow_liq_layer))  deallocate(noahmp(ie)%snow_liq_layer)
            if (allocated(noahmp(ie)%temperature_soil))  deallocate(noahmp(ie)%temperature_soil)
        End do 
        if (allocated(SNOFCS_at_stn))   DEALLOCATE(SNOFCS_at_stn)	
        if (allocated(SNOANL_atStn))   DEALLOCATE(SNOANL_atStn)
        if (allocated(Lat_stn))         DEALLOCATE(Lat_stn) 
        if (allocated(Lon_stn))         DEALLOCATE(Lon_stn) 
        if (allocated(SNOOBS_stn))      DEALLOCATE(SNOOBS_stn)
        if (allocated(index_back_atObs)) DEALLOCATE(index_back_atObs)

        if (allocated(Index_Obs_Excluded)) DEALLOCATE(Index_Obs_Excluded) 

        if (allocated(OROG_at_stn))  DEALLOCATE(OROG_at_stn)        
        if (allocated(Lat_atEvalPts))   DEALLOCATE(Lat_atEvalPts) 
        if (allocated(Lon_atEvalPts))   DEALLOCATE(Lon_atEvalPts)
        if (allocated(Obs_atEvalPts))   DEALLOCATE(Obs_atEvalPts)
        if (allocated(index_back_atEval)) DEALLOCATE(index_back_atEval)
        if (allocated(index_obs_atEval))   DEALLOCATE(index_obs_atEval)

        ! Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
        ! IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        if (myrank==PRINTRANK) PRINT*,'Finished EnSRF DA at datetime: ', y_str, m_str, d_str, h_str
    ! End do

999 CONTINUE
    !PRINT*,'Finished OI DA ON RANK: ', MYRANK
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    !STOP

    RETURN

END subroutine EnSRF_Snow_Analysis_NOAHMP


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
 
 
 subroutine map_outputs_toObs(MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, & 
	                          LENSFC, CURRENT_ANALYSIS_PATH)
							
	IMPLICIT NONE
	!
	include 'mpif.h'
	
	integer, parameter :: dp = kind(1.d0)

	INTEGER, intent(in) :: MAX_TASKS, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, LENSFC
	CHARACTER(LEN=*), Intent(In)   :: CURRENT_ANALYSIS_PATH
	CHARACTER(LEN=5)    :: TILE_NUM
	Character(LEN=3)    :: rank_str
	INTEGER			    :: IERR	
	REAL        :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC) !, OROG_UF(LENSFC)
	REAL                 :: SNOANL(LENSFC), SWDANL(LENSFC)	
	CHARACTER(len=250)   :: dim_name, anl_inp_path, da_out_file
	CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str, fvs_tile
	REAL, ALLOCATABLE    :: Lat_GHCND(:), Lon_GHCND(:), Lon_Obs(:)
	INTEGER              :: num_Eval, jndx
	Integer, Allocatable   :: index_back_atEval(:)
	REAL, ALLOCATABLE  :: SNOANL_atEvalPts(:)  !evalution points	
    ! for mpi par
	INTEGER            :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end

	IF (myrank ==0) PRINT*,"Total num proc ", NPROCS, " Num tiles /Max.tasks: ", MAX_TASKS
	Np_ext = MOD(NPROCS, MAX_TASKS)  ! extra/inactive procs
	if (MYRANK >  NPROCS - Np_ext - 1) goto 999
	Np_til = NPROCS / MAX_TASKS  ! num proc. per tile 
	p_tN = MOD(MYRANK, MAX_TASKS)  ! tile for proc.
	p_tRank = MYRANK / MAX_TASKS  ! proc. rank within tile
	N_sA = LENSFC / Np_til  ! sub array length per proc
	N_sA_Ext = LENSFC - N_sA * Np_til ! extra grid cells
	if(p_tRank == 0) then 
		mp_start = 1
	else
		mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
	endif
	mp_end = (p_tRank + 1) * N_sA + N_sA_Ext 		! end index of subarray for proc	
	If (myrank == 0 ) PRINT*,"sub array length ", N_sA, " extra sub array: ", N_sA_Ext
	! if (p_tN /= 2 ) goto 999
	
    ! READ THE OROGRAPHY AND GRID POINT LAT/LONS FOR THE CUBED-SPHERE TILE.
	CALL READ_LAT_LON_OROG_atRank(p_tN, RLA,RLO,OROG, TILE_NUM,IDIM,JDIM,LENSFC)	
	PRINT*,"Snow anl on ", MYRANK, " Tile group: ", p_tN, " Tile: ", TILE_NUM
 
	write(y_str, "(I4)") IY
	write(m_str, "(I0.2)") IM
	write(d_str, "(I0.2)") ID
	write(h_str, "(I0.2)") IH
	write(fvs_tile, "(I3)") IDIM
	!"/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Cur_Analysis/"
	anl_inp_path = TRIM(CURRENT_ANALYSIS_PATH)//"snow."// &
				TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
				TRIM(h_str)// "0000.sfcanl_data."//TILE_NUM//".nc"
	Call READ_Analysis_Data(anl_inp_path, LENSFC, SNOANL, SWDANL)
	Write(rank_str, '(I3.3)') (p_tN+1)	
	! "/scratch2/BMC/gsienkf/Tseganeh.Gichamo/SnowObs/Analysis2/
	da_out_file = "./SNOANL."// &  !
				  TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"18_tile"//rank_str//".nc"						
	Call Observation_Read_atEval(da_out_file, num_Eval, Lat_GHCND, Lon_GHCND, p_tN)					
	if ((p_tRank==0) ) then
		print*, "Tile ", p_tN, " num. Eval ", num_Eval
		! PRINT*, "Lat at Eval from rank: ", MYRANK
		! PRINT*, Lat_GHCND
		! PRINT*, "Lon at Eval from rank: ", MYRANK
		! PRINT*, Lon_GHCND
	endif
	ALLOCATE(index_back_atEval(num_Eval)) 
	ALLOCATE(SNOANL_atEvalPts(num_Eval)) 
	ALLOCATE(Lon_Obs(num_Eval)) 

    Lon_Obs = Lon_GHCND
    Where(Lon_Obs < 0) Lon_Obs = 360. + Lon_Obs
	! grid indices from coarse (eg C384) to high res (C768)
	Call Find_Nearest_GridIndices_Parallel(MYRANK, MAX_TASKS, p_tN, p_tRank, Np_til, &
					LENSFC, num_Eval, RLA, RLO, Lat_GHCND, Lon_Obs, index_back_atEval)	
	if ( (p_tRank==0) ) then 
		PRINT*, "Background Indices at eval points rank ",MYRANK
		PRINT*, index_back_atEval	
	endif

	if (MYRANK > MAX_TASKS - 1 ) goto 998   ! if(p_tRank /= 0 ) goto 998

	SNOANL_atEvalPts = IEEE_VALUE(SNOANL_atEvalPts, IEEE_QUIET_NAN)
	Do jndx = 1, num_Eval
		if (index_back_atEval(jndx) > 0) then
			SNOANL_atEvalPts(jndx) = SWDANL(index_back_atEval(jndx))
		endif
	End do
	dim_name = 'SNOANL_Cur_atEvalPts'
	call Append_DA_Outputs(da_out_file, MYRANK, num_Eval, dim_name, SNOANL_atEvalPts)  

998 CONTINUE
	DEALLOCATE(Lat_GHCND, Lon_GHCND, Lon_Obs, SNOANL_atEvalPts, index_back_atEval) 
999 CONTINUE
    PRINT*,'Finished ON RANK: ', MYRANK
	CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	STOP

	RETURN

End subroutine map_outputs_toObs

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
        If (p_tRank==0)  print*, "Tile ", p_tN, ' mean snow density', snodens
        tmp = SUM(SWEFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SWE', tmp
        tmp = SUM(SNDFCS,  Mask = (LANDMASK==1 .and. SNDFCS>0.01 )) &
                         / COUNT (LANDMASK==1 .and. SNDFCS> 0.01)
        If (p_tRank==0)  print*, "Tile ", p_tN,  ' mean SND', tmp

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

 ! the following code based on write_data() in read_write_data.f90
 Subroutine Write_DA_Outputs_ens(output_file, idim, jdim, lensfc, myrank,   &
                snwdforc, snwdanal, snoforc, snoanl, SNDANL_Cur, SWEANL_Cur, landmask,  &
                num_stn, Lat_atObs, Lon_atObs, OBS_stn, index_back_atObs,  &  ! FCS_at_stn,  SNOFCS_atObs_ens, OmB_innov_at_stn,  
                num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, index_back_atEval, &  !! SNOFCS_atEvalPts 
                SNCOV_IMS, IMS_Foot_Print, &   !SNOFCS_atEvalPts, incr_atEvalPts, SNOANL_atEvalPts, SNOANL_Cur_atEvalPts
                ens_size, SNDFCS_Ens, SWEFCS_Ens, incr_atGrid, anl_at_Grid_ens)  !, anl_fSCA) !updated snocov
        !------------------------------------------------------------------
        ! Write DA ouputs: 
        ! forecast SWE
        ! analysis SWE
        ! analysis Snow Depth
        ! innovation at grid
        !------------------------------------------------------------------
        implicit none

        CHARACTER(LEN=*), Intent(In)      :: output_file
        integer, intent(in)         :: idim, jdim, lensfc, num_stn, num_Eval, ens_size
        real, intent(in)            :: snoforc(lensfc), snwdforc(lensfc), snoanl(lensfc), snwdanal(lensfc)
        real, intent(in)            :: SWEANL_Cur(lensfc), SNDANL_Cur(lensfc)
        integer, intent(in)         :: landmask(lensfc)
        ! Real, intent(in)            ::   !, FCS_at_stn(num_stn), OmB_innov_at_stn(num_stn)
        Real, intent(in)            :: Lat_atObs(num_stn), Lon_atObs(num_stn), OBS_stn(num_stn)
        integer, intent(in)         :: index_back_atObs(num_stn), index_back_atEval(num_Eval)
        real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), Obs_atEvalPts(num_Eval)
        Real, intent(in)    :: SNCOV_IMS(lensfc), IMS_Foot_Print(lensfc)  !, anl_fSCA(lensfc) 
        Real, intent(in)    :: SNDFCS_Ens(ens_size, lensfc), SWEFCS_Ens(ens_size, lensfc)
        Real, intent(in)    :: incr_atGrid(lensfc), anl_at_Grid_ens(ens_size+1, lensfc)
        ! real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), SNOANL_atEvalPts(num_Eval) 
        ! real, intent(in)    :: incr_atEvalPts(num_Eval), SNOANL_Cur_atEvalPts(num_Eval)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_3dens(3), dims_strt(3), dims_end(3)
        integer                     :: error, i, ncid, ens_indx
        integer                     :: dim_x, dim_y, dim_time, dim_ens, dim_eval, dim_stn
        integer                     :: id_x, id_y, id_time, id_ens
        integer       :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_incr, id_imscov, id_imsftp   !, id_anlscov
        integer       :: id_swe_cur, id_snwd_cur, id_sndf_ens, id_swef_ens, id_snda_ens
        integer       :: id_latstn, id_lonstn, id_obsstn, id_landmask, id_indxobs!id_forcstn, id_innovstn, 
        integer       :: id_lateval, id_loneval, id_obseval, id_indxeval !id_forceval, id_anleval, id_increval, id_cur_anleval  !, , id_anlscov
        
        integer                     :: myrank 

        real(kind=4)                :: times
        real(kind=4), allocatable   :: x_data(:), y_data(:)
        integer, allocatable        :: ens_data(:)
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
        error = nf90_def_dim(ncid, 'Ensemble', ens_size, dim_ens)
        call netcdf_err(error, 'DEFINING Ensemble DIMENSION' )
        ! obs points
        error = nf90_def_dim(ncid, 'obs_points', num_stn, dim_stn)
        call netcdf_err(error, 'DEFINING obs_points' )
        ! eval obs points        
        if (num_Eval>0) then
                error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        else
                error = nf90_def_dim(ncid, 'eval_points', 1, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        endif

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

        error = nf90_def_var(ncid, 'Ensemble', NF90_INT, dim_ens, id_ens)
        call netcdf_err(error, 'DEFINING Ensemble FIELD' )
        error = nf90_put_att(ncid, id_ens, "long_name", "Ensemble members")
        call netcdf_err(error, 'DEFINING Ensemble LONG NAME' )
        error = nf90_put_att(ncid, id_ens, "units", "Ensemble rank")
        call netcdf_err(error, 'DEFINING Ensemble UNITS' )
        error = nf90_put_att(ncid, id_ens, "cartesian_axis", "Ens")
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

        error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_3d, id_swe_forc)
        call netcdf_err(error, 'DEFINING SWE_Forecast' )
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

        error = nf90_def_var(ncid, 'SWE_Analysis_Current', NF90_DOUBLE, dims_3d, id_swe_cur)
        call netcdf_err(error, 'DEFINING SWE_Analysis_Current' )
        error = nf90_put_att(ncid, id_swe_cur, "long_name", "Current Analysis Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING Cur SWE LONG NAME' )
        error = nf90_put_att(ncid, id_swe_cur, "units", "mm")
        call netcdf_err(error, 'DEFINING Cur SWE UNITS' )

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

        error = nf90_def_var(ncid, 'SND_Analysis_Current', NF90_DOUBLE, dims_3d, id_snwd_cur)
        call netcdf_err(error, 'DEFINING SND_Analysis_Current' )
        error = nf90_put_att(ncid, id_snwd_cur, "long_name", "Current Analysis Snow Depth")
        call netcdf_err(error, 'DEFINING Cur SND Analysis LONG NAME' )
        error = nf90_put_att(ncid, id_snwd_cur, "units", "mm")
        call netcdf_err(error, 'DEFINING Cur SND Analysis UNITS' )

        error = nf90_def_var(ncid, 'DA_Increment', NF90_DOUBLE, dims_3d, id_incr)
        call netcdf_err(error, 'DEFINING DA_Increment' )
        error = nf90_put_att(ncid, id_incr, "long_name", "DA_Increment at model grid")
        call netcdf_err(error, 'DEFINING DA_Increment LONG NAME' )
        error = nf90_put_att(ncid, id_incr, "units", "mm")
        call netcdf_err(error, 'DEFINING DA_Increment UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        error = nf90_def_var(ncid, 'ims_FootPrint', NF90_DOUBLE, dims_3d, id_imsftp)
        call netcdf_err(error, 'DEFINING ims_FootPrint' )
        error = nf90_put_att(ncid, id_imsftp, "long_name", "Assimilated IMS Data Foot Print")
        call netcdf_err(error, 'DEFINING ims_FootPrint LONG NAME' )
        error = nf90_put_att(ncid, id_imsftp, "units", "-")
        call netcdf_err(error, 'DEFINING ims_FootPrint UNITS' )

        !ens 
        dims_3dens(1) = dim_x
        dims_3dens(2) = dim_y
        dims_3dens(3) = dim_ens
        error = nf90_def_var(ncid, 'SND_Forecast_Ens', NF90_DOUBLE, dims_3dens, id_sndf_ens)
        call netcdf_err(error, 'DEFINING SND_Forecast_Ens' )
        error = nf90_put_att(ncid, id_sndf_ens, "long_name", "Ensemble Forecast Snow Depth")
        call netcdf_err(error, 'DEFINING Ensemble SND Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_sndf_ens, "units", "mm")
        call netcdf_err(error, 'DEFINING Ensemble SND Forecast UNITS' )
        
        error = nf90_def_var(ncid, 'SWE_Forecast_Ens', NF90_DOUBLE, dims_3dens, id_swef_ens)
        call netcdf_err(error, 'DEFINING SWE_Forecast_Ens' )
        error = nf90_put_att(ncid, id_swef_ens, "long_name", "Ensemble Forecast Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING Ensemble SWE Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_swef_ens, "units", "mm")
        call netcdf_err(error, 'DEFINING Ensemble SWE Forecast UNITS' )

        error = nf90_def_var(ncid, 'SND_Analysis_Ens', NF90_DOUBLE, dims_3dens, id_snda_ens)
        call netcdf_err(error, 'DEFINING SND_Analysis_Ens' )
        error = nf90_put_att(ncid, id_snda_ens, "long_name", "Ensemble Analysis Snow Depth")
        call netcdf_err(error, 'DEFINING Ensemble SND Analysis LONG NAME' )
        error = nf90_put_att(ncid, id_snda_ens, "units", "mm")
        call netcdf_err(error, 'DEFINING Ensemble SND Analysis UNITS' )

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
        
        error = nf90_def_var(ncid, 'snwdph@ObsValue', NF90_DOUBLE, dim_stn, id_obsstn)
        call netcdf_err(error, 'DEFINING snwdph@ObsValue' )
        error = nf90_put_att(ncid, id_obsstn, "long_name", "Observed at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue LONG NAME' )
        error = nf90_put_att(ncid, id_obsstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue UNITS' ) 
        
        ! error = nf90_def_var(ncid, 'snwdph@hofx', NF90_DOUBLE, dim_stn, id_forcstn) 
        ! call netcdf_err(error, 'DEFINING snwdph@hofx' )
        ! error = nf90_put_att(ncid, id_forcstn, "long_name", "Forecast at Observation Points")
        ! call netcdf_err(error, 'DEFINING snwdph@hofx LONG NAME' )
        ! error = nf90_put_att(ncid, id_forcstn, "units", "mm")
        ! call netcdf_err(error, 'DEFINING snwdph@hofx UNITS' )

        ! error = nf90_def_var(ncid, 'snwdph@omb', NF90_DOUBLE, dim_stn, id_innovstn)
        ! call netcdf_err(error, 'DEFINING snwdph@omb' )
        ! error = nf90_put_att(ncid, id_innovstn, "long_name", "Innovation(O-B) at Observation Points")
        ! call netcdf_err(error, 'DEFINING snwdph@omb LONG NAME' )
        ! error = nf90_put_att(ncid, id_innovstn, "units", "mm")
        ! call netcdf_err(error, 'DEFINING snwdph@omb UNITS' )

        error = nf90_def_var(ncid, 'index_back_at_obs', NF90_INT, dim_stn, id_indxobs) 
        call netcdf_err(error, 'DEFINING index_back_at_obs' )
        error = nf90_put_att(ncid, id_indxobs, "long_name", "Index of background at Observation Points")
        call netcdf_err(error, 'DEFINING index_back_at_obs LONG NAME' )
        error = nf90_put_att(ncid, id_indxobs, "units", "-")
        call netcdf_err(error, 'DEFINING index_back_at_obs UNITS' )
        
        ! error = nf90_def_var(ncid, 'anlfSCA', NF90_DOUBLE, dims_3d, id_anlscov)
        ! call netcdf_err(error, 'DEFINING anlfSCA' )
        ! error = nf90_put_att(ncid, id_anlscov, "long_name", "Analysis fractional Snow Covered Area")
        ! call netcdf_err(error, 'DEFINING anlfSCA LONG NAME' )
        ! error = nf90_put_att(ncid, id_anlscov, "units", "-")
        ! call netcdf_err(error, 'DEFINING anlfSCA UNITS' )

        ! eval points
        ! if (num_Eval>0) then 
                error = nf90_def_var(ncid, 'LatEvalPoints', NF90_DOUBLE, dim_eval, id_lateval)
                call netcdf_err(error, 'DEFINING LatEvalPoints' )
                error = nf90_put_att(ncid, id_lateval, "long_name", "Latitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LatEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_lateval, "units", "deg")
                call netcdf_err(error, 'DEFINING LatEvalPoints UNITS' )

                error = nf90_def_var(ncid, 'LonEvalPoints', NF90_DOUBLE, dim_eval, id_loneval)
                call netcdf_err(error, 'DEFINING LonEvalPoints' )
                error = nf90_put_att(ncid, id_loneval, "long_name", "Longitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LonEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_loneval, "units", "deg")
                call netcdf_err(error, 'DEFINING LonEvalPoints UNITS' )
                
                error = nf90_def_var(ncid, 'Obs_atEvalPts', NF90_DOUBLE, dim_eval, id_obseval)
                call netcdf_err(error, 'DEFINING Obs_atEvalPts' )
                error = nf90_put_att(ncid, id_obseval, "long_name", "Observed at Evaluation Points")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_obseval, "units", "mm")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts UNITS' )
                
                ! error = nf90_def_var(ncid, 'SNOFCS_atEvalPts', NF90_DOUBLE, dim_eval, id_forceval)
                ! call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts' )
                ! error = nf90_put_att(ncid, id_forceval, "long_name", "Forecast at Evaluation Points")
                ! call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts LONG NAME' )
                ! error = nf90_put_att(ncid, id_forceval, "units", "mm")
                ! call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts UNITS' )

                ! error = nf90_def_var(ncid, 'SNOANL_atEvalPts', NF90_DOUBLE, dim_eval, id_anleval)
                ! call netcdf_err(error, 'DEFINING SNOANL_atEvalPts' )
                ! error = nf90_put_att(ncid, id_anleval, "long_name", "Analysis at Evaluation Points")
                ! call netcdf_err(error, 'DEFINING SNOANL_atEvalPts LONG NAME' )
                ! error = nf90_put_att(ncid, id_anleval, "units", "mm")
                ! call netcdf_err(error, 'DEFINING SNOANL_atEvalPts UNITS' )

                ! error = nf90_def_var(ncid, 'SNOANL_Cur_atEvalPts', NF90_DOUBLE, dim_eval, id_cur_anleval)
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts' )
                ! error = nf90_put_att(ncid, id_cur_anleval, "long_name", "Current Analysis at Evaluation Points")
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts LONG NAME' )
                ! error = nf90_put_att(ncid, id_cur_anleval, "units", "mm")
                ! call netcdf_err(error, 'DEFINING SNOANL_Cur_atEvalPts UNITS' )
                
                ! error = nf90_def_var(ncid, 'incr_atEvalPts', NF90_DOUBLE, dim_eval, id_increval)
                ! call netcdf_err(error, 'DEFINING incr_atEvalPts' )
                ! error = nf90_put_att(ncid, id_increval, "long_name", "Increment at Evaluation Points")
                ! call netcdf_err(error, 'DEFINING incr_atEvalPts LONG NAME' )
                ! error = nf90_put_att(ncid, id_increval, "units", "mm")
                ! call netcdf_err(error, 'DEFINING incr_atEvalPts UNITS' )id_indxeval

                error = nf90_def_var(ncid, 'indx_atEvalPts', NF90_INT, dim_eval, id_indxeval)
                call netcdf_err(error, 'DEFINING indx_atEvalPts' )
                error = nf90_put_att(ncid, id_indxeval, "long_name", "Background index at Evaluation Points")
                call netcdf_err(error, 'DEFINING indx_atEvalPts LONG NAME')
                error = nf90_put_att(ncid, id_indxeval, "units", "-")
                call netcdf_err(error, 'DEFINING indx_atEvalPts UNITS')

        ! endif

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
        allocate(ens_data(ens_size))
        do i = 1, ens_size
            ens_data(i) = i
        enddo

        error = nf90_put_var( ncid, id_x, x_data)
        call netcdf_err(error, 'WRITING XAXIS RECORD' )
        error = nf90_put_var( ncid, id_y, y_data)
        call netcdf_err(error, 'WRITING YAXIS RECORD' )
        error = nf90_put_var( ncid, id_time, times)
        call netcdf_err(error, 'WRITING TIME RECORD')
        error = nf90_put_var( ncid, id_ens, ens_data)
        call netcdf_err(error, 'WRITING Ensemble RECORD')

        allocate(dum2d(idim,jdim))
        dims_strt(1:3) = 1
        dims_end(1) = idim
        dims_end(2) = jdim
        dims_end(3) = 1
        
        dum2d = reshape(landmask, (/idim,jdim/))
        error = nf90_put_var( ncid, id_landmask, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING LandMask RECORD' )

        dum2d = reshape(snoforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe_forc, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Forecast RECORD' )

        dum2d = reshape(snoanl, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Analysis RECORD' ) 

        dum2d = reshape(SWEANL_Cur, (/idim,jdim/))
        error = nf90_put_var( ncid, id_swe_cur, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING Cur SWE Analysis RECORD' )

        dum2d = reshape(snwdforc, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwdf, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND Forecast RECORD' )

        dum2d = reshape(snwdanal, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwd, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND analysis RECORD' )

        dum2d = reshape(SNDANL_Cur, (/idim,jdim/))
        error = nf90_put_var( ncid, id_snwd_cur, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING Cur SND analysis RECORD')

        dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )

        dum2d = reshape(IMS_Foot_Print, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imsftp, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD')
        ! dum2d = reshape(anl_fSCA, (/idim,jdim/))
        ! error = nf90_put_var( ncid, id_anlscov, dum2d, dims_strt, dims_end)
        ! call netcdf_err(error, 'WRITING anlfSCA RECORD')

        dum2d = reshape(incr_atGrid, (/idim,jdim/))
        error = nf90_put_var( ncid, id_incr, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING incr_atGrid RECORD')

        Do ens_indx = 1, ens_size            
            dims_strt(3) = ens_indx
            dum2d = reshape(SNDFCS_Ens(ens_indx, :), (/idim,jdim/))
            error = nf90_put_var( ncid, id_sndf_ens, dum2d, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING Ensemble SND Forecast RECORD' )
            dum2d = reshape(SWEFCS_Ens(ens_indx, :), (/idim,jdim/))
            error = nf90_put_var( ncid, id_swef_ens, dum2d, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING Ensemble SWE Forecast RECORD' )
            ! note ens mean is written as SND analysis earlier
            dum2d = reshape(anl_at_Grid_ens(ens_indx, :), (/idim,jdim/))
            error = nf90_put_var( ncid, id_snda_ens, dum2d, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING Ensemble SND Analysis RECORD' )
        End do
        
        ! obs points (obs, hofx, omb) 
        error = nf90_put_var( ncid, id_latstn, Lat_atObs)
        call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
        call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_obsstn, OBS_stn)
        call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

        ! error = nf90_put_var( ncid, id_forcstn, FCS_at_stn)
        ! call netcdf_err(error, 'WRITING SNOFCS_atObsPts RECORD' )

        ! error = nf90_put_var( ncid, id_innovstn, OmB_innov_at_stn)
        ! call netcdf_err(error, 'WRITING innov_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_indxobs, index_back_atObs)
        call netcdf_err(error, 'WRITING indx_atobs RECORD')

        ! eval points
        if (num_Eval>0) then 
                error = nf90_put_var( ncid, id_lateval, Lat_atEvalPts)
                call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_loneval, Lon_atEvalPts)
                call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_obseval, Obs_atEvalPts)
                call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                ! error = nf90_put_var( ncid, id_forceval, SNOFCS_atEvalPts)
                ! call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD')
                ! error = nf90_put_var( ncid, id_anleval, SNOANL_atEvalPts)
                ! call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_cur_anleval, SNOANL_Cur_atEvalPts)
                ! call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_increval, incr_atEvalPts)
                ! call netcdf_err(error, 'WRITING incr_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_indxeval, index_back_atEval)
                call netcdf_err(error, 'WRITING indx_atEval RECORD' )
        else
                error = nf90_put_var( ncid, id_lateval, 0.0)
                call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_loneval, 0.0)
                call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_obseval, 0.0)
                call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                ! error = nf90_put_var( ncid, id_forceval, 0.0)
                ! call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_anleval, 0.0)
                ! call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_cur_anleval, 0.0)
                ! call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )
                ! error = nf90_put_var( ncid, id_increval, 0.0)
                ! call netcdf_err(error, 'WRITING incr_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_indxeval, 0)
                call netcdf_err(error, 'WRITING indx_atEval RECORD' )
        endif

        deallocate(x_data, y_data, ens_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 End subroutine Write_DA_Outputs_ens

  Subroutine Write_DA_Outputs_ens_vect(da_output_file, vector_restart_prefix, &
        y_str, m_str, d_str, h_str, &
        MYRANK, NPROCS, LENSFC, LENSFC_proc, ens_size, num_stn, num_Eval, &
        max_num_nearStn,  &
        N_sA, N_sA_Ext, mpiReal_size, mpiInt_size,    &
        SNDnoDA, SWEnoDA,  noahmp,   &
        SNDFCS, SWEFCS, SNDANL, SWEANL, incr_at_Grid,   &  !
        index_obs_assmilated,   &
        SNCOV_IMS, IMS_Foot_Print, SNO_IMS_at_Grid,    &
        SNODENS_Grid, SCF_Grid, &
        noahmp_ensm_swe, noahmp_ensm_snow_depth, &            
        Lat_atObs, Lon_atObs, OROG_at_stn, OBS_stn, index_back_atObs, &! SNOFCS_at_stn, SNOFCS_atObs_ens, OmB_innov_at_stn, SNOANL_atStn  
        NEXC, Index_Obs_Excluded,  &
        Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, &
        index_back_atEval, index_obs_atEval)             
    ! SNDFCS_Ens, SWEFCS_Ens, incr_at_Grid_ensM, anl_at_Grid_ens)
    !------------------------------------------------------------------
    ! Write DA ouputs: 
    ! forecast SWE
    ! analysis SWE
    ! analysis Snow Depth
    ! innovation at grid
    !------------------------------------------------------------------
    implicit none
    include "mpif.h"
    
    CHARACTER(LEN=*), Intent(In)    :: da_output_file, vector_restart_prefix 
    CHARACTER(LEN=*), Intent(In)    :: y_str, m_str, d_str, h_str
    integer, intent(in)     :: myrank, NPROCS, lensfc, LENSFC_proc, ens_size
    integer, intent(in)     :: num_stn, num_Eval, max_num_nearStn  
    integer, intent(in)     :: N_sA, N_sA_Ext, mpiReal_size, mpiInt_size

    Real, intent(in)        :: SNDnoDA(LENSFC_proc), SWEnoDA(LENSFC_proc)
    type(noahmp_type)       :: noahmp(ens_size)    

    real, intent(in)  :: SNDFCS(ens_size+1,LENSFC_proc), SWEFCS(ens_size+1,LENSFC_proc), &
                        SNDANL(ens_size+1,LENSFC_proc), SWEANL(ens_size+1,LENSFC_proc), &
                        incr_at_Grid(ens_size+1, LENSFC_proc)

    real, intent(in) :: SNODENS_Grid(ens_size+1, LENSFC_proc), SCF_Grid(ens_size+1, LENSFC_proc)
    
    Real, intent(in) :: noahmp_ensm_swe(LENSFC_proc), noahmp_ensm_snow_depth(LENSFC_proc)

    Real, intent(in)    :: Lat_atObs(num_stn), Lon_atObs(num_stn) 
    Real, intent(in)    :: OROG_at_stn(num_stn), OBS_stn(num_stn)
    integer, intent(in) :: index_back_atObs(num_stn)
    real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), &
                           Obs_atEvalPts(num_Eval) 
    integer, intent(in) :: index_back_atEval(num_Eval), index_obs_atEval(num_Eval)
    !1.4 track obserations assimilated
    Integer          :: index_obs_assmilated(max_num_nearStn+1, LENSFC_proc)
    ! excluded obs
    Integer, intent(in) :: NEXC
    Integer, intent(in) :: Index_Obs_Excluded(NEXC)

    Real, intent(in)    :: SNCOV_IMS(LENSFC_proc), IMS_Foot_Print(LENSFC_proc), &
                           SNO_IMS_at_Grid(LENSFC_proc)  !, anl_fSCA(lensfc)

    integer            :: fsize=65536, inital=0
    integer            :: header_buffer_val = 16384
    integer            :: dims_2d(2), dims_strt(2), dims_end(2)
    integer            :: error, i, ncid
    integer    :: dim_x, dim_ens, dim_stn, dim_eval, dim_sn, dim_snso, dim_exc
    integer    :: dim_no, id_no
    integer    :: id_x, id_ens 
    integer    :: id_snwd_noda, id_swe_noda    !, id_snwd_noda_eval, id_swe_noda_eval
    integer    :: id_snwd_forc, id_swe_forc, id_imscov, id_imsftp, id_imssnd 
    integer    :: id_incr, id_snwd, id_swe
    integer    :: id_den, id_scf
    integer    :: id_swe_noah, id_snd_noah

    integer       :: id_latstn, id_lonstn, id_elestn, id_obsstn, id_bkobs  !, id_forcstn, id_innovstn, id_anlstn !, id_landmask
    integer       :: id_lateval, id_loneval, id_obseval, id_bkeval, id_indxobseval, id_exc 
    ! id_forceval, id_ele_evl, id_ele_fcs_evl, 
    
    real(kind=4)                :: times
    real(kind=4), allocatable   :: x_data(:)  !, y_data(:)
    integer                     :: dims_3d(3), iv, varid
    integer     :: id_weasd, id_sndpth_m, id_snowxy, id_sneqvoxy, id_stc, &
                    id_zsnsoxy, id_tsnoxy, id_snicexy, id_snliqxy

    Real  :: DUMMY1(lensfc), DUMMY2(lensfc), DUMMY3(lensfc), DUMMY4(lensfc), DUMMY45(lensfc)
    Real    :: DUMMY5(lensfc, 7), DUMMY6(lensfc, 3), &
               DUMMY7(lensfc, 3), DUMMY8(lensfc, 3), DUMMY9(lensfc), DUMMY10(lensfc)

    Real    :: DUMMY11(lensfc), DUMMY12(lensfc)

    integer :: DUMMY13(LENSFC)
    
    Integer     :: dest_Aoffset, dest_Aoffset_end, arLen, IERR, ixy, ie, jndx

    LOGICAL     :: file_exists
    CHARACTER(LEN=2)      :: ensCH 
    CHARACTER(LEN=500)    :: restart_filename 
    

    ! print*,"Process ", myrank, "writing output data to: ",trim(output_file)
    !--- create the file
    if (myrank ==0) then 
        error = NF90_CREATE(trim(da_output_file), IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
        call netcdf_err(error, 'CREATING FILE='//trim(da_output_file) )     
        ! print*,"Process ", myrank, "created file: ",trim(output_file)
        !--- define dimensions
        error = nf90_def_dim(ncid, 'location', lensfc, dim_x)
        call netcdf_err(error, 'DEFINING location DIMENSION' )
        ! error = nf90_def_dim(ncid, 'yaxis_1', jdim, dim_y)
        ! call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
        ! error = nf90_def_dim(ncid, 'Time', 1, dim_time)
        ! call netcdf_err(error, 'DEFINING TIME DIMENSION' )
        error = nf90_def_dim(ncid, 'nearest_obs', max_num_nearStn+1, dim_no)
        call netcdf_err(error, 'DEFINING nearest_obs DIMENSION' )
        !ens
        error = nf90_def_dim(ncid, 'Ensemble', ens_size+1, dim_ens)
        call netcdf_err(error, 'DEFINING Ensemble DIMENSION' )
        ! obs points
        error = nf90_def_dim(ncid, 'obs_points', num_stn, dim_stn)
        call netcdf_err(error, 'DEFINING obs_points' )
        ! eval obs points        
        if (num_Eval>0) then
                error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        else
                error = nf90_def_dim(ncid, 'eval_points', 1, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        endif
        !excluded obs indx 
        error = nf90_def_dim(ncid, 'obs_excluded', NEXC, dim_exc)
        call netcdf_err(error, 'DEFINING obs_excluded')

        !--- define fields
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

        error = nf90_def_var(ncid, 'elevation@MetaData', NF90_DOUBLE, dim_stn, id_elestn)
        call netcdf_err(error, 'DEFINING elevation@MetaData' )
        error = nf90_put_att(ncid, id_elestn, "long_name", "elevation at Observation Points")
        call netcdf_err(error, 'DEFINING elevation@MetaData LONG NAME' )
        error = nf90_put_att(ncid, id_elestn, "units", "m")
        call netcdf_err(error, 'DEFINING elevation@MetaData UNITS' )
        
        error = nf90_def_var(ncid, 'snwdph@ObsValue', NF90_DOUBLE, dim_stn, id_obsstn)
        call netcdf_err(error, 'DEFINING snwdph@ObsValue' )
        error = nf90_put_att(ncid, id_obsstn, "long_name", "Observed at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue LONG NAME' )
        error = nf90_put_att(ncid, id_obsstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@ObsValue UNITS' )

        error = nf90_def_var(ncid, 'index_back_at_obs', NF90_INT, dim_stn, id_bkobs) 
        call netcdf_err(error, 'DEFINING index_back_at_obs' )
        error = nf90_put_att(ncid, id_bkobs, "long_name", "Index of background at Observation Points")
        call netcdf_err(error, 'DEFINING index_back_at_obs LONG NAME' )
        error = nf90_put_att(ncid, id_bkobs, "units", "-")
        call netcdf_err(error, 'DEFINING index_back_at_obs UNITS' )
        ! eval points
! if (num_Eval>0) then 
        error = nf90_def_var(ncid, 'LatEvalPoints', NF90_DOUBLE, dim_eval, id_lateval)
        call netcdf_err(error, 'DEFINING LatEvalPoints' )
        error = nf90_put_att(ncid, id_lateval, "long_name", "Latitude at Evaluation Points")
        call netcdf_err(error, 'DEFINING LatEvalPoints LONG NAME' )
        error = nf90_put_att(ncid, id_lateval, "units", "deg")
        call netcdf_err(error, 'DEFINING LatEvalPoints UNITS' )

        error = nf90_def_var(ncid, 'LonEvalPoints', NF90_DOUBLE, dim_eval, id_loneval)
        call netcdf_err(error, 'DEFINING LonEvalPoints' )
        error = nf90_put_att(ncid, id_loneval, "long_name", "Longitude at Evaluation Points")
        call netcdf_err(error, 'DEFINING LonEvalPoints LONG NAME' )
        error = nf90_put_att(ncid, id_loneval, "units", "deg")
        call netcdf_err(error, 'DEFINING LonEvalPoints UNITS' ) 

        error = nf90_def_var(ncid, 'Obs_atEvalPts', NF90_DOUBLE, dim_eval, id_obseval)
        call netcdf_err(error, 'DEFINING Obs_atEvalPts' )
        error = nf90_put_att(ncid, id_obseval, "long_name", "Observed at Evaluation Points")
        call netcdf_err(error, 'DEFINING Obs_atEvalPts LONG NAME' )
        error = nf90_put_att(ncid, id_obseval, "units", "mm")
        call netcdf_err(error, 'DEFINING Obs_atEvalPts UNITS' )            
        ! index of background at eval points
        error = nf90_def_var(ncid, 'Index_Back_EvalPts', NF90_INT, dim_eval, id_bkeval)
        call netcdf_err(error, 'DEFINING Index_Back_EvalPts' )
        error = nf90_put_att(ncid, id_bkeval, "long_name", "Index of Background at Evaluation Points")
        call netcdf_err(error, 'DEFINING Index_Back_EvalPts LONG NAME' )
        error = nf90_put_att(ncid, id_bkeval, "units", "-")
        call netcdf_err(error, 'DEFINING Index_Back_EvalPts UNITS' )
        ! index of Obs at eval points 
        error = nf90_def_var(ncid, 'Index_Obs_EvalPts', NF90_INT, dim_eval, id_indxobseval)
        call netcdf_err(error, 'DEFINING Index_Obs_EvalPts' )
        error = nf90_put_att(ncid, id_indxobseval, "long_name", "Index of Observations at Evaluation Points")
        call netcdf_err(error, 'DEFINING Index_Obs_EvalPts LONG NAME' )
        error = nf90_put_att(ncid, id_indxobseval, "units", "-")
        call netcdf_err(error, 'DEFINING Index_Obs_EvalPts UNITS' )
    ! endif eval
        ! index of Obs excluded
        error = nf90_def_var(ncid, 'Index_Obs_Excluded', NF90_INT, dim_exc, id_exc)
        call netcdf_err(error, 'DEFINING Index_Obs_Excluded' )
        error = nf90_put_att(ncid, id_exc, "long_name", "Index of Observations Excluded")
        call netcdf_err(error, 'DEFINING Index_Obs_Excluded LONG NAME' )
        error = nf90_put_att(ncid, id_exc, "units", "-")
        call netcdf_err(error, 'DEFINING Index_Obs_Excluded UNITS' )

        ! No DA 
        error = nf90_def_var(ncid, 'SWE_Forecast_NoDA', NF90_DOUBLE, dim_x, id_swe_noda)
        call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA' )
        error = nf90_put_att(ncid, id_swe_noda, "long_name", "Forecast Snow Water Equivalent No DA")
        call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA LONG NAME' )
        error = nf90_put_att(ncid, id_swe_noda, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA UNITS' )

        error = nf90_def_var(ncid, 'SND_Forecast_NoDA', NF90_DOUBLE, dim_x, id_snwd_noda)
        call netcdf_err(error, 'DEFINING SND_Forecast_NoDA' )
        error = nf90_put_att(ncid, id_snwd_noda, "long_name", "Forecast Snow Depth No DA")
        call netcdf_err(error, 'DEFINING SND_Forecast_NoDA LONG NAME' )
        error = nf90_put_att(ncid, id_snwd_noda, "units", "mm")
        call netcdf_err(error, 'DEFINING SND_Forecast_NoDA UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dim_x, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        error = nf90_def_var(ncid, 'imsSND', NF90_DOUBLE, dim_x, id_imssnd)
        call netcdf_err(error, 'DEFINING imsSND' )
        error = nf90_put_att(ncid, id_imssnd, "long_name", "IMS SND")
        call netcdf_err(error, 'DEFINING imsSND LONG NAME' )
        error = nf90_put_att(ncid, id_imssnd, "units", "mm")
        call netcdf_err(error, 'DEFINING imsSND UNITS' )

        error = nf90_def_var(ncid, 'ims_FootPrint', NF90_DOUBLE, dim_x, id_imsftp)
        call netcdf_err(error, 'DEFINING ims_FootPrint' )
        error = nf90_put_att(ncid, id_imsftp, "long_name", "Assimilated IMS Data Foot Print")
        call netcdf_err(error, 'DEFINING ims_FootPrint LONG NAME' )
        error = nf90_put_att(ncid, id_imsftp, "units", "-")
        call netcdf_err(error, 'DEFINING ims_FootPrint UNITS' )

    ! noaa ensemble mean 
        error = nf90_def_var(ncid, 'SWE_EnsM_Noah', NF90_DOUBLE, dim_x, id_swe_noah)
        call netcdf_err(error, 'DEFINING SWE_EnsM_Noah' )
        error = nf90_put_att(ncid, id_swe_noah, "long_name", "Noah Ensemble Mean Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE_EnsM_Noah LONG NAME' )
        error = nf90_put_att(ncid, id_swe_noah, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE_EnsM_Noah UNITS' )
       
        error = nf90_def_var(ncid, 'SND_EnsM_Noah', NF90_DOUBLE, dim_x, id_snd_noah)
        call netcdf_err(error, 'DEFINING SND_EnsM_Noah' )
        error = nf90_put_att(ncid, id_snd_noah, "long_name", "Noah Ensemble Mean Snow Depth")
        call netcdf_err(error, 'DEFINING SND_EnsM_Noah LONG NAME' )
        error = nf90_put_att(ncid, id_snd_noah, "units", "mm")
        call netcdf_err(error, 'DEFINING SND_EnsM_Noah UNITS' )

        dims_2d(1) = dim_no
        dims_2d(2) = dim_x
        error = nf90_def_var(ncid, 'index_obs_assmilated', NF90_INT, dims_2d, id_no) 
        call netcdf_err(error, 'DEFINING index_obs_assmilated' )
        error = nf90_put_att(ncid, id_no, "long_name", "Indics of assimilated Observation Points")
        call netcdf_err(error, 'DEFINING index_obs_assmilated LONG NAME' )
        error = nf90_put_att(ncid, id_no, "units", "-")
        call netcdf_err(error, 'DEFINING index_obs_assmilated UNITS' )

        dims_2d(1) = dim_ens 
        ! dims_3d(2) = dim_ens
        dims_2d(2) = dim_x

        error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_2d, id_swe_forc)
        call netcdf_err(error, 'DEFINING SWE_Forecast' )
        error = nf90_put_att(ncid, id_swe_forc, "long_name", "Forecast Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_swe_forc, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE Forecast UNITS' )

        error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_2d, id_snwd_forc)
        call netcdf_err(error, 'DEFINING SND Forecast' )
        error = nf90_put_att(ncid, id_snwd_forc, "long_name", "Forecast Snow Depth")
        call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
        error = nf90_put_att(ncid, id_snwd_forc, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

        error = nf90_def_var(ncid, 'SWE_Analysis', NF90_DOUBLE, dims_2d, id_swe)
        call netcdf_err(error, 'DEFINING SWE_Analysis' )
        error = nf90_put_att(ncid, id_swe, "long_name", "Analysis Snow Water Equivalent")
        call netcdf_err(error, 'DEFINING SWE LONG NAME' )
        error = nf90_put_att(ncid, id_swe, "units", "mm")
        call netcdf_err(error, 'DEFINING SWE UNITS' )
        !
        error = nf90_def_var(ncid, 'SND_Analysis', NF90_DOUBLE, dims_2d, id_snwd)
        call netcdf_err(error, 'DEFINING SND Analyis' )
        error = nf90_put_att(ncid, id_snwd, "long_name", "Analysis Snow Depth")
        call netcdf_err(error, 'DEFINING SND Analysis LONG NAME' )
        error = nf90_put_att(ncid, id_snwd, "units", "mm")
        call netcdf_err(error, 'DEFINING SND Analysis UNITS' )
        !
        error = nf90_def_var(ncid, 'DA_Increment', NF90_DOUBLE, dims_2d, id_incr)
        call netcdf_err(error, 'DEFINING DA_Increment' )
        error = nf90_put_att(ncid, id_incr, "long_name", "DA_Increment at model grid")
        call netcdf_err(error, 'DEFINING DA_Increment LONG NAME' )
        error = nf90_put_att(ncid, id_incr, "units", "mm")
        call netcdf_err(error, 'DEFINING DA_Increment UNITS' )

        error = nf90_def_var(ncid, 'DA_Density', NF90_DOUBLE, dims_2d, id_den)
        call netcdf_err(error, 'DEFINING DA_Density' )
        error = nf90_put_att(ncid, id_den, "long_name", "DA_Density at model grid")
        call netcdf_err(error, 'DEFINING DA_Density LONG NAME' )
        error = nf90_put_att(ncid, id_den, "units", "kg/m3")
        call netcdf_err(error, 'DEFINING DA_Density UNITS' )
        
        error = nf90_def_var(ncid, 'DA_SCF', NF90_DOUBLE, dims_2d, id_scf)
        call netcdf_err(error, 'DEFINING DA_SCF' )
        error = nf90_put_att(ncid, id_scf, "long_name", "DA_SCF at model grid")
        call netcdf_err(error, 'DEFINING DA_SCF LONG NAME' )
        error = nf90_put_att(ncid, id_scf, "units", "-")
        call netcdf_err(error, 'DEFINING DA_SCF UNITS' )

        error = nf90_enddef(ncid, header_buffer_val,4,0,4)
        call netcdf_err(error, 'DEFINING HEADER' )

! Write out values            
        ! obs points (obs, hofx, omb) 
        error = nf90_put_var( ncid, id_latstn, Lat_atObs)
        call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
        call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' ) 

        error = nf90_put_var( ncid, id_elestn, OROG_at_stn)
        call netcdf_err(error, 'WRITING OROG_at_stn RECORD' )

        error = nf90_put_var( ncid, id_obsstn, OBS_stn)
        call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_bkobs, index_back_atObs)
        call netcdf_err(error, 'WRITING indx_atobs RECORD')

        ! eval points
        if (num_Eval>0) then 
            error = nf90_put_var( ncid, id_lateval, Lat_atEvalPts)
            call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

            error = nf90_put_var( ncid, id_loneval, Lon_atEvalPts)
            call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' ) 

            error = nf90_put_var( ncid, id_obseval, Obs_atEvalPts)
            call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )
            !eval back indx 
            error = nf90_put_var( ncid, id_bkeval, index_back_atEval)
            call netcdf_err(error, 'WRITING index_back_atEval RECORD' )
            !eval obs indx 
            error = nf90_put_var( ncid, id_indxobseval, index_obs_atEval)
            call netcdf_err(error, 'WRITING index_obs_atEval RECORD' )
        else
            error = nf90_put_var( ncid, id_lateval, 0.0)
            call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

            error = nf90_put_var( ncid, id_loneval, 0.0)
            call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )  

            error = nf90_put_var( ncid, id_obseval, 0.0)
            call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )
            !eval back indx 
            error = nf90_put_var( ncid, id_bkeval, 1)
            call netcdf_err(error, 'WRITING index_back_atEval RECORD' )
            !eval obs indx 
            error = nf90_put_var( ncid, id_indxobseval, -1)
            call netcdf_err(error, 'WRITING index_obs_atEval RECORD' )
        endif

        error = nf90_put_var(ncid, id_exc, Index_Obs_Excluded)
        call netcdf_err(error, 'WRITING Index_Obs_Excluded RECORD' )

        print*, "done writing 1d arrays"
    endif

    if (MYRANK > 0) then
        call MPI_SEND(SNDnoDA, LENSFC_proc, mpiReal_size, 0,   &
                                MYRANK+1, MPI_COMM_WORLD, IERR) 
        call MPI_SEND(SWEnoDA, LENSFC_proc, mpiReal_size, 0,   &
                                100*(MYRANK+1), MPI_COMM_WORLD, IERR)
        call MPI_SEND(SNCOV_IMS, LENSFC_proc, mpiReal_size, 0,   &
                                200*(MYRANK+1), MPI_COMM_WORLD, IERR)
        call MPI_SEND(IMS_Foot_Print, LENSFC_proc, mpiReal_size, 0,   &
                                300*(MYRANK+1), MPI_COMM_WORLD, IERR)
        call MPI_SEND(SNO_IMS_at_Grid, LENSFC_proc, mpiReal_size, 0,   &
                                500*(MYRANK+1), MPI_COMM_WORLD, IERR)

        call MPI_SEND(noahmp_ensm_swe, LENSFC_proc, mpiReal_size, 0,   &
                                900*(MYRANK+1), MPI_COMM_WORLD, IERR) 
        call MPI_SEND(noahmp_ensm_snow_depth, LENSFC_proc, mpiReal_size, 0,   &
                                1000*(MYRANK+1), MPI_COMM_WORLD, IERR)     

    Endif
    if (MYRANK == 0) then        
        DUMMY1(1:LENSFC_proc) = SNDnoDA
        DUMMY2(1:LENSFC_proc) = SWEnoDA
        DUMMY3(1:LENSFC_proc) = SNCOV_IMS
        DUMMY4(1:LENSFC_proc) = IMS_Foot_Print
        DUMMY45(1:LENSFC_proc) = SNO_IMS_at_Grid

        DUMMY9(1:LENSFC_proc) = noahmp_ensm_swe 
        DUMMY10(1:LENSFC_proc) = noahmp_ensm_snow_depth

        Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
            if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                arLen = N_sA + 1  
            else
                dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                arLen = N_sA
            endif
            call MPI_RECV(DUMMY1(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(DUMMY2(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(DUMMY3(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(DUMMY4(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(DUMMY45(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            500*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)

            call MPI_RECV(DUMMY9(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            900*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(DUMMY10(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                            1000*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        End do
        !no da 
        error = nf90_put_var(ncid, id_snwd_noda, DUMMY1)  !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SND Forecast NoDA RECORD' )
        ! dum2d = reshape(SWEFCS, (/idim,jdim/))
        error = nf90_put_var(ncid, id_swe_noda, DUMMY2)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING SWE Forecast NoDA RECORD' )

        error = nf90_put_var(ncid, id_imscov, DUMMY3)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )
        
        error = nf90_put_var(ncid, id_imssnd, DUMMY45)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imssnd RECORD' )

        error = nf90_put_var(ncid, id_imsftp, DUMMY4)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING ims foot print RECORD' )
        
        error = nf90_put_var(ncid, id_swe_noah, DUMMY9)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING id_swe_noah RECORD' )
        error = nf90_put_var(ncid, id_snd_noah, DUMMY10)   !, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING id_snd_noah RECORD' )

        print*, "done writing ens mean outs"
    endif
    
! allocate(dum2d(idim,jdim))
    Do ie=1, ens_size+1
        dims_strt(1) = ie
        dims_strt(2) = 1
        dims_end(1) = 1    !ens_size
        dims_end(2) = lensfc

        if (MYRANK > 0) then
            call MPI_SEND(SNDFCS(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                400*(MYRANK+1), MPI_COMM_WORLD, IERR)  
            call MPI_SEND(SWEFCS(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                    500*(MYRANK+1), MPI_COMM_WORLD, IERR)   
            call MPI_SEND(SNDANL(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                600*(MYRANK+1), MPI_COMM_WORLD, IERR)  
            call MPI_SEND(SWEANL(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                700*(MYRANK+1), MPI_COMM_WORLD, IERR)   
            call MPI_SEND(incr_at_Grid(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                800*(MYRANK+1), MPI_COMM_WORLD, IERR)      
                                
            call MPI_SEND(SNODENS_Grid(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                1100*(MYRANK+1), MPI_COMM_WORLD, IERR)      
            call MPI_SEND(SCF_Grid(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                1200*(MYRANK+1), MPI_COMM_WORLD, IERR)      
        Endif
        if (MYRANK == 0) then
            DUMMY1(1:LENSFC_proc) = SNDFCS(ie,:)
            DUMMY2(1:LENSFC_proc) = SWEFCS(ie,:)
            DUMMY3(1:LENSFC_proc) = SNDANL(ie,:)
            DUMMY4(1:LENSFC_proc) = SWEANL(ie,:)
            DUMMY9(1:LENSFC_proc) = incr_at_Grid(ie,:)

            DUMMY11(1:LENSFC_proc) = SNODENS_Grid(ie,:)
            DUMMY12(1:LENSFC_proc) = SCF_Grid(ie,:)

            ! DUMMY13(1:LENSFC_proc, :) = index_obs_assmilated(:, :)

            Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1  
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(DUMMY1(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                400*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY2(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                500*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY3(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                600*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY4(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                700*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY9(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                800*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)

                call MPI_RECV(DUMMY11(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                1100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY12(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                1200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
            error = nf90_put_var(ncid, id_snwd_forc, DUMMY1, dims_strt, dims_end)                                
            call netcdf_err(error, 'WRITING SND Forecast RECORD')
            error = nf90_put_var(ncid, id_swe_forc, DUMMY2, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SWE Forecast RECORD')
            error = nf90_put_var(ncid, id_snwd, DUMMY3, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SND Analysis RECORD')
            error = nf90_put_var(ncid, id_swe, DUMMY4, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SWE Analysis RECORD')
            error = nf90_put_var(ncid, id_incr, DUMMY9, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING incr_atGrid RECORD')

            error = nf90_put_var(ncid, id_den, DUMMY11, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_den RECORD')
            error = nf90_put_var(ncid, id_scf, DUMMY12, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_scf RECORD')

            print*, "done writing SWE SND Ens forc, ens ", ie
        endif
    End do

    Do ie=1, max_num_nearStn+1
        if (MYRANK > 0) then   
            call MPI_SEND(index_obs_assmilated(ie,:), LENSFC_proc, mpiInt_size, 0,   &
                                1300*(MYRANK+1), MPI_COMM_WORLD, IERR)      
        Endif
        if (MYRANK == 0) then
            DUMMY13(1:LENSFC_proc) = index_obs_assmilated(ie,:)
            Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1  
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(DUMMY13(dest_Aoffset:dest_Aoffset_end), arLen, mpiInt_size, ixy,      &
                            1300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do            
            dims_strt(1) = ie
            dims_strt(2) = 1
            dims_end(1) = 1    
            dims_end(2) = lensfc
            error = nf90_put_var(ncid, id_no, DUMMY13, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_no RECORD')
        endif
    End do
    if (MYRANK == 0) then
        error = nf90_close(ncid)
        call netcdf_err(error, 'closing DA file' )
        print*, "done writing DA file"
    endif

    Do ie=1, ens_size
        if (MYRANK > 0) then
            call MPI_SEND(noahmp(ie)%swe, LENSFC_proc, mpiReal_size, 0,   &
                                    MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(noahmp(ie)%snow_depth, LENSFC_proc, mpiReal_size, 0,   &
                                    100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(noahmp(ie)%active_snow_layers, LENSFC_proc, mpiReal_size, 0,   &
                                    200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(noahmp(ie)%swe_previous, LENSFC_proc, mpiReal_size, 0,   &
                                    300*(MYRANK+1), MPI_COMM_WORLD, IERR)
            Do iv=1, 7                        
                call MPI_SEND(noahmp(ie)%snow_soil_interface(:,iv), LENSFC_proc, mpiReal_size, 0,   &
                                400*(MYRANK+1) + 10*iv, MPI_COMM_WORLD, IERR)
            Enddo
            Do iv=1, 3
                call MPI_SEND(noahmp(ie)%temperature_snow(:,iv), LENSFC_proc, mpiReal_size, 0,   &
                                        500*(MYRANK+1) + 20*iv, MPI_COMM_WORLD, IERR)
                call MPI_SEND(noahmp(ie)%snow_ice_layer(:,iv), LENSFC_proc, mpiReal_size, 0,   &
                                        600*(MYRANK+1) + 30*iv, MPI_COMM_WORLD, IERR)
                call MPI_SEND(noahmp(ie)%snow_liq_layer(:,iv), LENSFC_proc, mpiReal_size, 0,   &
                                        700*(MYRANK+1) + 40*iv, MPI_COMM_WORLD, IERR)
            enddo
            call MPI_SEND(noahmp(ie)%temperature_soil, LENSFC_proc, mpiReal_size, 0,   &
                                        800*(MYRANK+1), MPI_COMM_WORLD, IERR)

            call MPI_SEND(SCF_Grid(ie,:), LENSFC_proc, mpiReal_size, 0,   &
                                1300*(MYRANK+1), MPI_COMM_WORLD, IERR)      
        Endif
        if (MYRANK == 0) then
            DUMMY1(1:LENSFC_proc) = noahmp(ie)%swe
            DUMMY2(1:LENSFC_proc) = noahmp(ie)%snow_depth
            DUMMY3(1:LENSFC_proc) = noahmp(ie)%active_snow_layers
            DUMMY4(1:LENSFC_proc) = noahmp(ie)%swe_previous

            DUMMY12(1:LENSFC_proc) = SCF_Grid(ie,:)

            Do iv=1, 7 
                DUMMY5(1:LENSFC_proc, iv) = noahmp(ie)%snow_soil_interface(:,iv)
            Enddo
            Do iv=1, 3
                DUMMY6(1:LENSFC_proc, iv) = noahmp(ie)%temperature_snow(:,iv)
                DUMMY7(1:LENSFC_proc, iv) = noahmp(ie)%snow_ice_layer(:,iv)
                DUMMY8(1:LENSFC_proc, iv) = noahmp(ie)%snow_liq_layer(:,iv)
            Enddo
            DUMMY9(1:LENSFC_proc) = noahmp(ie)%temperature_soil
            Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1  
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(DUMMY1(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY2(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY3(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY4(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Do iv=1, 7 
                    call MPI_RECV(DUMMY5(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                400*(ixy+1) + 10*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo
                Do iv=1, 3
                    call MPI_RECV(DUMMY6(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                500*(ixy+1) + 20*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY7(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                600*(ixy+1) + 30*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY8(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                700*(ixy+1) + 40*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo
                call MPI_RECV(DUMMY9(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                800*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)

                call MPI_RECV(DUMMY12(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                1300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do

        ! REstart file
            WRITE(ensCH, '(I2.2)') ie
            restart_filename=TRIM(vector_restart_prefix)//"/ens"//ensCH// &
            "/ufs_land_restart."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//&
            "-00-00.nc" 

            INQUIRE(FILE=trim(restart_filename), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(restart_filename) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif
            
            error = nf90_open(trim(restart_filename), NF90_WRITE, ncid)
            call netcdf_err(error, 'opening restart file' )
        
            ! Start writing restart file
            error = nf90_inq_varid(ncid, "weasd", varid)
            call netcdf_err(error, 'getting varid weasd' )
            error = nf90_put_var(ncid, varid , DUMMY1   , &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing weasd' )

            error = nf90_inq_varid(ncid, "snwdph", varid)
            call netcdf_err(error, 'getting varid snwdph' )
            error = nf90_put_var(ncid, varid , DUMMY2,  &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing snwdph' )

            error = nf90_inq_varid(ncid, "snowxy", varid)
            call netcdf_err(error, 'getting varid snowxy' )
            error = nf90_put_var(ncid, varid , DUMMY3  , &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing snowxy' )

            error = nf90_inq_varid(ncid, "sneqvoxy", varid)
            call netcdf_err(error, 'getting varid sneqvoxy' )
            error = nf90_put_var(ncid, varid , DUMMY4, &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing sneqvoxy' )

            error = nf90_inq_varid(ncid, "zsnsoxy", varid)
            call netcdf_err(error, 'getting varid zsnsoxy' )
            error = nf90_put_var(ncid, varid, DUMMY5, &
                start = (/1, 1, 1/), count = (/lensfc, 7, 1/))
            call netcdf_err(error, 'writing zsnsoxy' )

            error = nf90_inq_varid(ncid, "tsnoxy", varid)
            call netcdf_err(error, 'getting varid tsnoxy' )
            error = nf90_put_var(ncid, varid, DUMMY6, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing tsnoxy' )

            error = nf90_inq_varid(ncid, "snicexy", varid)
            call netcdf_err(error, 'getting varid snicexy' )
            error = nf90_put_var(ncid, varid, DUMMY7, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing snicexy' )

            error = nf90_inq_varid(ncid, "snliqxy", varid)
            call netcdf_err(error, 'getting varid snliqxy' )
            error = nf90_put_var(ncid, varid, DUMMY8, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing snliqxy' )

            error = nf90_inq_varid(ncid, "stc", varid)
            call netcdf_err(error, 'getting varid stc' )
            error = nf90_put_var(ncid, varid, DUMMY9, &
                start = (/1, 1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing stc')

            error = nf90_inq_varid(ncid, "snowc", varid)
            call netcdf_err(error, 'getting varid snowc' )
            error = nf90_put_var(ncid, varid, DUMMY12, &
                start = (/1, 1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing snowc')

            error = nf90_close(ncid)
            call netcdf_err(error, 'closing restart file' )

        endif
    enddo
             
 End subroutine Write_DA_Outputs_ens_vect

 Subroutine Write_DA_Outputs_restart_vector(myrank, NPROCS, da_output_file, restart_filename, &
            lensfc, LENSFC_land, num_stn, num_Eval, max_num_nearStn, &            
            N_sA, N_sA_Ext, mpiReal_size, mpiInt_size, noahmp,    &    !idim, jdim, 
            snwdforc, SWEFCS, incr_atGrid, snwdanal, SWEANL, &      !SWEANL_Cur, SNDANL_Cur, landmask,  &
            Lat_atObs, Lon_atObs, OROG_at_stn, OBS_stn, FCS_at_stn, OmB_innov_at_stn,  & 
            index_back_atObs, index_obs_assmilated, &
            NEXC, Index_Obs_Excluded,   &
            SNOANL_atStn, SNCOV_IMS, IMS_Foot_Print, SCF_Grid, SNO_IMS_at_Grid, &
            Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
            SNOFCS_atEvalPts, incr_atEvalPts, SNOANL_atEvalPts, &
            SNDnoDA, SWEnoDA, SNDnoDA_atEval, SWEnoDA_atEval, index_back_atEval, &
            Orog_obs_atEvalPts, Orog_fcs_atEvalPts) !SND_ml_noDA_atEval)    !, SNOANL_Cur_atEvalPts)  !, anl_fSCA) !updated snocov
        !------------------------------------------------------------------
        ! Write DA ouputs: 
        ! forecast SWE
        ! analysis SWE
        ! analysis Snow Depth
        ! innovation at grid
        !------------------------------------------------------------------
        implicit none
        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)    :: da_output_file, restart_filename
        integer, intent(in)     :: myrank, NPROCS, lensfc, LENSFC_land, num_stn, num_Eval    
        integer, intent(in)     :: max_num_nearStn, N_sA, N_sA_Ext, mpiReal_size, mpiInt_size
        type(noahmp_type)       :: noahmp          
        real, intent(in)    :: snwdforc(lensfc), SWEFCS(lensfc), incr_atGrid(lensfc), &
                               snwdanal(lensfc), SWEANL(lensfc)
        ! real, intent(in)    :: SWEANL_Cur(lensfc), SNDANL_Cur(lensfc)
        Real, intent(in)    :: Lat_atObs(num_stn), Lon_atObs(num_stn), OROG_at_stn(num_stn)
        Real, intent(in)    :: OBS_stn(num_stn), FCS_at_stn(num_stn), &
                               OmB_innov_at_stn(num_stn), SNOANL_atStn(num_stn)
        integer, intent(in) :: index_back_atObs(num_stn)
        !1.4 track obserations assimilated
        Integer          :: index_obs_assmilated(max_num_nearStn+1, LENSFC_land)
        ! excluded obs
        Integer, intent(in) :: NEXC
        Integer, intent(in) :: Index_Obs_Excluded(NEXC)

        Real, intent(in)    :: SNCOV_IMS(lensfc), IMS_Foot_Print(lensfc)  !, anl_fSCA(lensfc)
        Real, intent(in)    :: SCF_Grid(lensfc), SNO_IMS_at_Grid(lensfc)
        real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), &
                               Obs_atEvalPts(num_Eval), Orog_obs_atEvalPts(num_Eval), Orog_fcs_atEvalPts(num_Eval)
        real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), incr_atEvalPts(num_Eval), &
                               SNOANL_atEvalPts(num_Eval), &
                               SNDnoDA(lensfc), SWEnoDA(lensfc), &
                               SNDnoDA_atEval(num_Eval), SWEnoDA_atEval(num_Eval)
                               !, &  SND_ml_noDA_atEval(num_Eval)
        integer, intent(in)    :: index_back_atEval(num_Eval)
        real                   :: SNDANL_ml_atEval(num_Eval)

        integer            :: fsize=65536, inital=0
        integer            :: header_buffer_val = 16384
        integer            :: dims_2d(2), dims_strt(2), dims_end(2)
        integer            :: error, i, ncid
        integer            :: dim_x, dim_time, dim_stn, dim_eval, dim_sn, dim_snso, dim_exc 
        integer            :: id_time, id_x, id_sn, id_snso 
        integer            :: dim_no, id_no
        integer       :: id_snwd_forc, id_swe_forc, id_imscov, id_imsftp, id_fsca, id_imssnd 
        integer       :: id_snwd_noda, id_swe_noda, id_snwd_noda_eval, id_swe_noda_eval
        integer       :: id_snwd_ml_eval, id_indxval, id_exc
        integer       :: id_incr_1, id_snwd_1, id_swe_1, id_anleval_1, id_increval_1 
      
        integer       :: id_latstn, id_lonstn, id_elestn, id_obsstn, id_forcstn, id_innovstn, id_anlstn !, id_landmask
        integer       :: id_lateval, id_loneval, id_ele_evl, id_ele_fcs_evl, id_obseval, id_forceval, id_indxobs
       
        real(kind=4)                :: times
        real(kind=4), allocatable   :: x_data(:)  !, y_data(:)
        integer                     :: dims_3d(3), iv, varid
        integer     :: id_weasd, id_sndpth_m, id_snowxy, id_sneqvoxy, id_stc, &
                       id_zsnsoxy, id_tsnoxy, id_snicexy, id_snliqxy

        Real        :: DUMMY1(lensfc), DUMMY2(lensfc), DUMMY3(lensfc), &
                       DUMMY4(lensfc), DUMMY5(lensfc, 7), DUMMY6(lensfc, 3), &
                       DUMMY7(lensfc, 3), DUMMY8(lensfc, 3), DUMMY9(lensfc)

        integer :: DUMMY13(LENSFC)

        Integer     :: dest_Aoffset, dest_Aoffset_end, arLen, IERR, ixy, jndx, ie

        LOGICAL     :: file_exists

        if (MYRANK > 0) then
            call MPI_SEND(noahmp%swe, LENSFC_land, mpiReal_size, 0,   &
                                    MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(noahmp%snow_depth, LENSFC_land, mpiReal_size, 0,   &
                                    100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(noahmp%active_snow_layers, LENSFC_land, mpiReal_size, 0,   &
                                    200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(noahmp%swe_previous, LENSFC_land, mpiReal_size, 0,   &
                                    300*(MYRANK+1), MPI_COMM_WORLD, IERR)
            Do iv=1, 7                        
                call MPI_SEND(noahmp%snow_soil_interface(:,iv), LENSFC_land, mpiReal_size, 0,   &
                                400*(MYRANK+1) + 10*iv, MPI_COMM_WORLD, IERR)
            Enddo
            Do iv=1, 3
                call MPI_SEND(noahmp%temperature_snow(:,iv), LENSFC_land, mpiReal_size, 0,   &
                                        500*(MYRANK+1) + 20*iv, MPI_COMM_WORLD, IERR)
                call MPI_SEND(noahmp%snow_ice_layer(:,iv), LENSFC_land, mpiReal_size, 0,   &
                                        600*(MYRANK+1) + 30*iv, MPI_COMM_WORLD, IERR)
                call MPI_SEND(noahmp%snow_liq_layer(:,iv), LENSFC_land, mpiReal_size, 0,   &
                                        700*(MYRANK+1) + 40*iv, MPI_COMM_WORLD, IERR)
            enddo
            call MPI_SEND(noahmp%temperature_soil, LENSFC_land, mpiReal_size, 0,   &
                                        800*(MYRANK+1), MPI_COMM_WORLD, IERR)
        Endif
        if (MYRANK == 0) then
            DUMMY1(1:LENSFC_land) = noahmp%swe
            DUMMY2(1:LENSFC_land) = noahmp%snow_depth
            DUMMY3(1:LENSFC_land) = noahmp%active_snow_layers
            DUMMY4(1:LENSFC_land) = noahmp%swe_previous
            Do iv=1, 7 
                DUMMY5(1:LENSFC_land, iv) = noahmp%snow_soil_interface(:,iv)
            Enddo
            Do iv=1, 3
                DUMMY6(1:LENSFC_land, iv) = noahmp%temperature_snow(:,iv)
                DUMMY7(1:LENSFC_land, iv) = noahmp%snow_ice_layer(:,iv)
                DUMMY8(1:LENSFC_land, iv) = noahmp%snow_liq_layer(:,iv)
            Enddo
            DUMMY9(1:LENSFC_land) = noahmp%temperature_soil
            Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                    dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                    dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                    arLen = N_sA + 1  
                else
                    dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                    dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                    arLen = N_sA
                endif
                call MPI_RECV(DUMMY1(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY2(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY3(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(DUMMY4(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Do iv=1, 7 
                    call MPI_RECV(DUMMY5(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                400*(ixy+1) + 10*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo
                Do iv=1, 3
                    call MPI_RECV(DUMMY6(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                500*(ixy+1) + 20*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY7(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                600*(ixy+1) + 30*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY8(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                700*(ixy+1) + 40*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo
                call MPI_RECV(DUMMY9(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, ixy,      &
                                800*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        ! endif
            SNDANL_ml_atEval = IEEE_VALUE(SNDANL_ml_atEval, IEEE_QUIET_NAN)        
            Do jndx = 1, num_Eval
                if (index_back_atEval(jndx) > 0)  then   !                     
                    SNDANL_ml_atEval(jndx) = DUMMY2(index_back_atEval(jndx))
                endif
            enddo
        ! print*,"Process ", myrank, "writing output data to: ",trim(output_file)
        !--- create the file
        ! if (myrank ==) then 
            error = NF90_CREATE(trim(da_output_file), IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
            call netcdf_err(error, 'CREATING FILE='//trim(da_output_file) )     
            ! print*,"Process ", myrank, "created file: ",trim(output_file)
            !--- define dimensions
            error = nf90_def_dim(ncid, 'location', lensfc, dim_x)
            call netcdf_err(error, 'DEFINING location DIMENSION' )
            ! error = nf90_def_dim(ncid, 'yaxis_1', jdim, dim_y)
            ! call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
            error = nf90_def_dim(ncid, 'Time', 1, dim_time)
            call netcdf_err(error, 'DEFINING TIME DIMENSION' )
            ! obs points
            error = nf90_def_dim(ncid, 'obs_points', num_stn, dim_stn)
            call netcdf_err(error, 'DEFINING obs_points' )
            ! assimilated obs            
            error = nf90_def_dim(ncid, 'nearest_obs', max_num_nearStn+1, dim_no)
            call netcdf_err(error, 'DEFINING nearest_obs DIMENSION' )
            ! eval obs points        
            if (num_Eval>0) then
                    error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
                    call netcdf_err(error, 'DEFINING eval_points' )
            else
                    error = nf90_def_dim(ncid, 'eval_points', 1, dim_eval)
                    call netcdf_err(error, 'DEFINING eval_points' )
            endif
            !excluded obs indx 
            error = nf90_def_dim(ncid, 'obs_excluded', NEXC, dim_exc)
            call netcdf_err(error, 'DEFINING obs_excluded')

            error = nf90_def_dim(ncid, 'snow_levels', 3, dim_sn)
            call netcdf_err(error, 'DEFINING snow_levels DIMENSION' )
            error = nf90_def_dim(ncid, 'snso_levels', 7, dim_snso)
            call netcdf_err(error, 'DEFINING snso_levels DIMENSION' )

            !--- define fields
            error = nf90_def_var(ncid, 'location', NF90_FLOAT, dim_x, id_x)
            call netcdf_err(error, 'DEFINING location FIELD' )
            error = nf90_put_att(ncid, id_x, "long_name", "location")
            call netcdf_err(error, 'DEFINING location LONG NAME' )
            error = nf90_put_att(ncid, id_x, "units", "none")
            call netcdf_err(error, 'DEFINING location UNITS' )
            error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
            call netcdf_err(error, 'WRITING location FIELD' )

            error = nf90_def_var(ncid, 'Time', NF90_FLOAT, dim_time, id_time)
            call netcdf_err(error, 'DEFINING TIME FIELD' )
            error = nf90_put_att(ncid, id_time, "long_name", "Time")
            call netcdf_err(error, 'DEFINING TIME LONG NAME' )
            error = nf90_put_att(ncid, id_time, "units", "time level")
            call netcdf_err(error, 'DEFINING TIME UNITS' )
            error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
            call netcdf_err(error, 'WRITING TIME FIELD' )

            error = nf90_def_var(ncid, 'snow_levels', NF90_FLOAT, dim_sn, id_sn)
            call netcdf_err(error, 'DEFINING snow_levels FIELD' )
            error = nf90_put_att(ncid, id_sn, "long_name", "snow_levels")
            call netcdf_err(error, 'DEFINING snow_levels LONG NAME' )
            error = nf90_put_att(ncid, id_sn, "units", "none")
            call netcdf_err(error, 'DEFINING snow_levels UNITS' )

            error = nf90_def_var(ncid, 'snso_levels', NF90_FLOAT, dim_snso, id_snso)
            call netcdf_err(error, 'DEFINING snso_levels FIELD' )
            error = nf90_put_att(ncid, id_snso, "long_name", "snso_levels")
            call netcdf_err(error, 'DEFINING snso_levels LONG NAME' )
            error = nf90_put_att(ncid, id_snso, "units", "none")
            call netcdf_err(error, 'DEFINING snso_levels UNITS' )

            dims_2d(1) = dim_no
            dims_2d(2) = dim_x
            error = nf90_def_var(ncid, 'index_obs_assmilated', NF90_INT, dims_2d, id_no) 
            call netcdf_err(error, 'DEFINING index_obs_assmilated' )
            error = nf90_put_att(ncid, id_no, "long_name", "Indics of assimilated Observation Points")
            call netcdf_err(error, 'DEFINING index_obs_assmilated LONG NAME' )
            error = nf90_put_att(ncid, id_no, "units", "-")
            call netcdf_err(error, 'DEFINING index_obs_assmilated UNITS' )

            dims_2d(1) = dim_x
            ! dims_3d(2) = dim_y
            dims_2d(2) = dim_time

            error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_2d, id_swe_forc)
            call netcdf_err(error, 'DEFINING SWE_Forecast' )
            error = nf90_put_att(ncid, id_swe_forc, "long_name", "Forecast Snow Water Equivalent")
            call netcdf_err(error, 'DEFINING SWE Forecast LONG NAME' )
            error = nf90_put_att(ncid, id_swe_forc, "units", "mm")
            call netcdf_err(error, 'DEFINING SWE Forecast UNITS' )

            error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_2d, id_snwd_forc)
            call netcdf_err(error, 'DEFINING SND Forecast' )
            error = nf90_put_att(ncid, id_snwd_forc, "long_name", "Forecast Snow Depth")
            call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
            error = nf90_put_att(ncid, id_snwd_forc, "units", "mm")
            call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

            error = nf90_def_var(ncid, 'SWE_Analysis', NF90_DOUBLE, dims_2d, id_swe_1)
            call netcdf_err(error, 'DEFINING SWE_Analysis' )
            error = nf90_put_att(ncid, id_swe_1, "long_name", "Analysis Snow Water Equivalent")
            call netcdf_err(error, 'DEFINING SWE LONG NAME' )
            error = nf90_put_att(ncid, id_swe_1, "units", "mm")
            call netcdf_err(error, 'DEFINING SWE UNITS' )

            error = nf90_def_var(ncid, 'SND_Analysis', NF90_DOUBLE, dims_2d, id_snwd_1)
            call netcdf_err(error, 'DEFINING SND Analyis' )
            error = nf90_put_att(ncid, id_snwd_1, "long_name", "Analysis Snow Depth")
            call netcdf_err(error, 'DEFINING SND Analysis LONG NAME' )
            error = nf90_put_att(ncid, id_snwd_1, "units", "mm")
            call netcdf_err(error, 'DEFINING SND Analysis UNITS' )
            !
            error = nf90_def_var(ncid, 'DA_Increment', NF90_DOUBLE, dims_2d, id_incr_1)
            call netcdf_err(error, 'DEFINING DA_Increment' )
            error = nf90_put_att(ncid, id_incr_1, "long_name", "DA_Increment at model grid")
            call netcdf_err(error, 'DEFINING DA_Increment LONG NAME' )
            error = nf90_put_att(ncid, id_incr_1, "units", "mm")
            call netcdf_err(error, 'DEFINING DA_Increment UNITS' )
            !
            error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_2d, id_imscov)
            call netcdf_err(error, 'DEFINING imsfSCA' )
            error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
            call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
            error = nf90_put_att(ncid, id_imscov, "units", "-")
            call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

            error = nf90_def_var(ncid, 'ims_FootPrint', NF90_DOUBLE, dims_2d, id_imsftp)
            call netcdf_err(error, 'DEFINING ims_FootPrint' )
            error = nf90_put_att(ncid, id_imsftp, "long_name", "Assimilated IMS Data Foot Print")
            call netcdf_err(error, 'DEFINING ims_FootPrint LONG NAME' )
            error = nf90_put_att(ncid, id_imsftp, "units", "-")
            call netcdf_err(error, 'DEFINING ims_FootPrint UNITS' )

            error = nf90_def_var(ncid, 'fSCA', NF90_DOUBLE, dims_2d, id_fsca)
            call netcdf_err(error, 'DEFINING fSCA' )
            error = nf90_put_att(ncid, id_fsca, "long_name", "fractional Snow Covered Area")
            call netcdf_err(error, 'DEFINING fSCA LONG NAME' )
            error = nf90_put_att(ncid, id_fsca, "units", "-")
            call netcdf_err(error, 'DEFINING fSCA UNITS' )

            error = nf90_def_var(ncid, 'imsSND', NF90_DOUBLE, dims_2d, id_imssnd)
            call netcdf_err(error, 'DEFINING imsSND' )
            error = nf90_put_att(ncid, id_imssnd, "long_name", "IMS Snow Depth")
            call netcdf_err(error, 'DEFINING imsSND LONG NAME' )
            error = nf90_put_att(ncid, id_imssnd, "units", "mm")
            call netcdf_err(error, 'DEFINING imsSND UNITS' )
 
            ! No DA 
            error = nf90_def_var(ncid, 'SWE_Forecast_NoDA', NF90_DOUBLE, dims_2d, id_swe_noda)
            call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA' )
            error = nf90_put_att(ncid, id_swe_noda, "long_name", "Forecast Snow Water Equivalent No DA")
            call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA LONG NAME' )
            error = nf90_put_att(ncid, id_swe_noda, "units", "mm")
            call netcdf_err(error, 'DEFINING SWE_Forecast_NoDA UNITS' )

            error = nf90_def_var(ncid, 'SND_Forecast_NoDA', NF90_DOUBLE, dims_2d, id_snwd_noda)
            call netcdf_err(error, 'DEFINING SND_Forecast_NoDA' )
            error = nf90_put_att(ncid, id_snwd_noda, "long_name", "Forecast Snow Depth No DA")
            call netcdf_err(error, 'DEFINING SND_Forecast_NoDA LONG NAME' )
            error = nf90_put_att(ncid, id_snwd_noda, "units", "mm")
            call netcdf_err(error, 'DEFINING SND_Forecast_NoDA UNITS' )

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

            error = nf90_def_var(ncid, 'elevation@MetaData', NF90_DOUBLE, dim_stn, id_elestn)
            call netcdf_err(error, 'DEFINING elevation@MetaData' )
            error = nf90_put_att(ncid, id_elestn, "long_name", "elevation at Observation Points")
            call netcdf_err(error, 'DEFINING elevation@MetaData LONG NAME' )
            error = nf90_put_att(ncid, id_elestn, "units", "m")
            call netcdf_err(error, 'DEFINING elevation@MetaData UNITS' )
            
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

            error = nf90_def_var(ncid, 'snwdph@omb', NF90_DOUBLE, dim_stn, id_innovstn)
            call netcdf_err(error, 'DEFINING snwdph@omb' )
            error = nf90_put_att(ncid, id_innovstn, "long_name", "Innovation(O-B) at Observation Points")
            call netcdf_err(error, 'DEFINING snwdph@omb LONG NAME' )
            error = nf90_put_att(ncid, id_innovstn, "units", "mm")
            call netcdf_err(error, 'DEFINING snwdph@omb UNITS' )

            error = nf90_def_var(ncid, 'snwdph@anl', NF90_DOUBLE, dim_stn, id_anlstn)
            call netcdf_err(error, 'DEFINING snwdph@anl' )
            error = nf90_put_att(ncid, id_anlstn, "long_name", "Analysis at Observation Points")
            call netcdf_err(error, 'DEFINING snwdph@anl LONG NAME' )
            error = nf90_put_att(ncid, id_anlstn, "units", "mm")
            call netcdf_err(error, 'DEFINING snwdph@anl UNITS' )

            error = nf90_def_var(ncid, 'index_back_at_obs', NF90_INT, dim_stn, id_indxobs) 
            call netcdf_err(error, 'DEFINING index_back_at_obs' )
            error = nf90_put_att(ncid, id_indxobs, "long_name", "Index of background at Observation Points")
            call netcdf_err(error, 'DEFINING index_back_at_obs LONG NAME' )
            error = nf90_put_att(ncid, id_indxobs, "units", "-")
            call netcdf_err(error, 'DEFINING index_back_at_obs UNITS' )

            ! error = nf90_def_var(ncid, 'anlfSCA', NF90_DOUBLE, dims_3d, id_anlscov)
            ! call netcdf_err(error, 'DEFINING anlfSCA' )
            ! error = nf90_put_att(ncid, id_anlscov, "long_name", "Analysis fractional Snow Covered Area")
            ! call netcdf_err(error, 'DEFINING anlfSCA LONG NAME' )
            ! error = nf90_put_att(ncid, id_anlscov, "units", "-")
            ! call netcdf_err(error, 'DEFINING anlfSCA UNITS' )

            ! eval points
    ! if (num_Eval>0) then 
            error = nf90_def_var(ncid, 'LatEvalPoints', NF90_DOUBLE, dim_eval, id_lateval)
            call netcdf_err(error, 'DEFINING LatEvalPoints' )
            error = nf90_put_att(ncid, id_lateval, "long_name", "Latitude at Evaluation Points")
            call netcdf_err(error, 'DEFINING LatEvalPoints LONG NAME' )
            error = nf90_put_att(ncid, id_lateval, "units", "deg")
            call netcdf_err(error, 'DEFINING LatEvalPoints UNITS' )

            error = nf90_def_var(ncid, 'LonEvalPoints', NF90_DOUBLE, dim_eval, id_loneval)
            call netcdf_err(error, 'DEFINING LonEvalPoints' )
            error = nf90_put_att(ncid, id_loneval, "long_name", "Longitude at Evaluation Points")
            call netcdf_err(error, 'DEFINING LonEvalPoints LONG NAME' )
            error = nf90_put_att(ncid, id_loneval, "units", "deg")
            call netcdf_err(error, 'DEFINING LonEvalPoints UNITS' ) 

            error = nf90_def_var(ncid, 'ElevationEvalPoints', NF90_DOUBLE, dim_eval, id_ele_evl)
            call netcdf_err(error, 'DEFINING ElevationEvalPoints' )
            error = nf90_put_att(ncid, id_ele_evl, "long_name", "Elevation at Evaluation Points")
            call netcdf_err(error, 'DEFINING ElevationEvalPoints LONG NAME' )
            error = nf90_put_att(ncid, id_ele_evl, "units", "m")
            call netcdf_err(error, 'DEFINING ElevationEvalPoints UNITS' )   
            
            error = nf90_def_var(ncid, 'ElevationFCSEvalPoints', NF90_DOUBLE, dim_eval, id_ele_fcs_evl)
            call netcdf_err(error, 'DEFINING ElevationFCSEvalPoints' )
            error = nf90_put_att(ncid, id_ele_fcs_evl, "long_name", "Forecast Elevation at Evaluation Points")
            call netcdf_err(error, 'DEFINING ElevationFCSEvalPoints LONG NAME' )
            error = nf90_put_att(ncid, id_ele_fcs_evl, "units", "m")
            call netcdf_err(error, 'DEFINING ElevationFCSEvalPoints UNITS' ) 

            error = nf90_def_var(ncid, 'Obs_atEvalPts', NF90_DOUBLE, dim_eval, id_obseval)
            call netcdf_err(error, 'DEFINING Obs_atEvalPts' )
            error = nf90_put_att(ncid, id_obseval, "long_name", "Observed at Evaluation Points")
            call netcdf_err(error, 'DEFINING Obs_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_obseval, "units", "mm")
            call netcdf_err(error, 'DEFINING Obs_atEvalPts UNITS' )
            
            error = nf90_def_var(ncid, 'SNOFCS_atEvalPts', NF90_DOUBLE, dim_eval, id_forceval)
            call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts' )
            error = nf90_put_att(ncid, id_forceval, "long_name", "Forecast at Evaluation Points")
            call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_forceval, "units", "mm")
            call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts UNITS' )

            error = nf90_def_var(ncid, 'SNOANL_atEvalPts', NF90_DOUBLE, dim_eval, id_anleval_1)
            call netcdf_err(error, 'DEFINING SNOANL_atEvalPts' )
            error = nf90_put_att(ncid, id_anleval_1, "long_name", "Analysis at Evaluation Points")
            call netcdf_err(error, 'DEFINING SNOANL_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_anleval_1, "units", "mm")
            call netcdf_err(error, 'DEFINING SNOANL_atEvalPts UNITS' )
            !            
            error = nf90_def_var(ncid, 'incr_atEvalPts', NF90_DOUBLE, dim_eval, id_increval_1)
            call netcdf_err(error, 'DEFINING incr_atEvalPts' )
            error = nf90_put_att(ncid, id_increval_1, "long_name", "Increment at Evaluation Points")
            call netcdf_err(error, 'DEFINING incr_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_increval_1, "units", "mm")
            call netcdf_err(error, 'DEFINING incr_atEvalPts UNITS' )
            !
            ! multi layer
            error = nf90_def_var(ncid, 'SNDANL_ML_atEvalPts', NF90_DOUBLE, dim_eval, id_snwd_ml_eval)
            call netcdf_err(error, 'DEFINING SNDANL_ML_atEvalPts' )
            error = nf90_put_att(ncid, id_snwd_ml_eval, "long_name", "Multilayer Analysis at Evaluation Points")
            call netcdf_err(error, 'DEFINING SNDANL_ML_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_snwd_ml_eval, "units", "mm")
            call netcdf_err(error, 'DEFINING SNDANL_ML_atEvalPts UNITS' )
            ! no da
            error = nf90_def_var(ncid, 'SNDNoDA_atEvalPts', NF90_DOUBLE, dim_eval, id_snwd_noda_eval)
            call netcdf_err(error, 'DEFINING SNDNoDA_atEvalPts' )
            error = nf90_put_att(ncid, id_snwd_noda_eval, "long_name", "Forecast at Evaluation Points No DA")
            call netcdf_err(error, 'DEFINING SNDNoDA_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_snwd_noda_eval, "units", "mm")
            call netcdf_err(error, 'DEFINING SNDNoDA_atEvalPts UNITS' )

            error = nf90_def_var(ncid, 'SWENoDA_atEvalPts', NF90_DOUBLE, dim_eval, id_swe_noda_eval)
            call netcdf_err(error, 'DEFINING SWENoDA_atEvalPts' )
            error = nf90_put_att(ncid, id_swe_noda_eval, "long_name", "SWE Forecast at Evaluation Points No DA")
            call netcdf_err(error, 'DEFINING SWENoDA_atEvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_swe_noda_eval, "units", "mm")
            call netcdf_err(error, 'DEFINING SWENoDA_atEvalPts UNITS' )

            ! index of eval points
            error = nf90_def_var(ncid, 'Index_Back_EvalPts', NF90_INT, dim_eval, id_indxval)
            call netcdf_err(error, 'DEFINING Index_Back_EvalPts' )
            error = nf90_put_att(ncid, id_indxval, "long_name", "Index of Background at Evaluation Points")
            call netcdf_err(error, 'DEFINING Index_Back_EvalPts LONG NAME' )
            error = nf90_put_att(ncid, id_indxval, "units", "-")
            call netcdf_err(error, 'DEFINING Index_Back_EvalPts UNITS' )
        ! endif eval
            ! index of Obs excluded
            error = nf90_def_var(ncid, 'Index_Obs_Excluded', NF90_INT, dim_exc, id_exc)
            call netcdf_err(error, 'DEFINING Index_Obs_Excluded' )
            error = nf90_put_att(ncid, id_exc, "long_name", "Index of Observations Excluded")
            call netcdf_err(error, 'DEFINING Index_Obs_Excluded LONG NAME' )
            error = nf90_put_att(ncid, id_exc, "units", "-")
            call netcdf_err(error, 'DEFINING Index_Obs_Excluded UNITS' )
            
            error = nf90_def_var(ncid, 'weasd', NF90_DOUBLE, dims_2d, id_weasd)
            call netcdf_err(error, 'DEFINING weasd' )
            error = nf90_put_att(ncid, id_weasd, "long_name", "weasd")
            call netcdf_err(error, 'DEFINING weasd LONG NAME' )
            error = nf90_put_att(ncid, id_weasd, "units", "mm")
            call netcdf_err(error, 'DEFINING weasd UNITS' )

            error = nf90_def_var(ncid, 'snwdph', NF90_DOUBLE, dims_2d, id_sndpth_m)
            call netcdf_err(error, 'DEFINING sndpth' )
            error = nf90_put_att(ncid, id_sndpth_m, "long_name", "sndpth multi layer")
            call netcdf_err(error, 'DEFINING sndpth LONG NAME' )
            error = nf90_put_att(ncid, id_sndpth_m, "units", "mm")
            call netcdf_err(error, 'DEFINING sndpth UNITS' )
            
            error = nf90_def_var(ncid, 'snowxy', NF90_DOUBLE, dims_2d, id_snowxy)
            call netcdf_err(error, 'DEFINING snowxy' )
            error = nf90_put_att(ncid, id_snowxy, "long_name", "snowxy")
            call netcdf_err(error, 'DEFINING snowxy LONG NAME' )
            error = nf90_put_att(ncid, id_snowxy, "units", "mm")
            call netcdf_err(error, 'DEFINING snowxy UNITS' )

            error = nf90_def_var(ncid, 'sneqvoxy', NF90_DOUBLE, dims_2d, id_sneqvoxy)
            call netcdf_err(error, 'DEFINING sneqvoxy' )
            error = nf90_put_att(ncid, id_sneqvoxy, "long_name", "sneqvoxy")
            call netcdf_err(error, 'DEFINING sneqvoxy LONG NAME' )
            error = nf90_put_att(ncid, id_sneqvoxy, "units", "mm")
            call netcdf_err(error, 'DEFINING sneqvoxy UNITS' )

            dims_3d(1) = dim_x
            dims_3d(2) = dim_sn
            dims_3d(3) = dim_time
            error = nf90_def_var(ncid, 'tsnoxy', NF90_DOUBLE, dims_3d, id_tsnoxy)
            call netcdf_err(error, 'DEFINING tsnoxy' )
            error = nf90_put_att(ncid, id_tsnoxy, "long_name", "tsnoxy")
            call netcdf_err(error, 'DEFINING tsnoxy LONG NAME' )
            error = nf90_put_att(ncid, id_tsnoxy, "units", "mm")
            call netcdf_err(error, 'DEFINING tsnoxy UNITS' )
            
            error = nf90_def_var(ncid, 'snicexy', NF90_DOUBLE, dims_3d, id_snicexy)
            call netcdf_err(error, 'DEFINING snicexy' )
            error = nf90_put_att(ncid, id_snicexy, "long_name", "snicexy")
            call netcdf_err(error, 'DEFINING snicexy LONG NAME' )
            error = nf90_put_att(ncid, id_snicexy, "units", "mm")
            call netcdf_err(error, 'DEFINING snicexy UNITS' )
            
            error = nf90_def_var(ncid, 'snliqxy', NF90_DOUBLE, dims_3d, id_snliqxy)
            call netcdf_err(error, 'DEFINING snliqxy' )
            error = nf90_put_att(ncid, id_snliqxy, "long_name", "snliqxy")
            call netcdf_err(error, 'DEFINING snliqxy LONG NAME' )
            error = nf90_put_att(ncid, id_snliqxy, "units", "mm")
            call netcdf_err(error, 'DEFINING snliqxy UNITS' )

            dims_3d(1) = dim_x
            dims_3d(2) = dim_snso
            dims_3d(3) = dim_time
            error = nf90_def_var(ncid, 'zsnsoxy', NF90_DOUBLE, dims_3d, id_zsnsoxy)
            call netcdf_err(error, 'DEFINING zsnsoxy' )
            error = nf90_put_att(ncid, id_zsnsoxy, "long_name", "zsnsoxy")
            call netcdf_err(error, 'DEFINING zsnsoxy LONG NAME' )
            error = nf90_put_att(ncid, id_zsnsoxy, "units", "mm")
            call netcdf_err(error, 'DEFINING zsnsoxy UNITS' )
            
            error = nf90_def_var(ncid, 'stc', NF90_DOUBLE, dims_2d, id_stc)
            call netcdf_err(error, 'DEFINING stc' )
            error = nf90_put_att(ncid, id_stc, "long_name", "stc")
            call netcdf_err(error, 'DEFINING stc LONG NAME' )
            error = nf90_put_att(ncid, id_stc, "units", "-")
            call netcdf_err(error, 'DEFINING stc UNITS' )

            error = nf90_enddef(ncid, header_buffer_val,4,0,4)
            call netcdf_err(error, 'DEFINING HEADER' )

            allocate(x_data(lensfc))
            do i = 1, lensfc
            x_data(i) = float(i)
            enddo
            ! allocate(y_data(jdim))
            ! do i = 1, jdim
            ! y_data(i) = float(i)
            ! enddo
            times = 1.0

            error = nf90_put_var( ncid, id_x, x_data)
            call netcdf_err(error, 'WRITING location RECORD' )
            ! error = nf90_put_var( ncid, id_y, y_data)
            ! call netcdf_err(error, 'WRITING YAXIS RECORD' )
            error = nf90_put_var( ncid, id_time, times)
            call netcdf_err(error, 'WRITING TIME RECORD' )

            deallocate(x_data)
            allocate(x_data(3))
            do i = 1, 3
                x_data(i) = float(i)
            enddo
            error = nf90_put_var( ncid, id_sn, x_data)
            call netcdf_err(error, 'WRITING snow_levels RECORD' )
            deallocate(x_data)
            allocate(x_data(7))
            do i = 1, 7
                x_data(i) = float(i)
            enddo            
            error = nf90_put_var( ncid, id_snso, x_data)
            call netcdf_err(error, 'WRITING snow_soil_levels RECORD' )

            ! allocate(dum2d(idim,jdim))
            dims_strt(1:2) = 1
            dims_end(1) = lensfc
            dims_end(2) = 1
            
            error = nf90_put_var( ncid, id_snwd_forc, snwdforc, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SND Forecast RECORD' )
            ! dum2d = reshape(SWEFCS, (/idim,jdim/))
            error = nf90_put_var( ncid, id_swe_forc, SWEFCS, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SWE Forecast RECORD' )
            
            error = nf90_put_var( ncid, id_snwd_1, snwdanal, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SND analysis RECORD' )
            ! dum2d = reshape(SWEANL(1,:), (/idim,jdim/))
            error = nf90_put_var( ncid, id_swe_1, SWEANL, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SWE Analysis RECORD' ) 
            
            error = nf90_put_var(ncid, id_incr_1, incr_atGrid, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING incr_atGrid RECORD')

            error = nf90_put_var( ncid, id_imscov, SNCOV_IMS, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING imsfSCA RECORD' )
            error = nf90_put_var( ncid, id_imsftp, IMS_Foot_Print, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING ims foot print RECORD' )
            
            error = nf90_put_var( ncid, id_fsca, SCF_Grid, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING fSCA RECORD' )

            error = nf90_put_var( ncid, id_imssnd, SNO_IMS_at_Grid, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING imsSND RECORD' )
            ! dum2d = reshape(anl_fSCA, (/idim,jdim/))
            ! error = nf90_put_var( ncid, id_anlscov, dum2d, dims_strt, dims_end)
            ! call netcdf_err(error, 'WRITING anlfSCA RECORD' )

            !no da 
            error = nf90_put_var( ncid, id_snwd_noda, SNDnoDA, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SND Forecast NoDA RECORD' )
            ! dum2d = reshape(SWEFCS, (/idim,jdim/))
            error = nf90_put_var( ncid, id_swe_noda, SWEnoDA, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING SWE Forecast NoDA RECORD' )
            
            ! obs points (obs, hofx, omb) 
            error = nf90_put_var( ncid, id_latstn, Lat_atObs)
            call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

            error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
            call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' ) 

            error = nf90_put_var( ncid, id_elestn, OROG_at_stn)
            call netcdf_err(error, 'WRITING OROG_at_stn RECORD' )

            error = nf90_put_var( ncid, id_obsstn, OBS_stn)
            call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

            error = nf90_put_var( ncid, id_forcstn, FCS_at_stn)
            call netcdf_err(error, 'WRITING SNOFCS_atObsPts RECORD' )

            error = nf90_put_var( ncid, id_innovstn, OmB_innov_at_stn)
            call netcdf_err(error, 'WRITING innov_atObsPts RECORD' )

            error = nf90_put_var( ncid, id_anlstn, SNOANL_atStn)
            call netcdf_err(error, 'WRITING SNOANL_atStn RECORD' )

            error = nf90_put_var( ncid, id_indxobs, index_back_atObs)
            call netcdf_err(error, 'WRITING indx_atobs RECORD')

            ! eval points
            if (num_Eval>0) then 
                    error = nf90_put_var( ncid, id_lateval, Lat_atEvalPts)
                    call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_loneval, Lon_atEvalPts)
                    call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' ) 

                    error = nf90_put_var( ncid, id_ele_evl, Orog_obs_atEvalPts)
                    call netcdf_err(error, 'WRITING Orog_obs_atEvalPts RECORD' ) 
                   
                    error = nf90_put_var( ncid, id_ele_fcs_evl, Orog_fcs_atEvalPts)
                    call netcdf_err(error, 'WRITING Orog_fcs_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_obseval, Obs_atEvalPts)
                    call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_forceval, SNOFCS_atEvalPts)
                    call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_anleval_1, SNOANL_atEvalPts)
                    call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_increval_1, incr_atEvalPts)
                    call netcdf_err(error, 'WRITING incr_atEvalPts RECORD')
                    
                    error = nf90_put_var( ncid, id_snwd_ml_eval, SNDANL_ml_atEval)
                    call netcdf_err(error, 'WRITING SNDANL_ml_atEval RECORD' )

                    !no da 
                    error = nf90_put_var( ncid, id_snwd_noda_eval, SNDnoDA_atEval)
                    call netcdf_err(error, 'WRITING SNDnoDA_atEval RECORD' )
                    error = nf90_put_var( ncid, id_swe_noda_eval, SWEnoDA_atEval)
                    call netcdf_err(error, 'WRITING SWEnoDA_atEval RECORD' )

                    !eval indx 
                    error = nf90_put_var( ncid, id_indxval, index_back_atEval)
                    call netcdf_err(error, 'WRITING index_back_atEval RECORD' )

            else
                    error = nf90_put_var( ncid, id_lateval, 0.0)
                    call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_loneval, 0.0)
                    call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )  

                    error = nf90_put_var( ncid, id_ele_evl, 0.0)
                    call netcdf_err(error, 'WRITING Orog_obs_atEvalPts RECORD' )  

                    error = nf90_put_var( ncid, id_ele_fcs_evl, 0.0)
                    call netcdf_err(error, 'WRITING Orog_fcs_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_obseval, 0.0)
                    call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_forceval, 0.0)
                    call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )

                    error = nf90_put_var( ncid, id_anleval_1, 0.0)
                    call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
                   
                    error = nf90_put_var( ncid, id_increval_1, 0.0)
                    call netcdf_err(error, 'WRITING incr_atEvalPts RECORD' )
         
                    error = nf90_put_var( ncid, id_snwd_ml_eval, 0.0)
                    call netcdf_err(error, 'WRITING SNDANL_ml_atEval RECORD')
                    !no da 
                    error = nf90_put_var( ncid, id_snwd_noda_eval, 0.0)
                    call netcdf_err(error, 'WRITING SNDnoDA_atEval RECORD')
                    error = nf90_put_var( ncid, id_swe_noda_eval, 0.0)
                    call netcdf_err(error, 'WRITING SWEnoDA_atEval RECORD')
                    !eval indx 
                    error = nf90_put_var( ncid, id_indxval, -1)
                    call netcdf_err(error, 'WRITING index_back_atEval RECORD' )
            endif
            error = nf90_put_var(ncid, id_exc, Index_Obs_Excluded)
            call netcdf_err(error, 'WRITING Index_Obs_Excluded RECORD' )

            error = nf90_put_var( ncid, id_weasd, DUMMY1, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_weasd RECORD' )

            error = nf90_put_var( ncid, id_sndpth_m, DUMMY2, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_sndpth_m RECORD' )

            error = nf90_put_var( ncid, id_snowxy, DUMMY3, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_snowxy RECORD' )

            error = nf90_put_var( ncid, id_sneqvoxy, DUMMY4, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_sneqvoxy RECORD' )

            error = nf90_put_var( ncid, id_sneqvoxy, DUMMY9, dims_strt, dims_end)
            call netcdf_err(error, 'WRITING id_stc RECORD' )
            
            error = nf90_put_var( ncid, id_zsnsoxy, DUMMY5, &
                start = (/1, 1, 1/), count = (/lensfc, 7, 1/))
            call netcdf_err(error, 'WRITING id_zsnsoxy RECORD' )

            error = nf90_put_var( ncid, id_tsnoxy, DUMMY6, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'WRITING id_zsnsoxy RECORD' )

            error = nf90_put_var( ncid, id_snicexy, DUMMY7, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'WRITING id_snicexy RECORD' )

            error = nf90_put_var( ncid, id_snliqxy, DUMMY8, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'WRITING id_snicexy RECORD' )

            deallocate(x_data)
            ! deallocate(dum2d)

            error = nf90_close(ncid)
            call netcdf_err(error, 'closing DA file' )

    ! REstart file
            INQUIRE(FILE=trim(restart_filename), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(restart_filename) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif
            error = nf90_open(trim(restart_filename), NF90_WRITE, ncid)
            call netcdf_err(error, 'opening restart file' )

            ! Start writing restart file
            error = nf90_inq_varid(ncid, "weasd", varid)
            call netcdf_err(error, 'getting varid weasd' )
            error = nf90_put_var(ncid, varid , DUMMY1   , &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing weasd' )

            error = nf90_inq_varid(ncid, "snwdph", varid)
            call netcdf_err(error, 'getting varid snwdph' )
            error = nf90_put_var(ncid, varid , DUMMY2,  &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing snwdph' )

            error = nf90_inq_varid(ncid, "snowxy", varid)
            call netcdf_err(error, 'getting varid snowxy' )
            error = nf90_put_var(ncid, varid , DUMMY3  , &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing snowxy' )

            error = nf90_inq_varid(ncid, "sneqvoxy", varid)
            call netcdf_err(error, 'getting varid sneqvoxy' )
            error = nf90_put_var(ncid, varid , DUMMY4, &
                start = (/1,1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing sneqvoxy' )

            error = nf90_inq_varid(ncid, "zsnsoxy", varid)
            call netcdf_err(error, 'getting varid zsnsoxy' )
            error = nf90_put_var(ncid, varid, DUMMY5, &
                start = (/1, 1, 1/), count = (/lensfc, 7, 1/))
            call netcdf_err(error, 'writing zsnsoxy' )

            error = nf90_inq_varid(ncid, "tsnoxy", varid)
            call netcdf_err(error, 'getting varid tsnoxy' )
            error = nf90_put_var(ncid, varid, DUMMY6, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing tsnoxy' )

            error = nf90_inq_varid(ncid, "snicexy", varid)
            call netcdf_err(error, 'getting varid snicexy' )
            error = nf90_put_var(ncid, varid, DUMMY7, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing snicexy' )

            error = nf90_inq_varid(ncid, "snliqxy", varid)
            call netcdf_err(error, 'getting varid snliqxy' )
            error = nf90_put_var(ncid, varid, DUMMY8, &
                start = (/1, 1, 1/), count = (/lensfc, 3, 1/))
            call netcdf_err(error, 'writing snliqxy' )

            error = nf90_inq_varid(ncid, "stc", varid)
            call netcdf_err(error, 'getting varid stc' )
            error = nf90_put_var(ncid, varid, DUMMY9, &
                start = (/1, 1/), count = (/lensfc, 1/))
            call netcdf_err(error, 'writing stc')

            error = nf90_close(ncid)
            call netcdf_err(error, 'closing restart file' )

        endif

        ! obs within rad influence 
        if (MYRANK == 0) then
            INQUIRE(FILE=trim(da_output_file), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(da_output_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif
            error = nf90_open(trim(da_output_file), NF90_WRITE, ncid)
            call netcdf_err(error, 'opening DA file '//trim(da_output_file) )
        endif

        Do ie=1, max_num_nearStn+1
            if (MYRANK > 0) then   
                call MPI_SEND(index_obs_assmilated(ie,:), LENSFC_land, mpiInt_size, 0,   &
                                    1300*(MYRANK+1), MPI_COMM_WORLD, IERR)      
            Endif
            if (MYRANK == 0) then
                DUMMY13(1:LENSFC_land) = index_obs_assmilated(ie,:)
                Do ixy =  1, NPROCS - 1   ! sender proc index within tile group
                    if(ixy < N_sA_Ext) then   !p_tRank== 0) then 
                        dest_Aoffset = ixy * (N_sA + 1) + 1    ! 1
                        dest_Aoffset_end = dest_Aoffset + N_sA    !+ 1 
                        arLen = N_sA + 1  
                    else
                        dest_Aoffset = ixy * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
                        dest_Aoffset_end = dest_Aoffset + N_sA - 1 
                        arLen = N_sA
                    endif
                    call MPI_RECV(DUMMY13(dest_Aoffset:dest_Aoffset_end), arLen, mpiInt_size, ixy,      &
                                1300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                End do            
                error = nf90_put_var(ncid, id_no, DUMMY13, &
                        start = (/ie,1/), count = (/1, lensfc/))
                call netcdf_err(error, 'WRITING index_obs_assmilated')
            endif
        End do
        if (MYRANK == 0) then
            error = nf90_close(ncid)
            call netcdf_err(error, 'closing DA file'//trim(da_output_file) )
        endif
             
 End subroutine Write_DA_Outputs_restart_vector

 Subroutine Write_DA_Outputs(output_file, idim, jdim, lensfc, myrank,   &
            snoforc, snoanl, snwdforc, snwdanal, landmask,  &
            num_stn, Lat_atObs, Lon_atObs, OBS_stn, FCS_at_stn, OmB_innov_at_stn,  & 
                        incr_atGrid, SNCOV_IMS, IMS_Foot_Print, &
                    num_Eval, Lat_atEvalPts, Lon_atEvalPts, Obs_atEvalPts, & 
            SNOFCS_atEvalPts, incr_atEvalPts, SNOANL_atEvalPts)  !, anl_fSCA) !updated snocov
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
        real, intent(in)            :: snoforc(lensfc), snoanl(lensfc), &
                                       snwdforc(lensfc), snwdanal(lensfc)
        integer, intent(in)         :: landmask(lensfc)
        Real, intent(in)            :: OBS_stn(num_stn), FCS_at_stn(num_stn), &
                                       OmB_innov_at_stn(num_stn)
        Real, intent(in)            :: Lat_atObs(num_stn), Lon_atObs(num_stn)
        Real, intent(in)            :: incr_atGrid(lensfc), &
                                       SNCOV_IMS(lensfc), IMS_Foot_Print(lensfc)  !, anl_fSCA(lensfc)
        integer, intent(in)         :: num_Eval
        real, intent(in)    :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), Obs_atEvalPts(num_Eval)
        real, intent(in)    :: SNOFCS_atEvalPts(num_Eval), SNOANL_atEvalPts(num_Eval) 
        real, intent(in)    :: incr_atEvalPts(num_Eval)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        integer                     :: error, i, ncid
        integer                     :: dim_x, dim_y, dim_time, dim_eval, dim_stn
        integer                     :: id_x, id_y, id_time
        integer       :: id_swe_forc, id_swe, id_snwdf, id_snwd, id_incr, id_imscov, id_imsftp   !, id_anlscov
        integer       :: id_swe_cur, id_snwd_cur
      
        integer       :: id_latstn, id_lonstn, id_obsstn, id_forcstn, id_innovstn, id_landmask
        integer       :: id_lateval, id_loneval, id_obseval, id_forceval, id_anleval, id_increval, id_cur_anleval  !, , id_anlscov
        
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
        ! eval obs points        
        if (num_Eval>0) then
                error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        else
                error = nf90_def_dim(ncid, 'eval_points', 1, dim_eval)
                call netcdf_err(error, 'DEFINING eval_points' )
        endif

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

        error = nf90_def_var(ncid, 'SWE_Forecast', NF90_DOUBLE, dims_3d, id_swe_forc)
        call netcdf_err(error, 'DEFINING SWE_Forecast' )
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
        call netcdf_err(error, 'DEFINING SND Forecast' )
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


        error = nf90_def_var(ncid, 'DA_Increment', NF90_DOUBLE, dims_3d, id_incr)
        call netcdf_err(error, 'DEFINING DA_Increment' )
        error = nf90_put_att(ncid, id_incr, "long_name", "DA_Increment at model grid")
        call netcdf_err(error, 'DEFINING DA_Increment LONG NAME' )
        error = nf90_put_att(ncid, id_incr, "units", "mm")
        call netcdf_err(error, 'DEFINING DA_Increment UNITS' )

        error = nf90_def_var(ncid, 'imsfSCA', NF90_DOUBLE, dims_3d, id_imscov)
        call netcdf_err(error, 'DEFINING imsfSCA' )
        error = nf90_put_att(ncid, id_imscov, "long_name", "IMS fractional Snow Covered Area")
        call netcdf_err(error, 'DEFINING imsfSCA LONG NAME' )
        error = nf90_put_att(ncid, id_imscov, "units", "-")
        call netcdf_err(error, 'DEFINING imsfSCA UNITS' )

        error = nf90_def_var(ncid, 'ims_FootPrint', NF90_DOUBLE, dims_3d, id_imsftp)
        call netcdf_err(error, 'DEFINING ims_FootPrint' )
        error = nf90_put_att(ncid, id_imsftp, "long_name", "Assimilated IMS Data Foot Print")
        call netcdf_err(error, 'DEFINING ims_FootPrint LONG NAME' )
        error = nf90_put_att(ncid, id_imsftp, "units", "-")
        call netcdf_err(error, 'DEFINING ims_FootPrint UNITS' )

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

        error = nf90_def_var(ncid, 'snwdph@omb', NF90_DOUBLE, dim_stn, id_innovstn)
        call netcdf_err(error, 'DEFINING snwdph@omb' )
        error = nf90_put_att(ncid, id_innovstn, "long_name", "Innovation(O-B) at Observation Points")
        call netcdf_err(error, 'DEFINING snwdph@omb LONG NAME' )
        error = nf90_put_att(ncid, id_innovstn, "units", "mm")
        call netcdf_err(error, 'DEFINING snwdph@omb UNITS' )

        ! eval points
        ! if (num_Eval>0) then 
                error = nf90_def_var(ncid, 'LatEvalPoints', NF90_DOUBLE, dim_eval, id_lateval)
                call netcdf_err(error, 'DEFINING LatEvalPoints' )
                error = nf90_put_att(ncid, id_lateval, "long_name", "Latitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LatEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_lateval, "units", "deg")
                call netcdf_err(error, 'DEFINING LatEvalPoints UNITS' )

                error = nf90_def_var(ncid, 'LonEvalPoints', NF90_DOUBLE, dim_eval, id_loneval)
                call netcdf_err(error, 'DEFINING LonEvalPoints' )
                error = nf90_put_att(ncid, id_loneval, "long_name", "Longitude at Evaluation Points")
                call netcdf_err(error, 'DEFINING LonEvalPoints LONG NAME' )
                error = nf90_put_att(ncid, id_loneval, "units", "deg")
                call netcdf_err(error, 'DEFINING LonEvalPoints UNITS' )
                
                error = nf90_def_var(ncid, 'Obs_atEvalPts', NF90_DOUBLE, dim_eval, id_obseval)
                call netcdf_err(error, 'DEFINING Obs_atEvalPts' )
                error = nf90_put_att(ncid, id_obseval, "long_name", "Observed at Evaluation Points")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_obseval, "units", "mm")
                call netcdf_err(error, 'DEFINING Obs_atEvalPts UNITS' )
                
                error = nf90_def_var(ncid, 'SNOFCS_atEvalPts', NF90_DOUBLE, dim_eval, id_forceval)
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts' )
                error = nf90_put_att(ncid, id_forceval, "long_name", "Forecast at Evaluation Points")
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_forceval, "units", "mm")
                call netcdf_err(error, 'DEFINING SNOFCS_atEvalPts UNITS' )

                error = nf90_def_var(ncid, 'SNOANL_atEvalPts', NF90_DOUBLE, dim_eval, id_anleval)
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts' )
                error = nf90_put_att(ncid, id_anleval, "long_name", "Analysis at Evaluation Points")
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_anleval, "units", "mm")
                call netcdf_err(error, 'DEFINING SNOANL_atEvalPts UNITS' )
                
                error = nf90_def_var(ncid, 'incr_atEvalPts', NF90_DOUBLE, dim_eval, id_increval)
                call netcdf_err(error, 'DEFINING incr_atEvalPts' )
                error = nf90_put_att(ncid, id_increval, "long_name", "Increment at Evaluation Points")
                call netcdf_err(error, 'DEFINING incr_atEvalPts LONG NAME' )
                error = nf90_put_att(ncid, id_increval, "units", "mm")
                call netcdf_err(error, 'DEFINING incr_atEvalPts UNITS' )
        ! endif

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

        dum2d = reshape(incr_atGrid, (/idim,jdim/))
        error = nf90_put_var(ncid, id_incr, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING incr_atGrid RECORD' )

        dum2d = reshape(SNCOV_IMS, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imscov, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )

        dum2d = reshape(IMS_Foot_Print, (/idim,jdim/))
        error = nf90_put_var( ncid, id_imsftp, dum2d, dims_strt, dims_end)
        call netcdf_err(error, 'WRITING imsfSCA RECORD' )
        
        ! obs points (obs, hofx, omb) 
        error = nf90_put_var( ncid, id_latstn, Lat_atObs)
        call netcdf_err(error, 'WRITING Lat_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_lonstn, Lon_atObs)
        call netcdf_err(error, 'WRITING Lon_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_obsstn, OBS_stn)
        call netcdf_err(error, 'WRITING Obs_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_forcstn, FCS_at_stn)
        call netcdf_err(error, 'WRITING SNOFCS_atObsPts RECORD' )

        error = nf90_put_var( ncid, id_innovstn, OmB_innov_at_stn)
        call netcdf_err(error, 'WRITING innov_atObsPts RECORD' )

        ! eval points
        if (num_Eval>0) then 
                error = nf90_put_var( ncid, id_lateval, Lat_atEvalPts)
                call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_loneval, Lon_atEvalPts)
                call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_obseval, Obs_atEvalPts)
                call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_forceval, SNOFCS_atEvalPts)
                call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_anleval, SNOANL_atEvalPts)
                call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_increval, incr_atEvalPts)
                call netcdf_err(error, 'WRITING incr_atEvalPts RECORD' )
        else
                error = nf90_put_var( ncid, id_lateval, 0.0)
                call netcdf_err(error, 'WRITING Lat_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_loneval, 0.0)
                call netcdf_err(error, 'WRITING Lon_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_obseval, 0.0)
                call netcdf_err(error, 'WRITING Obs_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_forceval, 0.0)
                call netcdf_err(error, 'WRITING SNOFCS_atEvalPts RECORD' )

                error = nf90_put_var( ncid, id_anleval, 0.0)
                call netcdf_err(error, 'WRITING SNOANL_atEvalPts RECORD' )
               
                error = nf90_put_var( ncid, id_increval, 0.0)
                call netcdf_err(error, 'WRITING incr_atEvalPts RECORD' )
        endif

        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 End subroutine Write_DA_Outputs

 Subroutine Append_DA_Outputs(output_file, myrank,   &
                                 num_Eval, var_name, SNOANL_atEvalPts)  
	!------------------------------------------------------------------
	! Append var_name to existing nc file: 
	!------------------------------------------------------------------
	implicit none

	CHARACTER(LEN=*), Intent(In)      :: output_file, var_name
	integer, intent(in)         :: num_Eval
	real, intent(in)            :: SNOANL_atEvalPts(num_Eval)

	integer                     :: fsize=65536, inital=0
	integer                     :: header_buffer_val = 16384
	integer                     :: error, i, ncid, dim_eval, id_anleval
	integer                     :: myrank

	include "mpif.h"

	print*, "Process ", myrank, "writing output data to: ",trim(output_file)

	! error = NF90_CREATE(output_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), ncid, initialsize=inital, chunksize=fsize)
	error = NF90_OPEN(trim(output_file), NF90_WRITE, NCID)
	call netcdf_err(error, 'Opening FILE='//trim(output_file) )
	
	! obs points
	! error = nf90_def_dim(ncid, 'eval_points', num_Eval, dim_eval)
    error = NF90_INQ_DIMID(ncid, 'eval_points', dim_eval)
    call netcdf_err(error, 'ERROR READING Dimension' )
	
	error=NF90_REDEF(ncid) 
	call netcdf_err(error, 'ERROR Redefining header' )

	error = nf90_def_var(ncid, trim(var_name), NF90_DOUBLE, dim_eval, id_anleval)
	call netcdf_err(error, 'DEFINING '//trim(var_name) )
	error = nf90_put_att(ncid, id_anleval, "long_name", trim(var_name))
	call netcdf_err(error, 'DEFINING '//trim(var_name)//' LONG NAME' )
	error = nf90_put_att(ncid, id_anleval, "units", "mm")
	call netcdf_err(error, 'DEFINING '//trim(var_name)//' UNITS' )

	error = nf90_enddef(ncid, header_buffer_val,4,0,4)
	call netcdf_err(error, 'DEFINING HEADER' )

	
	! eval points
	error = nf90_put_var( ncid, id_anleval, SNOANL_atEvalPts)
	call netcdf_err(error, 'WRITING SNOANL_Cur_atEvalPts RECORD' )

	error = nf90_close(ncid)
    
 End subroutine Append_DA_Outputs

 SUBROUTINE Observation_Read_atEval(ghcnd_inp_file, &   !dim_name,                      &
                    NDIM, Lat_GHCND, Lon_GHCND, MYRANK)
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file     !, dim_name
        INTEGER                :: ERROR, NCID
        INTEGER                :: MYRANK
        INTEGER                :: ID_DIM, ID_VAR
        INTEGER, Intent(Out)   :: NDIM
    
        REAL, ALLOCATABLE, Intent(Out)     :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, 'eval_points', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
    
        ERROR=NF90_INQ_VARID(NCID, 'LatEvalPoints', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'LonEvalPoints', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
 End SUBROUTINE Observation_Read_atEval

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
     
 SUBROUTINE read_ensemble_forecast(TILE_NUM, ens_size, LENSFC, &
                IDIM_In, JDIM_In, ens_inp_path, &   !y_str, m_str, d_str, h_str, & ! hprev_str, &
                        SWEFCS_Ens, SNDFCS_Ens)    ! assim_SWE, 
    
        IMPLICIT NONE

        include "mpif.h"

        ! INTEGER, INTENT(IN)               :: MYRANK, NUM_TILES,    
        CHARACTER(LEN=5), INTENT(In)      :: TILE_NUM    
        INTEGER, Intent(In)               :: IDIM_In, JDIM_In, ens_size, LENSFC 
        ! Integer, Intent(In)               :: index_ens_atGrid(LENSFC)
        !REAL, INTENT(In)                  :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4), RLA(LENSFC), RLO(LENSFC)
        CHARACTER(LEN=*), Intent(In)      :: ens_inp_path  !, y_str, m_str, d_str, h_str  !, hprev_str        
        ! LOGICAL, Intent(In)               :: assim_SWE, save_Ens
        REAL, INTENT(OUT)         :: SWEFCS_Ens(ens_size, LENSFC), SNDFCS_Ens(ens_size, LENSFC) !, SNOCOV_Ens(ens_size, LENSFC)  !VEGFCS(LENSFC), 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)
        ! CHARACTER(LEN=*)          :: da_out_file

        CHARACTER(LEN=3)          :: ens_str
        CHARACTER(LEN=250)        :: ens_inp_path_full

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM, ens_indx       
        ! INTEGER                   :: ix, jy, dim_x, dim_y, dim_ens, id_x, id_y, id_ens, id_forc
        INTEGER                   :: ID_VAR

        REAL(KIND=8)              :: DUMMY(IDIM_In, JDIM_In)  ! DUMMY2D(IDIM_In, JDIM_In)  !, DUMMY(IDIM_In/2, JDIM_In/2), 
        ! REAL(KIND=8)              :: SNOFCS_Inp_Ens_In(LENSFC/4)
        ! real(kind=4), allocatable   :: x_data(:), y_data(:), ens_data(:)
        ! integer                     :: fsize=65536, inital=0
        ! integer                     :: header_buffer_val = 16384
        ! integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        LOGICAL                   :: file_exists
        
        Do ens_indx = 1, ens_size
                WRITE(ens_str, '(I3.3)') (ens_indx)
                ens_inp_path_full = TRIM(ens_inp_path)//"mem"//ens_str//"/INPUT/"// &
                      !TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // TRIM(h_str)//"0000."&
                                 "sfc_data."//TILE_NUM//".nc"

                INQUIRE(FILE=trim(ens_inp_path_full), EXIST=file_exists)
                if (.not. file_exists) then 
                        print *, 'read_ensemble_forecast error, file does not exist', &   
                                trim(ens_inp_path_full) , ' exiting'
                        call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
                endif
                ERROR=NF90_OPEN(TRIM(ens_inp_path_full),NF90_NOWRITE,NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ens_inp_path_full) )

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
                ! if(assim_SWE) then
                        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
                        ! Do ix = 1, IDIM_In
                        !       Do jy = 1, JDIM_In
                        !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                        !       End do
                        ! End do
                        ! SNOFCS_Inp_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))
                SWEFCS_Ens(ens_indx, :) = RESHAPE(DUMMY, (/LENSFC/))
                ! else
                        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
                SNDFCS_Ens(ens_indx, :) = RESHAPE(DUMMY, (/LENSFC/))
                ! endif
                ! SNOFCS_Inp_Ens_In = RESHAPE(DUMMY, (/(LENSFC/4)/))
                ! Do ix = 1, LENSFC
                !         SNOFCS_Inp_Ens(ens_indx, ix) = SNOFCS_Inp_Ens_In(index_ens_atGrid(ix))
                ! End do          
                ! ERROR=NF90_INQ_VARID(NCID, "sncovr", ID_VAR)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr ID' )
                ! ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr' )
                ! Do ix = 1, IDIM_In
                !       Do jy = 1, JDIM_In
                !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                !       End do
                ! End do
                ! SNOCOV_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))

                ERROR = NF90_CLOSE(NCID)
        End do

        !deallocate(DUMMY, DUMMY2D)
        
 END SUBROUTINE read_ensemble_forecast

 SUBROUTINE read_ensemble_forecast_temporal(TILE_NUM, ens_size, LENSFC, &
                IDIM_In, JDIM_In, ens_inp_path, &   !, y_str, m_str, d_str, h_str, & ! hprev_str, &
                IY, IM, ID, IH, dT_Asssim, &   
                SWEFCS_Ens, SNDFCS_Ens)    ! assim_SWE, 
    
        IMPLICIT NONE

        include "mpif.h"

        ! INTEGER, INTENT(IN)               :: MYRANK, NUM_TILES,    
        CHARACTER(LEN=5), INTENT(In)      :: TILE_NUM    
        INTEGER, Intent(In)               :: IDIM_In, JDIM_In, ens_size, LENSFC 
        INTEGER, Intent(In)               :: IY, IM, ID, IH
        REAL, INTENT(In)                  :: dT_Asssim
        ! Integer, Intent(In)               :: index_ens_atGrid(LENSFC)
        !REAL, INTENT(In)                  :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4), RLA(LENSFC), RLO(LENSFC)
        CHARACTER(LEN=*), Intent(In)      :: ens_inp_path  !, y_str, m_str, d_str, h_str  !, hprev_str        
        ! LOGICAL, Intent(In)               :: assim_SWE, save_Ens
        REAL, INTENT(OUT)         :: SWEFCS_Ens(ens_size, LENSFC), SNDFCS_Ens(ens_size, LENSFC) !, SNOCOV_Ens(ens_size, LENSFC)  !VEGFCS(LENSFC), 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)
        ! CHARACTER(LEN=*)          :: da_out_file

        CHARACTER(LEN=3)          :: ens_str
        CHARACTER(LEN=250)        :: ens_inp_path_full

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM, ens_indx       
        ! INTEGER                   :: ix, jy, dim_x, dim_y, dim_ens, id_x, id_y, id_ens, id_forc
        INTEGER                   :: ID_VAR

        REAL(KIND=8)              :: DUMMY(IDIM_In, JDIM_In)  ! DUMMY2D(IDIM_In, JDIM_In)  !, DUMMY(IDIM_In/2, JDIM_In/2), 
        ! REAL(KIND=8)              :: SNOFCS_Inp_Ens_In(LENSFC/4)
        ! real(kind=4), allocatable   :: x_data(:), y_data(:), ens_data(:)
        ! integer                     :: fsize=65536, inital=0
        ! integer                     :: header_buffer_val = 16384
        ! integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        LOGICAL                       :: file_exists
        CHARACTER(LEN=250)            :: y_str, m_str, d_str, h_str  
        INTEGER                       :: IY_loc, IM_loc, ID_loc, IH_loc  
        REAL                          :: IH_real, JulianDate
        
        IY_loc = IY; IM_loc=IM; ID_loc=ID; IH_loc=IH
        IH_real = Real(IH)

        Call JULDAT (IY_loc, IM_loc, ID_loc, IH_real, JulianDate)
        JulianDate = JulianDate - (0.5 * (dT_Asssim/24.0) * ens_size)
        Call CALDAT (JulianDate, IY_loc, IM_loc, ID_loc, IH_real)
        IH_loc = INT(IH_real)

        Do ens_indx = 1, ens_size

            write(y_str, "(I4)") IY_loc
            write(m_str, "(I0.2)") IM_loc
            write(d_str, "(I0.2)") ID_loc
            write(h_str, "(I0.2)") IH_loc

            ens_inp_path_full =  TRIM(ens_inp_path)//TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"."//&
                                    TRIM(h_str)//"0000.sfc_data."//TILE_NUM//".nc" 

            ! WRITE(ens_str, '(I3.3)') (ens_indx)              
            ! TRIM(ens_inp_path)//"mem"//ens_str//"/INPUT/"// &
            !       TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // TRIM(h_str)//"0000."&
            !                  "sfc_data."//TILE_NUM//".nc"

            INQUIRE(FILE=trim(ens_inp_path_full), EXIST=file_exists)
            if (.not. file_exists) then 
                    print *, 'read_ensemble_forecast error, file does not exist', &   
                            trim(ens_inp_path_full) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif
            ERROR=NF90_OPEN(TRIM(ens_inp_path_full),NF90_NOWRITE,NCID)
            CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ens_inp_path_full) )

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
            ! if(assim_SWE) then
                    ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
                    CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
                    ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                    CALL NETCDF_ERR(ERROR, 'READING sheleg' )
                    ! Do ix = 1, IDIM_In
                    !       Do jy = 1, JDIM_In
                    !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                    !       End do
                    ! End do
                    ! SNOFCS_Inp_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))
            SWEFCS_Ens(ens_indx, :) = RESHAPE(DUMMY, (/LENSFC/))
            ! else
                    ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
                    CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
                    ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                    CALL NETCDF_ERR(ERROR, 'READING snwdph' )
            SNDFCS_Ens(ens_indx, :) = RESHAPE(DUMMY, (/LENSFC/))
            ! endif
            ! SNOFCS_Inp_Ens_In = RESHAPE(DUMMY, (/(LENSFC/4)/))
            ! Do ix = 1, LENSFC
            !         SNOFCS_Inp_Ens(ens_indx, ix) = SNOFCS_Inp_Ens_In(index_ens_atGrid(ix))
            ! End do          
            ! ERROR=NF90_INQ_VARID(NCID, "sncovr", ID_VAR)
            ! CALL NETCDF_ERR(ERROR, 'READING sncovr ID' )
            ! ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
            ! CALL NETCDF_ERR(ERROR, 'READING sncovr' )
            ! Do ix = 1, IDIM_In
            !       Do jy = 1, JDIM_In
            !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
            !       End do
            ! End do
            ! SNOCOV_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))

            ERROR = NF90_CLOSE(NCID)
            Call UPDATEtime(IY_loc, IM_loc, ID_loc, IH_real, dT_Asssim)
            IH_loc = INT(IH_real)    ! (for the obs. available currently) this should be 18
        End do

        !deallocate(DUMMY, DUMMY2D)
        
 END SUBROUTINE read_ensemble_forecast_temporal
   
 SUBROUTINE read_ensemble_forecast_regrid(MYRANK, NUM_TILES, TILE_NUM, ens_size, LENSFC, &
        IDIM_In, JDIM_In, index_ens_atGrid,    &
                ens_inp_path, y_str, m_str, d_str, h_str, hprev_str, &
                        assim_SWE, save_Ens, SNOFCS_Inp_Ens, da_out_file)    ! assim_SWE, 
    
        IMPLICIT NONE

        include "mpif.h"

        INTEGER, INTENT(IN)               :: MYRANK, NUM_TILES, ens_size, LENSFC        
        CHARACTER(LEN=5), INTENT(In)      :: TILE_NUM
        INTEGER, Intent(In)               :: IDIM_In, JDIM_In
        Integer, Intent(In)               :: index_ens_atGrid(LENSFC)
        !REAL, INTENT(In)                  :: RLA_cg(LENSFC/4), RLO_cg(LENSFC/4), RLA(LENSFC), RLO(LENSFC)
        CHARACTER(LEN=*), Intent(In)      :: ens_inp_path, y_str, m_str, d_str, h_str, hprev_str        
        LOGICAL, Intent(In)                       :: assim_SWE, save_Ens
        REAL, INTENT(OUT)         :: SNOFCS_Inp_Ens(ens_size, LENSFC) !, SNOCOV_Ens(ens_size, LENSFC)  !VEGFCS(LENSFC), 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)
        CHARACTER(LEN=*)          :: da_out_file

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: ens_str
        CHARACTER(LEN=250)        :: ens_inp_path_full

        INTEGER                   :: ERROR, NCID
        INTEGER                   :: IDIM, JDIM, ID_DIM, ens_indx       
        INTEGER                   :: i, ix, jy, dim_x, dim_y, dim_ens, id_x, id_y, id_ens, id_forc
        INTEGER                   :: ID_VAR

        REAL(KIND=8)              :: DUMMY(IDIM_In/2, JDIM_In/2), DUMMY2D(IDIM_In, JDIM_In)
        REAL(KIND=8)              :: SNOFCS_Inp_Ens_In(LENSFC/4)
        real(kind=4), allocatable   :: x_data(:), y_data(:), ens_data(:)
        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dims_3d(3), dims_strt(3), dims_end(3)
        
        Do ens_indx = 1, ens_size
                WRITE(ens_str, '(I3.3)') (ens_indx)
                ens_inp_path_full = TRIM(ens_inp_path)//"mem"//ens_str//"/RESTART/"// &
                                                TRIM(y_str)//TRIM(m_str)//TRIM(d_str)//"." // &
                                                TRIM(h_str)// "0000.sfc_data."//TILE_NUM//".nc"

                ERROR=NF90_OPEN(TRIM(ens_inp_path_full),NF90_NOWRITE,NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ens_inp_path_full) )

                ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
                ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

                IF ((IDIM*JDIM) /= (LENSFC/4)) THEN
                PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
                CALL MPI_ABORT(MPI_COMM_WORLD, 88)
                ENDIF
                if(assim_SWE) then
                        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING sheleg' )
                else
                        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
                        ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                        CALL NETCDF_ERR(ERROR, 'READING snwdph' )
                endif
                ! Do ix = 1, IDIM_In
                !       Do jy = 1, JDIM_In
                !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                !       End do
                ! End do
                ! SNOFCS_Inp_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))
                SNOFCS_Inp_Ens_In = RESHAPE(DUMMY, (/(LENSFC/4)/))
                Do ix = 1, LENSFC
                        SNOFCS_Inp_Ens(ens_indx, ix) = SNOFCS_Inp_Ens_In(index_ens_atGrid(ix))
                End do          
                ! ERROR=NF90_INQ_VARID(NCID, "sncovr", ID_VAR)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr ID' )
                ! ERROR=NF90_GET_VAR(NCID, ID_VAR, DUMMY)
                ! CALL NETCDF_ERR(ERROR, 'READING sncovr' )
                ! Do ix = 1, IDIM_In
                !       Do jy = 1, JDIM_In
                !               DUMMY2D(ix, jy) = DUMMY((ix+1)/2, (jy+1)/2)
                !       End do
                ! End do
                ! SNOCOV_Ens(ens_indx, :) = RESHAPE(DUMMY2D, (/LENSFC/))

                ERROR = NF90_CLOSE(NCID)
        End do
        
        If (save_Ens .and. (Myrank < NUM_TILES)) then
        
                print*,"Process ", myrank, "writing ens snow forc data to: ",trim(da_out_file)
                
                !--- create the file
                error = NF90_CREATE(da_out_file, IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL), &
                                                                        ncid, initialsize=inital, chunksize=fsize)
                call netcdf_err(error, 'CREATING FILE='//trim(da_out_file) )

                !--- define dimensions
                error = nf90_def_dim(ncid, 'xaxis_1', IDIM_In, dim_x)
                call netcdf_err(error, 'DEFINING XAXIS DIMENSION' )
                error = nf90_def_dim(ncid, 'yaxis_1', JDIM_In, dim_y)
                call netcdf_err(error, 'DEFINING YAXIS DIMENSION' )
                error = nf90_def_dim(ncid, 'Ens_mem', ens_size, dim_ens)
                call netcdf_err(error, 'DEFINING ENS DIMENSION' )

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

                error = nf90_def_var(ncid, 'Ens_mem', NF90_FLOAT, dim_ens, id_ens)
                call netcdf_err(error, 'DEFINING ENS FIELD' )
                error = nf90_put_att(ncid, id_ens, "long_name", "Ensemble members")
                call netcdf_err(error, 'DEFINING ENS LONG NAME' )
                error = nf90_put_att(ncid, id_ens, "units", "-")
                call netcdf_err(error, 'DEFINING ENS UNITS' )
                error = nf90_put_att(ncid, id_ens, "cartesian_axis", "ENS")
                call netcdf_err(error, 'WRITING ENS FIELD' )

                dims_3d(1) = dim_x
                dims_3d(2) = dim_y
                dims_3d(3) = dim_ens

                error = nf90_def_var(ncid, 'SND_Forecast', NF90_DOUBLE, dims_3d, id_forc)
                call netcdf_err(error, 'DEFINING SND_Forecast' )
                error = nf90_put_att(ncid, id_forc, "long_name", "Forecast Snow Depth")
                call netcdf_err(error, 'DEFINING SND Forecast LONG NAME' )
                error = nf90_put_att(ncid, id_forc, "units", "mm")
                call netcdf_err(error, 'DEFINING SND Forecast UNITS' )

                error = nf90_enddef(ncid, header_buffer_val,4,0,4)
                call netcdf_err(error, 'DEFINING HEADER' )

                allocate(x_data(IDIM_In))
                do i = 1, IDIM_In
                x_data(i) = float(i)
                enddo
                allocate(y_data(JDIM_In))
                do i = 1, JDIM_In
                y_data(i) = float(i)
                enddo
                allocate(ens_data(ens_size))
                do i = 1, ens_size
                        ens_data(i) = float(i)
                enddo

                error = nf90_put_var( ncid, id_x, x_data)
                call netcdf_err(error, 'WRITING XAXIS RECORD' )
                error = nf90_put_var( ncid, id_y, y_data)
                call netcdf_err(error, 'WRITING YAXIS RECORD' )
                error = nf90_put_var( ncid, id_ens, ens_data)
                call netcdf_err(error, 'WRITING Ensmem RECORD' )

                Do ens_indx = 1, ens_size
                        dims_strt(1:2) = 1
                        dims_strt(3) = ens_indx
                        dims_end(1) = IDIM_In
                        dims_end(2) = JDIM_In
                        dims_end(3) = 1
                        
                        DUMMY2D = reshape(SNOFCS_Inp_Ens(ens_indx, :), (/IDIM_In,JDIM_In/))
                        error = nf90_put_var( ncid, id_forc, DUMMY2D, dims_strt, dims_end)
                        call netcdf_err(error, 'WRITING SWE Forecast RECORD' )
                End do
                error = nf90_close(ncid)
                deallocate(x_data, y_data, ens_data)
        Endif

        !deallocate(DUMMY, DUMMY2D)
        
 END SUBROUTINE read_ensemble_forecast_regrid


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

 SUBROUTINE READ_NoahMP_Restart(restart_path, LENSFC, noahmp) !VEGFCS, !SRFLAG
    

        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: restart_path
        INTEGER, INTENT(IN)               :: LENSFC
        type(noahmp_type)                 :: noahmp
    
        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID, i
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR, iv

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists

        ! Real :: swe                (LENSFC)
        ! Real :: snow_depth         (LENSFC)
        ! Real :: active_snow_layers (LENSFC)
        ! Real :: swe_previous       (LENSFC)
        ! Real :: snow_soil_interface(LENSFC,7)
        ! Real :: temperature_snow   (LENSFC,3)
        ! Real :: snow_ice_layer     (LENSFC,3)
        ! Real :: snow_liq_layer     (LENSFC,3)
        ! Real :: temperature_soil   (LENSFC)
        associate( &
            ! obs_snow_depth => obs%snow_depth            ,&
            swe => noahmp%swe                ,&
            snow_depth => noahmp%snow_depth         ,&
            active_snow_layers => noahmp%active_snow_layers ,&
            swe_previous => noahmp%swe_previous       ,&
            snow_soil_interface => noahmp%snow_soil_interface,&
            temperature_snow => noahmp%temperature_snow   ,&
            snow_ice_layer => noahmp%snow_ice_layer     ,&
            snow_liq_layer => noahmp%snow_liq_layer     ,&
            temperature_soil => noahmp%temperature_soil )

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! WRITE(RANKCH, '(I3.3)') (MYRANK+1)
        ! FNBGSI = "./fnbgsi." // RANKCH
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(forc_inp_path)
        
        INQUIRE(FILE=trim(restart_path), EXIST=file_exists)

        if (.not. file_exists) then 
                print *, 'READ_Forecast_Data_atPath error,file does not exist', &   
                        trim(restart_path) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif

        ERROR=NF90_OPEN(TRIM(restart_path), NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(restart_path) )
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
        print*, 'idim, jdim, lensfc = ', IDIM, JDIM, LENSFC
        CALL MPI_ABORT(MPI_COMM_WORLD, 88)
        ENDIF

        ALLOCATE(DUMMY(IDIM,JDIM))

        ! ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        ! CALL NETCDF_ERR(ERROR, 'READING sheleg' )
        ! SWEFCS = RESHAPE(DUMMY, (/LENSFC/))

        ! ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        ! CALL NETCDF_ERR(ERROR, 'READING snwdph' )
        ! SNDFCS = RESHAPE(DUMMY, (/LENSFC/))

        ! ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        ! CALL NETCDF_ERR(ERROR, 'READING vtype' )
        ! VETFCS = RESHAPE(DUMMY, (/LENSFC/))    

        ! ERROR=NF90_INQ_VARID(NCID, "slmsk", ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'READING slmsk ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        ! CALL NETCDF_ERR(ERROR, 'READING slmsk' )
        ! SLMASK = RESHAPE(DUMMY, (/LENSFC/))    
        
        ! ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        ! do i = 1, LENSFC 
        !    ! if land, but not land ice, set mask to 1.
        !    if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
        !         ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
        !         LANDMASK(i) = 1 
        !    else 
        !         LANDMASK(i) = 0
        !    endif 
        ! enddo

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)  !weasd
        CALL NETCDF_ERR(ERROR, 'READING weasd ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING weasd' )
        swe = RESHAPE(DUMMY, (/LENSFC/))   
        
        ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph_m ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph_m' )
        snow_depth = RESHAPE(DUMMY, (/LENSFC/))   
        
        ERROR=NF90_INQ_VARID(NCID, "snowxy", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snowxy ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snowxy' )
        active_snow_layers = RESHAPE(DUMMY, (/LENSFC/))   

        ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)           !sneqvoxy
        CALL NETCDF_ERR(ERROR, 'READING sneqvoxy ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING sneqvoxy' )
        swe_previous = RESHAPE(DUMMY, (/LENSFC/))   
            
        ERROR=NF90_INQ_VARID(NCID, "zsnsoxy", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING zsnsoxy ID' )        
        Do iv=1, 7 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING zsnsoxy' )     
            snow_soil_interface(:,iv) = RESHAPE(DUMMY, (/LENSFC/))
        Enddo

        Do iv=1, 3 
            ERROR=NF90_INQ_VARID(NCID, "tsnoxy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING tsnoxy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING tsnoxy' )     
            temperature_snow(:,iv) = RESHAPE(DUMMY, (/LENSFC/))

            ERROR=NF90_INQ_VARID(NCID, "snicexy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING snicexy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING snicexy' )     
            snow_ice_layer(:,iv) = RESHAPE(DUMMY, (/LENSFC/))

            ERROR=NF90_INQ_VARID(NCID, "snliqxy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING snliqxy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING snliqxy' )     
            snow_liq_layer(:,iv) = RESHAPE(DUMMY, (/LENSFC/))
        Enddo
!TBC This is tsc in ufs_land code
        ERROR=NF90_INQ_VARID(NCID, "tg3", ID_VAR)  
        CALL NETCDF_ERR(ERROR, 'READING tg3 ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING tg3' )
        temperature_soil = RESHAPE(DUMMY, (/LENSFC/))

        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)

        end associate
    
 END SUBROUTINE READ_NoahMP_Restart

 SUBROUTINE READ_Restart(restart_path, veg_type_landice, LENSFC, SWEFCS, SNDFCS, VETFCS, &
                         LANDMASK, swe, snow_depth, active_snow_layers, swe_previous, &
                         snow_soil_interface, temperature_snow, snow_ice_layer, &
                         snow_liq_layer, temperature_soil) !VEGFCS, !SRFLAG
    
        IMPLICIT NONE

        include "mpif.h"
        
        CHARACTER(LEN=*), Intent(In)      :: restart_path
        INTEGER, INTENT(IN)               :: LENSFC, veg_type_landice
        REAL, INTENT(OUT)                 :: SWEFCS(LENSFC), SNDFCS(LENSFC), VETFCS(LENSFC)  !VEGFCS(LENSFC), 
        INTEGER, INTENT(OUT)              :: LANDMASK(LENSFC) 
        !REAL, INTENT(OUT)        :: FMM(LENSFC), FHH(LENSFC), SRFLAG(LENSFC)
        Real :: swe                (LENSFC)
        Real :: snow_depth         (LENSFC)
        Real :: active_snow_layers (LENSFC)
        Real :: swe_previous       (LENSFC)
        Real :: snow_soil_interface(LENSFC,7)
        Real :: temperature_snow   (LENSFC,3)
        Real :: snow_ice_layer     (LENSFC,3)
        Real :: snow_liq_layer     (LENSFC,3)
        Real :: temperature_soil   (LENSFC)

        CHARACTER(LEN=50)         :: FNBGSI
        CHARACTER(LEN=3)          :: RANKCH

        INTEGER                   :: ERROR, NCID, i
        INTEGER                   :: IDIM, JDIM, ID_DIM
        INTEGER                   :: ID_VAR, iv

        REAL(KIND=8), ALLOCATABLE :: DUMMY(:,:)
        REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! WRITE(RANKCH, '(I3.3)') (MYRANK+1)
        ! FNBGSI = "./fnbgsi." // RANKCH
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(FNBGSI)

        !CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, ERROR)
        ! if (print_deb) PRINT*, "READ INPUT SFC DATA FROM: "//TRIM(forc_inp_path)
        
        INQUIRE(FILE=trim(restart_path), EXIST=file_exists)

        if (.not. file_exists) then 
                print *, 'READ_Forecast_Data_atPath error,file does not exist', &   
                        trim(restart_path) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif

        ERROR=NF90_OPEN(TRIM(restart_path), NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(restart_path) )
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
        print*, 'idim, jdim, lensfc = ', IDIM, JDIM, LENSFC
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
        
        ! slmsk is 0 - ocean, 1 - land, 2 -seaice 
        ! convert to integer  0 not land or glacier, 1 - non-glacier covered land
        do i = 1, LENSFC 
           ! if land, but not land ice, set mask to 1.
           if ( (NINT(SLMASK(i)) == 1 ) .and.   & 
                ( NINT(VETFCS(i)) /=  veg_type_landice  )) then 
                LANDMASK(i) = 1 
           else 
                LANDMASK(i) = 0
           endif 
        enddo

        ERROR=NF90_INQ_VARID(NCID, "weasd", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING weasd ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING weasd' )
        swe = RESHAPE(DUMMY, (/LENSFC/))   
        
        ERROR=NF90_INQ_VARID(NCID, "snwdph_m", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snwdph_m ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snwdph_m' )
        snow_depth = RESHAPE(DUMMY, (/LENSFC/))   
        
        ERROR=NF90_INQ_VARID(NCID, "snowxy", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING snowxy ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING snowxy' )
        active_snow_layers = RESHAPE(DUMMY, (/LENSFC/))        

        ERROR=NF90_INQ_VARID(NCID, "zsnsoxy", ID_VAR)
        CALL NETCDF_ERR(ERROR, 'READING zsnsoxy ID' )        
        Do iv=1, 7 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING zsnsoxy' )     
            snow_soil_interface(:,iv) = RESHAPE(DUMMY, (/LENSFC/))
        Enddo

        Do iv=1, 3 
            ERROR=NF90_INQ_VARID(NCID, "tsnoxy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING tsnoxy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING tsnoxy' )     
            temperature_snow(:,iv) = RESHAPE(DUMMY, (/LENSFC/))

            ERROR=NF90_INQ_VARID(NCID, "snicexy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING snicexy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING snicexy' )     
            snow_ice_layer(:,iv) = RESHAPE(DUMMY, (/LENSFC/))

            ERROR=NF90_INQ_VARID(NCID, "snliqxy", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING snliqxy ID' ) 
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy, &
            start = (/1, 1, iv, 1/), count =(/idim, jdim, 1, 1/))
            CALL NETCDF_ERR(ERROR, 'READING snliqxy' )     
            snow_liq_layer(:,iv) = RESHAPE(DUMMY, (/LENSFC/))
        Enddo
!TBC This is tsc in ufs_land code
        ERROR=NF90_INQ_VARID(NCID, "tg3", ID_VAR)  
        CALL NETCDF_ERR(ERROR, 'READING tg3 ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy)
        CALL NETCDF_ERR(ERROR, 'READING tg3' )
        temperature_soil = RESHAPE(DUMMY, (/LENSFC/))

        DEALLOCATE(DUMMY)

        ERROR = NF90_CLOSE(NCID)
    
 END SUBROUTINE READ_Restart

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
        print*, 'idim, jdim, lensfc = ', IDIM, JDIM, LENSFC
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
        print*, 'idim, jdim, lensfc = ', IDIM, JDIM, LENSFC
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
 
