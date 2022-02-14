 PROGRAM driver_snowda
 
! !2021/2022 adapted from ufs_utils/global_cycle/cycle.f90 
! the surface driver 
! modified to include OI and EnKF snow anlysis for NOAH and NOAH-MP 
!----------------------------------------------------------------------
!  ! NAMSNO NAMELIST VARIABLE DEFINITIONS:
!
!  IDIM,JDIM      i/j dimension of a cubed-sphere tile.
!  LSOIL          Number of soil layers.
!  IY,IM,ID,IH    Year, month, day, and hour of initial state.
!  FH             Forecast hour
!  IVEGSRC        Use igbp veg type when '1'.  Use sib when '2'.

! SNOW_DA_TYPE controls the OI of snow obs (in snow analysis)
! 0  - no DA of snow
! 1  - OI applied to NOAH snow depth on FV3 tiles
! 2  - OI applied to NOAH-MP multi-layer snow depth on land-only FV3 grids vector
! 3  - EnKF applied to NOAH-MP multi-layer snow depth on land-only vector
! 4  - EnSRF applied to NOAH-MP on land only 
!lsm_type
!  1 - NOAH
!  2 - NOAH-MP

!PERCENT_OBS_WITHHELD : % (e.g., 5%) of obs excluded from DA, for evaluation
!horz_len_scale, ver_len_scale: horizonal (in Km) and vertical (m) 
!                              correlation length Scales in Brassnet (1998) 
!                              correlation/localization function
!obs_tolerance: observations are rejected if obs-forecast >> obs_tolerance* standdev is 
!obs_srch_rad : look for snow depth observations within this distance from grid cell (Km)
!bkgst_srch_rad: for nearest neighbor interpolation of background state to observation 
! max_num_nearStn: maximum number of observations near a grid cell to be assimilated
! max_num_nearIMS: maximum number of IMS observations to be assimilated at grid cell (1)
! ims_max_ele: elvation cutoff for IMS data 
! num_subgrd_ims_cels: number of IMS grid cells within UFS/FV3 grid to calcluate snow cover fraction
! stdev_obsv_depth, stdev_obsv_sncov, stdev_back: snow depth, snow cover, background error
!                                                 standard deviations
! ENS_SIZE: number of ensemble members
! rcov_localize: boolean, localize the R error covariance matrix
! ens_inflate : boolean, inflate ensemble error covariance
! assim_SnowPack_obs: boolean, whether to assimilate snowpack obs or not
! assim_SnowCov_obs: whether to assimilate snow cover obs
! ims_correlated: are IMS obs assumed correlated (if max_num_nearIMS > 1)    
! STN_OBS_PREFIX: directory and name of snow depth obs (without the date) 
! STN_DIM_NAME,STN_VAR_NAME,STN_ELE_NAME: names of dimension, variable, and elevation within 
!                                         the snow depth obs netCDF file
! IMS_SNOWCOVER_PATH: IMS observation location 
! IMS_INDEXES_PATH:   location of IMS indices mapping IMS to FV3 grid
! restart_prefix, openloop_prefix: restart file and open loop file locations
! static_prefix: location of static file contaning land vector information 
                 ! e.g., lat/lon/land cover type....
! output_prefix: output location
! print_debg_info: debug information to be printed
! PRINTRANK: MPI rank to print debug info
! snowUpdateOpt: for NOAH-MP analysis snow depth partitioning
    ! DO_NOTHING = 0
    ! Add all increment to TOP_LAYER = 1
    ! Add all incrment to BOTTOM_LAYER = 2
    ! parition increment between ALL_LAYERS = 3
! begloc, endloc: start and end indexes of the land vector 
! exclude_obs_at_grid: exclude obsevation at this point (for evaluation of impact of
!                      assimilating observation away from the point)

! Note:SFCSUB coded for different snow obs. types, but OI calling program is not

! Clara S Draper and Tseganeh Z Gichamo
!-----------------------------------------------------------------------------------

    USE M_Snow_Analysis
    ! USE M_UFSLAND_SNOW_UPDATE
    ! USE MPI

    IMPLICIT NONE
    !
    ! include 'mpif.h'

    ! CHARACTER(LEN=3) :: DONST
    INTEGER :: IDIM, JDIM, LSOIL, IY, IM, ID, IH, NUM_TILES
    INTEGER :: IVEGSRC, ISOT, LENSFC, IERR     !, LUGB, IALB, ISOT, , ZSEA1_MM, ZSEA2_MM
    INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS    !, MAX_TASKS   !, SNOW_OI_TYPE
    INTEGER :: SNOW_DA_TYPE
    REAL    :: FH   !, DELTSFC, ZSEA1, ZSEA2
    ! LOGICAL :: USE_UFO, DO_NSST, ADJT_NST_ONLY
 
    REAL, ALLOCATABLE   :: SNDANL(:)   !, incr_at_Grid(:,:)  !, SWEANL(:), SNOFCS(:)    !, SNOANL_out(:), incr_at_Grid_out(:)
    REAL                :: PERCENT_OBS_WITHHELD
    REAL                :: horz_len_scale, ver_len_scale, obs_tolerance 
    REAL                :: obs_srch_rad, bkgst_srch_rad, ims_max_ele  !, dT_Asssim
    REAL                :: stdev_obsv_depth, stdev_obsv_sncov, stdev_back
    Integer             :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels
    Integer             :: ENS_SIZE   !num_assim_steps, ims_assm_hour, 
    LOGICAL             :: assim_SnowPack_obs, assim_SnowCov_obs, ims_correlated, &
                           rcov_localize, ens_inflate
    CHARACTER(LEN=500)  :: STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, & 
                           IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH,  &
                           !SFC_FORECAST_PREFIX, &  ! CURRENT_ANALYSIS_PREFIX, ENKFGDAS_TOP_DIR, &
                           restart_prefix, openloop_prefix, &
                           static_prefix, output_prefix    !, point_state_prefix                       
    ! CHARACTER(len=4)    :: stn_var 
    LOGICAL             :: print_debg_info   !STANDALONE_SNOWDA, , vector_inputs, fv3_index
    Integer             :: begloc, endloc, lsm_type
    LOGICAL             :: exclude_obs_at_grid
    !Integer	     :: s_assm_hour
! noah mp updte   
    integer             :: PRINTRANK, snowUpdateOpt
     ! for mpi par
    INTEGER      :: Np_ext, Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end
    Integer      :: LENSFC_proc

    ! NAMELIST/NAMCYC/ IDIM,JDIM,LSOIL,LUGB,IY,IM,ID,IH,FH,    &
    !                 DELTSFC,IALB,USE_UFO,DONST,             &
    !                 ADJT_NST_ONLY,ISOT,IVEGSRC,ZSEA1_MM,    &
    !                 ZSEA2_MM, MAX_TASKS

    NAMELIST/NAMSNO/ IDIM, JDIM, LSOIL, NUM_TILES,  &
        IY, IM, ID, IH, FH, IVEGSRC, ISOT,    &
    ! LUGB, DELTSFC,IALB,USE_UFO,DONST, ADJT_NST_ONLY,ISOT,ZSEA1_MM, ZSEA2_MM, MAX_TASKS
        SNOW_DA_TYPE, PERCENT_OBS_WITHHELD,    &  !, SNOW_OI_TYPE
        horz_len_scale, ver_len_scale, obs_tolerance, & 
        obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &
        ims_max_ele, num_subgrd_ims_cels, &
        stdev_obsv_depth, stdev_obsv_sncov, stdev_back, &   !num_assim_steps, dT_Asssim, ims_assm_hour, 
        ENS_SIZE, rcov_localize, ens_inflate, &
        assim_SnowPack_obs, assim_SnowCov_obs, ims_correlated,  &    !stn_var, &
        STN_OBS_PREFIX, STN_DIM_NAME,STN_VAR_NAME,STN_ELE_NAME, &
        IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, &      !SFC_FORECAST_PREFIX, &  !CURRENT_ANALYSIS_PREFIX, ENKFGDAS_TOP_DIR, &
        restart_prefix, openloop_prefix, static_prefix, output_prefix, &
        print_debg_info, & !STANDALONE_SNOWDA, , fv3_index, vector_inputs, point_state_prefix, &        
        PRINTRANK, snowUpdateOpt, begloc, endloc, lsm_type, &
        exclude_obs_at_grid
    !
    DATA IDIM,JDIM,LSOIL,NUM_TILES/96,96,4,6/ 
    DATA IY,IM,ID,IH,FH/1997,8,2,0,0./
    ! DATA LUGB/51/, DELTSFC/0.0/, IALB/1/, MAX_TASKS/99999/
    DATA ISOT/1/    !, ZSEA1_MM/0/, ZSEA2_MM/0/
    DATA IVEGSRC/2/
    DATA SNOW_DA_TYPE/0/
    ! DATA SNOW_OI_TYPE/0/
    DATA PERCENT_OBS_WITHHELD/0.0/
    DATA horz_len_scale/55.0/
    DATA ver_len_scale/800./
    DATA obs_tolerance/5./ 
    DATA obs_srch_rad/250.0/
    DATA bkgst_srch_rad/27.0/ 
    DATA max_num_nearStn/50/ 
    DATA max_num_nearIMS/5/
    DATA ims_max_ele/1500./
    DATA num_subgrd_ims_cels/30/
    DATA stdev_obsv_depth/40.0/
    DATA stdev_obsv_sncov/80.0/
    DATA stdev_back/30.0/
    ! DATA num_assim_steps/1/  ! For multiple time steps of assimilation
    ! DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
    ! DATA ims_assm_hour/18/
    DATA ENS_SIZE/20/
    DATA rcov_localize/.false./
    DATA ens_inflate/.false./
    DATA assim_SnowPack_obs/.false./
    DATA assim_SnowCov_obs/.false./
    DATA ims_correlated/.false./
    ! DATA stn_var/'SND'/
    DATA STN_OBS_PREFIX/'        '/ 
    DATA STN_DIM_NAME/"Site_Id"/
    DATA STN_VAR_NAME/"SNWD"/
    DATA STN_ELE_NAME/"elevation"/
    DATA IMS_SNOWCOVER_PATH/'        '/
    DATA IMS_INDEXES_PATH/'        '/
    ! DATA SFC_FORECAST_PREFIX/'        '/   ! leave this empty to use the default sfc_ files location
    ! DATA CURRENT_ANALYSIS_PREFIX/'        '/
    ! DATA ENKFGDAS_TOP_DIR/'        '/
    ! Data point_state_prefix/"./"/
    DATA restart_prefix/"./"/
    DATA openloop_prefix/"./"/
    DATA static_prefix/"./"/
    DATA output_prefix/"./"/
    ! DATA STANDALONE_SNOWDA/.false./
    DATA print_debg_info/.false./
    DATA PRINTRANK/0/
    DATA snowUpdateOpt/0/
    ! DATA vector_inputs/.false./
    ! DATA fv3_index/.false./
    DATA begloc/1/
    DATA endloc/1/ 
    DATA lsm_type/2/
    DATA exclude_obs_at_grid/.false./

    CALL MPI_INIT(IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

    if (myrank==0) call w3tagb('GLOBAL_CYCLE',2018,0179,0055,'NP20')

    NUM_THREADS = NUM_PARTHDS()

    PRINT*
    PRINT*,"STARTING CYCLE PROGRAM ON RANK ", MYRANK
    PRINT*,"RUNNING WITH ", NPROCS, "TASKS"
    PRINT*,"AND WITH ", NUM_THREADS, " THREADS."

    ! USE_UFO = .FALSE.
    ! DONST   = "NO"
    ! ADJT_NST_ONLY = .FALSE.

    ! IF (MYRANK==0) PRINT*, "READING NAMCYC NAMELIST."
    ! CALL BAOPENR(36, "fort.36", IERR)
    ! READ(36, NML=NAMCYC)
    ! IF (MYRANK==0) WRITE(6,NAMCYC)

    IF (MYRANK==0) PRINT*,"READING NAMSNO NAMELIST."

    CALL BAOPENR(360, "fort.360", IERR)             !"snowDA.nml"   !
    READ(360, NML=NAMSNO)
    IF (MYRANK==0) WRITE(6, NAMSNO)
    
    If(SNOW_DA_TYPE .eq. 1) then 
        ! Call ReadVectorLength(static_prefix, IDIM, begloc, endloc, LENSFC)
        LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE 
    else
        LENSFC = endloc - begloc + 1
    endif
    print*, "Vector length ", LENSFC

    ! Np_ext = MOD(NPROCS, NUM_TILES)  ! extra/inactive procs
        ! if (MYRANK >  NPROCS - Np_ext - 1) goto 999
    Np_til = NPROCS !/ NUM_TILES  ! num proc. per tile 
    p_tN = 0 !MOD(MYRANK, NUM_TILES)  ! tile for proc.
    p_tRank = MYRANK ! / NUM_TILES  ! proc. rank within tile
    N_sA = LENSFC/NPROCS    !/ Np_til  ! sub array length per proc
    N_sA_Ext = MOD(LENSFC, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
    if(MYRANK < N_sA_Ext) then   !p_tRank== 0) then 
        !  mp_start = 1
        mp_start = MYRANK * (N_sA + 1) + 1 + (begloc - 1)    
        mp_end = mp_start + N_sA    !+ (begloc - 1)           !+ 1 
        LENSFC_proc = N_sA + 1
    else
        !mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        mp_start = MYRANK * N_sA + N_sA_Ext + 1 + (begloc - 1)     !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
        mp_end = mp_start + N_sA - 1    ! + (begloc - 1)    
        LENSFC_proc = N_sA 
    endif
    
    !mp_end = (p_tRank + 1) * N_sA + N_sA_Ext  
    PRINT*,"Proc ", myrank, " snowDA: sub array length ", N_sA, &
           " extra sub array: ", N_sA_Ext
    PRINT*,"Proc ", myrank, " start ", mp_start, " end ", mp_end

    if (myrank == (NPROCS - 1)) then
        if (mp_end /= endloc) then
            print*, "location end ", endloc, " =/ subarry end for last proc ", mp_end
            stop
        endif
        print*, "location end ", endloc, " subarray end for last proc ", mp_end
    endif

    ! do snow DA, if requested
    ! IF (SNOW_DA_TYPE .gt. 0) then

    !PRINT*,"snowDA: calling DA on RANK" , MYRANK
    ALLOCATE(SNDANL(LENSFC))        
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    ! s_assm_hour =18
    ! if (IH == s_assm_hour) then
    IF (SNOW_DA_TYPE .eq. 1) then                                
        Call OI_Snow_Analysis_NOAH(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &  
                                ! num_assim_steps, dT_Asssim,  & 
                                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
                                horz_len_scale, ver_len_scale, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                        ims_max_ele, num_subgrd_ims_cels, &
                        STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, & 
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, restart_prefix, output_prefix, &  ! CURRENT_ANALYSIS_PREFIX, &
                                SNDANL, exclude_obs_at_grid)    !incr_at_Grid, SNOFCS,, SWEANL, LANDMASK)
    Else if(SNOW_DA_TYPE .eq. 2) then
        Call OI_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &
                IY, IM, ID, IH, &  
                ! num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, &
                horz_len_scale, ver_len_scale, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, & 
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, &
                restart_prefix, openloop_prefix, static_prefix, output_prefix, &
                snowUpdateOpt, PRINTRANK, print_debg_info, &  !, vector_inputs, , fv3_index
                SNDANL,  &   !SNOFCS, SWEANL, & incr_at_Grid, 
                Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)          
    Else if(SNOW_DA_TYPE .eq. 3) then
        Call EnKF_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &
                IY, IM, ID, IH, &
                ! num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, & 
                horz_len_scale, ver_len_scale, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                ens_size, rcov_localize, ens_inflate,  &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
                restart_prefix, openloop_prefix, static_prefix, &
                output_prefix, &
                snowUpdateOpt, PRINTRANK, print_debg_info, &  !fv3_index, vector_inputs, &
                SNDANL, &   !SNDFCS_out, SWEANL_out, & incr_at_Grid_out, 
                Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)

    Else if(SNOW_DA_TYPE .eq. 4) then    
        Call EnSRF_Snow_Analysis_NOAHMP(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, &
                IY, IM, ID, IH, &
                ! num_assim_steps, dT_Asssim,  & 
                LENSFC, IVEGSRC, PERCENT_OBS_WITHHELD, & 
                horz_len_scale, ver_len_scale, obs_tolerance, &
                obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &                                
                ims_max_ele, num_subgrd_ims_cels, &
                ens_size, rcov_localize, ens_inflate,  &
                assim_SnowPack_obs, assim_SnowCov_obs, &
                STN_OBS_PREFIX, STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, & 
                restart_prefix, openloop_prefix, static_prefix, &
                output_prefix, &
                snowUpdateOpt, PRINTRANK, print_debg_info, &  !fv3_index, vector_inputs, &
                SNDANL, &   !SNDFCS_out, SWEANL_out, & incr_at_Grid_out, 
                Np_til, p_tN, p_tRank, N_sA, N_sA_Ext, mp_start, mp_end, LENSFC_proc, &
                begloc, endloc, exclude_obs_at_grid)
    Endif 

    ! ENDIF
    ! IF (MAX_TASKS < 99999 .AND. MYRANK > (MAX_TASKS - 1)) THEN
    !     PRINT*,"USER SPECIFIED MAX NUMBER OF TASKS: ", MAX_TASKS
    !     PRINT*,"WILL NOT RUN CYCLE PROGRAM ON RANK: ", MYRANK
    !     GOTO 333
    ! ENDIF

! 333 CONTINUE
    
    ! IF (ALLOCATED(SNOFCS)) DEALLOCATE(SNOFCS)   
    IF (ALLOCATED(SNDANL)) DEALLOCATE(SNDANL)
    ! IF (ALLOCATED(SWEANL)) DEALLOCATE(SWEANL)
    ! IF (ALLOCATED(incr_at_Grid)) DEALLOCATE(incr_at_Grid)
    
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    PRINT*
    PRINT*,' SNOWDA PROGRAM COMPLETED NORMALLY ON RANK: ', MYRANK

    if (myrank==0) call w3tage('GLOBAL_CYCLE')

    CALL MPI_FINALIZE(IERR)

    STOP

 END PROGRAM driver_snowda