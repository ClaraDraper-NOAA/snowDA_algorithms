 PROGRAM driver_snowOI
! CSD fill in comments

 USE M_Snow_Analysis

 IMPLICIT NONE
!
 include 'mpif.h'

 INTEGER :: IDIM, JDIM, LSOIL, IY, IM, ID, IH, NUM_TILES
 INTEGER :: LENSFC, IVEGSRC,  IERR
 INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS 

! 
 REAL, ALLOCATABLE   :: SNOANL(:) 
 REAL                :: horz_len_scale, ver_len_scale, obs_tolerance 
 REAL                :: obs_srch_rad, bkgst_srch_rad, ims_max_ele, dT_Asssim
 Integer             :: max_num_nearStn, max_num_nearIMS, num_subgrd_ims_cels, num_assim_steps
 LOGICAL             :: assim_SnowPack_obs, assim_SnowCov_obs
 CHARACTER(LEN=500)  :: GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, &
                        IMS_INDEXES_PATH, SFC_FORECAST_PREFIX

 NAMELIST/NAMSNO/  IDIM,JDIM,LSOIL,IY,IM,ID,IH , IVEGSRC, & 
                  horz_len_scale, ver_len_scale, obs_tolerance, & 
                  obs_srch_rad, bkgst_srch_rad, max_num_nearStn, max_num_nearIMS, &
                  ims_max_ele, num_subgrd_ims_cels, & 
                  assim_SnowPack_obs, assim_SnowCov_obs,  &
                  GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, SFC_FORECAST_PREFIX

! model setup 
 DATA NUM_TILES/6/
! snow DA  defaults
 DATA horz_len_scale/55.0/
 DATA ver_len_scale/800./
 DATA obs_tolerance/5./ 
 DATA obs_srch_rad/250.0/
 DATA bkgst_srch_rad/1000.0/  ! set the default so as to ensure obs aren't discarded at course res. 
 DATA max_num_nearStn/50/ 
 DATA max_num_nearIMS/5/
 DATA ims_max_ele/1500./
 DATA num_subgrd_ims_cels/30/
 DATA assim_SnowPack_obs/.false./
 DATA assim_SnowCov_obs/.false./
 DATA GHCND_SNOWDEPTH_PATH/'        '/ 
 DATA IMS_SNOWCOVER_PATH/'        '/
 DATA IMS_INDEXES_PATH/'        '/
 DATA SFC_FORECAST_PREFIX/'        '/   ! leave this empty to use the default sfc_ files location

 !     If (IDIM == 96) then   
 !         num_subgrd_ims_cels = 627    ! (max) number of IMS subcells within a tile grid cell           
 !         bkgst_srch_rad = 240.     !Km distance from gridcell to search for corresponding background state
 !     elseif (IDIM == 128) then
 !         num_subgrd_ims_cels = 627               
 !        bkgst_srch_rad = 240.      
 !     elseif (IDIM == 768) then
 !         num_subgrd_ims_cels = 30
 !         bkgst_srch_rad = 27.                              !Km 
 
 CALL MPI_INIT(IERR)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

 NUM_THREADS = NUM_PARTHDS()

 PRINT*
 PRINT*,"STARTING CYCLE PROGRAM ON RANK ", MYRANK
 PRINT*,"RUNNING WITH ", NPROCS, "TASKS"
 PRINT*,"AND WITH ", NUM_THREADS, " THREADS."

 PRINT*
 PRINT*,"READ NAMSNO NAMELIST."

 CALL BAOPENR(360, "snowDA.nml", IERR)
 READ(360, NML=NAMSNO)
 IF (MYRANK==0) WRITE(6,NAMSNO)

 LENSFC = IDIM*JDIM ! TOTAL NUMBER OF POINTS FOR THE CUBED-SPHERE TILE

 PRINT*,"snowDA: calling OI on RANK" , MYRANK
 ALLOCATE(SNOANL(LENSFC))
 CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

 Call Snow_Analysis_OI(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &
                      LENSFC, IVEGSRC,   &
                      horz_len_scale, ver_len_scale, obs_tolerance, &
                      obs_srch_rad, bkgst_srch_rad, max_num_nearStn, &
                      ims_max_ele, num_subgrd_ims_cels, &
                      assim_SnowPack_obs, assim_SnowCov_obs, &
                      GHCND_SNOWDEPTH_PATH, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, SFC_FORECAST_PREFIX, &
                      SNOANL)
 PRINT*,"snowDA: returned from OI on RANK", MYRANK

 DEALLOCATE(SNOANL)

 CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

 if (myrank==0) call w3tage('GLOBAL_CYCLE')

 CALL MPI_FINALIZE(IERR)

 STOP

 END PROGRAM driver_snowOI
