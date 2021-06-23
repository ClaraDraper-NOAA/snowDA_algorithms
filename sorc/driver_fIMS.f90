 PROGRAM driver_fIMS

 USE M_Snow_Analysis

 IMPLICIT NONE
!
 include 'mpif.h'

 INTEGER :: IDIM, JDIM, LSOIL, IY, IM, ID, IH, NUM_TILES
 INTEGER :: LENSFC, IVEGSRC,  IERR
 INTEGER :: NPROCS, MYRANK, NUM_THREADS, NUM_PARTHDS 

! 
 REAL, ALLOCATABLE   :: SNOANL(:) 
 REAL                :: dT_Asssim
 Integer             :: max_num_nearIMS, num_subgrd_ims_cels, num_assim_steps
 CHARACTER(LEN=500)  :: IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH, SFC_FORECAST_PREFIX

 NAMELIST/NAMSNO/  IDIM,JDIM,LSOIL,IY,IM,ID,IH , IVEGSRC, & 
                  max_num_nearIMS, num_subgrd_ims_cels, & 
                  SFC_FORECAST_PREFIX, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH

! model setup 
 DATA NUM_TILES/6/
! snow DA  defaults
 DATA max_num_nearIMS/5/
 DATA num_subgrd_ims_cels/30/
 DATA num_assim_steps/1/  ! For multiple time steps of assimilation
 DATA dT_Asssim/24.0/     ! hrs. For multiple time steps of assimilation
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

 Call calculate_IMS_fSCA(NUM_TILES, MYRANK, NPROCS, IDIM, JDIM, IY, IM, ID, IH, &
                      LENSFC, IVEGSRC,   &
                      num_assim_steps, dT_Asssim, &
                      max_num_nearIMS, num_subgrd_ims_cels, &
                      SFC_FORECAST_PREFIX, IMS_SNOWCOVER_PATH, IMS_INDEXES_PATH)
 PRINT*,"snowDA: returned from OI on RANK", MYRANK

 DEALLOCATE(SNOANL)

 CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

 if (myrank==0) call w3tage('GLOBAL_CYCLE')

 CALL MPI_FINALIZE(IERR)

 STOP

 END PROGRAM driver_fIMS
