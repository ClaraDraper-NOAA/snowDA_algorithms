 program driver_fIMS

 use IMSaggregate_mod

 implicit none

 include 'mpif.h'

 integer :: idim, jdim,  num_tiles
 integer :: lensfc, ierr
 integer :: nprocs, myrank

! 
 integer             :: num_subgrd_IMS_cels, num_assim_steps
 character(len=10)   :: date_str ! yyyymmddhh
 character(len=500)  :: IMS_snowcover_path, IMS_indexes_path 

 namelist/fIMS/  idim,jdim,date_str, num_subgrd_IMS_cels, IMS_snowcover_path, IMS_indexes_path

! model setup 
 data num_tiles/6/
! snow da  defaults
 data num_subgrd_IMS_cels/30/

 call mpi_init(ierr)
 call mpi_comm_size(mpi_comm_world, nprocs, ierr)
 call mpi_comm_rank(mpi_comm_world, myrank, ierr)

 print*,"starting fIMS program on rank ", myrank

 if (nprocs .NE. num_tiles) then 
        print *, 'driver_fIMS error: nprocs must = num_tiles', nprocs, num_tiles 
        call mpi_abort(mpi_comm_world, 999)
 endif

 call baopenr(360, "fIMS.nml", ierr)
 read(360, nml=fIMS)

 lensfc = idim*jdim ! total number of points for the cubed-sphere tile

 call mpi_barrier(mpi_comm_world, ierr)

 call calculate_IMS_fsca(num_tiles, myrank, idim, jdim, &
                      lensfc, num_subgrd_IMS_cels, date_str,&
                      IMS_snowcover_path, IMS_indexes_path)
 print*,"calcfIMS returned on rank", myrank

 call mpi_barrier(mpi_comm_world, ierr)
 call mpi_finalize(ierr)

 end program driver_fIMS
