! subroutine set_num_parthds(num_threads)
! use omp_lib
! integer num_threads
! !$OMP PARALLEL
! call omp_set_num_threads(num_threads) 
! print*,"number of threads = ", omp_get_num_threads()
! !$OMP END PARALLEL
! ! return
! end subroutine set_num_parthds

 integer function num_parthds()
 use omp_lib
!$OMP PARALLEL
 num_parthds=omp_get_num_threads()
!$OMP END PARALLEL
 return
 end function num_parthds
