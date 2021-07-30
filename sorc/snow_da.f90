
integer(2) function compar(a1, a2) 

    IMPLICIT NONE
    Real :: a1, a2

    if (a1 .LT. a2) then 
        compar = -1
    else if (a1 .GT. a2) then
        compar = 1
    else
        compar = 0
    endif

end function compar

MODULE M_DA

    USE NETCDF
    Use, Intrinsic :: IEEE_ARITHMETIC
    use, Intrinsic :: iso_c_binding
    use IFPORT
    
    Logical, parameter :: print_debug = .False.

    Logical, parameter :: print_deb = .False.  
    
    Integer(2), External :: compar
    
    CONTAINS
    
    subroutine compute_covariances(RLA_jndx, RLO_jndx, Orog_jndx,  SNOforc_jndx,  &
        Lat_Obs, Lon_Obs, Ele_Obs, num_Obs,                             &
        Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,         &
        L_horz, h_ver,                                      &  
        assim_IMS,                                          &
        B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect)
!**
        
        IMPLICIT NONE
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In)        :: RLA_jndx, RLO_jndx, Orog_jndx, SNOforc_jndx
        Real, Intent(In)        :: Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims 
        Real, Intent(In)        :: Lon_Obs(num_Obs), Lat_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In)        :: L_horz, h_ver 
        Integer, Intent(In) :: num_Obs
        LOGICAL, Intent(In) :: assim_IMS

        Real(dp), Intent(Out)    :: B_cov_mat(num_obs,num_obs), b_cov_vect(num_obs)
        Real(dp), Intent(Out)    :: O_cov_mat(num_obs,num_obs), W_wght_vect(num_obs)
        
        Real(dp)    :: W_wght_vect_intr(1, num_obs)     
        
        Integer :: indx, jndx, zndx    
        Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: Innov_cov_mat(num_obs,num_obs), Innov_cov_mat_inv(num_obs,num_obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)    
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real, Parameter         :: snforc_tol = 0.001
        Real, Parameter         :: Bcov_scale_factor = 0.1
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2

        ! deg to rad    
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   

        !1. Observation covariance 
        ! O = stdev_o*stdev_o * I , I = Identitity matrix
        O_cov_mat = 0.
        Do indx = 1, num_Obs
            O_cov_mat(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) O_cov_mat(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims

        if (print_debug) then
            print*, "Obs cov"
            print*, O_cov_mat
        endif

        !2. Background states covariance (at observation space)
        ! B = stddev_back*stdev_back * Corr(j, k)
        ! Corr(j, k) = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        ! rjk = horizontal distance between j and k
        ! zjk = vertical distance between j and k
        ! L = horizontal correlation length (note in Km)
        ! h = vertical correlation length   (in m)
        
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
! 4.16.20 ToDO: This is a symmetric matrix: can revise later to doing only half of the computations
        Do jndx = 1, num_Obs 
            d_latArr = (Lat_Obs_rad(jndx) - Lat_Obs_rad) / 2.
            d_lonArr = (Lon_Obs_rad(jndx) - Lon_Obs_rad) / 2.
            haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(Lat_Obs_rad(jndx)) * sin(d_lonArr)**2
            Where (haversinArr > 1) haversinArr = 1
            ! Do indx = 1, num_Obs 
            !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
            ! end do
            rjk_distArr(jndx,:) = 2 * earth_rad * asin(sqrt(haversinArr))       ! rjk, k = 1, Num obs for a given j
            zjk_distArr(jndx,:) = Ele_Obs(jndx) - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        End do
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        if (print_debug) then
            print*, "Dist for Back corr at obs pts"
            print*, rjk_distArr
            print*, "Vertical dist for Back corr at obs pts"
            print*, zjk_distArr
        endif
        B_cov_mat = (1. + rjk_distArr/L_horz) * exp(-1. * rjk_distArr/L_horz) !L_horz in Km, h_ver in m
        B_cov_mat = B_cov_mat * exp(-1. * (zjk_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Backround corr at obs pts"
            print*, B_cov_mat
        endif   
        !B_cov_mat(num_Obs, num_Obs) = 1.
        B_cov_mat = B_cov_mat * Stdev_back * stdev_back   
        if (print_debug) then
            print*, "Backround cov at obs pts"
            print*, B_cov_mat
        endif   
        
        !3. b background covariance between model grid and obs points
        ! similar to the above (B_cov_mat) except just for one point againt N obs points 
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        ! Do indx = 1, num_Obs 
        !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
        ! end do
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))     ! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        b_cov_vect =  (1. + l_distArr/L_horz) * exp(-1. * l_distArr/L_horz)  !L_horz in Km, h_ver in m
        b_cov_vect = b_cov_vect * exp(-1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "b corr between model grid and obs pts"
            print*, b_cov_vect
        endif
        b_cov_vect = b_cov_vect * stdev_back * stdev_back 
        if (print_debug) then
            print*, "b cov between model grid and obs pts"
            print*, b_cov_vect
        endif

        ! 4. Weight vector
        ! W = (B_cov_mat + Obs_cov_mat)^-1 b
        Innov_cov_mat = B_cov_mat + O_cov_mat       
        Innov_cov_mat_inv = inv(Innov_cov_mat)        
        !W_wght_vect = matmul(Innov_cov_mat_inv, RESHAPE(b_cov_vect,(/num_obs,1/)))
        W_wght_vect_intr = matmul(RESHAPE(b_cov_vect,(/1, num_obs/)), Innov_cov_mat_inv) ! [1,m]x[m,m]=[1,m]
        W_wght_vect = RESHAPE(W_wght_vect_intr,(/num_obs/))
        if (print_debug) then
            print*, "Innov cov"
            print*, Innov_cov_mat
            print*, "inverse of Innov cov"
            print*, Innov_cov_mat_inv
            print*, "weights vector"
            print*, W_wght_vect
        endif
        
        RETURN

    END SUBROUTINE compute_covariances   

    SUBROUTINE Snow_DA_OI(back_at_Obs, obs_Array, num_Obs,  &
        W_wght_vect,            &
        back_at_Grid, anl_at_Grid, obs_Innov)
!**
    
        IMPLICIT NONE

        include 'mpif.h'

        integer, parameter :: dp = kind(1.d0)

        REAL, INTENT(In)    :: back_at_Obs(num_Obs), obs_Array(num_Obs)
        Integer, Intent(In)  :: num_Obs
        Real(dp), Intent(In)    :: W_wght_vect(num_Obs)    !, weighted_Innov
        REAL, INTENT(In)    :: back_at_Grid
        REAL, INTENT(Out)    :: anl_at_Grid
        Real, INTENT(Out)  :: obs_Innov(num_Obs)
        Real               :: anl_incr_interm(1,1)  ! this is the increment, not the innovation

        obs_Innov = obs_Array - back_at_Obs
        anl_incr_interm = matmul(RESHAPE(W_wght_vect,(/1,num_Obs/)), &
                                    RESHAPE(obs_Innov,(/num_Obs,1/))) ! weighted innovation
        
        anl_at_Grid = back_at_Grid + anl_incr_interm(1,1)

        RETURN

    END SUBROUTINE Snow_DA_OI
    
    !http://fortranwiki.org/fortran/show/Matrix+inversion
    ! ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
        
        IMPLICIT NONE
        
        integer, parameter :: dp = kind(1.d0)

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Ainv
    
        real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
    
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
    
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
    
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
    
        if (info /= 0) then
        stop 'Matrix is numerically singular!'
        end if
    
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
        if (info /= 0) then
        stop 'Matrix inversion failed!'
        end if
    end function inv

    subroutine nearest_Observations_Locations(RLA_jndx, RLO_jndx,    &
        Lat_Obs, Lon_Obs, num_Obs, max_distance, max_num_loc,        &
        Stdev_back, Stdev_Obs_depth,                                 &
        SNOFCS_atObs, OBS_atOBs,                                     &
        Loc_backSt_atObs,  num_loc) 
!**
        IMPLICIT NONE

        Real, Intent(In)        :: RLA_jndx, RLO_jndx  ! don't want to alter these
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs)
        Integer, Intent(In) :: num_Obs, max_num_loc
        Real, Intent(In)        :: max_distance   ! radius_of_influence
        Real, Intent(In)        :: Stdev_back, Stdev_Obs_depth
        Real, Intent(In)    :: SNOFCS_atObs(num_Obs)
        Real, Intent(In)    :: OBS_atOBs(num_Obs)
        Integer, Allocatable, Intent(Out)    :: Loc_backSt_atObs(:)
        Integer, Intent(Out) :: num_loc
        
        Integer :: indx, jndx, zndx, num_loc_counter        
        Real    :: distArr(num_Obs), haversinArr(num_Obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    :: Lon_Obs_2(num_Obs)  
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real, Allocatable   :: dist_atObs(:)
        Real, Allocatable   :: dist_atObs_dummy(:)
        Integer, Allocatable   :: Loc_backSt_atObs_dummy(:)
        Real                 :: max_value

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize 

        i_isize = SIZEOF(max_value)

        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2

        ! at each obs point compute its distance from RLA/RLO pairs 
        ! shortest distance over sphere using great circle distance     
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   
        
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
        ! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
        !Do jndx = 1, LENSFC 
        d_latArr = (Lat_Obs_rad - RLA_rad_jndx) / 2.
        d_lonArr = (Lon_Obs_rad - RLO_rad_jndx) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        Where (haversinArr < 0) haversinArr = 0.
        distArr = 2 * earth_rad * asin(sqrt(haversinArr))
        
        ! 4.15.20: can you do the following without loop?       
        num_loc_counter = 0
        Do indx = 1, num_Obs    
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. & 
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx))) ) then 
                num_loc_counter = num_loc_counter + 1
            endif
        End do
        num_loc = num_loc_counter
        Allocate(Loc_backSt_atObs_dummy(num_loc))
        Allocate(dist_atObs_dummy(num_loc))
        jndx = 1
        Do indx = 1, num_Obs    
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. &
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx))) ) then 
                Loc_backSt_atObs_dummy(jndx) = indx
                dist_atObs_dummy(jndx) = distArr(indx)
                jndx = jndx  + 1
            endif
        End do
        
        ! if num of obs > 50, choose the 50 closest obs
        if (num_loc > max_num_loc) then 
            Allocate(dist_atObs(num_loc))
            dist_atObs = dist_atObs_dummy
            i_len = num_loc
            Call QSORT (dist_atObs_dummy, i_len, i_isize, compar)
            max_value = dist_atObs_dummy(max_num_loc)
            Allocate(Loc_backSt_atObs(max_num_loc)) 
            indx = 1     
            Do jndx = 1, num_loc
                if(dist_atObs(jndx) <= max_value) then
                    Loc_backSt_atObs(indx) = Loc_backSt_atObs_dummy(jndx)
                    indx = indx + 1
                endif
            End do
            num_loc = max_num_loc                       
            Deallocate(dist_atObs)
        else
            Allocate(Loc_backSt_atObs(num_loc)) 
            Loc_backSt_atObs = Loc_backSt_atObs_dummy
        endif

        Deallocate(dist_atObs_dummy, Loc_backSt_atObs_dummy)
        
        RETURN

    End subroutine nearest_Observations_Locations

     ! Gets obs snow depth from IMS based on exponential/log 'depletion curve' 
    subroutine CalcSWEFromSnowCover(SNCOV_IMS, VETFCS_in, LENSFC, SWE_IMS_at_Grid, SNUP_Array)
!**

        
        IMPLICIT NONE
        !
        Real, Intent(In)        :: SNCOV_IMS(LENSFC),  VETFCS_in(LENSFC)
        INTEGER, Intent(In)     :: LENSFC
        !SNUP_array is the SWE (mm) at which SCF reaches 100% 
        Real, Intent(Out)       :: SWE_IMS_at_Grid(LENSFC), SNUP_Array(LENSFC)
        
        INTEGER            :: VETFCS(LENSFC)
        REAL               :: snupx(30), SNUP, SALP, RSNOW
        Integer                    :: indx, vtype_int


        !Fill background values to nan (to differentiate those that don't have value)
        SWE_IMS_at_Grid = IEEE_VALUE(SWE_IMS_at_Grid, IEEE_QUIET_NAN)
   
        ! NOTE: this is an empirical inversion of   snfrac rotuine in Noah 
        !  should really have a land model check in here. 

        !This is for the IGBP veg classification scheme.
        ! SWE at which snow cover reaches 100%, in m
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
    
        SALP = -4.0
        VETFCS = INT(VETFCS_in)
        Where(VETFCS==0) VETFCS = 7 
        
        Do indx = 1, LENSFC  
            if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
                SNUP = snupx(VETFCS(indx))
                if (SNUP == 0.) then
                    print*, " 0.0 snup value, check vegclasses", VETFCS(indx)
                    Stop
                endif

                SNUP_Array(indx) = SNUP * 1000. !  mm

                if (SNCOV_IMS(indx) >= 1.0) then
                    RSNOW = 1.
                elseif (SNCOV_IMS(indx) < 0.001) then
                    RSNOW = 0.0 
                else
                    RSNOW = min(LOG(1. - SNCOV_IMS(indx)) / SALP, 1.0) 
                endif  
                ! return SWE in mm (or SWE if snowdens = 1)
                SWE_IMS_at_Grid(indx) = RSNOW * SNUP * 1000. !  mm
            endif
        end do  
        RETURN
    
    END SUBROUTINE CalcSWEFromSnowCover
    
    subroutine resample_to_model_tiles_intrp(data_grid_ims, data_grid_ims_ind, &
                                                nlat_ims, nlon_ims, n_lat, n_lon, num_sub, & 
                                                grid_dat)
!**

                                                
        Use, Intrinsic :: IEEE_ARITHMETIC
    
        Implicit None
    
        Integer, Intent(In)     :: nlat_ims, nlon_ims, n_lat, n_lon, num_sub 
        Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, n_lon, n_lat) 
        Real, Intent(Out)       :: grid_dat(n_lon, n_lat)
    
        Integer   :: jc, jy, ix, num_loc_counter
        Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
        
        grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)
    
        Do jy=1, n_lat
    
            Do ix=1, n_lon
                
                num_loc_counter = data_grid_ims_ind(1, ix, jy)
                if (num_loc_counter < 1) then 
                    cycle
                end if

                grid_dat(ix, jy) = 0.
                Do jc = 2, num_loc_counter+1
                    lonlatcoord_ims = data_grid_ims_ind(jc, ix, jy) - 1 
                    latcoord_ims = lonlatcoord_ims / nlon_ims + 1
                    loncoord_ims = mod(lonlatcoord_ims, nlon_ims) + 1
                    if(latcoord_ims > nlat_ims) then
                        latcoord_ims = nlat_ims
                        print*, "Warning! lat coordinate outside domain boundary"
                    endif
                    if(loncoord_ims > nlon_ims) then
                        loncoord_ims = nlon_ims
                        print*, "Warning! lon coordinate outside domain boundary"
                    endif
                    grid_dat(ix, jy) =  grid_dat(ix, jy) + data_grid_ims(loncoord_ims, latcoord_ims)              
                End do
    
                grid_dat(ix, jy) =  grid_dat(ix, jy) / num_loc_counter ! first location, num obs
    
            End do
    
        End do
    
        return !grid_dat
    
    End subroutine resample_to_model_tiles_intrp
    
     !This reads the whole IMS file and uses a-priori prepared indices to sample those wihin the grid cel
     SUBROUTINE Observation_Read_IMS_Full(inp_file, inp_file_indices, &
                    MYRANK, n_lat, n_lon, num_sub, &
                    SNCOV_IMS)
!**
                    
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)   :: inp_file, inp_file_indices 
        INTEGER, Intent(In)            :: MYRANK, n_lat, n_lon, num_sub
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)  
        Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon)     
    
        INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)   
        Integer                 :: data_grid_ims_ind(num_sub, n_lon, n_lat) 
        Real                    :: grid_dat(n_lon, n_lat)
        
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN, DIM_LEN_lat, DIM_LEN_lon
        LOGICAL                :: file_exists

        INQUIRE(FILE=trim(inp_file), EXIST=file_exists)

        if (.not. file_exists) then
                print *, 'iObservation_Read_IMS_Full error,file does not exist', &
                        trim(inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
    
        ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, 'lat', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lat' )
        
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lat)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lat' )
        
        ERROR=NF90_INQ_DIMID(NCID, 'lon', ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension lon' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=DIM_LEN_lon)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension Lon' )
    
        ALLOCATE(SNCOV_IMS_2D_full(DIM_LEN_lon, DIM_LEN_lat))   

        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
        
        ERROR = NF90_CLOSE(NCID)
! read index file for mapping IMS to model grid 

        INQUIRE(FILE=trim(inp_file_indices), EXIST=file_exists)

        if (.not. file_exists) then
                print *, 'iObservation_Read_IMS_Full error, index file does not exist', &
                        trim(inp_file_indices) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) 
        endif
    
        ERROR=NF90_OPEN(TRIM(inp_file_indices),NF90_NOWRITE, NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )
    
        ERROR=NF90_INQ_VARID(NCID, 'IMS_Indices', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, data_grid_ims_ind, start = (/ 1, 1, 1 /), &
                                count = (/ num_sub, n_lon, n_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices' )
    
        ERROR = NF90_CLOSE(NCID)
    
        Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
        Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1
        
        call resample_to_model_tiles_intrp(SNCOV_IMS_2D_full, data_grid_ims_ind, &
                                           DIM_LEN_lat, DIM_LEN_lon, n_lat, n_lon, num_sub, &  !myrank, &
                                           grid_dat)
    
        SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))
    
        DEALLOCATE(SNCOV_IMS_2D_full)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_IMS_Full
    
     SUBROUTINE Observation_Operator_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
                            RLA, RLO, Lat_Obs, Lon_Obs, stn_obs,              &  !OROG, 
                            LENSFC, num_Obs, max_distance, SNOFCS_back, LANDMASK,  &
                            gross_thold,SNOFCS_atObs, index_back_atObs ) !OROGFCS_atObs, 
! Draper, edited to make generic

        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"
    
        Real, Intent(In)        :: RLA(LENSFC), RLO(LENSFC) !, OROG(LENSFC)
        integer, intent(in)     ::  LANDMASK(LENSFC)
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs), stn_obs(num_obs)  ! don't want to alter these
        INTEGER             :: Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, LENSFC, num_Obs 
        Real                :: max_distance   ! extent from center of grid cell to search for obs
        Real, Intent(In)        :: SNOFCS_back(LENSFC)
        Real, Intent(Out)       :: SNOFCS_atObs(num_Obs)  !, OROGFCS_atObs(num_Obs)
        Integer, Intent(Out)    :: index_back_atObs(num_Obs)   ! model index of grid cell nearst the obs.
        
        Real    :: Lon_Obs_2(num_Obs)          
        Real    :: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        INTEGER :: indx, jndx, jzndx, zndx, min_indx
        Real    :: distArr(LENSFC), haversinArr(LENSFC)
        Real    :: d_latArr(LENSFC), d_lonArr(LENSFC)
        real    :: gross_thold
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR
    
        !Np_til ! num proc. per tile p_tRank ! proc. rank within tile !p_tN  ! tile for proc.
        N_sA = num_Obs / Np_til  ! sub array length per proc
        N_sA_Ext = num_Obs - N_sA * Np_til ! extra grid cells
        if(p_tRank == 0) then 
            mp_start = 1
        else
            mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
        endif
        mp_end = (p_tRank + 1) * N_sA + N_sA_Ext                ! end index of subarray for proc
    
        !Fill background values to nan (to differentiate those htat don't have value)
        SNOFCS_atObs = IEEE_VALUE(SNOFCS_atObs, IEEE_QUIET_NAN) 
        ! OROGFCS_atObs = IEEE_VALUE(OROGFCS_atObs, IEEE_QUIET_NAN)       
        index_back_atObs = -1   ! when corresponding value doesn't exit 
    
        ! RLO from 0 to 360 (no -ve lon)
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! RLO_2 = RLO
        ! Where(RLO_2 > 180) RLO_2 = RLO_2 - 360
    
        ! shortest distance over sphere using great circle distance     
        RLA_rad =  pi_div_180 * RLA
        RLO_rad =  pi_div_180 * RLO
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2           
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
        ! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
        Do indx = mp_start, mp_end   !1, num_Obs 
            d_latArr = (Lat_Obs_rad(indx) - RLA_rad) / 2.
            d_lonArr = (Lon_Obs_rad(indx) - RLO_rad) / 2.
            haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad(indx)) * cos(RLA_rad) * sin(d_lonArr)**2
            WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
            
            distArr = 2 * earth_rad * asin(sqrt(haversinArr))           
            !distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
            min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
    
            if(distArr(min_indx) < max_distance) then ! if too far away, don't use
                if ( (LANDMASK(min_indx) == 1)   .and. &   ! if nearest cell is no land, fcs value remains NaN
                    ( abs( SNOFCS_back(min_indx) - stn_obs(indx)   ) < gross_thold) ) then 
                    SNOFCS_atObs(indx) = SNOFCS_back(min_indx) 
                    ! OROGFCS_atObs(indx) = OROG(min_indx)
                    index_back_atObs(indx) = min_indx
                endif
            endif
        end do
    
    ! ToDO: Better way to handle this?
    ! Real data type size corresponding to mpi
        rsize = SIZEOF(max_distance) 
        Call MPI_TYPE_SIZE(MPI_REAL, mpiReal_size, IERR) 
        If (rsize == 4 ) then 
            mpiReal_size = MPI_REAL4
        elseif (rsize == 8 ) then 
            mpiReal_size = MPI_REAL8
        elseif (rsize == 16 ) then 
            mpiReal_size = MPI_REAL16
            PRINT*," Possible mismatch between Fortran Real ", rsize," and Mpi Real", mpiReal_size
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
    
        if (MYRANK > (MAX_TASKS - 1) ) then
            call MPI_SEND(SNOFCS_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                          MYRANK, MPI_COMM_WORLD, IERR) 
            ! call MPI_SEND(OROGFCS_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
            !               MYRANK*100, MPI_COMM_WORLD, IERR)
            call MPI_SEND(index_back_atObs(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
                          MYRANK*1000, MPI_COMM_WORLD, IERR)
        else !if (MYRANK == p_tN ) then  
            Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                send_proc = MYRANK +  pindex * MAX_TASKS
                call MPI_RECV(SNOFCS_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,  &
                          send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                ! call MPI_RECV(OROGFCS_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,   &
                !           send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
                          send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            enddo
        endif
    !ToDO: better way to do this?
        ! now share the whole array
        if (MYRANK < MAX_TASKS ) then   !if (MYRANK == p_tN ) then      
            Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
                rec_proc = MYRANK +  pindex * MAX_TASKS
                call MPI_SEND(SNOFCS_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK, MPI_COMM_WORLD, IERR) 
                ! call MPI_SEND(OROGFCS_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
                call MPI_SEND(index_back_atObs, num_Obs, mpiInt_size, rec_proc, MYRANK*1000, MPI_COMM_WORLD, IERR)
            enddo
        else 
            call MPI_RECV(SNOFCS_atObs, num_Obs, mpiReal_size, p_tN, p_tN, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            ! call MPI_RECV(OROGFCS_atObs, num_Obs, mpiReal_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(index_back_atObs, num_Obs, mpiInt_size, p_tN, p_tN*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        endif
        
        RETURN
        
     END SUBROUTINE Observation_Operator_Parallel
   
    SUBROUTINE Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name,       &
                    lat_min, lat_max, lon_min, lon_max, &
                    NDIM,                       &
                    SND_GHCND, Ele_GHCND,         &
                    Lat_GHCND,      &
                    Lon_GHCND,          &
                    MYRANK)
! ** 
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        INTEGER, Intent(In)               :: p_tN
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
        REAL, Intent(In)       :: lat_min, lat_max, lon_min, lon_max 
        INTEGER, Intent(Out)   :: NDIM
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)    :: Lat_GHCND(:), Lon_GHCND(:), Ele_GHCND(:)
    
        INTEGER                :: MYRANK, ERROR, NCID, ID_DIM, ID_VAR, NDIM_In
        REAL, ALLOCATABLE      :: Lat_GHCND_In(:), Lon_GHCND_In(:), &
                                  SND_GHCND_In(:), Ele_GHCND_In(:)
        INTEGER, ALLOCATABLE   :: index_Array(:)
        INTEGER                :: jndx, jcounter
        LOGICAL                :: file_exists

        INQUIRE(FILE=trim(ghcnd_inp_file), EXIST=file_exists)

        if (.not. file_exists) then
                print *, 'Observation_Read_GHCND_Tile_excNaN erro,,file does not exist', &
                        trim(ghcnd_inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(SND_GHCND_In(NDIM_In))
        ALLOCATE(Lat_GHCND_In(NDIM_In))
        ALLOCATE(Lon_GHCND_In(NDIM_In))
        ALLOCATE(Ele_GHCND_In(NDIM_In))
    
        ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SND_GHCND_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
        ! need to read corresponding elevation values 
        ERROR=NF90_INQ_VARID(NCID, 'elevation', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        ALLOCATE(index_Array(NDIM_In))
        NDIM = 0
        jcounter = 1
        ! exclude all data if not within lat/lon range, or negative snow depth
        ! separate criteria for tile4 (ptn/root 3) because -180/180 lon falls inside its region
        If(p_tN == 3) then
            Do jndx = 1, NDIM_In
                ! Print*, "jndx = ", jndx
                If((Lat_GHCND_In(jndx) >= lat_min) .and. (Lat_GHCND_In(jndx) <= lat_max) .and. &
                   (((Lon_GHCND_In(jndx) >= lon_min) .and. (Lon_GHCND_In(jndx) <= 180.)) .or. &
                    ((Lon_GHCND_In(jndx) >= -180.) .and. (Lon_GHCND_In(jndx) <= lon_max))) .and. &
                    (SND_GHCND_In(jndx) >= 0 )) then   !(.NOT. IEEE_IS_NAN(SND_GHCND_In(jndx)))) then
                        NDIM = NDIM + 1
                        index_Array(jcounter) = jndx
                        jcounter = jcounter + 1
                Endif
            End do
        Else
            Do jndx = 1, NDIM_In
                If((Lat_GHCND_In(jndx) >= lat_min) .and. (Lat_GHCND_In(jndx) <= lat_max) .and. &
                    (Lon_GHCND_In(jndx) >= lon_min) .and. (Lon_GHCND_In(jndx) <= lon_max) .and. &
                    (SND_GHCND_In(jndx) >= 0 )) then  !(.NOT. IEEE_IS_NAN(SND_GHCND_In(jndx)))) then                        
                        NDIM = NDIM + 1
                        index_Array(jcounter) = jndx
                        jcounter = jcounter + 1
                Endif
            End do
        Endif
        ALLOCATE(SND_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        ALLOCATE(Ele_GHCND(NDIM))
        If(NDIM > 0) then
            Do jndx = 1, NDIM
                SND_GHCND(jndx) = SND_GHCND_In(index_Array(jndx))
                Lat_GHCND(jndx) = Lat_GHCND_In(index_Array(jndx))
                Lon_GHCND(jndx) = Lon_GHCND_In(index_Array(jndx))
                Ele_GHCND(jndx) = Ele_GHCND_In(index_Array(jndx))
            End do
        Endif
        DEALLOCATE(index_Array)
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND_Tile_excNaN

     SUBROUTINE NETCDF_ERR( ERR, STRING )
    
        !--------------------------------------------------------------
        ! IF AT NETCDF CALL RETURNS AN ERROR, PRINT OUT A MESSAGE
        ! AND STOP PROCESSING.
        !--------------------------------------------------------------
        
            IMPLICIT NONE
        
            include 'mpif.h'
        
            INTEGER, INTENT(IN) :: ERR
            CHARACTER(LEN=*), INTENT(IN) :: STRING
            CHARACTER(LEN=80) :: ERRMSG
        
            IF( ERR == NF90_NOERR )RETURN
            ERRMSG = NF90_STRERROR(ERR)
            PRINT*,''
            PRINT*,'FATAL ERROR: ', TRIM(STRING), ': ', TRIM(ERRMSG)
            PRINT*,'STOP.'
            CALL MPI_ABORT(MPI_COMM_WORLD, 999)
        
            RETURN
         END SUBROUTINE NETCDF_ERR
        
         SUBROUTINE debug_print(STRING, num_val )
            
            !--------------------------------------------------------------
            ! prints ERROR  MESSAGE
            !--------------------------------------------------------------
            
             IMPLICIT NONE
            
             CHARACTER(LEN=*), INTENT(IN) :: STRING
             real, Intent(in)                     :: num_val
             CHARACTER(LEN=20)                    :: numval_Str
        
             write(numval_Str, "(F18.3)"),  num_val
                 
             IF(print_deb) PRINT*, TRIM(STRING), " ", numval_Str
            
             RETURN
         END SUBROUTINE debug_print

END MODULE M_DA
     
    
