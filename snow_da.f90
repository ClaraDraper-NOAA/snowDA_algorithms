
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


    
    subroutine compute_covariances_multiIMS(RLA_jndx, RLO_jndx, Orog_jndx,  &  !SNOforc_jndx,  &
        Lat_Obs, Lon_Obs, Ele_Obs, num_Obs, num_loc_1, num_loc_2,                 &
        Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,         &
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS, ims_correlated,                                  &
        B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect) !,      &LENSFC,
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In)        :: RLA_jndx, RLO_jndx, Orog_jndx  !, SNOforc_jndx
        Real, Intent(In)        :: Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims 
        Real, Intent(In)        :: Lon_Obs(num_Obs), Lat_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In)        :: L_horz, h_ver  !L_horz in Km, h_ver in m
        Integer, Intent(In)     :: num_Obs, num_loc_1, num_loc_2
        LOGICAL, Intent(In)     :: assim_IMS, ims_correlated

        Real(dp), Intent(Out)    :: B_cov_mat(num_obs,num_obs), b_cov_vect(num_obs)
        Real(dp), Intent(Out)    :: O_cov_mat(num_obs,num_obs), W_wght_vect(num_obs) !6.10.20: W_wght_vect(num_obs, 1)
        
        Real(dp)    :: W_wght_vect_intr(1, num_obs)     
        
        Integer :: indx, jndx, zndx    
        Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: Innov_cov_mat(num_obs,num_obs), Innov_cov_mat_inv(num_obs,num_obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real, Parameter         :: snforc_tol = 0.001
        Real, Parameter         :: Bcov_scale_factor = 0.1
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180      
        
        !Lon between -180 and 180 for some inputs; RLO 0 - 360
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2

        ! deg to rad    
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   

! 1. Spatial correlation between obs. locations
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

!2. Observation covariance 
 ! O = stdev_o*stdev_o * I , I = Identitity matrix
        O_cov_mat = 0.
        if (num_loc_1 > 0) then
            Do indx = 1, num_loc_1
                O_cov_mat(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
            end do
        end if
!CSDCSD - change to have >1 IMS obs.
        if (assim_IMS) then 
            Do indx = 1, num_loc_2
                O_cov_mat(num_loc_1+indx, num_loc_1+indx) = Stdev_Obs_ims * Stdev_Obs_ims
            end do
            if (ims_correlated) then
                O_cov_mat(num_loc_1+1:num_Obs, num_loc_1+1:num_Obs) = &
                B_cov_mat(num_loc_1+1:num_Obs, num_loc_1+1:num_Obs) * Stdev_Obs_ims * Stdev_Obs_ims
            endif
        end if        
        if (print_debug) then
            print*, "Obs cov"
            print*, O_cov_mat
        endif

!3. Background states covariance (at observation space)
! B = stddev_back*stdev_back * Corr(j, k)
        !B_cov_mat(num_Obs, num_Obs) = 1.
        B_cov_mat = B_cov_mat * Stdev_back * stdev_back   
        if (print_debug) then
            print*, "Backround cov at obs pts"
            print*, B_cov_mat
        endif   
        
!4. b background covariance between model grid and obs points
! similar to B_cov_mat except just for one point againt N obs points 
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
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
! 5. Weight vector
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

    END SUBROUTINE compute_covariances_multiIMS
    

    subroutine compute_covariances(RLA_jndx, RLO_jndx, Orog_jndx, &   ! SNOforc_jndx,  &
        Lat_Obs, Lon_Obs, Ele_Obs, num_Obs,                             &
        Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims,         &
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS,                                          &
        B_cov_mat, b_cov_vect, O_cov_mat, W_wght_vect) !,      &LENSFC,
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In)        :: RLA_jndx, RLO_jndx, Orog_jndx   !, SNOforc_jndx
        Real, Intent(In)        :: Stdev_back, Stdev_Obs_depth, Stdev_Obs_ims 
        Real, Intent(In)        :: Lon_Obs(num_Obs), Lat_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In)        :: L_horz, h_ver  !L_horz in Km, h_ver in m
        Integer, Intent(In) :: num_Obs
        LOGICAL, Intent(In) :: assim_IMS

        Real(dp), Intent(Out)    :: B_cov_mat(num_obs,num_obs), b_cov_vect(num_obs)
        Real(dp), Intent(Out)    :: O_cov_mat(num_obs,num_obs), W_wght_vect(num_obs) !6.10.20: W_wght_vect(num_obs, 1)
        
        Real(dp)    :: W_wght_vect_intr(1, num_obs)     
        
        Integer :: indx, jndx, zndx    
        Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: Innov_cov_mat(num_obs,num_obs), Innov_cov_mat_inv(num_obs,num_obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real, Parameter         :: snforc_tol = 0.001
        Real, Parameter         :: Bcov_scale_factor = 0.1
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180      
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! Do zndx = 1, num_Obs 
        !     if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
        ! end do   
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
        back_at_Grid, anl_incr_at_Grid, anl_at_Grid, obs_Innov)
    
        IMPLICIT NONE

        include 'mpif.h'

        integer, parameter :: dp = kind(1.d0)

        REAL, INTENT(In)    :: back_at_Obs(num_Obs), obs_Array(num_Obs)
        Integer, Intent(In)  :: num_Obs
        Real(dp), Intent(In)    :: W_wght_vect(num_Obs)    !, weighted_Innov
        REAL, INTENT(In)    :: back_at_Grid
        REAL, INTENT(Out)    :: anl_incr_at_Grid, anl_at_Grid
        Real, INTENT(Out)  :: obs_Innov(num_Obs)
        Real               :: anl_incr_interm(1,1)  ! this is the increment, not the innovation

        obs_Innov = obs_Array - back_at_Obs
        anl_incr_interm = matmul(RESHAPE(W_wght_vect,(/1,num_Obs/)), &
                                    RESHAPE(obs_Innov,(/num_Obs,1/))) ! weighted innovation
        
        anl_incr_at_Grid = anl_incr_interm(1,1)
        anl_at_Grid = back_at_Grid + anl_incr_at_Grid

        RETURN

    END SUBROUTINE Snow_DA_OI

    
    Subroutine snow_DA_EnKF_LocCov(RLA_jndx, RLO_jndx, Orog_jndx,   &
        num_Obs_1, num_Obs, num_stn, jindx, ens_size, &
        Lat_Obs, Lon_Obs, Ele_Obs,                              &
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        assim_IMS,  rcov_localize,                  &
        SNOFCS_Inp_Ens,          &       !LENSFC, 
        loc_nearest_Obs, SNOFCS_atObs_ens,   &
        Stdev_Obs_depth, Stdev_Obs_ims,     &   !stdev_back,     &
        obs_Array,                          &
        obs_Innov, incr_atGrid_ens, anl_at_Grid_ens)  !obs_Innov_ens
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        Real, Intent(In)        :: RLA_jndx, RLO_jndx, Orog_jndx
        Integer, Intent(In)     :: num_Obs_1, num_Obs, num_stn, jindx, ens_size  !, LENSFC 
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs), Ele_Obs(num_Obs)
        Real, Intent(In)        :: L_horz, h_ver  !L_horz in Km, h_ver in m
        LOGICAL, Intent(In)     :: assim_IMS, rcov_localize    !, ens_inflate
        
        Real, Intent(In)        :: SNOFCS_Inp_Ens(ens_size)   !, LENSFC)
        Integer, Intent(In)     :: loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In)        :: SNOFCS_atObs_ens(ens_size, num_stn)
        Real, Intent(In)        :: Stdev_Obs_depth, Stdev_Obs_ims
        REAL, INTENT(In)        :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov(num_Obs)  !, ens_size)
        Real, INTENT(Out)   :: incr_atGrid_ens(ens_size+1)
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_h_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_cov(num_Obs, num_Obs), R_cov(num_obs, num_obs)
        Real(dp)            :: Pzz_cov_inv(num_obs, num_obs), K_gain(1, num_obs)
        Real(dp)            :: obs_Innov_ens(num_Obs, ens_size)
        Real(dp)            :: innov_at_Grid_ens_interm(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size) 

        !Integer             :: indx, zndx 
        Integer :: indx, jndx, zndx    
        ! Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)    :: R_cov_loc(num_obs)  !,num_obs)
        
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2          
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! Do jndx = 1, num_Obs 
        !     d_latArr = (Lat_Obs_rad(jndx) - Lat_Obs_rad) / 2.
        !     d_lonArr = (Lon_Obs_rad(jndx) - Lon_Obs_rad) / 2.
        !     haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(Lat_Obs_rad(jndx)) * sin(d_lonArr)**2
        !     Where (haversinArr > 1) haversinArr = 1
        !     rjk_distArr(jndx,:) = 2 * earth_rad * asin(sqrt(haversinArr))     ! rjk, k = 1, Num obs for a given j
        !     zjk_distArr(jndx,:) = Ele_Obs(jndx) - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        ! End do
        !Corr_j_k = (1+rjk/L)exp(-rjk/L)exp(-(zjk/h)^2)
        ! if (print_debug) then
        !     print*, "Dist for Back corr at obs pts"
        !     print*, rjk_distArr
        !     print*, "Vertical dist for Back corr at obs pts"
        !     print*, zjk_distArr
        ! endif
        ! factor for localization of R covariance 
        ! R_cov_loc = exp(1. * (rjk_distArr/L_horz)**2) !L_horz in Km, h_ver in m (1. + rjk_distArr/L_horz) * 
        ! R_cov_loc = R_cov_loc * exp(1. * (zjk_distArr/h_ver)**2)
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))     ! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif      
        ! factor for localization of R covariance 
        R_cov_loc = exp(1. * (l_distArr/L_horz)**2) !L_horz in Km, h_ver in m (1. + rjk_distArr/L_horz) * 
        R_cov_loc = R_cov_loc * exp(1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Obs cov localization factor"
            print*, R_cov_loc
        endif   
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:)       !, jindx)
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
        End do
        if(assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:)     !, jindx)
        X_ave_State = SUM(X_state(1,:)) / ens_size
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  
        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do
        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)
        ! State at obs Cov: Pzz_h = Xh'*Xh'_T [m,m] 
        Pzz_h_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_h_cov = Pzz_h_cov / (ens_size - 1)
        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov before localization"
            print*, R_cov
        endif
        ! R cov localization 
        ! R_cov = R_cov * R_cov_loc
        if(rcov_localize) then 
            Do indx = 1, num_Obs
                R_cov(indx, indx) = R_cov(indx, indx) * R_cov_loc(indx)
            end do
            if (print_debug) then
                print*, "Obs cov after localization"
                print*, R_cov
            endif
        endif
        ! Innovvation Cov: Pzz = Pzz_h + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_h_cov + R_cov
        Pzz_cov_inv = inv(Pzz_cov)                       
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m]
        K_gain = matmul(Pxz_cov, Pzz_cov_inv)        
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, Pzz_cov_inv
            print*, "Gain vector"
            print*, K_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov_ens(zndx, :) = obs_Array(zndx) - Xh_state_atObs(zndx, :)
            obs_Innov(zndx) = SUM(obs_Innov_ens(zndx, :)) / ens_size
        End do
        ! Innovation at model dX = K*d [1,m][m,E] = [1,E]
        innov_at_Grid_ens_interm = matmul(K_gain, obs_Innov_ens)
        ! Update (1,E)
        anl_at_Grid_ens_interm = X_state + innov_at_Grid_ens_interm
        !
        incr_atGrid_ens(1:ens_size) = RESHAPE(innov_at_Grid_ens_interm,(/ens_size/))
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        ! ens mean
        anl_at_Grid_ens(ens_size + 1) = SUM(anl_at_Grid_ens(1:ens_size)) / ens_size
        incr_atGrid_ens(ens_size + 1) = SUM(incr_atGrid_ens(1:ens_size)) / ens_size
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov_ens
            print*, "Increment at grid"
            print*, incr_atGrid_ens
            print*, "Analysis"
            print*, anl_at_Grid_ens
        endif

        RETURN

    End Subroutine snow_DA_EnKF_LocCov

    Subroutine snow_DA_EnKF(assim_IMS,                    &
        jindx, ens_size, LENSFC, SNOFCS_Inp_Ens,          &
        num_Obs_1, num_Obs, loc_nearest_Obs, SNOFCS_atObs_ens,   &
        Stdev_Obs_depth, Stdev_Obs_ims,                  &
        obs_Array,                          &
        obs_Innov_ens, innov_at_Grid_ens, anl_at_Grid_ens)
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)

        LOGICAL, Intent(In) :: assim_IMS
        Integer, Intent(In) :: jindx, ens_size, LENSFC 
        Real, Intent(In)        :: SNOFCS_Inp_Ens(ens_size, LENSFC)
        Integer, Intent(In) :: num_Obs_1, num_Obs, loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In)        :: SNOFCS_atObs_ens(ens_size, num_Obs_1)
        Real, Intent(In)        :: Stdev_Obs_depth, Stdev_Obs_ims
        REAL, INTENT(In)    :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov_ens(num_Obs, ens_size)
        Real, INTENT(Out)   :: innov_at_Grid_ens(ens_size+1)
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_h_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_cov(num_Obs, num_Obs), R_cov(num_obs, num_obs)
        Real(dp)            :: Pzz_cov_inv(num_obs, num_obs), K_gain(1, num_obs)
        Real(dp)            :: innov_at_Grid_ens_interm(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size) 

        Integer             :: indx, zndx 
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:, jindx)
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
        End do
        if(assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:, jindx)
        X_ave_State = SUM(X_state(1,:)) / ens_size
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  
        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do
        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)
        ! State at obs Cov: Pzz_h = Xh'*Xh'_T [m,m] 
        Pzz_h_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_h_cov = Pzz_h_cov / (ens_size - 1)
        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov"
            print*, R_cov
        endif
        ! Innovvation Cov: Pzz = Pzz_h + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_h_cov + R_cov
        Pzz_cov_inv = inv(Pzz_cov)                       
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m]
        K_gain = matmul(Pxz_cov, Pzz_cov_inv)        
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, Pzz_cov_inv
            print*, "Gain vector"
            print*, K_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov_ens(zndx, :) = obs_Array(zndx) - Xh_state_atObs(zndx, :)
        End do
        ! Innovation at model dX = K*d [1,m][m,E] = [1,E]
        innov_at_Grid_ens_interm = matmul(K_gain, obs_Innov_ens)
        ! Update (1,E)
        anl_at_Grid_ens_interm = X_state + innov_at_Grid_ens_interm
        !
        innov_at_Grid_ens(1:ens_size) = RESHAPE(innov_at_Grid_ens_interm,(/ens_size/))
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        ! ens mean
        anl_at_Grid_ens(ens_size + 1) = SUM(anl_at_Grid_ens(1:ens_size)) / ens_size
        innov_at_Grid_ens(ens_size + 1) = SUM(innov_at_Grid_ens(1:ens_size)) / ens_size
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov_ens
            print*, "Innov at grid"
            print*, innov_at_Grid_ens
            print*, "Analysis"
            print*, anl_at_Grid_ens
        endif

        RETURN

    End Subroutine snow_DA_EnKF
    
    Subroutine snow_DA_EnSRF(RLA_jndx, RLO_jndx, Orog_jndx,   &
        num_Obs_1, num_Obs, num_stn, jindx, ens_size,   &
        L_horz, h_ver,                                      &   !L_horz in Km, h_ver in m
        Lat_Obs, Lon_Obs, Ele_Obs,                              &
        assim_IMS, rcov_localize, ens_inflate,                    &
        SNOFCS_Inp_Ens,          &     !LENSFC, 
        loc_nearest_Obs, SNOFCS_atObs_ens,   &
        Stdev_Obs_depth, Stdev_Obs_ims,  stdev_back,                &
        obs_Array,                          &
        obs_Innov, incr_atGrid_ens, anl_at_Grid_ens)
        
        IMPLICIT NONE
        !USE intrinsic::ieee_arithmetic
        integer, parameter :: dp = kind(1.d0)
       
        
        Real, Intent(In)        :: RLA_jndx, RLO_jndx, Orog_jndx
        Integer, Intent(In)     :: num_Obs_1, num_Obs, num_stn, jindx, ens_size
        Real, Intent(In)        :: L_horz, h_ver  !L_horz in Km, h_ver in m
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs), Ele_Obs(num_Obs)

        LOGICAL, Intent(In)     :: assim_IMS, rcov_localize, ens_inflate
        Real, Intent(In)        :: SNOFCS_Inp_Ens(ens_size)  !, LENSFC)
        Integer, Intent(In)     :: loc_nearest_Obs(num_Obs_1)        
        Real, Intent(In)        :: SNOFCS_atObs_ens(ens_size, num_stn)
        Real, Intent(In)        :: Stdev_Obs_depth, Stdev_Obs_ims, stdev_back
        REAL, INTENT(In)        :: obs_Array(num_Obs)

        Real, INTENT(Out)   :: obs_Innov(num_Obs)
        Real, INTENT(Out)   :: incr_atGrid_ens(ens_size+1)
        Real, INTENT(Out)   :: anl_at_Grid_ens(ens_size+1)          

        Real                :: X_state(1, ens_size), Xh_state_atObs(num_Obs, ens_size)
        Real(dp)            :: X_ave_State, Xh_ave_State(num_obs)
        Real(dp)            :: X_ens_Anomaly(1, ens_size), Xh_ens_Anomaly(num_Obs, ens_size)
        Real(dp)            :: Pxz_cov(1, num_obs), Pzz_cov(num_Obs, num_Obs)
        Real(dp)            :: Pzz_s(num_Obs, num_Obs), Pzz_s_i(num_Obs, num_Obs)
        Real(dp)            :: R_cov(num_obs, num_obs), R_cov_s(num_obs, num_obs)
        Real(dp)            :: K_gain(1, num_obs), K_anom_gain(1, num_obs)   !, Pzz_cov_inv(num_obs, num_obs)
        Real(dp)            :: innov_ens_Anomaly(1, ens_size)
        Real(dp)            :: innov_atGrid_ensM_interm(1, 1) 
        Real(dp)            :: Xu_ens_Anomaly_upd(1, ens_size)
        Real(dp)            :: anl_at_Grid_ens_interm(1, ens_size)

        !Integer             :: indx, zndx 
        Integer :: indx, jndx, zndx    
        ! Real    :: rjk_distArr(num_Obs, num_Obs), zjk_distArr(num_Obs, num_Obs)    
        Real    :: l_distArr(num_Obs), h_distArr(num_Obs), haversinArr(num_Obs)
        Real(dp)            :: R_cov_loc(num_obs)  !,num_obs)
        ! this is only used in 'inflation of ensembles'
        Real(dp)            :: Pxx_cov(1, 1) 
        ! Real(dp)            :: R_div_Pzz(num_obs, num_obs)
        Real(dp)            :: alpha_anom_gain(1, num_obs)
        
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    ::  Lon_Obs_2(num_Obs)      !RLO_2_jndx,
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real    :: ens_inflation_fact
        Integer :: info
        Integer :: print_i

        print_i = 0
        
        !Lon between -180 and 180 for some inputs
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2          
        ! https://en.wikipedia.org/wiki/Haversine_formula
        d_latArr = (RLA_rad_jndx - Lat_Obs_rad) / 2.
        d_lonArr = (RLO_rad_jndx - Lon_Obs_rad) / 2.
        haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad) * cos(RLA_rad_jndx) * sin(d_lonArr)**2
        Where (haversinArr > 1) haversinArr = 1.
        l_distArr = 2 * earth_rad * asin(sqrt(haversinArr))     ! rjk, k = 1, Num obs for a given j
        h_distArr = Orog_jndx - Ele_Obs       ! zjk, k = 1, Num obs for a given j
        if (print_debug) then
            print*, "Horz Dist for Back corr at obs pts and model grid"
            print*, l_distArr
            print*, "Vertical dist for Back corr at obs pts and model grid"
            print*, h_distArr
        endif      
        ! factor for localization of R covariance 
        R_cov_loc = exp(1. * (l_distArr/L_horz)**2) !L_horz in Km, h_ver in m 
        R_cov_loc = R_cov_loc * exp(1. * (h_distArr/h_ver)**2)
        if (print_debug) then
            print*, "Obs cov localization factor"
            print*, R_cov_loc
        endif   
        
        ! Background states X and Xh (Dim [E] [m, E], E = ensemble size, m = num obs)
        X_state(1, :) = SNOFCS_Inp_Ens(:)  !, jindx)
        
        Do zndx = 1, num_Obs_1
            Xh_state_atObs(zndx, :) = SNOFCS_atObs_ens(:, loc_nearest_Obs(zndx))
            ! Xh_ave_State(zndx) = SNOFCS_atObs_ens(ens_size+1, loc_nearest_Obs(zndx))
        End do
        if (assim_IMS) Xh_state_atObs(num_Obs, :) = SNOFCS_Inp_Ens(:)  !, jindx)
        ! Xh_ave_State(num_Obs) = SNOFCS_Inp_Ens(ens_size+1)
        
        !ens ave
        X_ave_State = SUM(X_state(1,:)) / ens_size   !SNOFCS_Inp_Ens(ens_size+1)   !
        Do zndx = 1, num_Obs
            Xh_ave_State(zndx) = SUM(Xh_state_atObs(zndx, :)) / ens_size
        End do  

        ! ens anomaly X' and Xh'    
        X_ens_Anomaly(1,:) = X_state(1,:) - X_ave_State
        Do zndx = 1, num_Obs
            Xh_ens_Anomaly(zndx, :) = Xh_state_atObs(zndx, :) - Xh_ave_State(zndx)
        End do

        ! infl
        ! because the ensemble spread from GFS forecast snow depth is narrow
        ! we use simple std.dev ratio to expand the ensemble anomalies 
        ! expect this to lead to background std.dev comparable to that used in the OI DA
        ens_inflation_fact = 1.0
        if (ens_inflate) then             
            Pxx_cov = matmul(X_ens_Anomaly, TRANSPOSE(X_ens_Anomaly)) ![1,E][E,1] = [1,1]
            Pxx_cov = Pxx_cov / (ens_size - 1)
            if (Pxx_cov(1,1) > 0.1) then 
                ens_inflation_fact = stdev_back / sqrt(Pxx_cov(1,1))
                ! if (ens_inflation_fact > 1000) then
                !     print*, "ensemble inflation factor: ", ens_inflation_fact
                !     print*, " background cov: " , Pxx_cov(1,1) 
                !     print*, "X_state: " , X_state
                !     print*, "Xh_state_atObs: " , Xh_state_atObs
                ! endif
            else    !if (Pxx_cov(1,1) .eq. 0) 
                ens_inflation_fact = 1.0
                !print*, "zero background cov: " , Pxx_cov(1,1)                
            endif   
            if (ens_inflation_fact > 1) then ! .and. (ens_inflation_fact < 100)) then
                X_ens_Anomaly = ens_inflation_fact * X_ens_Anomaly
                Xh_ens_Anomaly = ens_inflation_fact * Xh_ens_Anomaly
            endif
        end if

        ! State-obs xCov: Pxz = X'*Xh'_T [1,m] (_T = transpose)
        Pxz_cov = matmul(X_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![1,E][E,m] = [1,m]
        Pxz_cov = Pxz_cov / (ens_size - 1)

        ! Observation covariance 
        ! R = stdev_o*stdev_o * I , I = Identitity matrix
        R_cov = 0.
        Do indx = 1, num_Obs
            R_cov(indx, indx) = Stdev_Obs_depth * Stdev_Obs_depth
        end do
        if (assim_IMS) R_cov(num_Obs, num_Obs) = Stdev_Obs_ims * Stdev_Obs_ims
        if (print_debug) then
            print*, "Obs cov before localization"
            print*, R_cov
        endif
        ! R cov localization 
        ! R_cov = R_cov * R_cov_loc
        if(rcov_localize) then 
            Do indx = 1, num_Obs
                R_cov(indx, indx) = R_cov(indx, indx) * R_cov_loc(indx)
            end do
            if (print_debug) then
                print*, "Obs cov after localization"
                print*, R_cov
            endif
        endif  
        ! State at obs Cov: Pzz_cov = Xh'*Xh'_T [m,m] 
        Pzz_cov = matmul(Xh_ens_Anomaly, TRANSPOSE(Xh_ens_Anomaly)) ![m,E][E,m] = [m,m]
        Pzz_cov = Pzz_cov / (ens_size - 1)        
        ! Innovvation Cov: Pzz = Pzz + R = Xh'*Xh'_T + R [m,m] 
        Pzz_cov = Pzz_cov + R_cov
        ! Kalman gain K = Pxz*Pzz_i [1,m][m,m] = [1,m] !𝐾=  𝑃_𝑥𝑧 𝑃_𝑧𝑧^(−1)
        K_gain = matmul(Pxz_cov, inv(Pzz_cov))   
        if (print_debug) then
            print*, "Innov cov"
            print*, Pzz_cov
            print*, "inverse of Innov cov"
            print*, inv(Pzz_cov)
            print*, "Ens mean Gain"
            print*, K_gain
        endif
        ! Innovation at obs d = Z - Xh [m,E]
        Do zndx = 1, num_Obs
            obs_Innov(zndx) = obs_Array(zndx) - Xh_ave_State(zndx)
        End do
        ! Ens mean Innovation at model dX = K*d [1,m][m,1] = [1,1]
        innov_atGrid_ensM_interm = matmul(K_gain, RESHAPE(obs_Innov, (/num_Obs, 1/))) 
        incr_atGrid_ens(ens_size+1) = innov_atGrid_ensM_interm(1, 1) 
        ! Ens mean Update 
        anl_at_Grid_ens(ens_size + 1) = X_ave_State + incr_atGrid_ens(ens_size+1)
        if (print_debug) then
            print*, "Obs Innov"
            print*, obs_Innov
            print*, "grid ens mean increment"
            print*, incr_atGrid_ens(ens_size+1)
            print*, "Analysis ensemble mean"
            print*, anl_at_Grid_ens(ens_size + 1)
        endif
!        K' = PH_T [{sqrt(HPH^T + R)}^-1]^T * [sqrt(HPH^T + R) + sqrt(R)]^-1
!        K' = Pxz [{sqrt(Pzz)}^-1]^T * [sqrt(Pzz) + sqrt(R)]^-1
!        K' = Pxz [{Pzz_s}^-1]^T * [Pzz_s + R_s]^-1
!        K' = Pxz [Pzz_s_i]^T * [Pzz_s_Pl_R_s]^-1
!        K' = Pxz Pzz_s_i_T * Pzz_s_Pl_R_s_i
        !sqrt of R
        print_i = -99
        call matsqrt_sub(print_i, R_cov, R_cov_s, info)  !Rs_cov = matsqrt(R_cov)
        if (info /= 0) then
            stop 'Matrix factorization failed at R cov!'
        end if
        ! sqrt pzz
        call matsqrt_sub(print_i, Pzz_cov, Pzz_s, info)    !Pzz_s = matsqrt(Pzz_cov)
        if (info /= 0) then
            print*, "grid index: ", jindx
            print*, "ensemble inflation factor: ", ens_inflation_fact
            print*, "esemble mean: ", Xh_ave_State
            print*, "Xh_ens_anomaly: ", Xh_ens_Anomaly
            print*, "Pzz_cov: ", Pzz_cov
            stop "Matrix factorization failed at Pzz_cov!"
        end if
        ! ! inv
        Pzz_s_i = inv(Pzz_s)                
        ! 𝑲^′=  𝑃_𝑥𝑧 ∗ (𝑃_𝐳𝐳_𝒔_𝒊)^𝑻 ∗ (𝑃_(𝐳𝑧_𝒔)+𝑹_𝒔 )^(−𝟏)
        K_anom_gain = matmul(Pxz_cov, TRANSPOSE(Pzz_s_i))
        K_anom_gain = matmul(K_anom_gain, inv(Pzz_s + R_cov_s))       
        if (print_debug) then
            print*, "Sqrt Innov cov"
            print*, Pzz_s
            print*, "inverse of Sqrt Innov cov"
            print*, Pzz_s_i
            print*, "Sqrt Obs cov"
            print*, R_cov_s
            print*, "Ens anomaly Gain"
            print*, K_anom_gain
        endif

        if (ens_inflate .and. ens_inflation_fact > 1) then ! .and. (ens_inflation_fact < 100)) then
 !##TBC is correct??
            X_ens_Anomaly = X_ens_Anomaly/ens_inflation_fact
            Xh_ens_Anomaly = Xh_ens_Anomaly/ens_inflation_fact 
        endif
        ! Ens anomaly Innovation = 𝑲′ * 𝑿_𝒉′ [1,m][m,E] = [1,E]
        innov_ens_Anomaly = matmul(K_anom_gain, Xh_ens_Anomaly)
        ! Ens anomaly Update (1,E) 𝑿′𝒂= 𝑿′ − 𝑲′*𝑿_𝒉′ 
        Xu_ens_Anomaly_upd = X_ens_Anomaly - innov_ens_Anomaly
        ! non-ens_mean, ensemble members update
        anl_at_Grid_ens_interm(1,:) = Xu_ens_Anomaly_upd(1, :) + anl_at_Grid_ens(ens_size + 1)
        anl_at_Grid_ens(1:ens_size) = RESHAPE(anl_at_Grid_ens_interm,(/ens_size/))
        !10.26.21 ensemble increments
        incr_atGrid_ens(1:ens_size) = anl_at_Grid_ens(1:ens_size) - X_state(1, :)
        if (print_debug) then
            print*, "Innov at grid ens anomaly"
            print*, innov_ens_Anomaly
            print*, "Analysis ensemble "
            print*, anl_at_Grid_ens
            print*, "Ensemble increments"
            print*, incr_atGrid_ens
        endif
        
        RETURN

    End Subroutine snow_DA_EnSRF

    ! ! LAPACK.
    function matsqrt(A) result(Asqrt)
        
        IMPLICIT NONE
        
        integer, parameter :: dp = kind(1.d0)

        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)) :: Asqrt
        integer :: N, info
     
        ! External procedures defined in LAPACK
        external DPOTRF
    
        ! Backup
        Asqrt = A
        N = size(A,1)
    
        call DPOTRF ('L', N, Asqrt, N, info)
    
        if (info /= 0) then
            stop 'Matrix factorization failed!'
        end if

    end function matsqrt

    ! ! LAPACK.
    subroutine matsqrt_sub(print_i, A, Asqrt, info) !result()
        
        IMPLICIT NONE
        
        integer, parameter :: dp = kind(1.d0)
        
        integer, intent(in) :: print_i
        real(dp), dimension(:,:), intent(in) :: A
        real(dp), dimension(size(A,1),size(A,2)), intent(out) :: Asqrt
        integer, intent(out) :: info
        integer :: N, ir, jc
     
        ! External procedures defined in LAPACK
        external DPOTRF
    
        ! Backup
        Asqrt = A
        N = size(A,1)
    
        call DPOTRF ('L', N, Asqrt, N, info)
        
        if (print_i == 1) print*, "Asqrt = ", Asqrt

        Do ir = 1, N - 1 
            Do jc = ir + 1, N
                Asqrt(ir, jc) = 0.0
            Enddo
        End do

        if (print_i == 1) print*, "Asqrt after filling upper diag 0 = ", Asqrt
    
        if (info /= 0) then
            stop 'Matrix factorization failed!'
        end if

    end subroutine matsqrt_sub

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

    subroutine nearest_IMS_Locations(num_Obs, max_num_loc, max_distance, ims_max_ele, &
        RLA_jndx, RLO_jndx, RLA, RLO, OROG, SNOFCS_atObs, SWEFCS, SNUP_Array, & 
        SNO_IMS_at_Grid, num_loc, Loc_backSt_atObs)
        !intp_mode) 
        !  SWE_atObs  SWE_back
        IMPLICIT NONE
        
        Integer, Intent(In)     :: num_Obs, max_num_loc
        Real, Intent(In)        :: max_distance, ims_max_ele   ! radius_of_influence
        Real, Intent(In)        :: RLA_jndx, RLO_jndx  !, OROG_jndx  ! don't want to alter these
        Real, Intent(In)        :: RLA(num_Obs), RLO(num_Obs), OROG(num_Obs)
        Real, Intent(In)        :: SNOFCS_atObs(num_Obs), SNO_IMS_at_Grid(num_Obs) 
        Real, Intent(In)        :: SWEFCS(num_Obs) , SNUP_Array(num_Obs) 
              
        Integer, Intent(Out)                 :: num_loc
        Integer, Allocatable, Intent(Out)    :: Loc_backSt_atObs(:)

        Integer :: indx, jndx, zndx, num_loc_counter        
        Real    :: distArr(num_Obs), haversinArr(num_Obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        !Real    :: Lon_Obs_2(num_Obs)   !RLO_2_jndx, 
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
    
        Real, Allocatable   :: dist_atObs(:)
        Real, Allocatable   :: dist_atObs_dummy(:)
        Integer, Allocatable   :: Loc_backSt_atObs_dummy(:)
        Real                 :: max_value

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize  !C_SIZE_T
        i_isize = SIZEOF(max_value)
        !print*, "size of real: ", i_
    
        !Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2

        ! shortest distance over sphere using great circle distance     
        RLA_rad_jndx =  pi_div_180 * RLA_jndx
        RLO_rad_jndx =  pi_div_180 * RLO_jndx
        Lat_Obs_rad =  pi_div_180 * RLA
        Lon_Obs_rad =  pi_div_180 * RLO   
        
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
               (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(indx))) .AND. &
               (OROG(indx) <= ims_max_ele) .AND. &
               ( .not.(SWEFCS(indx) >= SNUP_Array(indx) .AND. & 
                      SNO_IMS_at_Grid(indx) >= SNUP_Array(indx)) ) ) then                
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
               (.NOT. IEEE_IS_NAN(SNO_IMS_at_Grid(indx))) .AND. &
               (OROG(indx) <= ims_max_ele) .AND. &
               ( .NOT.(SWEFCS(indx) >= SNUP_Array(indx) .AND. & 
                      SNO_IMS_at_Grid(indx) >= SNUP_Array(indx)) ) ) then
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
            ! print*, "unsorted dist array ",  dist_atObs
            ! print*
            ! print*, "Sorted dist array ",  dist_atObs_dummy
            max_value = dist_atObs_dummy(max_num_loc)
            !loc_indices = findloc(dist_atObs <= max_value) 
            Allocate(Loc_backSt_atObs(max_num_loc)) 
            indx = 1     
            Do jndx = 1, num_loc
                if(dist_atObs(jndx) <= max_value) then
                    Loc_backSt_atObs(indx) = Loc_backSt_atObs_dummy(jndx)
                    indx = indx + 1
                endif
            End do
            ! print*, "unsorted locations ",  Loc_backSt_atObs_dummy
            ! print*
            ! print*, "Sorted locations ",  Loc_backSt_atObs
            num_loc = max_num_loc                       
            Deallocate(dist_atObs)
        else
            Allocate(Loc_backSt_atObs(num_loc)) 
            Loc_backSt_atObs = Loc_backSt_atObs_dummy
        endif

        Deallocate(dist_atObs_dummy, Loc_backSt_atObs_dummy)
        
        RETURN

    End subroutine nearest_IMS_Locations

    subroutine nearest_Observations_Locations(RLA_jndx, RLO_jndx,    &
        Lat_Obs, Lon_Obs, num_Obs, max_distance, max_num_loc,                           &
        Stdev_back, Stdev_Obs_depth, obs_tolerance,                 &
        SNOFCS_atObs, OBS_atOBs,                                 &
        Loc_backSt_atObs,  num_loc) !,      &LENSFC,
        !intp_mode) 
        !  SWE_atObs  SWE_back
        IMPLICIT NONE

        Real, Intent(In)        :: RLA_jndx, RLO_jndx  ! don't want to alter these
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs)
        Integer, Intent(In) :: num_Obs, max_num_loc
        Real, Intent(In)        :: max_distance   ! radius_of_influence
        Real, Intent(In)        :: Stdev_back, Stdev_Obs_depth, obs_tolerance
        Real, Intent(In)    :: SNOFCS_atObs(num_Obs)
        Real, Intent(In)    :: OBS_atOBs(num_Obs)
        Integer, Allocatable, Intent(Out)    :: Loc_backSt_atObs(:)
        !Real, Allocatable, Intent(Out)    :: dist_atObs(:)
        Integer, Intent(Out) :: num_loc
        
        Integer :: indx, jndx, zndx, num_loc_counter        
        Real    :: distArr(num_Obs), haversinArr(num_Obs)
        Real    :: d_latArr(num_Obs), d_lonArr(num_Obs)
        Real    :: Lon_Obs_2(num_Obs)   !RLO_2_jndx, 
        Real    :: RLA_rad_jndx, RLO_rad_jndx
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        Real                :: innov_criteria
        Real, Allocatable   :: dist_atObs_sorted(:)
        Real, Allocatable   :: dist_atObs_dummy(:)
        Integer, Allocatable   :: Loc_backSt_atObs_dummy(:), Loc_Sorted(:)
        Real                 :: max_value

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize  !C_SIZE_T
        i_isize = SIZEOF(max_value)
        !print*, "size of real: ", i_isize

        !Integer              :: loc_indices(max_num_loc)
        
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180      
    
        !Do jndx = 1, LENSFC 
        !if (RLO_jndx > 180) RLO_2_jndx = RLO_jndx - 360.0 
        !end do
        !print*, " here 1"
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! Do zndx = 1, num_Obs 
        !     if (Lon_Obs(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs(zndx)
        ! end do

        ! at each obs point compute its distance from RLA/RLO pairs 
        ! then find the position of the minimum
        !print*, " here 2"
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
        ! Do indx = 1, num_Obs 
        !     if (haversinArr(indx) > 1) haversinArr(indx) = 1 ! ensure numerical errors don't make h>1
        ! end do
        Where (haversinArr < 0) haversinArr = 0.
        ! Do indx = 1, num_Obs 
                !       if (haversinArr(indx) < 0) haversinArr(indx) = 0 ! ensure <0
                ! end do
        distArr = 2 * earth_rad * asin(sqrt(haversinArr))
        
        innov_criteria =  obs_tolerance * sqrt(Stdev_back**2 + Stdev_Obs_depth**2)
        ! 4.15.20: can you do the following without loop?       
        num_loc_counter = 0
        Do indx = 1, num_Obs    
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. & 
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx))) .AND. &    ! /OBS_atOBs(indx)
               (abs(OBS_atOBs(indx) - SNOFCS_atObs(indx)) < innov_criteria) ) then                
                num_loc_counter = num_loc_counter + 1
            endif
        End do
        num_loc = num_loc_counter
        Allocate(Loc_backSt_atObs_dummy(num_loc))
        Allocate(dist_atObs_dummy(num_loc))
        Allocate(Loc_Sorted(num_loc))
        jndx = 1
        Do indx = 1, num_Obs    
            if((distArr(indx) < max_distance) .AND. &
               (.NOT. IEEE_IS_NAN(SNOFCS_atObs(indx))) .AND. &
               (.NOT. IEEE_IS_NAN(OBS_atOBs(indx))) .AND. &    !/OBS_atOBs(indx)
               (abs(OBS_atOBs(indx) - SNOFCS_atObs(indx)) < innov_criteria) ) then 
                Loc_backSt_atObs_dummy(jndx) = indx
                dist_atObs_dummy(jndx) = distArr(indx)
                jndx = jndx  + 1
            endif
        End do
        
        Allocate(dist_atObs_sorted(num_loc))
        dist_atObs_sorted = dist_atObs_dummy
        i_len = num_loc
        Call QSORT(dist_atObs_sorted, i_len, i_isize, compar)
        ! print*, "unsorted dist array ",  dist_atObs
        ! print*
        ! print*, "Sorted dist array ",  dist_atObs_dummy
        Do indx = 1, num_loc
            Do jndx = 1, num_loc
                if(dist_atObs_sorted(indx) == dist_atObs_dummy(jndx)) then
                    Loc_Sorted(indx) = Loc_backSt_atObs_dummy(jndx)
                    exit
                endif
                if(jndx == num_loc) print*, "Warning, index not found for nearest Obs"
            End do            
        Enddo
        ! if num of obs > 50, choose the 50 closest obs
        if (num_loc > max_num_loc) then             
            Allocate(Loc_backSt_atObs(max_num_loc)) 
            Loc_backSt_atObs = Loc_Sorted(1:max_num_loc)
            num_loc = max_num_loc                  
        else
            Allocate(Loc_backSt_atObs(num_loc)) 
            Loc_backSt_atObs = Loc_Sorted
        endif
        Deallocate(dist_atObs_sorted)
        Deallocate(dist_atObs_dummy, Loc_backSt_atObs_dummy, Loc_Sorted)
        
        ! if (num_loc > max_num_loc) then 
        !     Allocate(dist_atObs(num_loc))
        !     dist_atObs = dist_atObs_dummy
        !     i_len = num_loc
        !     Call QSORT(dist_atObs_dummy, i_len, i_isize, compar)
        !     print*, "unsorted dist array ",  dist_atObs
        !     print*
        !     print*, "Sorted dist array ",  dist_atObs_dummy            
        !     max_value = dist_atObs_dummy(max_num_loc)
        !     !loc_indices = findloc(dist_atObs <= max_value) 
        !     Allocate(Loc_backSt_atObs(max_num_loc)) 
        !     indx = 1     
        !     Do jndx = 1, num_loc
        !         if(dist_atObs(jndx) <= max_value) then  !if(dist_atObs(jndx) <= max_value) then
        !             Loc_backSt_atObs(indx) = Loc_backSt_atObs_dummy(jndx)
        !             indx = indx + 1
        !         endif
        !     End do
        !       Deallocate(dist_atObs)
        !     ! print*, "unsorted locations ",  Loc_backSt_atObs_dummy
        !     ! print*
        !     ! print*, "Sorted locations ",  Loc_backSt_atObs
        !     num_loc = max_num_loc                  
        ! else
        !     Allocate(Loc_backSt_atObs(num_loc)) 
        !     Loc_backSt_atObs = Loc_backSt_atObs_dummy
        ! endif      
        ! Deallocate(dist_atObs_dummy, Loc_backSt_atObs_dummy)
        
        RETURN

    End subroutine nearest_Observations_Locations

    subroutine select_max_num_obs(dist_atObs, SNOFCS_atObs, OBS_atOBs, &
        num_loc, max_num_loc,                           &
        out_SNOFCS_atObs, out_OBS_atOBs) 

        IMPLICIT NONE
    
        !Integer(2), External :: compar
        !USE intrinsic::ieee_arithmetic
        
        Real, Intent(In)    :: SNOFCS_atObs(num_loc), OBS_atOBs(num_loc), dist_atObs(num_loc)
        Integer, Intent(In) :: num_loc, max_num_loc
        Real, Intent(Out)    :: out_SNOFCS_atObs(max_num_loc), out_OBS_atOBs(max_num_loc)
        Real                 :: dist_atObs_dummy(num_loc)
        Real                 :: max_value
        !Integer              :: loc_indices(max_num_loc)
        Integer :: indx, jndx

        INTEGER(SIZEOF_SIZE_T)   :: i_len, i_isize    !C_SIZE_T
        i_isize = SIZEOF(max_value)
        i_len = num_loc
        
        dist_atObs_dummy = dist_atObs
        Call QSORT (dist_atObs_dummy, i_len, i_isize, compar)
        max_value = dist_atObs_dummy(max_num_loc)
        !loc_indices = findloc(dist_atObs <= max_value) 
        indx = 1       
        Do jndx = 1, max_num_loc
            if(dist_atObs(jndx) <= max_value) then
                out_SNOFCS_atObs(indx) = SNOFCS_atObs(jndx)
                out_OBS_atOBs(indx) = OBS_atOBs(jndx)
                indx = indx + 1
            endif
        End do
        
    End subroutine select_max_num_obs

     ! Gets obs SND from IMS based on NoahMP 'depletion curve' 
    subroutine CalcSWEFromSnowCover_NoahMP(SNCOV_IMS, VETFCS_in, LENSFC, &
                                               SND_IMS_at_Grid, SN_Ratio)
        
        IMPLICIT NONE
        !
        Real, Intent(In)        :: SNCOV_IMS(LENSFC), VETFCS_in(LENSFC), &
                                   SN_Ratio(LENSFC)
        INTEGER, Intent(In)     :: LENSFC
        Real, Intent(Out)       :: SND_IMS_at_Grid(LENSFC)
        
        INTEGER                 :: VETFCS(LENSFC)
        REAL                    :: mfsno_table(20), scffac_table(20)
        REAL                    :: mfsno, scffac, fmelt

        Integer                 :: indx, vtype_int

        !Fill background values to nan (to differentiate those that don't have value)
        SND_IMS_at_Grid = IEEE_VALUE(SND_IMS_at_Grid, IEEE_QUIET_NAN)
   
        ! NOTE: this is an empirical function from NoahMP
        mfsno_table = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00, &
                         2.00, 2.00, 2.00, 3.00, 3.00, 4.00, 4.00, &
                         2.50, 3.00, 3.00, 3.50, 3.50, 3.50 /)

        scffac_table = (/ 0.005, 0.005, 0.005, 0.005, 0.005, 0.008, &
                          0.008, 0.010, 0.010, 0.010, 0.010, 0.007, 0.021, &  
                          0.013, 0.015, 0.008, 0.015, 0.015, 0.015, 0.015 /)

        VETFCS = INT(VETFCS_in)
        
        Do indx = 1, LENSFC  
            vtype_int = VETFCS(indx)
            mfsno  =  mfsno_table(vtype_int)
            scffac = scffac_table(vtype_int)
            if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
                if (SNCOV_IMS(indx) < 0.001) then
                    SND_IMS_at_Grid(indx) = 0.0
                else
                    !SN_Ratio = SWE/SND   !bdsno = sneqv/snowh = (SWE*1000/SND)                    
                    fmelt = (10*SN_Ratio(indx))**mfsno      ! (bdsno/100.)**mfsno
                    SND_IMS_at_Grid(indx) = scffac * fmelt * ATANH(SNCOV_IMS(indx)) 
                endif  
            endif
        end do  
        RETURN
    
    END SUBROUTINE CalcSWEFromSnowCover_NoahMP

     ! Gets obs SWE from IMS based on exponential/log 'depletion curve' 
    subroutine CalcSWEFromSnowCover(SNCOV_IMS, VETFCS_in, LENSFC, &
                                    SWE_IMS_at_Grid, SNUP_Array)

! CSD - even though it's labelled as snow depth in the Noah code, this
! routine is performed on SWE, not SND. So I think your output should be SWE.
! Also, this is not an observation operator, hence the rename
! do my changeg affect your inversion?
        
        IMPLICIT NONE
        !
        Real, Intent(In)        :: SNCOV_IMS(LENSFC), VETFCS_in(LENSFC)
        INTEGER, Intent(In)     :: LENSFC
        Real, Intent(Out)       :: SWE_IMS_at_Grid(LENSFC), SNUP_Array(LENSFC)
        
        INTEGER                 :: VETFCS(LENSFC)
        REAL                    :: snupx(30), SNUP, SALP, RSNOW
        Integer                 :: indx, vtype_int

        ! Real, parameter    :: ims_threshold=0.5 ! threshold for converting IMS fSCA to binary 1, 0 values 

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
        Where(VETFCS==0) VETFCS = 7  !vtype_tile[vtype_tile==0] = 7
        
        Do indx = 1, LENSFC  
            if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
                if (VETFCS(indx) > 30) then
                    print*, " >30  vegclasses", VETFCS(indx)
                    Stop
                endif
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
    
     ! Gets obs snow depth from IMS based on threshold fSCA 
    SUBROUTINE Observation_Operator_IMS_fSCA_Threshold(SNCOV_IMS, SNOFCS, SNWDEN, assim_SWE,       &
                                LENSFC, ims_threshold,          &
                                SND_IMS_at_Grid) !,      &
    
        IMPLICIT NONE
        !
        Real, Intent(In)        :: SNCOV_IMS(LENSFC), SNOFCS(LENSFC), SNWDEN(LENSFC)
        Logical                         :: assim_SWE
        INTEGER                         :: LENSFC
        Real, Intent(In)    :: ims_threshold
    
        Real, Intent(Out)       :: SND_IMS_at_Grid(LENSFC)
        
        INTEGER :: indx
        Real, Parameter         :: SWE_Tolerance = 0.001    ! smallest swe value
    
        !Fill background values to nan (to differentiate those that don't have value)
        SND_IMS_at_Grid = IEEE_VALUE(SND_IMS_at_Grid, IEEE_QUIET_NAN)
    
        ! call debug_print("Here ", 1.)
    
        Do indx = 1, LENSFC  
            if (.NOT. IEEE_IS_NAN(SNCOV_IMS(indx))) then
                if (SNCOV_IMS(indx) >= ims_threshold) then
                    ! ims snow, model no snow => obs=50 mm
                    if (SNOFCS(indx) < SWE_Tolerance) SND_IMS_at_Grid(indx) = 50.
                    ! ims snow, model snow => no assimilation
                else  ! IMS fSCA < thresh => Ims obs = 0
                    SND_IMS_at_Grid(indx) = 0.
                endif   ! all others nan        
            endif
        end do
        
        ! print*, "IMS Sndpth at model grids"
        ! print*, SND_IMS_at_Grid
        ! print*
        ! assim_SWE = True >> swe assimilated
        ! d = swe/snd
        if (assim_SWE) SND_IMS_at_Grid = SND_IMS_at_Grid * SNWDEN
    
        ! print*, "IMS SWE at model grids"
        ! print*, SND_IMS_at_Grid
    
        RETURN
    
    END SUBROUTINE Observation_Operator_IMS_fSCA_Threshold
    
    subroutine resample_to_model_tiles_intrp(data_grid_ims, data_grid_ims_ind, &
                                                nlat_ims, nlon_ims, n_lat, n_lon, num_sub, &  !myrank, &
                                                grid_dat)
                                                
        Use, Intrinsic :: IEEE_ARITHMETIC
    
        Implicit None
    
        Integer, Intent(In)     :: nlat_ims, nlon_ims, n_lat, n_lon, num_sub 
        Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, n_lon, n_lat) 
        Real, Intent(Out)       :: grid_dat(n_lon, n_lat)
    
        Integer   :: jc, jy, ix, num_loc_counter
        Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
        
        grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)
    
        Do jy=1, n_lat
        !print*, "process: ", myrank, "loop ", indx
    
            Do ix=1, n_lon
                
                num_loc_counter = data_grid_ims_ind(1, ix, jy)
                if (num_loc_counter < 1) then 
                    !print*, "no matching values!"
                    cycle
                end if
                !print*, "jy ", jy, " ix ", ix
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

    ! This uses G. Gayno's 'inside polygon function'
    ! subroutine resample_to_fv3_grid_inside(data_grid_ims, Lat_IMS_1D, Lon_IMS_1D, &
    !                                             nlat_ims, nlon_ims, n_lat, n_lon, &  !myrank, &
    !                                             grid_dat)
                                                
    !     Use, Intrinsic :: IEEE_ARITHMETIC
    
    !     Implicit None
    
    !     Integer, Intent(In)     :: nlat_ims, nlon_ims, n_lat, n_lon
    !     Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims)
    !     Real, Intent(In)        :: Lat_IMS_1D(nlat_ims), Lon_IMS_1D(nlon_ims)
    !     Real, Intent(Out)       :: grid_dat(n_lon, n_lat)
    
    !     Integer   :: jc, jy, ix, num_loc_counter
    !     Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
        
    !     grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)
    
    !     Do jy=1, n_lat
    !     !print*, "process: ", myrank, "loop ", indx
    
    !         Do ix=1, n_lon

    !             ! if(inside_a_polygon(LONI*D2R, LATI*D2R, 6, &
    !             !    LONO*D2R,LATO*D2R)) then
                    
    !             ! endif
                
               
    !             !print*, "jy ", jy, " ix ", ix
    !             grid_dat(ix, jy) = 0.
    !             Do jc = 2, num_loc_counter+1
    !                 lonlatcoord_ims = data_grid_ims_ind(jc, ix, jy) - 1 
    !                 latcoord_ims = lonlatcoord_ims / nlon_ims + 1
    !                 loncoord_ims = mod(lonlatcoord_ims, nlon_ims) + 1
    !                 if(latcoord_ims > nlat_ims) then
    !                     latcoord_ims = nlat_ims
    !                     print*, "Warning! lat coordinate outside domain boundary"
    !                 endif
    !                 if(loncoord_ims > nlon_ims) then
    !                     loncoord_ims = nlon_ims
    !                     print*, "Warning! lon coordinate outside domain boundary"
    !                 endif
    !                 grid_dat(ix, jy) =  grid_dat(ix, jy) + data_grid_ims(loncoord_ims, latcoord_ims)              
    !             End do
    
    !             grid_dat(ix, jy) =  grid_dat(ix, jy) / num_loc_counter ! first location, num obs
    
    !         End do
    
    !     End do
    
    !     return !grid_dat
    
    ! End subroutine resample_to_fv3_grid_inside

    subroutine resample_interpolate_toTarget(data_grid_ims, data_grid_ims_ind, &
                                            nlat_ims, nlon_ims, LENSFC, num_sub, &  !myrank, &
                                                grid_dat)
                                                
        Use, Intrinsic :: IEEE_ARITHMETIC
    
        Implicit None
    
        Integer, Intent(In)     :: nlat_ims, nlon_ims, LENSFC, num_sub 
        Integer, Intent(In)     :: data_grid_ims(nlon_ims, nlat_ims), data_grid_ims_ind(num_sub, LENSFC) 
        Real, Intent(Out)       :: grid_dat(LENSFC)
    
        Integer   :: jc, num_loc_counter, ixy !jy, ix, 
        Integer   :: lonlatcoord_ims, loncoord_ims, latcoord_ims
        
        grid_dat = IEEE_VALUE(grid_dat, IEEE_QUIET_NAN)
    
        ! Do jy=1, n_lat
        ! !print*, "process: ", myrank, "loop ", indx
    
        !     Do ix=1, n_lon
        Do ixy = 1, LENSFC

            num_loc_counter = data_grid_ims_ind(1, ixy)
            if (num_loc_counter < 1) then 
                !print*, "no matching values!"
                cycle
            end if
            !print*, "jy ", jy, " ix ", ix
            grid_dat(ixy) = 0.
            Do jc = 2, num_loc_counter+1
                lonlatcoord_ims = data_grid_ims_ind(jc, ixy) - 1 
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
                grid_dat(ixy) =  grid_dat(ixy) + data_grid_ims(loncoord_ims, latcoord_ims)              
            End do

            grid_dat(ixy) =  grid_dat(ixy) / num_loc_counter ! first location, num obs

            ! End do
    
        End do
    
        return !grid_dat
    
    End subroutine resample_interpolate_toTarget

    !This reads the whole IMS file and uses a-priori prepared indices to sample those wihin the grid cel
     SUBROUTINE Read_IMS_and_Resample_toTarget(inp_file, IMS_INDEXES_PATH, &  !fv3_index,
                    LENSFC, num_sub, IDIM, JDIM, tile_xy, Idim_xy, Jdim_xy, SNCOV_IMS)
                    ! Ele               &
                    
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        ! logical                        :: fv3_index
        CHARACTER(LEN=*), Intent(In)   :: inp_file, IMS_INDEXES_PATH !, dim_name
        INTEGER, Intent(In)            :: LENSFC, num_sub, IDIM, JDIM
        Integer, Intent(In)      ::  tile_xy(LENSFC), Idim_xy(LENSFC), Jdim_xy(LENSFC)
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)  
        Real, Intent(Out)       :: SNCOV_IMS(LENSFC)     
    
        CHARACTER(LEN=250)      :: inp_file_indices !, dim_name
        CHARACTER(LEN=3)        :: fvs_tile
        INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)    !SNCOV_IMS_1D(:), 
        Integer                 :: data_grid_ims_ind(num_sub, LENSFC) 
        ! Integer                 :: data_grid_ims_ind_in(num_sub, JDIM, IDIM) 
        ! Real                    :: grid_dat(LENSFC)        
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN_lat, DIM_LEN_lon
        LOGICAL                :: file_exists
        integer                :: myindx, ixy  
        integer                :: DUMMY1(6, num_sub, JDIM, IDIM)
        CHARACTER(LEN=1)       :: RANKCH
        
        write(fvs_tile, "(I3)") IDIM
        ! If(fv3_index) then
!index files from fv3 type files
            Do myindx = 1, 6 
                WRITE(RANKCH, '(I1.1)') myindx
                inp_file_indices = TRIM(IMS_INDEXES_PATH)//"/C"//TRIM(ADJUSTL(fvs_tile)) &
                                            //".IMS.Indices.tile"//RANKCH//".nc"
                INQUIRE(FILE=trim(inp_file_indices), EXIST=file_exists)
                if (.not. file_exists) then 
                    print *, 'error,file does not exist', &   
                            trim(inp_file_indices) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
                endif           

                ERROR=NF90_OPEN(TRIM(inp_file_indices), NF90_NOWRITE,NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )

                ERROR=NF90_INQ_VARID(NCID, "IMS_Indices", ID_VAR)
                CALL NETCDF_ERR(ERROR, 'READING IMS_Indices ID' )
                ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy1(myindx,:, :,:),  &
                        start = (/ 1, 1, 1 /), count = (/ num_sub, IDIM, JDIM/))
                CALL NETCDF_ERR(ERROR, 'READING snwdph' )
                
                ERROR = NF90_CLOSE(NCID)    
                CALL NETCDF_ERR(ERROR, 'closing file' )
            End do 
            ! 7.20.21 map fv3 to land
            Do ixy=1, LENSFC 
                data_grid_ims_ind(:,ixy)=dummy1(tile_xy(ixy),:,Idim_xy(ixy),Jdim_xy(ixy))
            end do
!         else
! ! read index file for grid/vector array           
!             inp_file_indices = TRIM(IMS_INDEXES_PATH)//"C"//TRIM(ADJUSTL(fvs_tile))//&
!                                                     ".IMS.Indices"//".nc" 
!             INQUIRE(FILE=trim(inp_file_indices), EXIST=file_exists)

!             if (.not. file_exists) then
!                     print *, 'Observation_Read_IMS_Full error, index file does not exist', &
!                             trim(inp_file_indices) , ' exiting'
!                     call MPI_ABORT(MPI_COMM_WORLD, 10) 
!             endif
        
!             ERROR=NF90_OPEN(TRIM(inp_file_indices),NF90_NOWRITE, NCID)
!             CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file_indices) )
        
!             ERROR=NF90_INQ_VARID(NCID, 'IMS_Indices', ID_VAR)
!             CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices ID' )
!             ERROR=NF90_GET_VAR(NCID, ID_VAR, data_grid_ims_ind, start = (/ 1, 1 /), &
!                                     count = (/ num_sub, LENSFC/))
!             CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV Indices' )
        
!             ERROR = NF90_CLOSE(NCID)
!         Endif
        ! print*, "IMS Indices 3D"
        ! print*, data_grid_ims_ind
        ! print*

! IMS raw data
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
        ! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        
        ERROR = NF90_CLOSE(NCID)
        ! print*, "IMS 2D"
        ! print*, SNCOV_IMS_2D_full
        ! print*
        
        ! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
        ! ALLOCATE(SNCOV_IMS_1D(DIM_LEN))
        ! SNCOV_IMS_1D = Reshape(SNCOV_IMS_2D_full, (/DIM_LEN/))
    
        ! print*, "IMS 1D"
        ! print*, SNCOV_IMS_1D
        ! print*
        Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
        Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1

        call resample_interpolate_toTarget(SNCOV_IMS_2D_full, data_grid_ims_ind, &
                                           DIM_LEN_lat, DIM_LEN_lon, LENSFC, num_sub, &  !myrank, &
                                           SNCOV_IMS)
    
        ! SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))
    
        DEALLOCATE(SNCOV_IMS_2D_full)
                  
        RETURN
        
     End SUBROUTINE Read_IMS_and_Resample_toTarget

     SUBROUTINE Observation_Read_IMS_Full(inp_file, inp_file_indices, &
                    MYRANK, n_lat, n_lon, num_sub, &
                    SNCOV_IMS)
                    ! Ele               &
                    
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)   :: inp_file, inp_file_indices !, dim_name
        INTEGER, Intent(In)            :: MYRANK, n_lat, n_lon, num_sub
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)  
        Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon)     
    
        INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D_full(:,:)    !SNCOV_IMS_1D(:), 
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
        ! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
        ! ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        
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
    
        ! print*, "IMS 2D"
        ! print*, SNCOV_IMS_2D_full
        ! print*
        
        ! DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
        ! ALLOCATE(SNCOV_IMS_1D(DIM_LEN))
        ! SNCOV_IMS_1D = Reshape(SNCOV_IMS_2D_full, (/DIM_LEN/))
    
        ! print*, "IMS 1D"
        ! print*, SNCOV_IMS_1D
        ! print*
        Where(SNCOV_IMS_2D_full /= 4) SNCOV_IMS_2D_full = 0
        Where(SNCOV_IMS_2D_full == 4) SNCOV_IMS_2D_full = 1
        
    
        ! print*, "IMS binary"
        ! print*, SNCOV_IMS_1D
        ! print*
    
        ! print*, "IMS Indices 3D"
        ! print*, data_grid_ims_ind
        ! print*
    
        call resample_to_model_tiles_intrp(SNCOV_IMS_2D_full, data_grid_ims_ind, &
                                           DIM_LEN_lat, DIM_LEN_lon, n_lat, n_lon, num_sub, &  !myrank, &
                                           grid_dat)
    
        SNCOV_IMS = Reshape(grid_dat, (/n_lat * n_lon/))
    
        DEALLOCATE(SNCOV_IMS_2D_full)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_IMS_Full


     SUBROUTINE Observation_Read_IMS_tile(inp_file, MYRANK, n_lat, n_lon, &
                                          SNCOV_IMS)
                                                  
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)   :: inp_file
        INTEGER, Intent(In)            :: MYRANK, n_lat, n_lon
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)  
        Real, Intent(Out)       :: SNCOV_IMS(n_lat * n_lon) 
        Real                    :: SNCOV_IMS_2D_full( n_lon, n_lat)
        
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN, DIM_LEN_lat, DIM_LEN_lon
    
        ERROR=NF90_OPEN(TRIM(inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(inp_file) )
    
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
        ERROR=NF90_INQ_VARID(NCID, 'imsfSCA', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/n_lon, n_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
    
        ! print*, "IMS 2D"
        ! print*, SNCOV_IMS_2D_full
        ! print*
        
        SNCOV_IMS = Reshape(SNCOV_IMS_2D_full, (/n_lat * n_lon/))
        ! print*, "IMS 1D"
        ! print*, SNCOV_IMS
        ! print*
                  
        RETURN
        
     End SUBROUTINE Observation_Read_IMS_tile
    
     SUBROUTINE Observation_Read_IMS(inp_file, MYRANK,  &
                    lat_min, lat_max, lon_min, lon_max, &
                    DIM_LEN,       & !_lat, DIM_LEN_lon,  &
                    SNCOV_IMS, Lat_IMS, Lon_IMS)
                    ! Ele               &
                    
        IMPLICIT NONE
    
        include 'mpif.h'                  
    
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)   :: inp_file !, dim_name
        INTEGER, Intent(In)            :: MYRANK
        REAL, Intent(In)                   :: lat_min, lat_max, lon_min, lon_max    
        ! ToDO: ims snow cover is of 'byte' type (Chcek the right one)
        INTEGER, Intent(Out)    :: DIM_LEN      
        INTEGER, ALLOCATABLE, Intent(Out)    :: SNCOV_IMS(:)    
        REAL, ALLOCATABLE, Intent(Out)     :: Lat_IMS(:), Lon_IMS(:)    !, Ele(:,:)
    
        INTEGER, ALLOCATABLE    :: SNCOV_IMS_2D(:,:), SNCOV_IMS_2D_full(:,:)
        REAL, ALLOCATABLE          :: Lat_IMS_1D(:), Lon_IMS_1D(:)      
        REAL, ALLOCATABLE          :: Lat_IMS_2D(:,:), Lon_IMS_2D(:,:)
        REAL, ALLOCATABLE          :: Lat_minmax_Diff(:), Lon_minmax_Diff(:)    
        
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, DIM_LEN_lat, DIM_LEN_lon
        INTEGER                :: indx, jndx, minlat_indx, maxlat_indx, minlon_indx, maxlon_indx
        INTEGER                :: icounter, jcounter, iincr, jincr
    
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
        ALLOCATE(Lat_IMS_1D(DIM_LEN_lat))
        ALLOCATE(Lon_IMS_1D(DIM_LEN_lon))
    
        ! print*, "initial IMS array size (lon, lat)= ", DIM_LEN_lon, " ",DIM_LEN_lat
        !ALLOCATE(Ele(DIM_LEN_lat, DIM_LEN_lon))
    
        ERROR=NF90_INQ_VARID(NCID, 'Band1', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SNCOV_IMS_2D_full, start = (/ 1, 1 /), &
                                count = (/ DIM_LEN_lon, DIM_LEN_lat/))
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNCOV RECORD' )
        
        ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat var ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_IMS_1D)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon var ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_IMS_1D)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! find indices within grid cell boundaries
        ALLOCATE(Lat_minmax_Diff(DIM_LEN_lat))
        ALLOCATE(Lon_minmax_Diff(DIM_LEN_lon))
        ! find index with min abs value difference
        Lat_minmax_Diff = abs(Lat_IMS_1D - lat_min)
        minlat_indx = MINLOC(Lat_minmax_Diff, dim = 1)
        Lat_minmax_Diff = abs(Lat_IMS_1D - lat_max)
        maxlat_indx = MINLOC(Lat_minmax_Diff, dim = 1)
        Lon_minmax_Diff = abs(Lon_IMS_1D - lon_min)
        minlon_indx = MINLOC(Lon_minmax_Diff, dim =1 )
        Lon_minmax_Diff = abs(Lon_IMS_1D - lon_max)
        maxlon_indx = MINLOC(Lon_minmax_Diff, dim =1 )
        
        DIM_LEN_lat = 1 + abs(maxlat_indx - minlat_indx)
        DIM_LEN_lon = 1 + abs(maxlon_indx - minlon_indx)
        ! print*, "New (Tile-specific) IMS array size = ", DIM_LEN_lon, " ",DIM_LEN_lat
        
        DIM_LEN = DIM_LEN_lat * DIM_LEN_lon
        ALLOCATE(SNCOV_IMS(DIM_LEN))
        ALLOCATE(Lat_IMS(DIM_LEN))
        ALLOCATE(Lon_IMS(DIM_LEN))
        ALLOCATE(Lat_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))
        ALLOCATE(Lon_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))
        ALLOCATE(SNCOV_IMS_2D(DIM_LEN_lon, DIM_LEN_lat))        
        
        ! print*," Dims of SNCOV_IMS_2D_full ", size(SNCOV_IMS_2D_full,1), size(SNCOV_IMS_2D_full,2)
        iincr = 1; jincr = 1;
        if (maxlat_indx < minlat_indx) jincr = -1
        if (maxlon_indx < minlon_indx) iincr = -1
        ! print*, "lat subarray indices = ", minlat_indx, " ", maxlat_indx
        ! print*, "lon subarray indices = ", minlon_indx, " ", maxlon_indx
        ! print*, "jincr, iincr = ", jincr, " ", iincr
        jcounter = 1
        Do jndx=minlat_indx, maxlat_indx, jincr
            !print*, "jcounter ", jcounter
            icounter = 1
            Do indx=minlon_indx, maxlon_indx, iincr
                !print*, "icounter ", icounter
                Lat_IMS_2D(icounter, jcounter) = Lat_IMS_1D(jndx)
                Lon_IMS_2D(icounter, jcounter) = Lon_IMS_1D(indx)
                SNCOV_IMS_2D(icounter, jcounter) = SNCOV_IMS_2D_full(indx, jndx)
                icounter = icounter + 1
            End do
            jcounter = jcounter + 1
        End do
        Lat_IMS = Reshape(Lat_IMS_2D, (/DIM_LEN/))
        Lon_IMS = Reshape(Lon_IMS_2D, (/DIM_LEN/))
        SNCOV_IMS = Reshape(SNCOV_IMS_2D, (/DIM_LEN/))
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        
        ERROR = NF90_CLOSE(NCID)
    
        DEALLOCATE(SNCOV_IMS_2D, SNCOV_IMS_2D_full)
        DEALLOCATE(Lat_IMS_1D, Lon_IMS_1D, Lat_IMS_2D, Lon_IMS_2D)
        DEALLOCATE(Lat_minmax_Diff, Lon_minmax_Diff)    
                  
        RETURN
        
     End SUBROUTINE Observation_Read_IMS
    
     ! Get model states at obs points
     ! Warning: This assumes all distance coordinates are valid; 
     ! do quality control of coordinates beforehand
     SUBROUTINE Observation_Operator(RLA, RLO, OROG, Lat_Obs, Lon_Obs,   &
                            LENSFC, num_Obs, max_distance,              &
                            SND_back,                                  &
                            SND_atObs, Ele_atObs, index_back_atObs) !,      &
                            !intp_mode) 
                            !  SWE_atObs  SWE_back
    
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        Real, Intent(In)        :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
        INTEGER :: LENSFC, num_Obs
        Real    :: max_distance   ! radius_of_influence
        Real, Intent(In)        :: SND_back(LENSFC)
    
        Real, Intent(Out)       :: SND_atObs(num_Obs), Ele_atObs(num_Obs)
        Integer, Intent(Out)    :: index_back_atObs(num_Obs)   ! the location of the corresponding obs
        
        Real    ::  Lon_Obs_2(num_Obs)          !RLO_2(LENSFC),         
        Real    :: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        INTEGER :: indx, jndx, zndx, min_indx
        Real    :: distArr(LENSFC), haversinArr(LENSFC)
        Real    :: d_latArr(LENSFC), d_lonArr(LENSFC)
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        ! PRINT*, "PI: ", PI_16
        ! PRINT*, "PI / 180: ", pi_div_180
    
        !Fill background values to nan (to differentiate those htat don't have value)
        SND_atObs = IEEE_VALUE(SND_atObs, IEEE_QUIET_NAN)     
        Ele_atObs = IEEE_VALUE(Ele_atObs, IEEE_QUIET_NAN)       
        index_back_atObs = -1   ! when corresponding value doesn't exit
        
        !if intp_mode == 'near'         ! [bilinear, customInterpol])
    
        ! RLO from 0 to 360 (no -ve lon)
        Lon_Obs_2 = Lon_Obs
        Where(Lon_Obs_2 < 0) Lon_Obs_2 = 360. + Lon_Obs_2
        ! Do zndx = 1, num_Obs 
        !       if (Lon_Obs_2(zndx) < 0) Lon_Obs_2(zndx) = 360. + Lon_Obs_2(zndx)
        ! end do
    
        ! at each obs point compute its distance from RLA/RLO pairs 
        ! then find the position of the minimum
    
        ! shortest distance over sphere using great circle distance     
        RLA_rad =  pi_div_180 * RLA
        RLO_rad =  pi_div_180 * RLO
        Lat_Obs_rad =  pi_div_180 * Lat_Obs
        Lon_Obs_rad =  pi_div_180 * Lon_Obs_2   
        
        ! https://en.wikipedia.org/wiki/Haversine_formula
        ! https://www.geeksforgeeks.org/program-distance-two-points-earth/
        ! Distance, d = R * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
        ! dist = 2 * R * asin { sqrt [sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2]}
        Do indx = 1, num_Obs 
            d_latArr = (Lat_Obs_rad(indx) - RLA_rad) / 2.
            d_lonArr = (Lon_Obs_rad(indx) - RLO_rad) / 2.
            haversinArr = sin(d_latArr)**2 + cos(Lat_Obs_rad(indx)) * cos(RLA_rad) * sin(d_lonArr)**2
            WHERE(haversinArr > 1) haversinArr = 1.   ! ensure numerical errors don't make h>1
            
            distArr = 2 * earth_rad * asin(sqrt(haversinArr))           
            !distArr = (Lat_Obs(indx) - RLA)**2 + (Lon_Obs_2(indx) - RLO)**2 
            min_indx = MINLOC(distArr, dim = 1)  !, MASK=ieee_is_nan(distArr))
    
            if(distArr(min_indx) < max_distance) then
                SND_atObs(indx) = SND_back(min_indx) 
                Ele_atObs(indx) = OROG(min_indx)
                index_back_atObs(indx) = min_indx
            ! else
            !     Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
            endif
        end do
        
        RETURN
        
     END SUBROUTINE Observation_Operator
 
    SUBROUTINE Observation_Operator_Parallel_vect(Myrank, NPROCS, LENSFC, LENSFC_land, &
            num_Obs, num_Eval, RLA_in, RLO_in, Lat_Obs, Lon_Obs, OROG_at_stn, max_distance, &
            SNOFCS_back_in, stn_obs, LANDMASK_in, SNOFCS_atObs, &
            index_back_atObs, index_back_atEval, index_obs_atEval, & !OROGFCS_atObs,
            Obs_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)    !, Orog_obs_atEvalPts) !SNOFCS_atEvalPts, 
! Draper, edited to make generic

! CSD -  split out the eval into a separate call. (later)
                            
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"

        INTEGER                 :: Myrank, NPROCS, LENSFC, LENSFC_land, num_Obs, num_Eval
        Real, Intent(In)        :: RLA_in(LENSFC_land), RLO_in(LENSFC_land), &
                                   SNOFCS_back_in(LENSFC_land)
        integer, intent(in)     :: LANDMASK_in(LENSFC_land)
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs), OROG_at_stn(num_Obs)  ! don't want to alter these
       
        Real                    :: max_distance   ! extent from center of grid cell to search for obs
       
        Real, Intent(InOut)     :: stn_obs(num_Obs) 
        Real, Intent(Out)       :: SNOFCS_atObs(num_Obs), Obs_atEvalPts(num_Eval) !, OROGFCS_atObs(num_Obs)
        Integer, Intent(Out)    :: index_back_atEval(num_Eval), &
                                   index_obs_atEval(num_Eval)   ! the location of evaluation points
        Integer, Intent(Out)    :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real, Intent(Out)       :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval)
        ! Real, Intent(Out)       :: Orog_obs_atEvalPts(num_Eval)
                                !    SNOFCS_atEvalPts(num_Eval), &
        
        Real                    :: RLA(LENSFC), RLO(LENSFC), SNOFCS_back(LENSFC)
        Integer                 :: LANDMASK(LENSFC)
        !Integer        :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real                    :: Lon_Obs_2(num_Obs)           ! RLO_2(LENSFC)  !         
        Real                    :: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real                    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        INTEGER                 :: indx, jndx, jzndx, zndx, min_indx
        Real                    :: distArr(LENSFC), haversinArr(LENSFC)
        Real                    :: d_latArr(LENSFC), d_lonArr(LENSFC)
        Real(16), Parameter     :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter     :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end, ixy 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), &
                              dest_Aoffset, dest_Aoffset_end
        INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR, arLen, irank
        Real               :: rand_nextVal  ! randomly select evalution points to exclude from DA
        Integer            :: rand_evalPoint(num_Eval)    ! randomly select evalution points to exclude from DA

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
        else
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

    ! background subarray for proc
        N_sA = LENSFC/NPROCS    !/ Np_til  ! sub array length per proc
        N_sA_Ext = MOD(LENSFC, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
        !global background array (sub) indices for a iproc
        if (MYRANK > 0) then  
            call MPI_SEND(RLA_in, LENSFC_land, mpiReal_size, 0, &
                            MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(RLO_in, LENSFC_land, mpiReal_size, 0, &
                            100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SNOFCS_back_in, LENSFC_land, mpiReal_size, 0, &
                            200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(LANDMASK_in, LENSFC_land, mpiInt_size, 0, &
                            300*(MYRANK+1), MPI_COMM_WORLD, IERR)
        endif
        if (MYRANK == 0) then 
            RLA(1:LENSFC_land) = RLA_in
            RLO(1:LENSFC_land) = RLO_in
            SNOFCS_back(1:1+N_sA) = SNOFCS_back_in
            LANDMASK(1:1+N_sA) = LANDMASK_in
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
                call MPI_RECV(RLA(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(RLO(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SNOFCS_back(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(LANDMASK(dest_Aoffset:dest_Aoffset_end), arLen, mpiInt_size, &
                            ixy, 300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        endif
        call MPI_BCAST(RLA, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(RLO, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(SNOFCS_back, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(LANDMASK, LENSFC, mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        print*, " proc ", myrank, " done broadcast RLO RLA SNFC_Back data"     

    ! end index of subarray for proc
        N_sA = num_Obs/NPROCS    !/ Np_til  ! sub array length per proc
        N_sA_Ext = MOD(num_Obs, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
        if(MYRANK < N_sA_Ext) then   !p_tRank== 0) then 
            mp_start = MYRANK * (N_sA + 1) + 1    ! 1
            mp_end = mp_start + N_sA        !+ 1 
            arLen = N_sA + 1
        else
            mp_start = MYRANK * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
            !mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
            mp_end = mp_start + N_sA - 1 
            arLen = N_sA 
        endif
        !Fill background values to nan (to differentiate those htat don't have value)
        SNOFCS_atObs = IEEE_VALUE(SNOFCS_atObs, IEEE_QUIET_NAN) 
        !OROGFCS_atObs = IEEE_VALUE(OROGFCS_atObs, IEEE_QUIET_NAN)       
        index_back_atObs = -1   ! when corresponding value doesn't exit 
        rand_evalPoint = -1
        index_back_atEval = -1
    
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
    
            if ( (distArr(min_indx) < max_distance) .and. &    ! if loc further than max dist
                 (LANDMASK(min_indx) == 1) ) then              ! & if nearest cell is no land, fcs value remains NaN
                SNOFCS_atObs(indx) = SNOFCS_back(min_indx) 
                !OROGFCS_atObs(indx) = OROG(min_indx)
                index_back_atObs(indx) = min_indx
            ! else 
                ! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
            endif
        end do

        if (MYRANK > 0) then  
            call MPI_SEND(SNOFCS_atObs(mp_start:mp_end), arLen, &
                        mpiReal_size, 0, 300*(MYRANK+1), MPI_COMM_WORLD, IERR) 
            call MPI_SEND(index_back_atObs(mp_start:mp_end), arLen, &
                        mpiInt_size, 0, 400*(MYRANK+1), MPI_COMM_WORLD, IERR)
        endif
        if (MYRANK == 0) then 
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
                call MPI_RECV(SNOFCS_atObs(dest_Aoffset:dest_Aoffset_end), arLen, &
                mpiReal_size, ixy, 300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset_end), arLen, &
                mpiInt_size, ixy, 400*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
            print*, " done receiving background index and data" 
        endif 
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_BCAST(SNOFCS_atObs, num_Obs, mpiReal_size, 0, MPI_COMM_WORLD, IERR)  
        print*, " proc ", myrank, " done broadcast background data"          
        call MPI_BCAST(index_back_atObs, num_Obs, mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        print*, " proc ", myrank, " done broadcast background index"        
        if (MYRANK == 0) then   !if (MYRANK == p_tN ) then      
            ! Print*, "Started selecting obs points"
            ! Select obs points to exclude from DA 
            if(num_Eval > 0) then
                Call random_number(rand_nextVal)
                rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                if(.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(1)))) then   !/   & stn_obs(rand_evalPoint(1))
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                else   ! try one more time
                    Call random_number(rand_nextVal)
                    rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                endif
            end if
            jndx = 2
            jzndx = 1
            Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
                Call random_number(rand_nextVal)
                ! Print*, rand_nextVal
                rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
                if((rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND.      &
                   (.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(jndx)))) )  then  !/stn_obs(rand_evalPoint(jndx))
                    index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
                    jndx = jndx + 1
                else
                    cycle
                endif
                jzndx = jzndx + 1
                If (jzndx >= 2*num_Eval) exit 
            Enddo 
        endif
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_BCAST(index_back_atEval, num_Eval, & 
                mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(rand_evalPoint, num_Eval, & 
                mpiInt_size, 0, MPI_COMM_WORLD, IERR)             
           
        Obs_atEvalPts = IEEE_VALUE(Obs_atEvalPts, IEEE_QUIET_NAN)
        ! SNOFCS_atEvalPts = IEEE_VALUE(SNOFCS_atEvalPts, IEEE_QUIET_NAN)
        Lat_atEvalPts = IEEE_VALUE(Lat_atEvalPts, IEEE_QUIET_NAN)
        Lon_atEvalPts = IEEE_VALUE(Lon_atEvalPts, IEEE_QUIET_NAN)
        ! Orog_obs_atEvalPts = IEEE_VALUE(Orog_obs_atEvalPts, IEEE_QUIET_NAN)
        Do  jndx = 1, num_Eval
            if (index_back_atEval(jndx) > 0) then
                Obs_atEvalPts(jndx) = stn_obs(rand_evalPoint(jndx))
                ! SNOFCS_atEvalPts(jndx) = SNOFCS_atObs(rand_evalPoint(jndx))
                Lat_atEvalPts(jndx) = Lat_Obs(rand_evalPoint(jndx)) 
                Lon_atEvalPts(jndx) = Lon_Obs(rand_evalPoint(jndx))
                ! Orog_obs_atEvalPts(jndx) = OROG_at_stn(rand_evalPoint(jndx))
                stn_obs(rand_evalPoint(jndx)) = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN) ! exclude point from DA     
            endif       
        Enddo   
        index_obs_atEval = rand_evalPoint
        
        RETURN
        
     END SUBROUTINE Observation_Operator_Parallel_vect 

    SUBROUTINE Observation_Operator_Parallel(Myrank, NPROCS, LENSFC, LENSFC_land, &
            num_Obs, num_Eval, RLA_in, RLO_in, Lat_Obs, Lon_Obs, OROG_at_stn, max_distance, &
            SNOFCS_back_in, stn_obs, LANDMASK_in, SNOFCS_atObs, &
            index_back_atObs, index_back_atEval,     &  !index_obs_atEval, & !OROGFCS_atObs,
            Obs_atEvalPts,SNOFCS_atEvalPts,Lat_atEvalPts,Lon_atEvalPts,Orog_obs_atEvalPts) ! 
! Draper, edited to make generic

! CSD -  split out the eval into a separate call. (later)
                            
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"

        INTEGER                 :: Myrank, NPROCS, LENSFC, LENSFC_land, num_Obs, num_Eval
        Real, Intent(In)        :: RLA_in(LENSFC_land), RLO_in(LENSFC_land), &
                                   SNOFCS_back_in(LENSFC_land)
        integer, intent(in)     :: LANDMASK_in(LENSFC_land)
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs), OROG_at_stn(num_Obs)  ! don't want to alter these
       
        Real                    :: max_distance   ! extent from center of grid cell to search for obs
       
        Real, Intent(InOut)     :: stn_obs(num_Obs) 
        Real, Intent(Out)       :: SNOFCS_atObs(num_Obs), Obs_atEvalPts(num_Eval) !, OROGFCS_atObs(num_Obs)
        Integer, Intent(Out)    :: index_back_atEval(num_Eval) !, &
                                !    index_obs_atEval(num_Eval)   ! the location of evaluation points
        Integer, Intent(Out)    :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real, Intent(Out)       :: Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval), &
                        Orog_obs_atEvalPts(num_Eval), SNOFCS_atEvalPts(num_Eval)
        
        Real                    :: RLA(LENSFC), RLO(LENSFC), SNOFCS_back(LENSFC)
        Integer                 :: LANDMASK(LENSFC)
        !Integer        :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real                    :: Lon_Obs_2(num_Obs)           ! RLO_2(LENSFC)  !         
        Real                    :: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real                    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        INTEGER                 :: indx, jndx, jzndx, zndx, min_indx
        Real                    :: distArr(LENSFC), haversinArr(LENSFC)
        Real                    :: d_latArr(LENSFC), d_lonArr(LENSFC)
        Real(16), Parameter     :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter     :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end, ixy 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), &
                              dest_Aoffset, dest_Aoffset_end
        INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR, arLen, irank
        Real               :: rand_nextVal  ! randomly select evalution points to exclude from DA
        Integer            :: rand_evalPoint(num_Eval)    ! randomly select evalution points to exclude from DA

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
        else
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

    ! background subarray for proc
        N_sA = LENSFC/NPROCS    !/ Np_til  ! sub array length per proc
        N_sA_Ext = MOD(LENSFC, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
        !global background array (sub) indices for a iproc
        if (MYRANK > 0) then  
            call MPI_SEND(RLA_in, LENSFC_land, mpiReal_size, 0, &
                            MYRANK+1, MPI_COMM_WORLD, IERR) 
            call MPI_SEND(RLO_in, LENSFC_land, mpiReal_size, 0, &
                            100*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(SNOFCS_back_in, LENSFC_land, mpiReal_size, 0, &
                            200*(MYRANK+1), MPI_COMM_WORLD, IERR)
            call MPI_SEND(LANDMASK_in, LENSFC_land, mpiInt_size, 0, &
                            300*(MYRANK+1), MPI_COMM_WORLD, IERR)
        endif
        if (MYRANK == 0) then 
            RLA(1:LENSFC_land) = RLA_in
            RLO(1:LENSFC_land) = RLO_in
            SNOFCS_back(1:1+N_sA) = SNOFCS_back_in
            LANDMASK(1:1+N_sA) = LANDMASK_in
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
                call MPI_RECV(RLA(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, ixy+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(RLO(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 100*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(SNOFCS_back(dest_Aoffset:dest_Aoffset_end), arLen, mpiReal_size, &
                            ixy, 200*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(LANDMASK(dest_Aoffset:dest_Aoffset_end), arLen, mpiInt_size, &
                            ixy, 300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
        endif
        call MPI_BCAST(RLA, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(RLO, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(SNOFCS_back, LENSFC, mpiReal_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(LANDMASK, LENSFC, mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        print*, " proc ", myrank, " done broadcast RLO RLA SNFC_Back data"     

    ! end index of subarray for proc
        N_sA = num_Obs/NPROCS    !/ Np_til  ! sub array length per proc
        N_sA_Ext = MOD(num_Obs, NPROCS)  !LENSFC - N_sA * Np_til ! extra grid cells
        if(MYRANK < N_sA_Ext) then   !p_tRank== 0) then 
            mp_start = MYRANK * (N_sA + 1) + 1    ! 1
            mp_end = mp_start + N_sA        !+ 1 
            arLen = N_sA + 1
        else
            mp_start = MYRANK * N_sA + N_sA_Ext + 1  !N_sA_Ext*(N_sA+1) + (MYRANK-N_sA_Ext)*N_sA     
            !mp_start = p_tRank * N_sA + N_sA_Ext + 1   ! start index of subarray for proc
            mp_end = mp_start + N_sA - 1 
            arLen = N_sA 
        endif
        !Fill background values to nan (to differentiate those htat don't have value)
        SNOFCS_atObs = IEEE_VALUE(SNOFCS_atObs, IEEE_QUIET_NAN) 
        !OROGFCS_atObs = IEEE_VALUE(OROGFCS_atObs, IEEE_QUIET_NAN)       
        index_back_atObs = -1   ! when corresponding value doesn't exit 
        rand_evalPoint = -1
        index_back_atEval = -1
    
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
    
            if ( (distArr(min_indx) < max_distance) .and. &    ! if loc further than max dist
                 (LANDMASK(min_indx) == 1) ) then              ! & if nearest cell is no land, fcs value remains NaN
                SNOFCS_atObs(indx) = SNOFCS_back(min_indx) 
                !OROGFCS_atObs(indx) = OROG(min_indx)
                index_back_atObs(indx) = min_indx
            ! else 
                ! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
            endif
        end do

        if (MYRANK > 0) then  
            call MPI_SEND(SNOFCS_atObs(mp_start:mp_end), arLen, &
                        mpiReal_size, 0, 300*(MYRANK+1), MPI_COMM_WORLD, IERR) 
            call MPI_SEND(index_back_atObs(mp_start:mp_end), arLen, &
                        mpiInt_size, 0, 400*(MYRANK+1), MPI_COMM_WORLD, IERR)
        endif
        if (MYRANK == 0) then 
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
                call MPI_RECV(SNOFCS_atObs(dest_Aoffset:dest_Aoffset_end), arLen, &
                mpiReal_size, ixy, 300*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset_end), arLen, &
                mpiInt_size, ixy, 400*(ixy+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            End do
            print*, " done receiving background index and data" 
        endif 
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_BCAST(SNOFCS_atObs, num_Obs, mpiReal_size, 0, MPI_COMM_WORLD, IERR)  
        print*, " proc ", myrank, " done broadcast background data"          
        call MPI_BCAST(index_back_atObs, num_Obs, mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        print*, " proc ", myrank, " done broadcast background index"        
        if (MYRANK == 0) then   !if (MYRANK == p_tN ) then      
            ! Print*, "Started selecting obs points"
            ! Select obs points to exclude from DA 
            if(num_Eval > 0) then
                Call random_number(rand_nextVal)
                rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                if(.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(1)))) then   !/   & stn_obs(rand_evalPoint(1))
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                else   ! try one more time
                    Call random_number(rand_nextVal)
                    rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                endif
            end if
            jndx = 2
            jzndx = 1
            Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
                Call random_number(rand_nextVal)
                ! Print*, rand_nextVal
                rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
                if((rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND.      &
                   (.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(jndx)))) )  then  !/stn_obs(rand_evalPoint(jndx))
                    index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
                    jndx = jndx + 1
                else
                    cycle
                endif
                jzndx = jzndx + 1
                If (jzndx >= 2*num_Eval) exit 
            Enddo 
        endif
        ! MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        call MPI_BCAST(index_back_atEval, num_Eval, & 
                mpiInt_size, 0, MPI_COMM_WORLD, IERR) 
        call MPI_BCAST(rand_evalPoint, num_Eval, & 
                mpiInt_size, 0, MPI_COMM_WORLD, IERR)             
           
        Obs_atEvalPts = IEEE_VALUE(Obs_atEvalPts, IEEE_QUIET_NAN)
        SNOFCS_atEvalPts = IEEE_VALUE(SNOFCS_atEvalPts, IEEE_QUIET_NAN)
        Lat_atEvalPts = IEEE_VALUE(Lat_atEvalPts, IEEE_QUIET_NAN)
        Lon_atEvalPts = IEEE_VALUE(Lon_atEvalPts, IEEE_QUIET_NAN)
        Orog_obs_atEvalPts = IEEE_VALUE(Orog_obs_atEvalPts, IEEE_QUIET_NAN)
        Do  jndx = 1, num_Eval
            if (index_back_atEval(jndx) > 0) then
                Obs_atEvalPts(jndx) = stn_obs(rand_evalPoint(jndx))
                SNOFCS_atEvalPts(jndx) = SNOFCS_atObs(rand_evalPoint(jndx))
                Lat_atEvalPts(jndx) = Lat_Obs(rand_evalPoint(jndx)) 
                Lon_atEvalPts(jndx) = Lon_Obs(rand_evalPoint(jndx))
                Orog_obs_atEvalPts(jndx) = OROG_at_stn(rand_evalPoint(jndx))
                stn_obs(rand_evalPoint(jndx)) = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN) ! exclude point from DA     
            endif       
        Enddo   
        ! index_obs_atEval = rand_evalPoint
        
        RETURN
        
     END SUBROUTINE Observation_Operator_Parallel

     SUBROUTINE Observation_Operator_tiles_Parallel(Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, & 
                            RLA, RLO, OROG, Lat_Obs, Lon_Obs,               &
                            LENSFC, num_Obs, num_Eval, max_distance, SNOFCS_back, stn_obs, LANDMASK,  &
                            SNOFCS_atObs, OROGFCS_atObs, index_back_atObs, index_back_atEval,       &
                            Obs_atEvalPts, SNOFCS_atEvalPts, Lat_atEvalPts, Lon_atEvalPts)
! Draper, edited to make generic

! CSD -  split out the eval into a separate call. (later)
                            
        IMPLICIT NONE
        !
        !USE intrinsic::ieee_arithmetic
        include "mpif.h"
    
        Real, Intent(In)        :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC)
        integer, intent(in)     ::  LANDMASK(LENSFC)
        Real, Intent(In)        :: Lat_Obs(num_Obs), Lon_Obs(num_Obs)  ! don't want to alter these
        INTEGER             :: Myrank, MAX_TASKS, p_tN, p_tRank, Np_til, LENSFC, num_Obs, num_Eval
        Real                :: max_distance   ! extent from center of grid cell to search for obs
        Real, Intent(In)        :: SNOFCS_back(LENSFC)
        Real, Intent(InOut)     :: stn_obs(num_Obs) 
        Real, Intent(Out)       :: SNOFCS_atObs(num_Obs), OROGFCS_atObs(num_Obs), Obs_atEvalPts(num_Eval)
        Integer, Intent(Out)    :: index_back_atEval(num_Eval)   ! the location of evaluation points
        Integer, Intent(Out)    :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real, Intent(Out)           :: SNOFCS_atEvalPts(num_Eval), Lat_atEvalPts(num_Eval), Lon_atEvalPts(num_Eval)
        
        !Integer        :: index_back_atObs(num_Obs)   ! the location of background corresponding obs
        Real    :: Lon_Obs_2(num_Obs)           ! RLO_2(LENSFC)  !         
        Real    :: RLA_rad(LENSFC), RLO_rad(LENSFC)
        Real    :: Lat_Obs_rad(num_Obs), Lon_Obs_rad(num_Obs)   
        INTEGER :: indx, jndx, jzndx, zndx, min_indx
        Real    :: distArr(LENSFC), haversinArr(LENSFC)
        Real    :: d_latArr(LENSFC), d_lonArr(LENSFC)
        Real(16), Parameter :: PI_16 = 4 * atan (1.0_16)        
        Real(16), Parameter :: pi_div_180 = PI_16/180.0
        Real, Parameter         :: earth_rad = 6371.
        
        ! for mpi par
        INTEGER            :: N_sA, N_sA_Ext, mp_start, mp_end 
        INTEGER            :: send_proc, rec_proc, rec_stat(MPI_STATUS_SIZE), dest_Aoffset, pindex
        INTEGER            :: mpiReal_size, rsize, mpiInt_size, isize, IERR
        Real               :: rand_nextVal  ! randomly select evalution points to exclude from DA
        Integer            :: rand_evalPoint(num_Eval)    ! randomly select evalution points to exclude from DA
    
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
        OROGFCS_atObs = IEEE_VALUE(OROGFCS_atObs, IEEE_QUIET_NAN)       
        index_back_atObs = -1   ! when corresponding value doesn't exit 
        rand_evalPoint = -1
        index_back_atEval = -1
    
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
    
            if ( (distArr(min_indx) < max_distance) .and.   &  ! if too far away, don't use
	             (LANDMASK(min_indx) == 1) ) then           ! if nearest cell is no land, fcs value remains NaN
                SNOFCS_atObs(indx) = SNOFCS_back(min_indx) 
                OROGFCS_atObs(indx) = OROG(min_indx)
                index_back_atObs(indx) = min_indx
            ! else 
                ! Print*, " Warning! distance greater than ",max_distance," km ", distArr(min_indx)
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
        else
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
            call MPI_SEND(OROGFCS_atObs(mp_start:mp_end), N_sA, mpiReal_size, p_tN,   &
                          MYRANK*100, MPI_COMM_WORLD, IERR)
            call MPI_SEND(index_back_atObs(mp_start:mp_end), N_sA, mpiInt_size, p_tN,   &
                          MYRANK*1000, MPI_COMM_WORLD, IERR)
        else !if (MYRANK == p_tN ) then  
            Do pindex =  1, (Np_til - 1)   ! sender proc index within tile group
                dest_Aoffset = pindex * N_sA + N_sA_Ext + 1   ! dest array offset
                send_proc = MYRANK +  pindex * MAX_TASKS
                call MPI_RECV(SNOFCS_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,  &
                          send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(OROGFCS_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiReal_size, send_proc,   &
                          send_proc*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                call MPI_RECV(index_back_atObs(dest_Aoffset:dest_Aoffset+N_sA-1), N_sA, mpiInt_size, send_proc, &
                          send_proc*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            enddo
        endif
    !ToDO: better way to do this?
        ! now share the whole array
        if (MYRANK < MAX_TASKS ) then   !if (MYRANK == p_tN ) then      
            ! Print*, "Started selecting obs points"
            ! Select obs points to exclude from DA 
            if(num_Eval > 0) then
                Call random_number(rand_nextVal)
                rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                if(.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(1)))) then   !/   & stn_obs(rand_evalPoint(1))
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                else   ! try one more time
                    Call random_number(rand_nextVal)
                    rand_evalPoint(1) = floor(rand_nextVal * num_Obs) + 1
                    index_back_atEval(1) = index_back_atObs(rand_evalPoint(1))
                endif
            end if
            jndx = 2
            jzndx = 1
            Do  While (jndx <= num_Eval)        !jndx = 2, num_Eval
                Call random_number(rand_nextVal)
                ! Print*, rand_nextVal
                rand_evalPoint(jndx) = floor(rand_nextVal * num_Obs) + 1
                if((rand_evalPoint(jndx) /= rand_evalPoint(jndx-1)) .AND.      &
                   (.NOT. IEEE_IS_NAN(stn_obs(rand_evalPoint(jndx)))) )  then  !/stn_obs(rand_evalPoint(jndx))
                    index_back_atEval(jndx) = index_back_atObs(rand_evalPoint(jndx))
                    jndx = jndx + 1
                else
                    cycle
                endif
                jzndx = jzndx + 1
                If (jzndx >= 2*num_Eval) exit 
            Enddo 
            ! Print*, "Finished selecting obs points"   
            Do pindex =  1, (Np_til - 1)   ! receiving proc index within tile group
                rec_proc = MYRANK +  pindex * MAX_TASKS
                call MPI_SEND(SNOFCS_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK, MPI_COMM_WORLD, IERR) 
                call MPI_SEND(OROGFCS_atObs, num_Obs, mpiReal_size, rec_proc, MYRANK*100, MPI_COMM_WORLD, IERR)
                call MPI_SEND(index_back_atObs, num_Obs, mpiInt_size, rec_proc, MYRANK*1000, MPI_COMM_WORLD, IERR)
                call MPI_SEND(index_back_atEval, num_Eval, mpiInt_size, rec_proc, MYRANK*10000, MPI_COMM_WORLD, IERR)
                call MPI_SEND(rand_evalPoint, num_Eval, mpiInt_size, rec_proc, MYRANK*100000, MPI_COMM_WORLD, IERR)
            enddo
        else 
            call MPI_RECV(SNOFCS_atObs, num_Obs, mpiReal_size, p_tN, p_tN, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(OROGFCS_atObs, num_Obs, mpiReal_size, p_tN, p_tN*100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(index_back_atObs, num_Obs, mpiInt_size, p_tN, p_tN*1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(index_back_atEval, num_Eval, mpiInt_size, p_tN, p_tN*10000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
            call MPI_RECV(rand_evalPoint, num_Eval, mpiInt_size, p_tN, p_tN*100000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
        endif
        Obs_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
        SNOFCS_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
        Lat_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
        Lon_atEvalPts = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN)
        Do  jndx = 1, num_Eval
            if (index_back_atEval(jndx) > 0) then
                Obs_atEvalPts(jndx) = stn_obs(rand_evalPoint(jndx))
                SNOFCS_atEvalPts(jndx) = SNOFCS_atObs(rand_evalPoint(jndx))
                Lat_atEvalPts(jndx) = Lat_Obs(rand_evalPoint(jndx)) 
                Lon_atEvalPts(jndx) = Lon_Obs(rand_evalPoint(jndx))
                stn_obs(rand_evalPoint(jndx)) = IEEE_VALUE(rand_nextVal, IEEE_QUIET_NAN) ! exclude point from DA     
            endif       
        Enddo   
        
        RETURN
        
     END SUBROUTINE Observation_Operator_tiles_Parallel

    SUBROUTINE Observation_Read_GHCND_IODA(ghcnd_inp_file, &
                    STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                    NDIM, SND_GHCND, Ele_GHCND, Lat_GHCND, Lon_GHCND, &
                    NEXC, Index_Obs_Excluded)
        
        IMPLICIT NONE    
        include 'mpif.h'

        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, STN_DIM_NAME
        CHARACTER(LEN=*), Intent(In)      :: STN_VAR_NAME, STN_ELE_NAME
        INTEGER, Intent(Out)              :: NDIM, NEXC
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:), Ele_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)    :: Lat_GHCND(:), Lon_GHCND(:)
        INTEGER, ALLOCATABLE, Intent(Out) :: Index_Obs_Excluded(:)
    
        INTEGER                :: ERROR, NCID, grp_ncid, ID_DIM, ID_VAR
        LOGICAL                :: file_exists

        INQUIRE(FILE=trim(ghcnd_inp_file), EXIST=file_exists)
        if (.not. file_exists) then
                print *, 'Observation_Read_GHCND_IODA error,file does not exist', &
                        trim(ghcnd_inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
      
        ERROR=NF90_INQ_DIMID(NCID, TRIM(STN_DIM_NAME), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(SND_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        ALLOCATE(Ele_GHCND(NDIM))

        ERROR=nf90_inq_ncid(ncid, "MetaData", grp_ncid)
        CALL NETCDF_ERR(ERROR, 'ERROR MetaData ID' )
         ! need to read corresponding elevation values 
        ERROR=NF90_INQ_VARID(grp_ncid, TRIM(STN_ELE_NAME), ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, Ele_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )

        ERROR=NF90_INQ_VARID(grp_ncid, 'latitude', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING latitude ID' )
        ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, Lat_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
        ERROR=NF90_INQ_VARID(grp_ncid, 'longitude', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, Lon_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
        
        ERROR=nf90_inq_ncid(ncid, "ObsValue", grp_ncid)
        CALL NETCDF_ERR(ERROR, 'ERROR ObsValue ID' )
        ERROR=NF90_INQ_VARID(grp_ncid, TRIM(STN_VAR_NAME), ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
        ERROR=NF90_GET_VAR(grp_ncid, ID_VAR, SND_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )

        ! read excluded indices
        ERROR=NF90_INQ_DIMID(NCID, "obs_excluded", ID_DIM)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
        If (ERROR == NF90_NOERR) then                
            ERROR=NF90_INQUIRE_DIMENSION(NCID, ID_DIM, LEN=NEXC)
            CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )        
            ALLOCATE(Index_Obs_Excluded(NEXC))            
            ERROR=NF90_INQ_VARID(ncid, "Index_Obs_Excluded", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'ERROR Index_Obs_Excluded ID' )
            ERROR=NF90_GET_VAR(ncid, ID_VAR, Index_Obs_Excluded)
            CALL NETCDF_ERR(ERROR, 'ERROR READING Index_Obs_Excluded RECORD' )
        Else
            ! print*, " obs_excluded not found"
            NEXC = 1
            ALLOCATE(Index_Obs_Excluded(NEXC))
            Index_Obs_Excluded(NEXC) = -999
        Endif
    
        ERROR = NF90_CLOSE(NCID)
        CALL NETCDF_ERR(ERROR, 'ERROR closing obs file' )
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND_IODA

    SUBROUTINE Observation_Read_GHCND_All_excNaN(ghcnd_inp_file, &
                    STN_DIM_NAME, STN_VAR_NAME, STN_ELE_NAME, &
                    NDIM,                       &
                    SND_GHCND, Ele_GHCND,  Lat_GHCND,  Lon_GHCND)
                    ! ,          &
                    ! MYRANK)
                    !Ele_GHCND,         & lat_min, lat_max, lon_min, lon_max, &
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        ! INTEGER, Intent(In)               :: p_tN
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, STN_DIM_NAME
        CHARACTER(LEN=*), Intent(In)      :: STN_VAR_NAME, STN_ELE_NAME
        ! REAL, Intent(In)       :: lat_min, lat_max, lon_min, lon_max 
        INTEGER, Intent(Out)   :: NDIM
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:), Ele_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)    :: Lat_GHCND(:), Lon_GHCND(:)
    
        INTEGER                :: ERROR, NCID, ID_DIM, ID_VAR, NDIM_In !MYRANK, 
        REAL, ALLOCATABLE      :: SND_GHCND_In(:), Lat_GHCND_In(:), Lon_GHCND_In(:), &
                                  Ele_GHCND_In(:) 
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
      
        ERROR=NF90_INQ_DIMID(NCID, TRIM(STN_DIM_NAME), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(SND_GHCND_In(NDIM_In))
        ALLOCATE(Lat_GHCND_In(NDIM_In))
        ALLOCATE(Lon_GHCND_In(NDIM_In))
        ALLOCATE(Ele_GHCND_In(NDIM_In))
    
        ERROR=NF90_INQ_VARID(NCID, TRIM(STN_VAR_NAME), ID_VAR)
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
        ERROR=NF90_INQ_VARID(NCID, TRIM(STN_ELE_NAME), ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        ALLOCATE(index_Array(NDIM_In))
        NDIM = 0
        jcounter = 1
    ! exclude all data if not within lat/lon range, or negative snow depth
        Do jndx = 1, NDIM_In
            ! Print*, "jndx = ", jndx
            ! If((Lat_GHCND_In(jndx) >= lat_min) .and. (Lat_GHCND_In(jndx) <= lat_max) .and. &
            !     (Lon_GHCND_In(jndx) >= lon_min) .and. (Lon_GHCND_In(jndx) <= lon_max) .and. &
            If (SND_GHCND_In(jndx) >= 0 .and. Ele_GHCND_In(jndx) >= 0 ) then  !(.NOT. IEEE_IS_NAN(SND_GHCND_In(jndx)))) then                        
                NDIM = NDIM + 1
                index_Array(jcounter) = jndx
                jcounter = jcounter + 1
            Endif
        End do
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
        DEALLOCATE(SND_GHCND_In)
        DEALLOCATE(Lat_GHCND_In)
        DEALLOCATE(Lon_GHCND_In)
        DEALLOCATE(Ele_GHCND_In)
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND_All_excNaN
     
    SUBROUTINE Observation_Read_GHCND_Tile_excNaN(p_tN, ghcnd_inp_file, dim_name,       &
                    lat_min, lat_max, lon_min, lon_max, &
                    NDIM,                       &
                    SND_GHCND,         &
                    Lat_GHCND,      &
                    Lon_GHCND,          &
                    !Ele_GHCND,         &
                    MYRANK)
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        INTEGER, Intent(In)               :: p_tN
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
        REAL, Intent(In)       :: lat_min, lat_max, lon_min, lon_max 
        INTEGER, Intent(Out)   :: NDIM
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)    :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        INTEGER                :: MYRANK, ERROR, NCID, ID_DIM, ID_VAR, NDIM_In
        REAL, ALLOCATABLE      :: SND_GHCND_In(:), Lat_GHCND_In(:), Lon_GHCND_In(:) !, Ele_GHCND(:)
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
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND_In)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        ALLOCATE(index_Array(NDIM_In))
        NDIM = 0
        jcounter = 1
    ! exclude all data if not within lat/lon range, or negative snow depth
        Do jndx = 1, NDIM_In
            ! Print*, "jndx = ", jndx
            If((Lat_GHCND_In(jndx) >= lat_min) .and. (Lat_GHCND_In(jndx) <= lat_max) .and. &
                (Lon_GHCND_In(jndx) >= lon_min) .and. (Lon_GHCND_In(jndx) <= lon_max) .and. &
                (SND_GHCND_In(jndx) >= 0 )) then  !(.NOT. IEEE_IS_NAN(SND_GHCND_In(jndx)))) then                        
                    NDIM = NDIM + 1
                    index_Array(jcounter) = jndx
                    jcounter = jcounter + 1
            Endif
        End do
        ALLOCATE(SND_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
        If(NDIM > 0) then
            Do jndx = 1, NDIM
                SND_GHCND(jndx) = SND_GHCND_In(index_Array(jndx))
                Lat_GHCND(jndx) = Lat_GHCND_In(index_Array(jndx))
                Lon_GHCND(jndx) = Lon_GHCND_In(index_Array(jndx))
            End do
        Endif
        DEALLOCATE(index_Array)
        DEALLOCATE(SND_GHCND_In)
        DEALLOCATE(Lat_GHCND_In)
        DEALLOCATE(Lon_GHCND_In)
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND_Tile_excNaN

     SUBROUTINE Observation_Read_GHCND_Tile(ghcnd_inp_file, dim_name,                   &
                    lat_min, lat_max, lon_min, lon_max, &
                    NDIM,                       &
                    SND_GHCND,         &
                    Lat_GHCND,      &
                    Lon_GHCND,          &
                    !Ele_GHCND,         &
                    MYRANK)
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
        REAL, Intent(In)       :: lat_min, lat_max, lon_min, lon_max 
        INTEGER, Intent(Out)   :: NDIM
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)    :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        INTEGER                :: MYRANK, ERROR, NCID, ID_DIM, ID_VAR, NDIM_In
        REAL, ALLOCATABLE      :: SND_GHCND_In(:), Lat_GHCND_In(:), Lon_GHCND_In(:) !, Ele_GHCND(:)
        INTEGER, ALLOCATABLE   :: index_Array(:)
        INTEGER                :: jndx, jcounter
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM_In)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(SND_GHCND_In(NDIM_In))
        ALLOCATE(Lat_GHCND_In(NDIM_In))
        ALLOCATE(Lon_GHCND_In(NDIM_In))
        !ALLOCATE(Ele_GHCND(NDIM_In))
    
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
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND_In)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
        ALLOCATE(index_Array(NDIM_In))
        NDIM = 0
        jcounter = 1
        Do jndx = 1, NDIM_In
            ! Print*, "jndx = ", jndx
            If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max) .and. &
                (Lon_GHCND_In(jndx) > lon_min) .and. (Lon_GHCND_In(jndx) < lon_max)) then
                    NDIM = NDIM + 1
                    index_Array(jcounter) = jndx
                    jcounter = jcounter + 1
            Endif
        End do
        ALLOCATE(SND_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
        If(NDIM > 0) then
            Do jndx = 1, NDIM
                SND_GHCND(jndx) = SND_GHCND_In(index_Array(jndx))
                Lat_GHCND(jndx) = Lat_GHCND_In(index_Array(jndx))
                Lon_GHCND(jndx) = Lon_GHCND_In(index_Array(jndx))
            End do
        Endif
        DEALLOCATE(index_Array)
        ! jcounter = 1
        ! If(NDIM > 0) then
        !       Do jndx = 1, NDIM_In
        !               If((Lat_GHCND_In(jndx) > lat_min) .and. (Lat_GHCND_In(jndx) < lat_max)) then
        !                       If((Lon_GHCND_In(jndx) > lon_min) .and. (Lon_GHCND_In(jndx) < lon_max)) then
        !                               SND_GHCND(jcounter) = SND_GHCND_In(jndx)
        !                               Lat_GHCND(jcounter) = Lat_GHCND_In(jndx)
        !                               Lon_GHCND(jcounter) = Lon_GHCND_In(jndx)
        !                               jcounter = jcounter + 1
        !                       Endif
        !               Endif
        !       End do
        ! Endif 
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND_Tile
    
     SUBROUTINE Observation_Read_GHCND(ghcnd_inp_file, dim_name,                        &
                    NDIM,                       &
                    SND_GHCND,         &
                    Lat_GHCND,      &
                    Lon_GHCND,          &
                    !Ele_GHCND,         &
                    MYRANK)
        
        IMPLICIT NONE
    
        include 'mpif.h'
        !Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: ghcnd_inp_file, dim_name
        INTEGER                :: ERROR, NCID
        INTEGER                :: MYRANK
        INTEGER                :: ID_DIM, ID_VAR
        INTEGER, Intent(Out)   :: NDIM
    
        REAL, ALLOCATABLE, Intent(Out)    :: SND_GHCND(:)
        REAL, ALLOCATABLE, Intent(Out)     :: Lat_GHCND(:), Lon_GHCND(:) !, Ele_GHCND(:)
    
        ERROR=NF90_OPEN(TRIM(ghcnd_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(ghcnd_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
        ALLOCATE(SND_GHCND(NDIM))
        ALLOCATE(Lat_GHCND(NDIM))
        ALLOCATE(Lon_GHCND(NDIM))
        !ALLOCATE(Ele_GHCND(NDIM))
    
        ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SND_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_GHCND)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_GHCND)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_GHCND

     SUBROUTINE Observation_Read_SNOTEL(snotel_inp_file,  &
                    dim_name,                   &
                    NDIM, &
                    SWE_SNOTEL,      &
                    SND_SNOTEL,                &
                    Lat_SNOTEL,      &
                    Lon_SNOTEL,         &
                    !  Ele_SNOTEL               &
                    MYRANK)
        IMPLICIT NONE
    
        include 'mpif.h'
    
        !Open netCDF for a snotel and read the SWE, SnowDepth,..., Lat, Lon, at a given datetime
        !ToDO: Can you use variable length char array ?
        CHARACTER(LEN=*), Intent(In)      :: snotel_inp_file, dim_name
    
        INTEGER                :: ERROR, NCID
        INTEGER                :: MYRANK
        INTEGER                :: ID_DIM, ID_VAR 
        INTEGER, Intent(out)   :: NDIM
    
        REAL, ALLOCATABLE, Intent(Out)    :: SWE_SNOTEL(:), SND_SNOTEL(:)
        REAL, ALLOCATABLE, Intent(Out)     :: Lat_SNOTEL(:), Lon_SNOTEL(:) !, Ele_SNOTEL(:)
        REAL, ALLOCATABLE          :: SWE_Ratio(:) !, Ele_SNOTEL(:)
        Integer                :: num_invalid
    
        ERROR=NF90_OPEN(TRIM(snotel_inp_file),NF90_NOWRITE,NCID)
        CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(snotel_inp_file) )
    
        ERROR=NF90_INQ_DIMID(NCID, TRIM(dim_name), ID_DIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Dimension' )
    
        ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=NDIM)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Size of Dimension' )
    
    
        ALLOCATE(SWE_SNOTEL(NDIM))
        ALLOCATE(SND_SNOTEL(NDIM))
        ALLOCATE(Lat_SNOTEL(NDIM))
        ALLOCATE(Lon_SNOTEL(NDIM))
        ALLOCATE(SWE_Ratio(NDIM))       
        !ALLOCATE(Ele_SNOTEL(NDIM))
    
        ERROR=NF90_INQ_VARID(NCID, 'SWE', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SWE ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SWE_SNOTEL)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SWE RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'SNWD', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, SND_SNOTEL)
        CALL NETCDF_ERR(ERROR, 'ERROR READING SNWD RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lat', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lat_SNOTEL)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lat RECORD' )
    
        ERROR=NF90_INQ_VARID(NCID, 'lon', ID_VAR)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon ID' )
        ERROR=NF90_GET_VAR(NCID, ID_VAR, Lon_SNOTEL)
        CALL NETCDF_ERR(ERROR, 'ERROR READING Lon RECORD' )
    
        SWE_Ratio = SWE_SNOTEL / SWE_SNOTEL
        num_invalid = COUNT (IEEE_IS_NAN(SWE_Ratio))
        call debug_print("number of invalid values ", float(num_invalid))
    
        ! need to read corresponding elevation values 
        ! ERROR=NF90_INQ_VARID(NCID, 'Elevation', ID_VAR)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation ID' )
        ! ERROR=NF90_GET_VAR(NCID, ID_VAR, Ele_SNOTEL)
        ! CALL NETCDF_ERR(ERROR, 'ERROR READING Elevation RECORD' )
    
        ERROR = NF90_CLOSE(NCID)
    
        DEALLOCATE(SWE_Ratio)
                  
        RETURN
        
     End SUBROUTINE Observation_Read_SNOTEL

    ! want to run point model at obs and compute correlations based on point
     subroutine Read_state_atObs(point_state_prefix, var_name, &
        y_str, m_str, d_str, h_str, ens_size, vector_length, State_obsPt)
        
        ! use time_utilities
    
        IMPLICIT NONE

        include 'mpif.h'
             
        integer           :: ens_size, vector_length !, LENSFC
    
        CHARACTER(LEN=*)  :: point_state_prefix, var_name
        CHARACTER(LEN=*)  :: y_str, m_str, d_str, h_str
        Real              :: State_obsPt(ens_size)
    
        
        integer           :: ncid, dimid, varid, status, ie
        
        CHARACTER(LEN=500)  :: state_inp_file
        LOGICAL             :: file_exists
        
        ! vector_length =  endIndx - startIndx + 1
        
        state_inp_file=TRIM(point_state_prefix)//    &    !"/noda_ensSWE."//
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//".nc" 
    
        INQUIRE(FILE=trim(state_inp_file), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(state_inp_file) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif        
        
        status = nf90_open(trim(state_inp_file), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'OPENING FILE: '//TRIM(state_inp_file) )
    
        status = nf90_inq_varid(ncid, trim(var_name), varid)
        CALL NETCDF_ERR(status, 'getting varid:'//trim(var_name) )
        status = nf90_get_var(ncid, varid , State_obsPt, &
            start = (/1,1, 1, 1/), count = (/ens_size, 1, 1, 1/)) !vector_length/))
        CALL NETCDF_ERR(status, 'reading:'//trim(var_name) )
    
        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'Closing FILE: '//TRIM(state_inp_file) )
    
       
    end subroutine Read_state_atObs

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

! codes to resample high-res (IMS, elevation..) data to fv3 (polygon) grid
! from George Gayno et al.
! https://raw.githubusercontent.com/NOAA-EMC/UFS_UTILS/develop/sorc/orog_mask_tools.fd/orog.fd/mtnlm7_oclsm.f
          
!> Convert from latitude and longitude to x,y,z coordinates.
!!
!! @param[in] siz Number of points to convert.
!! @param[in] lon Longitude of points to convert.
!! @param[in] lat Latitude of points to convert.
!! @param[out] x 'x' coordinate of the converted points.
!! @param[out] y 'y' coordinate of the converted points.
!! @param[out] z 'z' coordinate of the converted points.
!! @author GFDL programmer
      subroutine latlon2xyz(siz,lon, lat, x, y, z)
      implicit none
      integer, intent(in) :: siz
      real, intent(in) :: lon(siz), lat(siz)
      real, intent(out) :: x(siz), y(siz), z(siz)
      
      integer n

      do n = 1, siz
        x(n) = cos(lat(n))*cos(lon(n))
        y(n) = cos(lat(n))*sin(lon(n))
        z(n) = sin(lat(n))
      enddo
      end subroutine latlon2xyz

!> Compute spherical angle.
!!
!! @param[in] v1 Vector 1.
!! @param[in] v2 Vector 2.
!! @param[in] v3 Vector 3.
!! @return spherical_angle Spherical Angle.
!! @author GFDL programmer
    FUNCTION spherical_angle(v1, v2, v3)
        implicit none
        real, parameter :: EPSLN30 = 1.e-30
        real, parameter :: PI=3.1415926535897931
        real v1(3), v2(3), v3(3)
        real  spherical_angle
 
        real px, py, pz, qx, qy, qz, ddd;
  
        ! vector product between v1 and v2 
        px = v1(2)*v2(3) - v1(3)*v2(2)
        py = v1(3)*v2(1) - v1(1)*v2(3)
        pz = v1(1)*v2(2) - v1(2)*v2(1)
        ! vector product between v1 and v3 
        qx = v1(2)*v3(3) - v1(3)*v3(2);
        qy = v1(3)*v3(1) - v1(1)*v3(3);
        qz = v1(1)*v3(2) - v1(2)*v3(1);

        ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
        if ( ddd <= 0.0 ) then
          spherical_angle = 0. 
        else 
          ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
          if( abs(ddd-1) < EPSLN30 ) ddd = 1;
          if( abs(ddd+1) < EPSLN30 ) ddd = -1;
          if ( ddd>1. .or. ddd<-1. ) then
            !FIX to correctly handle co-linear points (angle near pi or 0) */
            if (ddd < 0.) then
              spherical_angle = PI
            else
              spherical_angle = 0.
            endif
          else
            spherical_angle = acos( ddd )
          endif
        endif  

        return
    END FUNCTION spherical_angle
      
!> Check if a point is inside a polygon.
!!
!! @param[in] lon1 Longitude of the point to check.
!! @param[in] lat1 Latitude of the point to check.
!! @param[in] npts Number of polygon vertices.
!! @param[in] lon2 Longitude of the polygon vertices.
!! @param[in] lat2 Latitude of the polygon vertices.
!! @return inside_a_polygon When true, point is within
!! the polygon.
!! @author GFDL programmer
      FUNCTION inside_a_polygon(lon1, lat1, npts, lon2, lat2)
        implicit none
        real, parameter :: EPSLN10 = 1.e-10
        real, parameter :: EPSLN8 = 1.e-8
        real, parameter :: PI=3.1415926535897931
        real, parameter :: RANGE_CHECK_CRITERIA=0.05
        real :: anglesum, angle !, spherical_angle
        integer i, ip1
        real lon1, lat1
        integer npts
        real lon2(npts), lat2(npts)
        real x2(npts), y2(npts), z2(npts)
        real lon1_1d(1), lat1_1d(1)
        real x1(1), y1(1), z1(1)
        real pnt0(3),pnt1(3),pnt2(3)
        logical inside_a_polygon
        real max_x2,min_x2,max_y2,min_y2,max_z2,min_z2
        !first convert to cartesian grid */
        call latlon2xyz(npts,lon2, lat2, x2, y2, z2);
        lon1_1d(1) = lon1
        lat1_1d(1) = lat1
        call latlon2xyz(1,lon1_1d, lat1_1d, x1, y1, z1);
        inside_a_polygon = .false.
        max_x2 = maxval(x2)
        if( x1(1) > max_x2+RANGE_CHECK_CRITERIA ) return
        min_x2 = minval(x2)
        if( x1(1)+RANGE_CHECK_CRITERIA < min_x2 ) return
        max_y2 = maxval(y2)
        if( y1(1) > max_y2+RANGE_CHECK_CRITERIA ) return
        min_y2 = minval(y2)
        if( y1(1)+RANGE_CHECK_CRITERIA < min_y2 ) return
        max_z2 = maxval(z2)
        if( z1(1) > max_z2+RANGE_CHECK_CRITERIA ) return
        min_z2 = minval(z2)
        if( z1(1)+RANGE_CHECK_CRITERIA < min_z2 ) return

        pnt0(1) = x1(1)
        pnt0(2) = y1(1)
        pnt0(3) = z1(1)
        
        anglesum = 0;
        do i = 1, npts
           if(abs(x1(1)-x2(i)) < EPSLN10 .and. &
               abs(y1(1)-y2(i)) < EPSLN10 .and. &
              abs(z1(1)-z2(i)) < EPSLN10 ) then ! same as the corner point
              inside_a_polygon = .true.
              return
           endif
           ip1 = i+1
           if(ip1>npts) ip1 = 1
           pnt1(1) = x2(i)
           pnt1(2) = y2(i)
           pnt1(3) = z2(i)
           pnt2(1) = x2(ip1)
           pnt2(2) = y2(ip1)
           pnt2(3) = z2(ip1)

           angle = spherical_angle(pnt0, pnt2, pnt1);
!           anglesum = anglesum + spherical_angle(pnt0, pnt2, pnt1);
           anglesum = anglesum + angle
        enddo

        if(abs(anglesum-2*PI) < EPSLN8) then
           inside_a_polygon = .true.
        else
           inside_a_polygon = .false.
        endif

        return
        
      end FUNCTION inside_a_polygon

END MODULE M_DA
     
    
