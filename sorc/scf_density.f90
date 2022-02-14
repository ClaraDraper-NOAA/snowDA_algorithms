module M_SCF_Den
    
    ! modfified from IMSaggregate_mod.f90
    ! the calculations here exclude file read 
    ! (to use in the DA where arrays of of model inputs are divided in parallel to different 
    ! mpi processes)
    
    Use, Intrinsic :: IEEE_ARITHMETIC

    real, parameter    ::  nodata_real = -999. 
    integer, parameter ::  nodata_int = -999
    real, parameter    ::  nodata_tol = 0.1

    ! IMS noah-MP snow depth retrieval parameters
    real, parameter :: qc_limit = 0.90     ! QC limit for IMS obs (data will be removed if both model and IMS 
                                        ! are above this limit)
    real, parameter :: trunc_scf = 0.95 ! SCF asymptotes to 1. as SD increases 
                    ! use this value when calculating SD to represent "full" coverage


    ! snow depletion curve parameters for IGBP snow depletion curve.
    real, dimension(20), parameter ::  & 
        mfsno_table = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00, &
                        2.00, 2.00, 2.00, 3.00, 3.00, 4.00, 4.00, &
                        2.50, 3.00, 3.00, 3.50, 3.50, 3.50 /)

    real, dimension(20), parameter ::  & 
        scffac_table = (/ 0.005, 0.005, 0.005, 0.005, 0.005, 0.008, &
                        0.008, 0.010, 0.010, 0.010, 0.010, 0.007, 0.021, &  
                        0.013, 0.015, 0.008, 0.015, 0.015, 0.015, 0.015 /)

    contains

! calling pattern 
    ! call calc_density(lensfc, lsm, landmask, swefcs, sndfcs, stcfcs, denfcs)
    ! call calcSWE_noah(lensfc, scfIMS, vetfcs_in, swefcs, sweIMS)
    ! ! calculate SD from IMS SCF
    ! call calcSD_noahmp(lensfc, scfIMS, vetfcs_in, denfcs, sndIMS)
    ! ! calculate SCF from model forecast SD and SWE (since not always in restart)
    ! call calcSCF_noahmp(lensfc, vetfcs_in, denfcs, sndfcs, scffcs)

!====================================
! calculate snow density from forecast fields.
! density = SWE/SND where snow present. 
!         = average from snow forecasts over land, where snow not present

 subroutine calc_density(lensfc, lsm, landmask, swe, snd, stc, density)

    implicit none 

    integer, intent(in) :: lensfc, lsm
    integer, intent(in) :: landmask(lensfc)
    real, intent(in)    :: swe(lensfc), snd(lensfc), stc(lensfc)
    real, intent(out)   :: density(lensfc)

    real :: dens_mean
    integer :: i

    ! density = swe/snd
    do i =1, lensfc
        if (snd(i) > 0.01 ) then
            density(i) = swe(i)/snd(i)
        elseif (snd(i) <= 0.01 .and. lsm==2) then ! for noah-MP, calculate from stc
        ! divide by 1000 as noah-MP has snd in m
            density(i) = max(80.0, min(120., &
                                       67.92+51.25*exp((stc(i)-273.15)/2.59)))/1000.
        endif
    enddo

    where (density < 0.0001) density = 0.08

    if (lsm==1) then ! for noah, use mean density 
        ! calculate mean density over land
        if (count (landmask==1 .and. snd> 0.01) > 0) then
            ! mean density over snow-covered land
            dens_mean = sum(density, mask = (landmask==1 .and. snd>0.01 )) &
                            / count (landmask==1 .and. snd> 0.01)
            ! print *, 'mean density: ', dens_mean
        else
            dens_mean = 0.1  ! default value if have no snow 
            print *, 'no snow, using default density: ', dens_mean
        endif
        ! for grid cells with no valid density, fill in the average snodens
        where( snd <= 0.01 ) density = dens_mean
    endif

 end subroutine calc_density

!====================================
! calculate SWE from fractional snow cover, using the noah model relationship 
! uses empirical inversion of snow depletion curve in the model 
 subroutine calcSWE_noah(lensfc, scfIMS, vetfcs_in, swefcs, sweIMS)
    
    implicit none
    !
    integer, intent(in)     :: lensfc
    real, intent(in)        :: vetfcs_in(lensfc)
    real, intent(in)        :: swefcs(lensfc)
    real, intent(inout)     :: scfIMS(lensfc) 
    real, intent(out)       :: sweIMS(lensfc)
    
    integer            :: vetfcs(lensfc)
    !snup_array is the swe (mm) at which scf reaches 100% 
    real               :: snupx(30), snup, salp, rsnow
    integer            :: i   !, vtype_int

    ! fill background values to nan
    ! sweIMS = nodata_real
    sweIMS = IEEE_VALUE(sweIMS, IEEE_QUIET_NAN)

    ! note: this is an empirical inversion of   snfrac rotuine in noah 
    !  should really have a land model check in here. 

    !this is for the igbp veg classification scheme.
    ! swe at which snow cover reaches 100%, in m
    snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
            0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
            0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
            0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)

    salp = -4.0
    vetfcs = nint(vetfcs_in)
    ! this is done in the noaa code, but we don't want to create a value outside land
    ! where(vetfcs==0) vetfcs = 7 
    do i = 1, lensfc 
        if ( abs( scfIMS(i) - nodata_real ) > nodata_tol ) then  ! is have IMS data
            if  (vetfcs(i)>0)  then ! if model has land
                snup = snupx(vetfcs(i))*1000. ! convert to mm
                if (snup == 0.) then
                    print*, " 0.0 snup value, check vegclasses", vetfcs(i)
                    stop
                endif

                ! if model and IMS both have 100% snow cover, don't convert IMS to a snow depth
                if ( (swefcs(i) >= snup)  .and. (scfIMS(i) >= 1.0 ) ) cycle 

                if (scfIMS(i) >= 1.0) then
                    rsnow = 1.
                elseif (scfIMS(i) < 0.001) then
                    rsnow = 0.0 
                else
                    rsnow = min(log(1. - scfIMS(i)) / salp, 1.0) 
                endif  
                ! return swe in mm 
                sweIMS(i) = rsnow * snup  !  mm
            else  ! if model is not land, remove the IMS data
                scfIMS(i) = nodata_real 
            endif
        endif
    enddo  
    return
    
 end subroutine calcSWE_noah

!====================================
! calculate IMS SD from fractional IMS snow cover, using the noah-MP model relationship
! (this is the inverse of calcSCF_noahmp). 
 subroutine calcSD_noahmp(lensfc, scfIMS, vetfcs_in, denfcs, sndIMS)

    implicit none
    !
    integer, intent(in)     :: lensfc 
    real, intent(in)        :: vetfcs_in(lensfc)
    real, intent(in)        :: denfcs(lensfc)
    real, intent(inout)     :: scfIMS(lensfc)
    real, intent(out)       :: sndIMS(lensfc)

    integer            :: vetfcs
    real               :: mfsno, scffac, bdsno, fmelt
    integer            :: i !, vtype_int

    ! fill background values to nan
    ! sndIMS = nodata_real
    sndIMS = IEEE_VALUE(sndIMS, IEEE_QUIET_NAN)

    do i = 1, lensfc
        if (abs(scfIMS(i) - nodata_real) > nodata_tol) then  ! if have IMS data
            vetfcs = int(vetfcs_in(i))
            if  (vetfcs>0)  then ! if model has land
                if (scfIMS(i) < 0.01 ) then 
                    sndIMS(i) = 0. 
                else 
                    ! calculate snow depth
                    mfsno  =  mfsno_table(vetfcs)
                    scffac = scffac_table(vetfcs)
                    bdsno  = max(50., min(650.,denfcs(i)*1000.) )   ! x1000, as noah-mp has SND in m.
                    fmelt  = (bdsno/100.)**mfsno
                    sndIMS(i) = (scffac * fmelt)*atanh(min(scfIMS(i), trunc_scf))*1000. ! x1000 into mm 
                endif
            endif
        endif
    enddo

 end subroutine calcSD_noahmp

!====================================
! calculate fractional snow cover from SD and density for Noah-MP
! copied from module_sf_noahmplsm.f90 
 subroutine calcSCF_noahmp(lensfc, vetfcs_in, denfcs, sndfcs, scffcs) 

    implicit none
    !
    integer, intent(in)     :: lensfc
    real, intent(in)        :: vetfcs_in(lensfc)
    real, intent(in)        :: denfcs(lensfc)
    real, intent(in)        :: sndfcs(lensfc)
    real, intent(out)       :: scffcs(lensfc)

    integer            :: vetfcs
    real               :: mfsno, scffac,  bdsno, fmelt, snowh
    integer            :: i !, vtype_int

    do i = 1, lensfc
        vetfcs = int(vetfcs_in(i))
        if (vetfcs>0) then ! if model has land
            mfsno  =  mfsno_table(vetfcs)
            scffac = scffac_table(vetfcs)
            snowh  = sndfcs(i)*0.001 ! into m
            if(sndfcs(i) .gt.0.)  then
                bdsno = denfcs(i)*1000. ! 1000, as noah-mp has snd in m 
                fmelt = (bdsno/100.)**mfsno
                scffcs(i) = tanh( snowh /(scffac * fmelt))
            else 
                scffcs(i) = 0.
            endif 
        endif 
    enddo 

 end subroutine calcSCF_noahmp

 

 end module M_SCF_Den
 
