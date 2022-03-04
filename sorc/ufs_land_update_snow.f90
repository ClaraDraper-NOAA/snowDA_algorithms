MODULE M_UFSLAND_SNOW_UPDATE

    use netcdf
    Use, Intrinsic :: IEEE_ARITHMETIC
    implicit none
    !
    include 'mpif.h'

    type noahmp_type
    Real, allocatable :: swe                (:)
    Real, allocatable :: snow_depth         (:)
    Real, allocatable :: active_snow_layers (:)
    Real, allocatable :: swe_previous       (:)
    Real, allocatable :: snow_soil_interface(:,:)
    Real, allocatable :: temperature_snow   (:,:)
    Real, allocatable :: snow_ice_layer     (:,:)
    Real, allocatable :: snow_liq_layer     (:,:)
    Real, allocatable :: temperature_soil   (:)
    end type noahmp_type    

    type observation_type
    Real, allocatable :: snow_depth (:)
    end type observation_type    

  contains

    subroutine UpdateTopLayer(myrank, LENSFC, noahmp, increment, analysis_1d, LANDMASK) !,  &                       !, vector_size, LENSFC_land,
                    ! & IDIM, JDIM, mp_start, tile_xy, Idim_xy, Jdim_xy, )

        integer                :: myrank, LENSFC   !, NUM_TILES, vector_size, LENSFC 
        ! INTEGER                :: IDIM, JDIM, mp_start
        type(noahmp_type)      :: noahmp
        ! type(observation_type) :: obs
        Real       :: increment(LENSFC), analysis_1d(LENSFC)
        Integer    :: LANDMASK(LENSFC)
        ! INTEGER    :: tile_xy(vector_size), Idim_xy(vector_size), Jdim_xy(vector_size), &
        !               LANDMASK_DUMY(NUM_TILES, LENSFC)
        Real       :: to_remove, layer_density, swe_increment, liq_ratio, remove_ratio
        integer    :: iloc, ilayer, iinter, active_layers, vector_loc, pathway
        ! integer    :: arr_indx, land_mask
        Real       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)

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

        !   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
        print*, "proc ", myrank, " starting update top layer"  !. vector length ", LENSFC
        ! IF (MYRANK==4) then
        !     print*, LANDMASK
        ! Endif
        do iloc = 1, LENSFC

            if( LANDMASK(iloc) == 1 ) then 
                ! print*, "iloc", iloc 
                ! IF (MYRANK==4) print*, land_mask                
                pathway = 0
                ! if(obs_snow_depth(iloc) == 0.0) then
                !     swe                (iloc)   = 0.0
                !     snow_depth         (iloc)   = 0.0
                !     active_snow_layers (iloc)   = 0.0
                !     swe_previous       (iloc)   = 0.0
                !     snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                !     temperature_snow   (iloc,:) = 0.0
                !     snow_ice_layer     (iloc,:) = 0.0
                !     snow_liq_layer     (iloc,:) = 0.0
                ! else

                active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)
               
                ! print*, "1D OI snowdepth", analysis_1d(iloc)
                ! print*,"active layers ", active_layers
                ! print*,"increment ", increment(iloc)
                ! print*,"snow_depth ", snow_depth(iloc)
                ! print*,"swe ", swe(iloc)
                ! print*,"swe_previous ", swe_previous(iloc)
                ! print*,"snow_soil_interface ", snow_soil_interface(iloc,:)
                ! print*,"temperature_snow ", temperature_snow(iloc,:)
                ! print*,"snow_ice_layer ", snow_ice_layer(iloc,:)
                ! print*,"snow_liq_layer ", snow_liq_layer(iloc,:)
                ! print*,"temperature_soil ", temperature_soil(iloc)
                if(active_layers < 0) then  ! in multi-layer mode

                    if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

                        pathway = 1

                        vector_loc = 4 + active_layers

                        layer_density = (snow_ice_layer(iloc,vector_loc)+snow_liq_layer(iloc,vector_loc)) / &
                                        (-snow_soil_interface(iloc,vector_loc))
                                 
                        swe_increment = increment(iloc) * layer_density / 1000.d0
                        
                        liq_ratio = snow_liq_layer(iloc,vector_loc) / &
                                    ( snow_ice_layer(iloc,vector_loc) + snow_liq_layer(iloc,vector_loc) )
                        snow_ice_layer(iloc,vector_loc) = snow_ice_layer(iloc,vector_loc) + &
                                                            (1.0 - liq_ratio) * swe_increment
                        snow_liq_layer(iloc,vector_loc) = snow_liq_layer(iloc,vector_loc) + &
                                                            liq_ratio * swe_increment
                        do ilayer = vector_loc, 3
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - increment(iloc)/1000.d0
                        end do

                        ! print*,"layer_density ", layer_density  
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio  
                        ! print*,"snow_ice_layer(iloc,vector_loc) ", snow_ice_layer(iloc,vector_loc)
                        ! print*,"snow_liq_layer(iloc,vector_loc) ", snow_liq_layer(iloc,vector_loc) 
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)

                    elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

                        pathway = 2

                        to_remove = increment(iloc)  ! depth[mm] to be removed

                        vector_loc = 4 + active_layers  ! location in vector of top layer
                        
                        layerloop: do ilayer = vector_loc, 3
                        
                        if(to_remove < 1000.d0*snow_soil_interface(iloc,ilayer)) then  ! this entire layer will be removed
                            
                            to_remove = to_remove - 1000.d0*snow_soil_interface(iloc,ilayer)
                            snow_ice_layer(iloc,ilayer)      = 0.0
                            snow_liq_layer(iloc,ilayer)      = 0.0
                            temperature_snow(iloc,ilayer)    = 0.0
                            do iinter = 3,ilayer, -1  ! remove snow from each snow layer
                            snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - snow_soil_interface(iloc,ilayer)
                            end do
                            ! print*,"layer_density ", layer_density
                            ! print*,"swe_increment ", swe_increment
                            ! print*,"liq_ratio ", liq_ratio  
                            ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:) 
                            ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                            ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)
                                    
                        else  ! this layer will be partially removed
                            
                            active_snow_layers(iloc) = ilayer - 4  ! new number of layers
                            remove_ratio = 1.d0 - (to_remove/1000.d0/snow_soil_interface(iloc,ilayer)) ! fraction to remove from layer
                            snow_ice_layer(iloc,ilayer) = remove_ratio * snow_ice_layer(iloc,ilayer)
                            snow_liq_layer(iloc,ilayer) = remove_ratio * snow_liq_layer(iloc,ilayer)
                            do iinter = ilayer, 3  ! remove snow from each snow layer
                            snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - to_remove/1000.d0
                            end do
                           
                            ! print*,"layer_density ", layer_density 
                            ! print*,"swe_increment ", swe_increment
                            ! print*,"liq_ratio ", liq_ratio
                            ! print*,"snow_ice_layer(iloc,ilayer) ", snow_ice_layer(iloc,ilayer)
                            ! print*,"snow_liq_layer(iloc,ilayer) ", snow_ice_layer(iloc,ilayer)
                            ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)
                                
                            ! print*,"exiting layerloop"
                        
                            exit layerloop

                        end if

                        end do layerloop
                        
                    end if  ! increment
                    
                    ! For multi-layer mode, recalculate interfaces and sum depth/swe

                    do ilayer = 4, 7
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
                    end do

                    snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

                    swe(iloc) = 0.0

                    do ilayer = 1, 3
                        swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
                    end do

                    swe_previous(iloc) = swe(iloc)

                    ! print*,"recalculated"
                    ! print*,"layer_density ", layer_density
                    ! print*,"swe_increment ", swe_increment
                    ! print*,"liq_ratio ", liq_ratio 
                    ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:)
                    ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                    ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)
        

                    if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
                        active_snow_layers (iloc) = 0.d0
                        snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                        temperature_snow   (iloc,:) = 0.0
                        snow_ice_layer     (iloc,:) = 0.0
                        snow_liq_layer     (iloc,:) = 0.0

                        ! print*,"out of multi layer mode ", layer_density 
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio
                        ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:)
                        ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)
                                
                    end if

                    elseif(active_layers == 0) then  ! snow starts in zero-layer mode

                    if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode

                        if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
                        pathway = 3
                        layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
                        else   ! use existing density
                        pathway = 4
                        layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        end if
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                        ! print*,"here 1 "
                        ! print*,"layer_density ", layer_density
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio  
                        ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:)  
                        ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)

                        if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
                        pathway = 5
                        active_snow_layers(iloc) = -1.0
                        snow_ice_layer(iloc,3)   = swe(iloc)
                        temperature_snow(iloc,3) = temperature_soil(iloc)
                        do ilayer = 3, 7
                            snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
                        end do

                        ! print*,"here 2"
                        ! print*,"layer_density ", layer_density  
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio  
                        ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:)
                        ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)
                            
                        end if
                        
                    elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

                        ! print*,"remove snow in zero-layer mode "
                        ! print*,"layer_density", layer_density  
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio 
                        ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:) 
                        ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:) 
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)

                        pathway = 6

                        if(snow_depth(iloc) <= 0.0) then 
                            print*, "proc ", myrank, " index ", iloc, &
                            " inconsistency in snow_depth ",snow_depth(iloc), "and increment", increment(iloc), &
                            " 1d analysis ", analysis_1d(iloc)
                            ! stop
                        endif
                        layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                        ! print*,"remove snow in zero-layer mode "
                        ! print*,"layer_density", layer_density  
                        ! print*,"swe_increment ", swe_increment
                        ! print*,"liq_ratio ", liq_ratio  
                        ! print*,"snow_ice_layer(iloc,:) ", snow_ice_layer(iloc,:)
                        ! print*,"snow_liq_layer(iloc,:) ", snow_liq_layer(iloc,:)  
                        ! print*,"snow_soil_interface(iloc:) ", snow_soil_interface(iloc,:)


                    end if  ! increment

                end if  ! active_layers
                    
                ! end if  ! obs_snow_depth == 0

                ! do some gross checks

                if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
                    print*, "proc ", myrank, " index ", iloc
                    print*, "Depth of soil not 2m pathway ", pathway
                    print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
                !      stop
                end if

                if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
                    print*, "proc ", myrank, " index ", iloc
                    print*, "snow_depth and snow_soil_interface inconsistent pathway ", pathway
                    print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if

                if(abs(analysis_1d(iloc) - snow_depth(iloc)) > 0.0000001) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "1D update snow and updated model snow inconsistent"
                    print*, pathway
                    print*, analysis_1d(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if

                if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
                    print*, "proc ", myrank, " index ", iloc
                    print*, "observed snow and updated model snow inconsistent pathway ", pathway
                    print*, snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if
                ! print*, "final"
                ! print*,"active layers ", active_layers
                ! print*,"increment ", increment(iloc)
                ! print*,"snow_depth ", snow_depth(iloc)
                ! print*,"swe ", swe(iloc)
                ! print*,"swe_previous ", swe_previous(iloc)
                ! print*,"snow_soil_interface ", snow_soil_interface(iloc,:)
                ! print*,"temperature_snow ", temperature_snow(iloc,:)
                ! print*,"snow_ice_layer ", snow_ice_layer(iloc,:)
                ! print*,"snow_liq_layer ", snow_liq_layer(iloc,:)
                ! print*,"temperature_soil ", temperature_soil(iloc)
            endif   ! land mask

        end do

        end associate

    end subroutine UpdateTopLayer

    subroutine UpdateTopLayer_single_cell(myrank, iloc, increment, analysis_1d, & !vector_length, noahmp, obs, obs_snow_depth, &
                            swe, snow_depth, active_snow_layers, swe_previous, &
                            snow_soil_interface, temperature_snow, snow_ice_layer, &
                            snow_liq_layer, temperature_soil)

        integer                :: myrank  !, vector_length 
        ! type(noahmp_type)      :: noahmp
        ! type(observation_type) :: obs
        Real       :: increment , analysis_1d 
        ! Real       :: obs_snow_depth 
        Real       :: swe, snow_depth, active_snow_layers, swe_previous, temperature_soil 
        Real       :: snow_soil_interface(7), temperature_snow(3), &
                      snow_ice_layer(3), snow_liq_layer(3)

        Real       :: to_remove, layer_density, swe_increment, liq_ratio, remove_ratio
        integer    :: iloc, ilayer, iinter, active_layers, pathway, vector_loc 
        ! integer    :: arr_indx, land_mask
        Real       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
        
        !   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
        ! print*, "proc ", myrank, " starting update top layer. vector length ", LENSFC_land
       
        ! do iloc = 1, vector_length 
            ! IF (MYRANK==4) print*, land_mask            
            pathway = 0
            ! if(obs_snow_depth(iloc) == 0.0) then

            !     swe                (iloc)   = 0.0
            !     snow_depth         (iloc)   = 0.0
            !     active_snow_layers (iloc)   = 0.0
            !     swe_previous       (iloc)   = 0.0
            !     snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
            !     temperature_snow   (iloc,:) = 0.0
            !     snow_ice_layer     (iloc,:) = 0.0
            !     snow_liq_layer     (iloc,:) = 0.0

            ! else

                active_layers = nint(active_snow_layers)  ! number of active layers (0,-1,-2,-3)

                if(active_layers < 0) then  ! in multi-layer mode

                if(increment > 0.0) then  ! add snow in multi-layer mode

                    pathway = 1

                    vector_loc = 4 + active_layers

                    layer_density = (snow_ice_layer(vector_loc)+snow_liq_layer(vector_loc)) / &
                                    (-snow_soil_interface(vector_loc))
                    swe_increment = increment * layer_density / 1000.d0
                    liq_ratio = snow_liq_layer(vector_loc) / &
                                ( snow_ice_layer(vector_loc) + snow_liq_layer(vector_loc) )
                    snow_ice_layer(vector_loc) = snow_ice_layer(vector_loc) + &
                                                        (1.0 - liq_ratio) * swe_increment
                    snow_liq_layer(vector_loc) = snow_liq_layer(vector_loc) + &
                                                        liq_ratio * swe_increment
                    do ilayer = vector_loc, 3
                        snow_soil_interface(ilayer) = snow_soil_interface(ilayer) - increment/1000.d0
                    end do
                    
                elseif(increment < 0.0) then  ! remove snow in multi-layer mode

                    pathway = 2

                    to_remove = increment  ! depth[mm] to be removed

                    vector_loc = 4 + active_layers  ! location in vector of top layer
                    
                    layerloop: do ilayer = vector_loc, 3
                    
                    if(to_remove < 1000.d0*snow_soil_interface(ilayer)) then  ! this entire layer will be removed
                        
                        to_remove = to_remove - 1000.d0*snow_soil_interface(ilayer)
                        snow_ice_layer(ilayer)      = 0.0
                        snow_liq_layer(ilayer)      = 0.0
                        temperature_snow(ilayer)    = 0.0
                        do iinter = 3,ilayer, -1  ! remove snow from each snow layer
                            snow_soil_interface(iinter) = snow_soil_interface(iinter) - snow_soil_interface(ilayer)
                        end do
                        
                    else  ! this layer will be partially removed
                        
                        active_snow_layers = ilayer - 4  ! new number of layers
                        remove_ratio = 1.d0 - (to_remove/1000.d0/snow_soil_interface(ilayer)) ! fraction to remove from layer
                        snow_ice_layer(ilayer) = remove_ratio * snow_ice_layer(ilayer)
                        snow_liq_layer(ilayer) = remove_ratio * snow_liq_layer(ilayer)
                        do iinter = ilayer, 3  ! remove snow from each snow layer
                        snow_soil_interface(iinter) = snow_soil_interface(iinter) - to_remove/1000.d0
                        end do
                    
                        exit layerloop

                    end if

                    end do layerloop
                    
                end if  ! increment
                
                ! For multi-layer mode, recalculate interfaces and sum depth/swe

                do ilayer = 4, 7
                    snow_soil_interface(ilayer) = snow_soil_interface(3) - soil_interfaces(ilayer)
                end do

                snow_depth = -snow_soil_interface(3) * 1000.d0

                swe = 0.0

                do ilayer = 1, 3
                    swe = swe + snow_ice_layer(ilayer) + snow_liq_layer(ilayer)
                end do

                swe_previous = swe

                if(snow_depth < 25.d0) then  ! go out of multi-layer mode
                    active_snow_layers = 0.d0
                    snow_soil_interface(:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                    temperature_snow   (:) = 0.0
                    snow_ice_layer     (:) = 0.0
                    snow_liq_layer     (:) = 0.0
                end if

                elseif(active_layers == 0) then  ! snow starts in zero-layer mode

                if(increment > 0.0) then  ! add snow in zero-layer mode

                    if(snow_depth == 0) then   ! no snow present, so assume density based on soil temperature
                    pathway = 3
                    layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil-273.15)/2.59)))
                    else   ! use existing density
                    pathway = 4
                    layer_density = swe / snow_depth * 1000.d0
                    end if
                    snow_depth = snow_depth + increment
                    swe = swe + increment * layer_density / 1000.d0
                    swe_previous = swe

                    active_snow_layers      = 0.0
                    snow_ice_layer(:)        = 0.0
                    snow_liq_layer(:)        = 0.0
                    temperature_snow(:)      = 0.0
                    snow_soil_interface(1:3) = 0.0

                    if(snow_depth > 25.0) then  ! snow depth is > 25mm so put in a layer
                    pathway = 5
                    active_snow_layers = -1.0
                    snow_ice_layer(3)   = swe
                    temperature_snow(3) = temperature_soil
                    do ilayer = 3, 7
                        snow_soil_interface(ilayer) = snow_soil_interface(ilayer) - snow_depth/1000.d0
                    end do
                    end if
                    
                elseif(increment < 0.0) then  ! remove snow in zero-layer mode

                    pathway = 6

                    if(snow_depth <= 0.0) then 
                        print*, "proc ", myrank, " index ", iloc, &
                        " inconsistency in snow_depth ",snow_depth, "and increment", increment, &
                        " 1d analysis ", analysis_1d
                        ! stop
                    endif
                    layer_density = swe / snow_depth * 1000.d0
                    snow_depth = snow_depth + increment
                    swe = swe + increment * layer_density / 1000.d0
                    swe_previous = swe

                    active_snow_layers     = 0.0
                    snow_ice_layer(:)        = 0.0
                    snow_liq_layer(:)        = 0.0
                    temperature_snow(:)      = 0.0
                    snow_soil_interface(1:3) = 0.0

                end if  ! increment

                end if  ! active_layers
                
            ! end if  ! obs_snow_depth == 0

            ! do some gross checks

            if(abs(snow_soil_interface(7) - snow_soil_interface(3) + 2.d0) > 0.0000001) then
                print*, "proc ", myrank, " index ", iloc
                print*, "Depth of soil not 2m pathway ", pathway
                print*, snow_soil_interface(7), snow_soil_interface(3)
            !      stop
            end if

            if(active_snow_layers < 0.0 .and. abs(snow_depth + 1000.d0*snow_soil_interface(3)) > 0.0000001) then
                print*, "proc ", myrank, " index ", iloc
                print*, "snow_depth and snow_soil_interface inconsistent pathway ", pathway
                print*, active_snow_layers, snow_depth, snow_soil_interface(3)
            !      stop
            end if

            if(abs(analysis_1d - snow_depth) > 0.0000001) then
                print*, " proc ", myrank, " loc ", iloc
                print*, "1D snow and updated model snow inconsistent"
                print*, "proc ", myrank, " index ", iloc
                print*, analysis_1d, snow_depth, snow_soil_interface(3)
            !      stop
            end if

            if(snow_depth < 0.0 .or. snow_soil_interface(3) > 0.0 ) then
                print*, "proc ", myrank, " index ", iloc
                print*, "observed snow and updated model snow inconsistent pathway ", pathway
                print*, snow_depth, snow_soil_interface(3)
            !      stop
            end if

        ! end do

    end subroutine UpdateTopLayer_single_cell

    subroutine UpdateBottomLayer(myrank, vector_length, noahmp, increment, analysis_1d, LANDMASK)

        type(noahmp_type)      :: noahmp
        ! type(observation_type) :: obs
        integer                :: myrank, vector_length
        Real       :: increment(vector_length), analysis_1d(vector_length)
        Integer    :: LANDMASK(vector_length)
        Real       :: to_remove
        Real       :: layer_density, swe_increment, liq_ratio, remove_ratio
        integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
        Real       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
        Real       :: temp_vector(3), layer_depths(3)

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

        !   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]

        do iloc = 1, vector_length
            if(LANDMASK(iloc) == 1) then

                pathway = 0

                ! if(obs_snow_depth(iloc) == 0.0) then

                !     swe                (iloc)   = 0.0
                !     snow_depth         (iloc)   = 0.0
                !     active_snow_layers (iloc)   = 0.0
                !     swe_previous       (iloc)   = 0.0
                !     snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                !     temperature_snow   (iloc,:) = 0.0
                !     snow_ice_layer     (iloc,:) = 0.0
                !     snow_liq_layer     (iloc,:) = 0.0

                ! else

                active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

                if(active_layers < 0) then  ! in multi-layer mode

                if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

                    pathway = 1

                    vector_loc = 3

                    layer_density = (snow_ice_layer(iloc,vector_loc)+snow_liq_layer(iloc,vector_loc)) / &
                                    (-snow_soil_interface(iloc,vector_loc))
                    swe_increment = increment(iloc) * layer_density / 1000.d0
                    liq_ratio = snow_liq_layer(iloc,vector_loc) / &
                                ( snow_ice_layer(iloc,vector_loc) + snow_liq_layer(iloc,vector_loc) )
                    snow_ice_layer(iloc,vector_loc) = snow_ice_layer(iloc,vector_loc) + &
                                                        (1.0 - liq_ratio) * swe_increment
                    snow_liq_layer(iloc,vector_loc) = snow_liq_layer(iloc,vector_loc) + &
                                                        liq_ratio * swe_increment
                    snow_soil_interface(iloc,3) = snow_soil_interface(iloc,3) - increment(iloc)/1000.d0
                    
                elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

                    pathway = 2
                    
                    layer_depths(1) = snow_soil_interface(iloc,1)  ! layer depth [m] (negative)
                    layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
                    layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

                    removed = 0

                    to_remove = increment(iloc)  ! depth[mm] to be removed

                    vector_loc = 4 + active_layers  ! location in vector of top layer
                    
                    layerloop: do ilayer = 3, vector_loc, -1
                    
                    if(to_remove < 1000.d0*layer_depths(ilayer)) then  ! this entire layer will be removed
                        
                        to_remove = to_remove - 1000.d0*layer_depths(ilayer)
                        removed = removed + 1

                        snow_ice_layer(iloc,ilayer)      = 0.0
                        snow_liq_layer(iloc,ilayer)      = 0.0
                        temperature_snow(iloc,ilayer)    = 0.0
                        snow_soil_interface(iloc,ilayer) = 0.0

                    else  ! this layer will be partially removed
                        
                        remove_ratio = 1.d0 - (to_remove/1000.d0/snow_soil_interface(iloc,ilayer)) ! fraction to remove from layer
                        snow_ice_layer(iloc,ilayer) = remove_ratio * snow_ice_layer(iloc,ilayer)
                        snow_liq_layer(iloc,ilayer) = remove_ratio * snow_liq_layer(iloc,ilayer)
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - to_remove/1000.d0

                        exit layerloop

                    end if

                    end do layerloop
                    
                    active_snow_layers(iloc) = active_snow_layers(iloc) + removed  ! new number of layers
                    
                    temp_vector = temperature_snow(iloc,:)
                    temperature_snow(iloc,:) = 0.0
                    temperature_snow(iloc,removed+1:3) = temp_vector(1:3-removed)

                    temp_vector = snow_ice_layer(iloc,:)
                    snow_ice_layer(iloc,:) = 0.0
                    snow_ice_layer(iloc,removed+1:3) = temp_vector(1:3-removed)

                    temp_vector = snow_liq_layer(iloc,:)
                    snow_liq_layer(iloc,:) = 0.0
                    snow_liq_layer(iloc,removed+1:3) = temp_vector(1:3-removed)

                    temp_vector = snow_soil_interface(iloc,:)
                    snow_soil_interface(iloc,:) = 0.0
                    snow_soil_interface(iloc,removed+1:3) = temp_vector(1:3-removed)
                    
                end if  ! increment
                
                ! For multi-layer mode, recalculate interfaces and sum depth/swe

                do ilayer = 4, 7
                    snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
                end do

                snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

                swe(iloc) = 0.0

                do ilayer = 1, 3
                    swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
                end do

                swe_previous(iloc) = swe(iloc)
                
                if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
                    active_snow_layers (iloc) = 0.d0
                    snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                    temperature_snow   (iloc,:) = 0.0
                    snow_ice_layer     (iloc,:) = 0.0
                    snow_liq_layer     (iloc,:) = 0.0
                end if

                elseif(active_layers == 0) then  ! snow starts in zero-layer mode

                if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode

                    if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
                    pathway = 3
                    layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
                    else   ! use existing density
                    pathway = 4
                    layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                    end if
                    snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                    swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                    swe_previous(iloc) = swe(iloc)

                    active_snow_layers(iloc)      = 0.0
                    snow_ice_layer(iloc,:)        = 0.0
                    snow_liq_layer(iloc,:)        = 0.0
                    temperature_snow(iloc,:)      = 0.0
                    snow_soil_interface(iloc,1:3) = 0.0

                    if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
                    pathway = 5
                    active_snow_layers(iloc) = -1.0
                    snow_ice_layer(iloc,3)   = swe(iloc)
                    temperature_snow(iloc,3) = temperature_soil(iloc)
                    do ilayer = 3, 7
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
                    end do
                    end if
                    
                elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

                    pathway = 6

                    if(snow_depth(iloc) <= 0.0) then 
                        print*, " proc ", myrank, " loc ", iloc
                        print*, "inconsistency in snow_depth and increment"
                         stop
                    endif
                    layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                    snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                    swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                    swe_previous(iloc) = swe(iloc)

                    active_snow_layers(iloc)      = 0.0
                    snow_ice_layer(iloc,:)        = 0.0
                    snow_liq_layer(iloc,:)        = 0.0
                    temperature_snow(iloc,:)      = 0.0
                    snow_soil_interface(iloc,1:3) = 0.0

                end if  ! increment

                end if  ! active_layers
                
                ! end if  ! obs_snow_depth == 0

                ! do some gross checks

                if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "Depth of soil not 2m"
                    print*, pathway
                    print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
                !      stop
                end if

                if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "snow_depth and snow_soil_interface inconsistent"
                    print*, pathway
                    print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if

                ! if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
                !     print*, "observed snow and updated model snow inconsistent"
                !     print*, pathway
                !     print*, obs_snow_depth(iloc), snow_depth(iloc)
                ! !      stop
                ! end if

                if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "observed snow and updated model snow inconsistent"
                    print*, pathway
                    print*, snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if
            endif   ! land
        end do

        end associate

    end subroutine UpdateBottomLayer

    subroutine UpdateAllLayers_Ens(ens_mem, myrank, vector_length, noahmp, increment, &
                               analysis_1d, SNODENS_Background, LANDMASK)

        type(noahmp_type)      :: noahmp
        ! type(observation_type) :: obs
        integer                :: myrank, vector_length, ens_mem
        Real       :: increment(vector_length), analysis_1d(vector_length)
        Real       :: SNODENS_Background(vector_length)
        Integer    :: LANDMASK(vector_length)
        Real       :: to_remove
        Real       :: layer_density, swe_increment, liq_ratio, remove_ratio
        integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
        Real       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
        Real       :: partition_ratio, layer_depths(3)

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


        !   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]

        do iloc = 1, vector_length
            if(LANDMASK(iloc) == 1) then

                pathway = 0

                ! if(obs_snow_depth(iloc) == 0.0) then

                !     swe                (iloc)   = 0.0
                !     snow_depth         (iloc)   = 0.0
                !     active_snow_layers (iloc)   = 0.0
                !     swe_previous       (iloc)   = 0.0
                !     snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                !     temperature_snow   (iloc,:) = 0.0
                !     snow_ice_layer     (iloc,:) = 0.0
                !     snow_liq_layer     (iloc,:) = 0.0

                ! else

                active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

                if(active_layers < 0) then  ! in multi-layer mode
                
                    layer_depths(1) = snow_soil_interface(iloc,1)
                    layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
                    layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

                    !if(iloc ==  10962)then
                    !print*, 'proc ', myrank, ' depths',  iloc,obs_snow_depth(iloc)  ,snow_depth(iloc)
                    !print*, 'increment',  increment(iloc)  ,layer_depths
                    !print*, 'interfaces',    snow_soil_interface(iloc,:)
                    !print*, 'ice',    snow_ice_layer(iloc,:), swe(iloc)
                    !print*, 'liq',    snow_liq_layer(iloc,:)
                    !end if
                    if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

                        pathway = 1

                        vector_loc = 4 + active_layers  ! location in vector of top layer
                        
                        layerloop: do ilayer = vector_loc, 3
                        
                            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
        !10.26.21 CSD changed the denom to -layer_depths(ilayer)
                            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                                                (-layer_depths(ilayer))    !(-snow_soil_interface(iloc,ilayer))
                            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
                            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                                            ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
                            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                                                (1.0 - liq_ratio) * swe_increment
                            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                                                liq_ratio * swe_increment
                            do iinter = ilayer, 3  ! remove snow from each snow layer
                                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                                    partition_ratio * increment(iloc)/1000.d0
                            end do

                        end do layerloop
                        
                    elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

                        pathway = 2
                        
                        vector_loc = 4 + active_layers  ! location in vector of top layer
                        
                        layerloop1: do ilayer = vector_loc, 3
                        
                            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
                            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                                                (-layer_depths(ilayer))
                            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
                            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                                            ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
                            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                                                (1.0 - liq_ratio) * swe_increment
                            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                                                liq_ratio * swe_increment
                            do iinter = ilayer, 3  ! remove snow from each snow layer
                                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                                    partition_ratio * increment(iloc)/1000.d0
                            end do

                        end do layerloop1
                        
                    end if  ! increment
                    
                    ! For multi-layer mode, recalculate interfaces and sum depth/swe

                    do ilayer = 4, 7
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
                    end do

                    snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

                    swe(iloc) = 0.0

                    do ilayer = 1, 3
                        swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
                    end do

                    swe_previous(iloc) = swe(iloc)

                    if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
                        active_snow_layers (iloc) = 0.d0
                        snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                        temperature_snow   (iloc,:) = 0.0
                        snow_ice_layer     (iloc,:) = 0.0
                        snow_liq_layer     (iloc,:) = 0.0
                    end if

                    elseif(active_layers == 0) then  ! snow starts in zero-layer mode

                    if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode

                        if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
                        pathway = 3
                        layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
                        else   ! use existing density
                        pathway = 4
                        layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        end if
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                        if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
                        pathway = 5
                        active_snow_layers(iloc) = -1.0
                        snow_ice_layer(iloc,3)   = swe(iloc)
                        temperature_snow(iloc,3) = temperature_soil(iloc)
                        do ilayer = 3, 7
                            snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
                        end do
                        end if
                        
                    elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

                        pathway = 6

                        if(snow_depth(iloc) <= 0.0) then 
                            
                            ! print*, " proc ", myrank, " loc ", iloc, "ens ", ens_mem
                            ! print*, " inconsistency in snow_depth and increment"

                            !print*, " exiting"
                            ! stop
                            layer_density = SNODENS_Background(iloc)
                        else
                            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        endif
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                    end if  ! increment

                end if  ! active_layers
                
                ! end if  ! obs_snow_depth == 0

                ! do some gross checks

                if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
                    
                    ! print*, " proc ", myrank, " loc ", iloc
                    ! print*, "Depth of soil not 2m"
                    ! print*, pathway
                    ! print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
                
                    !      stop
                end if

                if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
                !     print*, " proc ", myrank, " loc ", iloc
                !     print*, "snow_depth and snow_soil_interface inconsistent"
                !     print*, pathway
                !     print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
                ! !      stop
                end if

                ! if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
               
                !     print*, "observed snow and updated model snow inconsistent"
                !     print*, pathway
                !     print*, obs_snow_depth(iloc), snow_depth(iloc)
               
                ! !      stop
                ! end if

                if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
                   
                    ! if (print_noah_debug) then 
                    !     print*, " proc ", myrank, " loc ", iloc, "ens ", ens_mem
                    !     print*, "observed snow and updated model snow inconsistent"
                    !     print*, pathway
                    !     print*, snow_depth(iloc), snow_soil_interface(iloc,3)
                    ! endif
               
                    !      stop
                end if

            endif !land

        end do

        end associate

    end subroutine UpdateAllLayers_Ens

   subroutine UpdateAllLayers(myrank, vector_length, noahmp, increment, &
                               analysis_1d, SNODENS_Background, LANDMASK)

        type(noahmp_type)      :: noahmp
        ! type(observation_type) :: obs
        integer                :: myrank, vector_length
        Real       :: increment(vector_length), analysis_1d(vector_length)
        Real       :: SNODENS_Background(vector_length)
        Integer    :: LANDMASK(vector_length)
        Real       :: to_remove
        Real       :: layer_density, swe_increment, liq_ratio, remove_ratio
        integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
        Real       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
        Real       :: partition_ratio, layer_depths(3)

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


        !   increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]

        do iloc = 1, vector_length
            if(LANDMASK(iloc) == 1) then

                pathway = 0

                ! if(obs_snow_depth(iloc) == 0.0) then

                !     swe                (iloc)   = 0.0
                !     snow_depth         (iloc)   = 0.0
                !     active_snow_layers (iloc)   = 0.0
                !     swe_previous       (iloc)   = 0.0
                !     snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                !     temperature_snow   (iloc,:) = 0.0
                !     snow_ice_layer     (iloc,:) = 0.0
                !     snow_liq_layer     (iloc,:) = 0.0

                ! else

                active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

                if(active_layers < 0) then  ! in multi-layer mode
                
                    layer_depths(1) = snow_soil_interface(iloc,1)
                    layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
                    layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

                    !if(iloc ==  10962)then
                    !print*, 'proc ', myrank, ' depths',  iloc,obs_snow_depth(iloc)  ,snow_depth(iloc)
                    !print*, 'increment',  increment(iloc)  ,layer_depths
                    !print*, 'interfaces',    snow_soil_interface(iloc,:)
                    !print*, 'ice',    snow_ice_layer(iloc,:), swe(iloc)
                    !print*, 'liq',    snow_liq_layer(iloc,:)
                    !end if
                    if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

                        pathway = 1

                        vector_loc = 4 + active_layers  ! location in vector of top layer
                        
                        layerloop: do ilayer = vector_loc, 3
                        
                            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
        !10.26.21 CSD changed the denom to -layer_depths(ilayer)
                            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                                                (-layer_depths(ilayer))    !(-snow_soil_interface(iloc,ilayer))
                            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
                            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                                            ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
                            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                                                (1.0 - liq_ratio) * swe_increment
                            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                                                liq_ratio * swe_increment
                            do iinter = ilayer, 3  ! remove snow from each snow layer
                                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                                    partition_ratio * increment(iloc)/1000.d0
                            end do

                        end do layerloop
                        
                    elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

                        pathway = 2
                        
                        vector_loc = 4 + active_layers  ! location in vector of top layer
                        
                        layerloop1: do ilayer = vector_loc, 3
                        
                            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
                            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                                                (-layer_depths(ilayer))
                            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
                            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                                            ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
                            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                                                (1.0 - liq_ratio) * swe_increment
                            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                                                liq_ratio * swe_increment
                            do iinter = ilayer, 3  ! remove snow from each snow layer
                                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                                    partition_ratio * increment(iloc)/1000.d0
                            end do

                        end do layerloop1
                        
                    end if  ! increment
                    
                    ! For multi-layer mode, recalculate interfaces and sum depth/swe

                    do ilayer = 4, 7
                        snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
                    end do

                    snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

                    swe(iloc) = 0.0

                    do ilayer = 1, 3
                        swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
                    end do

                    swe_previous(iloc) = swe(iloc)

                    if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
                        active_snow_layers (iloc) = 0.d0
                        snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
                        temperature_snow   (iloc,:) = 0.0
                        snow_ice_layer     (iloc,:) = 0.0
                        snow_liq_layer     (iloc,:) = 0.0
                    end if

                    elseif(active_layers == 0) then  ! snow starts in zero-layer mode

                    if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode

                        if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
                        pathway = 3
                        layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
                        else   ! use existing density
                        pathway = 4
                        layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        end if
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                        if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
                        pathway = 5
                        active_snow_layers(iloc) = -1.0
                        snow_ice_layer(iloc,3)   = swe(iloc)
                        temperature_snow(iloc,3) = temperature_soil(iloc)
                        do ilayer = 3, 7
                            snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
                        end do
                        end if
                        
                    elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

                        pathway = 6

                        if(snow_depth(iloc) <= 0.0) then 
                            print*, " proc ", myrank, " loc ", iloc
                            print*, " inconsistency in snow_depth and increment"
                            !print*, " exiting"
                            ! stop
                            layer_density = SNODENS_Background(iloc)
                        else
                            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
                        endif
                        snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
                        swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
                        swe_previous(iloc) = swe(iloc)

                        active_snow_layers(iloc)      = 0.0
                        snow_ice_layer(iloc,:)        = 0.0
                        snow_liq_layer(iloc,:)        = 0.0
                        temperature_snow(iloc,:)      = 0.0
                        snow_soil_interface(iloc,1:3) = 0.0

                    end if  ! increment

                end if  ! active_layers
                
                ! end if  ! obs_snow_depth == 0

                ! do some gross checks

                if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "Depth of soil not 2m"
                    print*, pathway
                    print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
                !      stop
                end if

                if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "snow_depth and snow_soil_interface inconsistent"
                    print*, pathway
                    print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if

                ! if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
                !     print*, "observed snow and updated model snow inconsistent"
                !     print*, pathway
                !     print*, obs_snow_depth(iloc), snow_depth(iloc)
                ! !      stop
                ! end if

                if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
                    print*, " proc ", myrank, " loc ", iloc
                    print*, "observed snow and updated model snow inconsistent"
                    print*, pathway
                    print*, snow_depth(iloc), snow_soil_interface(iloc,3)
                !      stop
                end if

            endif !land

        end do

        end associate

    end subroutine UpdateAllLayers


    subroutine WriteRestartNoahMP(myrank, vector_restart_prefix, IY, IM, ID, IH, &
                                 N_sA, N_sA_Ext, mpiReal_size, LENSFC_land, &
                                 vector_length, NUM_TILES, noahmp)
                                    ! IDIM, JDIM, tile_xy, Idim_xy, Jdim_xy)
        use netcdf

        type(noahmp_type)  :: noahmp
        CHARACTER(LEN=*)  :: vector_restart_prefix
        CHARACTER(LEN=250)   :: filename
        integer           :: IY, IM, ID, IH
        integer            :: myrank, N_sA, N_sA_Ext, mpiReal_size, &
                              LENSFC_land, vector_length, NUM_TILES
        integer            :: ncid, varid, status, ixy, iv
        ! Integer            :: tile_xy , Idim_xy , &
        !                       Jdim_xy , IDIM, JDIM
        integer            :: IERR, dest_Aoffset, dest_Aoffset_end, arLen
        Real       :: DUMMY1(vector_length), DUMMY2(vector_length), DUMMY3(vector_length), &
                      DUMMY4(vector_length), DUMMY5(vector_length, 7), &
        DUMMY6(vector_length, 3), DUMMY7(vector_length, 3), DUMMY8(vector_length, 3)

        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str
        LOGICAL              :: file_exists

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
        Endif
        if (MYRANK == 0) then
            DUMMY1(1:1+N_sA) = noahmp%swe
            DUMMY2(1:1+N_sA) = noahmp%snow_depth
            DUMMY3(1:1+N_sA) = noahmp%active_snow_layers
            DUMMY4(1:1+N_sA) = noahmp%swe_previous
            Do iv=1, 7 
                DUMMY5(1:1+N_sA, iv) = noahmp%snow_soil_interface(:,iv)
            Enddo
            Do iv=1, 3
                DUMMY6(1:1+N_sA, iv) = noahmp%temperature_snow(:,iv)
                DUMMY7(1:1+N_sA, iv) = noahmp%snow_ice_layer(:,iv)
                DUMMY8(1:1+N_sA, iv) = noahmp%snow_liq_layer(:,iv)
            Enddo
            Do ixy =  1, NUM_TILES-1   ! sender proc index within tile group
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
            End do
             
            write(y_str, "(I4)") IY
            write(m_str, "(I0.2)") IM
            write(d_str, "(I0.2)") ID
            write(h_str, "(I0.2)") IH
            filename = trim(vector_restart_prefix)//"/ufs_land_restart."// &
                TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc"
            INQUIRE(FILE=trim(filename), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(filename) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif
            status = nf90_open(trim(filename), NF90_WRITE, ncid)
            ! Start writing restart file
            status = nf90_inq_varid(ncid, "weasd", varid)
            status = nf90_put_var(ncid, varid , DUMMY1   , &
                start = (/1,1/), count = (/vector_length, 1/))

            status = nf90_inq_varid(ncid, "snwdph", varid)
            status = nf90_put_var(ncid, varid , DUMMY2,  &
                start = (/1,1/), count = (/vector_length, 1/))

            status = nf90_inq_varid(ncid, "snowxy", varid)
            status = nf90_put_var(ncid, varid , DUMMY3  , &
                start = (/1,1/), count = (/vector_length, 1/))

            status = nf90_inq_varid(ncid, "sneqvoxy", varid)
            status = nf90_put_var(ncid, varid , DUMMY4, &
                start = (/1,1/), count = (/vector_length, 1/))

            status = nf90_inq_varid(ncid, "zsnsoxy", varid)
            status = nf90_put_var(ncid, varid , DUMMY5 , &
                start = (/1            , 1, 1/), &
                count = (/vector_length, 7, 1/))

            status = nf90_inq_varid(ncid, "tsnoxy", varid)
            status = nf90_put_var(ncid, varid , DUMMY6  , &
                start = (/1            , 1, 1/), &
                count = (/vector_length, 3, 1/))

            status = nf90_inq_varid(ncid, "snicexy", varid)
            status = nf90_put_var(ncid, varid , DUMMY7 , &
                start = (/1            , 1, 1/), &
                count = (/vector_length, 3, 1/))

            status = nf90_inq_varid(ncid, "snliqxy", varid)
            status = nf90_put_var(ncid, varid , DUMMY8 , &
                start = (/1            , 1, 1/), &
                count = (/vector_length, 3, 1/))

            status = nf90_close(ncid)

        endif

    end subroutine WriteRestartNoahMP   
       
    subroutine ReadRestartNoahMP_Ens(myrank, vector_length, &
                    vector_restart_prefix, noda_inp_path, &
                    y_str, m_str, d_str, h_str, &
                    ens_size, startIndx, endIndx, &
                    noahmp, SNDnoDA, SWEnoDA, SCF_Grid) !, SNDFCS, SWEFCS) LENSFC, 

        use netcdf
        
        implicit none 
        integer           :: myrank, ens_size, startIndx, endIndx, vector_length !, LENSFC
        type(noahmp_type) :: noahmp(ens_size)
        CHARACTER(LEN=*)  :: vector_restart_prefix, noda_inp_path
        CHARACTER(LEN=*)  :: y_str, m_str, d_str, h_str
        Real              :: SNDnoDA(vector_length), SWEnoDA(vector_length)
        Real              :: SCF_Grid(ens_size+1, vector_length)
        ! Real              :: SNDFCS(vector_length), SWEFCS(vector_length)
        ! CHARACTER(LEN=250)   :: filename
        
        integer           :: ncid, dimid, varid, status, ie
        
        CHARACTER(LEN=500)  :: forc_inp_file
        LOGICAL           :: file_exists
        CHARACTER(LEN=2)       :: ensCH 
        
        ! vector_length =  endIndx - startIndx + 1
        
        Do ie = 1, ens_size
            WRITE(ensCH, '(I2.2)') ie

            forc_inp_file=TRIM(vector_restart_prefix)//"/ens"//ensCH//"/ufs_land_restart."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc" 

            INQUIRE(FILE=trim(forc_inp_file), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(forc_inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif        
        
            status = nf90_open(trim(forc_inp_file), NF90_NOWRITE, ncid)
            CALL NETCDF_ERR(status, 'OPENING FILE: '//TRIM(forc_inp_file) )

            status = nf90_inq_varid(ncid, "weasd", varid)
            CALL NETCDF_ERR(status, 'getting varid: weasd' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%swe   , &
                start = (/startIndx,1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: weasd' )

            status = nf90_inq_varid(ncid, "snwdph", varid)
            CALL NETCDF_ERR(status, 'getting varid: snwdph' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%snow_depth  , &
                start = (/startIndx,1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: snwdph' )

            status = nf90_inq_varid(ncid, "snowxy", varid)
            CALL NETCDF_ERR(status, 'getting varid: snowxy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%active_snow_layers  , &
                start = (/startIndx,1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: snowxy' )

            status = nf90_inq_varid(ncid, "sneqvoxy", varid)
            CALL NETCDF_ERR(status, 'getting varid: sneqvoxy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%swe_previous, &
                start = (/startIndx,1/), count = (/vector_length, 1/)) 
            CALL NETCDF_ERR(status, 'reading: sneqvoxy' )

            status = nf90_inq_varid(ncid, "tsnoxy", varid)
            CALL NETCDF_ERR(status, 'getting varid: tsnoxy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%temperature_snow  , &
                start = (/startIndx            , 1, 1/)                , &
                count = (/vector_length, 3, 1/))
            CALL NETCDF_ERR(status, 'reading: tsnoxy' )

            status = nf90_inq_varid(ncid, "zsnsoxy", varid)
            CALL NETCDF_ERR(status, 'getting varid: zsnsoxy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%snow_soil_interface , &
                start = (/startIndx            , 1, 1/)                , &
                count = (/vector_length, 7, 1/))
            CALL NETCDF_ERR(status, 'reading: zsnsoxy' )

            status = nf90_inq_varid(ncid, "snicexy", varid)
            CALL NETCDF_ERR(status, 'getting varid: snicexy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%snow_ice_layer , &
                start = (/startIndx            , 1, 1/)                , &
                count = (/vector_length, 3, 1/))
            CALL NETCDF_ERR(status, 'reading: snicexy' )

            status = nf90_inq_varid(ncid, "snliqxy", varid)
            CALL NETCDF_ERR(status, 'getting varid: snliqxy' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%snow_liq_layer , &
                start = (/startIndx            , 1, 1/)                , &
                count = (/vector_length, 3, 1/))
            CALL NETCDF_ERR(status, 'reading: snliqxy' )

            status = nf90_inq_varid(ncid, "stc", varid)
            CALL NETCDF_ERR(status, 'getting varid: stc' )
            status = nf90_get_var(ncid, varid , noahmp(ie)%temperature_soil, &
                start = (/startIndx            , 1, 1/)                , &
                count = (/vector_length, 1, 1/))
            CALL NETCDF_ERR(status, 'reading: stc' )

! double sncovr1(time, location) ;
! sncovr1:long_name = "snow cover over land" ;
! double snowc(time, location) ;
! snowc:long_name = "fractional snow cover" ;
            status = nf90_inq_varid(ncid, "snowc", varid)
            CALL NETCDF_ERR(status, 'getting varid: snowc' )
            status = nf90_get_var(ncid, varid, SCF_Grid(ie, :), &
                start = (/startIndx, 1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: snowc' )

            status = nf90_close(ncid)
            CALL NETCDF_ERR(status, 'Closing FILE: '//TRIM(forc_inp_file) )
        End do
        
        ! if (myrank == 0) then
        INQUIRE(FILE=trim(noda_inp_path), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(noda_inp_path) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif   
        ! Allocate(SNDnoDA(LENSFC))
        ! Allocate(SWEnoDA(LENSFC))     
    
        status = nf90_open(trim(noda_inp_path), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'Opening FILE: '//TRIM(noda_inp_path) )

        status = nf90_inq_varid(ncid, "weasd", varid)
        CALL NETCDF_ERR(status, 'getting varid: weasd'//TRIM(noda_inp_path) )
        status = nf90_get_var(ncid, varid , SWEnoDA   , &
            start = (/startIndx,1/), count = (/vector_length, 1/))
        CALL NETCDF_ERR(status, 'reading: weasd'//TRIM(noda_inp_path) )

        status = nf90_inq_varid(ncid, "snwdph", varid)
        CALL NETCDF_ERR(status, 'getting varid: snwdph'//TRIM(noda_inp_path) )
        status = nf90_get_var(ncid, varid, SNDnoDA, &
            start = (/startIndx,1/), count = (/vector_length, 1/))
        CALL NETCDF_ERR(status, 'reading: weasd'//TRIM(noda_inp_path) )

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'Closing FILE: '//TRIM(noda_inp_path) )
        ! endif

    end subroutine ReadRestartNoahMP_Ens

    subroutine ReadRestartNoah_Ens(myrank, vector_length, &
                    vector_restart_prefix, noda_inp_path, &
                    y_str, m_str, d_str, h_str, &
                    ens_size, startIndx, endIndx, &
                    SNDnoDA, SWEnoDA, SNDFCS, SWEFCS, &
                    SCF_Grid, SCF_Grid_Land) !, SNDFCS, SWEFCS) !noahmp, LENSFC, 

        use netcdf
        
        implicit none 
        integer           :: myrank, ens_size, startIndx, endIndx, vector_length !, LENSFC
        ! type(noahmp_type) :: noahmp(ens_size)
        CHARACTER(LEN=*)  :: vector_restart_prefix, noda_inp_path
        CHARACTER(LEN=*)  :: y_str, m_str, d_str, h_str
        Real              :: SNDnoDA(vector_length), SWEnoDA(vector_length)
        Real     :: SCF_Grid(ens_size+1, vector_length), SCF_Grid_Land(ens_size+1, vector_length)
        Real     :: SNDFCS(ens_size+1, vector_length), SWEFCS(ens_size+1, vector_length)
        ! CHARACTER(LEN=250)   :: filename
        
        integer           :: ncid, dimid, varid, status, ie
        
        CHARACTER(LEN=500)  :: forc_inp_file
        LOGICAL           :: file_exists
        CHARACTER(LEN=2)       :: ensCH 
        
        ! vector_length =  endIndx - startIndx + 1
        Do ie = 1, ens_size
            WRITE(ensCH, '(I2.2)') ie

            forc_inp_file=TRIM(vector_restart_prefix)//"/ens"//ensCH//"/ufs_land_restart."// &
            TRIM(y_str)//"-"//TRIM(m_str)//"-"//TRIM(d_str)//"_"//TRIM(h_str)//"-00-00.nc" 

            INQUIRE(FILE=trim(forc_inp_file), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(forc_inp_file) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif        
        
            status = nf90_open(trim(forc_inp_file), NF90_NOWRITE, ncid)
            CALL NETCDF_ERR(status, 'OPENING FILE: '//TRIM(forc_inp_file) )

            status = nf90_inq_varid(ncid, "weasd", varid)
            CALL NETCDF_ERR(status, 'getting varid: weasd' )
            status = nf90_get_var(ncid, varid , SWEFCS(ie, :),     &
                start = (/startIndx,1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: weasd' )

            status = nf90_inq_varid(ncid, "snwdph", varid)
            CALL NETCDF_ERR(status, 'getting varid: snwdph' )
            status = nf90_get_var(ncid, varid , SNDFCS(ie, :), &
                start = (/startIndx,1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: snwdph' )
! double sncovr1(time, location) ;
! sncovr1:long_name = "snow cover over land" ;
! double snowc(time, location) ;
! snowc:long_name = "fractional snow cover" ;
            status = nf90_inq_varid(ncid, "snowc", varid)
            CALL NETCDF_ERR(status, 'getting varid: snowc' )
            status = nf90_get_var(ncid, varid, SCF_Grid(ie, :), &
                start = (/startIndx, 1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: snowc' )

            status = nf90_inq_varid(ncid, "sncovr1", varid)
            CALL NETCDF_ERR(status, 'getting varid: sncovr1' )
            status = nf90_get_var(ncid, varid, SCF_Grid_Land(ie, :), &
                start = (/startIndx, 1/), count = (/vector_length, 1/))
            CALL NETCDF_ERR(status, 'reading: sncovr1' )

            ! status = nf90_inq_varid(ncid, "snowxy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: snowxy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%active_snow_layers  , &
            !     start = (/startIndx,1/), count = (/vector_length, 1/))
            ! CALL NETCDF_ERR(status, 'reading: snowxy' )

            ! status = nf90_inq_varid(ncid, "sneqvoxy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: sneqvoxy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%swe_previous, &
            !     start = (/startIndx,1/), count = (/vector_length, 1/)) 
            ! CALL NETCDF_ERR(status, 'reading: sneqvoxy' )

            ! status = nf90_inq_varid(ncid, "tsnoxy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: tsnoxy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%temperature_snow  , &
            !     start = (/startIndx            , 1, 1/)                , &
            !     count = (/vector_length, 3, 1/))
            ! CALL NETCDF_ERR(status, 'reading: tsnoxy' )

            ! status = nf90_inq_varid(ncid, "zsnsoxy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: zsnsoxy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%snow_soil_interface , &
            !     start = (/startIndx            , 1, 1/)                , &
            !     count = (/vector_length, 7, 1/))
            ! CALL NETCDF_ERR(status, 'reading: zsnsoxy' )

            ! status = nf90_inq_varid(ncid, "snicexy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: snicexy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%snow_ice_layer , &
            !     start = (/startIndx            , 1, 1/)                , &
            !     count = (/vector_length, 3, 1/))
            ! CALL NETCDF_ERR(status, 'reading: snicexy' )

            ! status = nf90_inq_varid(ncid, "snliqxy", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: snliqxy' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%snow_liq_layer , &
            !     start = (/startIndx            , 1, 1/)                , &
            !     count = (/vector_length, 3, 1/))
            ! CALL NETCDF_ERR(status, 'reading: snliqxy' )

            ! status = nf90_inq_varid(ncid, "stc", varid)
            ! CALL NETCDF_ERR(status, 'getting varid: stc' )
            ! status = nf90_get_var(ncid, varid , noahmp(ie)%temperature_soil, &
            !     start = (/startIndx            , 1, 1/)                , &
            !     count = (/vector_length, 1, 1/))
            ! CALL NETCDF_ERR(status, 'reading: stc' )

            status = nf90_close(ncid)
            CALL NETCDF_ERR(status, 'Closing FILE: '//TRIM(forc_inp_file) )
        End do
        
        ! if (myrank == 0) then
        INQUIRE(FILE=trim(noda_inp_path), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(noda_inp_path) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif   
        ! Allocate(SNDnoDA(LENSFC))
        ! Allocate(SWEnoDA(LENSFC))     
    
        status = nf90_open(trim(noda_inp_path), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'Opening FILE: '//TRIM(noda_inp_path) )

        status = nf90_inq_varid(ncid, "weasd", varid)
        CALL NETCDF_ERR(status, 'getting varid: weasd'//TRIM(noda_inp_path) )
        status = nf90_get_var(ncid, varid , SWEnoDA   , &
            start = (/startIndx,1/), count = (/vector_length, 1/))
        CALL NETCDF_ERR(status, 'reading: weasd'//TRIM(noda_inp_path) )

        status = nf90_inq_varid(ncid, "snwdph", varid)
        CALL NETCDF_ERR(status, 'getting varid: snwdph'//TRIM(noda_inp_path) )
        status = nf90_get_var(ncid, varid, SNDnoDA, &
            start = (/startIndx,1/), count = (/vector_length, 1/))
        CALL NETCDF_ERR(status, 'reading: weasd'//TRIM(noda_inp_path) )

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'Closing FILE: '//TRIM(noda_inp_path) )
        ! endif

    end subroutine ReadRestartNoah_Ens

    subroutine ReadRestartNoahMP(myrank, forc_inp_file, noda_inp_path, &
                    startIndx, endIndx, vector_length, &
                    LENSFC, noahmp, SNDnoDA, SWEnoDA, SCF_Grid) !, SNDFCS, SWEFCS)

        use netcdf

        type(noahmp_type) :: noahmp
        integer           :: myrank, startIndx, endIndx, vector_length, LENSFC
        CHARACTER(LEN=*)  :: forc_inp_file, noda_inp_path
        Real              :: SNDnoDA(LENSFC), SWEnoDA(LENSFC)
        Real              :: SCF_Grid(vector_length)
        ! Real              :: SNDFCS(vector_length), SWEFCS(vector_length)
        ! CHARACTER(LEN=250)   :: filename        
        integer           :: ncid, dimid, varid, status

        LOGICAL           :: file_exists
        
        ! vector_length =  endIndx - startIndx + 1

        INQUIRE(FILE=trim(forc_inp_file), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(forc_inp_file) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif        
    
        status = nf90_open(trim(forc_inp_file), NF90_NOWRITE, ncid)

        status = nf90_inq_varid(ncid, "weasd", varid)
        status = nf90_get_var(ncid, varid , noahmp%swe   , &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "snwdph", varid)
        status = nf90_get_var(ncid, varid , noahmp%snow_depth  , &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "snowxy", varid)
        status = nf90_get_var(ncid, varid , noahmp%active_snow_layers  , &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "sneqvoxy", varid)
        status = nf90_get_var(ncid, varid , noahmp%swe_previous, &
            start = (/startIndx,1/), count = (/vector_length, 1/))

        status = nf90_inq_varid(ncid, "tsnoxy", varid)
        status = nf90_get_var(ncid, varid , noahmp%temperature_snow  , &
            start = (/startIndx, 1, 1/)                , &
            count = (/vector_length, 3, 1/))

        status = nf90_inq_varid(ncid, "zsnsoxy", varid)
        status = nf90_get_var(ncid, varid , noahmp%snow_soil_interface , &
            start = (/startIndx, 1, 1/)                , &
            count = (/vector_length, 7, 1/))

        status = nf90_inq_varid(ncid, "snicexy", varid)
        status = nf90_get_var(ncid, varid , noahmp%snow_ice_layer , &
            start = (/startIndx            , 1, 1/)                , &
            count = (/vector_length, 3, 1/))

        status = nf90_inq_varid(ncid, "snliqxy", varid)
        status = nf90_get_var(ncid, varid , noahmp%snow_liq_layer , &
            start = (/startIndx            , 1, 1/)                , &
            count = (/vector_length, 3, 1/))

        status = nf90_inq_varid(ncid, "stc", varid)
        status = nf90_get_var(ncid, varid , noahmp%temperature_soil , &
            start = (/startIndx, 1, 1/)                , &
            count = (/vector_length, 1, 1/))

        status = nf90_inq_varid(ncid, "snowc", varid)
        CALL NETCDF_ERR(status, 'getting varid: snowc' )
        status = nf90_get_var(ncid, varid, SCF_Grid, &
            start = (/startIndx, 1/),                &
            count = (/vector_length, 1/))
        CALL NETCDF_ERR(status, 'reading: snowc' )

        status = nf90_close(ncid)
        
        if (myrank == 0) then
            INQUIRE(FILE=trim(noda_inp_path), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(noda_inp_path) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif   
            ! Allocate(SNDnoDA(LENSFC))
            ! Allocate(SWEnoDA(LENSFC))     
        
            status = nf90_open(trim(noda_inp_path), NF90_NOWRITE, ncid)

            status = nf90_inq_varid(ncid, "weasd", varid)
            status = nf90_get_var(ncid, varid , SWEnoDA   , &
                start = (/1,1/), count = (/LENSFC, 1/))

            status = nf90_inq_varid(ncid, "snwdph", varid)
            status = nf90_get_var(ncid, varid, SNDnoDA, &
                start = (/1,1/), count = (/LENSFC, 1/))

            status = nf90_close(ncid)
        endif

    end subroutine ReadRestartNoahMP

    subroutine ReadObservation(filename, vector_length, observation)

        use netcdf

        type(observation_type) :: observation
        CHARACTER(LEN=250)     :: filename
        integer                :: vector_length
        integer                :: ncid, dimid, varid, status

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)

        status = nf90_inq_varid(ncid, "snwdph", varid)
        status = nf90_get_var(ncid, varid , observation%snow_depth  , &
            start = (/1,1/), count = (/vector_length, 1/))

        status = nf90_close(ncid)

    end subroutine ReadObservation

    subroutine ReadVectorLength(static_prefix, IDIM, vector_length)

        use netcdf

        CHARACTER(LEN=*)    :: static_prefix
        integer             :: IDIM, vector_length
        integer             :: ncid, dimid, status
        LOGICAL             :: file_exists
        CHARACTER(LEN=500)  :: filename
        Character(LEN=4)    :: grid_str

        write(grid_str, "(I0.2)") IDIM

        filename = trim(static_prefix)//"/ufs-land_C"//trim(grid_str)//"_static_fields.nc"
        

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'OPENING FILE: '//TRIM(filename) )

        status = nf90_inq_dimid(ncid, "location", dimid)
        CALL NETCDF_ERR(status, 'reading location dim id')
        status = nf90_inquire_dimension(ncid, dimid, len = vector_length)
        CALL NETCDF_ERR(status, 'reading location dim')

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing FILE: '//TRIM(filename) )

    end subroutine ReadVectorLength

    subroutine ReadNamelist(restart_dir,restart_date, static_file, fv3_restart_dir, &
        fv3_analysis_dir, IDIM, JDIM, snowUpdateOpt, READ_FV3_RESTART, WRITE_FV3_RESTART, &
        fv3_index)

        character(LEN=128)  :: static_file
        character(LEN=128)  :: init_file
        character(LEN=128)  :: forcing_dir
        character(LEN=128)  :: output_dir
        character(LEN=256)  :: fv3_restart_dir, fv3_analysis_dir
        LOGICAL             :: READ_FV3_RESTART, WRITE_FV3_RESTART, fv3_index

        logical        :: separate_output = .false.

        integer        :: timestep_seconds = -999
        integer        :: IDIM, JDIM, snowUpdateOpt

        integer        :: restart_frequency_s = 0
        logical        :: restart_simulation = .false.
        character*19   :: restart_date
        character*128  :: restart_dir

        character*19   :: simulation_start = ""
        character*19   :: simulation_end = ""

        integer        :: run_days = -999
        integer        :: run_hours = -999
        integer        :: run_minutes = -999
        integer        :: run_seconds = -999
        integer        :: run_timesteps = -999

        integer        :: begloc = 1
        integer        :: endloc = 1


        namelist / run_setup  / static_file, init_file, forcing_dir, output_dir, timestep_seconds, &
                                simulation_start, simulation_end, run_days, run_hours, run_minutes, &
                                run_seconds, run_timesteps, separate_output, begloc, endloc, &
                                restart_dir, restart_frequency_s, restart_simulation, restart_date, &
                                fv3_restart_dir, fv3_analysis_dir, IDIM, JDIM, snowUpdateOpt, &
                                READ_FV3_RESTART, WRITE_FV3_RESTART, fv3_index 

        open(30, file="ufs-land.namelist", form="formatted")
            read(30, run_setup)
        close(30)

    end subroutine ReadNamelist

    subroutine restart_from_fv3(myrank, fv3_restart_prefix, vector_length, IDIM, JDIM, &
                        tile_xy, Idim_xy, Jdim_xy, &   !
                        SWEFCS, SNDFCS)

        use netcdf
        implicit none
        include "mpif.h"
        ! type(noahmp_type) :: noahmp
        integer           :: myindx, vector_length, myrank
        real              :: SWEFCS(vector_length) , SNDFCS(vector_length)
        !Real              :: VETFCS    !, VETFCS2 
        CHARACTER(LEN=*)  :: fv3_restart_prefix
        CHARACTER(LEN=250)   :: filename
        integer           :: ncid, dimid, varid, status
        Integer           :: tile_xy , Idim_xy , &
                                Jdim_xy 
        integer           :: ID_VAR, IDIM, JDIM, ERROR, ixy   !, ID_DIM,       !LENSFC
        real        :: DUMMY1(6, IDIM, JDIM), DUMMY2(6, IDIM, JDIM), DUMMY3(6, IDIM, JDIM) !, &
                    !    DUMMY4(6, IDIM, JDIM), DUMMY5(6, IDIM, JDIM), DUMMY6(6, IDIM, JDIM)
        ! REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists
        CHARACTER(LEN=1)          :: RANKCH

        Do myindx = 1, 6 
            WRITE(RANKCH, '(I1.1)') myindx
            filename = trim(fv3_restart_prefix)//"tile"//RANKCH//".nc"
            INQUIRE(FILE=trim(filename), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(filename) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif           

            ERROR=NF90_OPEN(TRIM(filename), NF90_NOWRITE,NCID)
            CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(filename) )

            ! ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
            ! CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
            ! ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
            ! CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

            ! ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
            ! CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
            ! ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
            ! CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

            ! IF ((IDIM*JDIM) /= LENSFC) THEN
            ! PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
            ! CALL MPI_ABORT(MPI_COMM_WORLD, 88)
            ! ENDIF
            ! ALLOCATE(DUMMY(IDIM,JDIM))

            ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING sheleg ID' )
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy1(myindx,:,:))
            CALL NETCDF_ERR(ERROR, 'READING sheleg' )

            ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING snwdph ID' )
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy2(myindx,:,:))
            CALL NETCDF_ERR(ERROR, 'READING snwdph' )

            ERROR=NF90_INQ_VARID(NCID, "vtype", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'READING vtype ID' )
            ERROR=NF90_GET_VAR(NCID, ID_VAR, dummy3(myindx,:,:))
            CALL NETCDF_ERR(ERROR, 'READING vtype' )
            ! VETFCS = RESHAPE(DUMMY, (/LENSFC/))  
            
            ERROR = NF90_CLOSE(NCID)    
            CALL NETCDF_ERR(ERROR, 'closing file' )
        End do 

        ! 7.20.21 map fv3 to land
        Do ixy=1, vector_length 
            SWEFCS(ixy) = dummy1(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
            SNDFCS(ixy) = dummy2(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
            ! VETFCS(ixy) = dummy3(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy))
        end do
        ! if (myrank==0) then
        !     print*, "Vet type 1 "
        !     print*, VETFCS
        !     ! print*, "Vet type 2 "
        !     ! print*, VETFCS2
        ! endif  

    end subroutine restart_from_fv3

    subroutine restart_to_fv3(fv3_restart_prefix, IY, IM, ID, IH, IDIM, JDIM, &
                        vector_length, tile_xy, Idim_xy, Jdim_xy, SWEFCS, SNDFCS, &
                        NUM_TILES, LENSFC, SNDFfull) 

        use netcdf
        implicit none
        include "mpif.h"
        ! type(noahmp_type) :: noahmp
        integer           :: IY, IM, ID, IH, IDIM, JDIM, vector_length,NUM_TILES, LENSFC   !,myrank
        real              :: SWEFCS , SNDFCS 
        CHARACTER(LEN=*)  :: fv3_restart_prefix
        CHARACTER(LEN=250)   :: filename    
        integer           :: ncid, dimid, varid, status, myindx 
        Integer           :: tile_xy , Idim_xy , Jdim_xy 
        integer           :: ID_VAR,  ixy, ERROR   !, ID_DIM,       !LENSFC
        real              :: DUMMY1(NUM_TILES, IDIM, JDIM), DUMMY2(NUM_TILES, IDIM, JDIM)
        real              :: SNDFfull(NUM_TILES, LENSFC)
        ! REAL                      :: SLMASK(LENSFC)
        LOGICAL                   :: file_exists
        CHARACTER(LEN=1)          :: RANKCH
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str
        integer              :: dims_strt(3), dims_end(3)

        write(y_str, "(I4)") IY
        write(m_str, "(I0.2)") IM
        write(d_str, "(I0.2)") ID
        write(h_str, "(I0.2)") IH
            ! 7.20.21 map land to fv3
        DUMMY1 = IEEE_VALUE(DUMMY1, IEEE_QUIET_NAN)
        DUMMY2 = IEEE_VALUE(DUMMY2, IEEE_QUIET_NAN)
        Do myindx = 1, NUM_TILES
            dummy2(myindx,:,:) = reshape(SNDFfull(myindx,:), (/IDIM, JDIM/))
        enddo
        Do ixy=1, vector_length 
            dummy1(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy)) = SWEFCS(ixy)
            dummy2(tile_xy(ixy), Jdim_xy(ixy), Idim_xy(ixy)) = SNDFCS(ixy)
        end do
        ! print*, "snd"
        ! print*, dummy2(3,:,:)

        dims_strt(1:3) = 1
        dims_end(1) = idim
        dims_end(2) = jdim
        dims_end(3) = 1

        Do myindx = 1, NUM_TILES
            WRITE(RANKCH, '(I1.1)') myindx
            ! filename = trim(fv3_restart_prefix)//"tile"//RANKCH//".nc"
            filename = trim(fv3_restart_prefix)//"/"//TRIM(y_str)//TRIM(m_str)// &
                TRIM(d_str)//"." //TRIM(h_str)//"0000.sfcanl_data.tile"//RANKCH//".nc"
            INQUIRE(FILE=trim(filename), EXIST=file_exists)
            if (.not. file_exists) then 
                print *, 'error,file does not exist', &   
                        trim(filename) , ' exiting'
                call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
            endif           

            ERROR=NF90_OPEN(TRIM(filename), NF90_WRITE, NCID)
            CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(filename) )

            ! ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', ID_DIM)
            ! CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
            ! ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=IDIM)
            ! CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

            ! ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', ID_DIM)
            ! CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
            ! ERROR=NF90_INQUIRE_DIMENSION(NCID,ID_DIM,LEN=JDIM)
            ! CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

            ! IF ((IDIM*JDIM) /= LENSFC) THEN
            ! PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
            ! CALL MPI_ABORT(MPI_COMM_WORLD, 88)
            ! ENDIF
            ! ALLOCATE(DUMMY(IDIM,JDIM))

            ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'VarID sheleg ID' )
            ERROR=NF90_PUT_VAR(NCID, ID_VAR, dummy1(myindx,:,:), dims_strt, dims_end)
            CALL NETCDF_ERR(ERROR, 'Writing sheleg' )

            ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
            CALL NETCDF_ERR(ERROR, 'VarID snwdph ID' )
            ERROR=NF90_PUT_VAR(NCID, ID_VAR, dummy2(myindx,:,:), dims_strt, dims_end)
            CALL NETCDF_ERR(ERROR, 'Writing snwdph' )
            
            ERROR = NF90_CLOSE(NCID)    
            CALL NETCDF_ERR(ERROR, 'closing file' )
        End do 

    end subroutine restart_to_fv3    

    subroutine restart_noah_to_fv3(fv3_restart_prefix, myrank, N_sA, N_sA_Ext, mpiReal_size, &
            IY, IM, ID, IH, NUM_TILES, IDIM, JDIM, LENSFC, vector_length, LENSFC_land,  &
            tile_xy, Idim_xy, Jdim_xy, SWEANL, SNDANL, noahmp) 

        use netcdf
        implicit none
        include "mpif.h"

        
        CHARACTER(LEN=*)  :: fv3_restart_prefix 
        integer           :: myrank, N_sA, N_sA_Ext, mpiReal_size
        integer   :: IY, IM, ID, IH, NUM_TILES, IDIM, JDIM, LENSFC, vector_length, LENSFC_land  
        Integer   :: tile_xy , Idim_xy , Jdim_xy  
        real      :: SWEANL(NUM_TILES, LENSFC), SNDANL(NUM_TILES, LENSFC)
        type(noahmp_type)  :: noahmp
        ! REAL                      :: SLMASK(LENSFC)
        integer            :: myindx, ncid, ixy, iv, i
        integer            :: dest_Aoffset, dest_Aoffset_end, arLen
        Real    :: DUMMY1(vector_length) , DUMMY2(vector_length), DUMMY3(vector_length), &
                DUMMY4(vector_length) , DUMMY5(vector_length, 7), &
            DUMMY6(vector_length, 3), DUMMY7(vector_length, 3), DUMMY8(vector_length, 3)
        LOGICAL              :: file_exists
        CHARACTER(LEN=250)   :: filename   
        CHARACTER(LEN=1)     :: RANKCH
        CHARACTER(len=5)     :: y_str, m_str, d_Str, h_str
        integer              :: dims_strt(3), dims_end(3)
        integer           :: dim_x, dim_y, dim_time, IDIM_In, JDIM_In
        real              :: x_data(3), y_data(7)
        integer           :: ID_VAR, IERR, ERROR   !, ID_DIM,       !LENSFC
        real              :: DUMMY_out(IDIM, JDIM)  !, DUMMY2(NUM_TILES, IDIM, JDIM)
        real              :: DUMMY_out2(IDIM, JDIM, 7), DUMMY_out3(IDIM, JDIM, 3)
        integer           :: dims_3d(3), dims_4d(4)
        integer   :: dim_sn, dim_snso, id_sn, id_snso, id_weasd, id_sndpth_m, id_snowxy, &
                     id_sneqvoxy, id_tsnoxy, id_snicexy, id_snliqxy, id_zsnsoxy
        integer           :: header_buffer_val = 16384

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
        Endif
        if (MYRANK == 0) then
            DUMMY1(1:1+N_sA) = noahmp%swe
            DUMMY2(1:1+N_sA) = noahmp%snow_depth
            DUMMY3(1:1+N_sA) = noahmp%active_snow_layers
            DUMMY4(1:1+N_sA) = noahmp%swe_previous
            Do iv=1, 7 
                DUMMY5(1:1+N_sA, iv) = noahmp%snow_soil_interface(:,iv)
            Enddo
            Do iv=1, 3
                DUMMY6(1:1+N_sA, iv) = noahmp%temperature_snow(:,iv)
                DUMMY7(1:1+N_sA, iv) = noahmp%snow_ice_layer(:,iv)
                DUMMY8(1:1+N_sA, iv) = noahmp%snow_liq_layer(:,iv)
            Enddo
            Do ixy =  1, NUM_TILES-1   ! sender proc index within tile group
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
                
                    call MPI_RECV(DUMMY5(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                400*(ixy+1) + 10*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Do iv=1, 3
                    call MPI_RECV(DUMMY6(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                500*(ixy+1) + 20*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY7(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                600*(ixy+1) + 30*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                    call MPI_RECV(DUMMY8(dest_Aoffset:dest_Aoffset_end, iv), arLen, mpiReal_size, ixy,      &
                                700*(ixy+1) + 40*iv, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)
                Enddo
            End do

            write(y_str, "(I4)") IY
            write(m_str, "(I0.2)") IM
            write(d_str, "(I0.2)") ID
            write(h_str, "(I0.2)") IH           

            Do myindx = 1, NUM_TILES
                WRITE(RANKCH, '(I1.1)') myindx
                ! filename = trim(fv3_restart_prefix)//"tile"//RANKCH//".nc"
                filename = trim(fv3_restart_prefix)//"/"//TRIM(y_str)//TRIM(m_str)// &
                    TRIM(d_str)//"." //TRIM(h_str)//"0000.sfcanl_data.tile"//RANKCH//".nc"
                INQUIRE(FILE=trim(filename), EXIST=file_exists)
                if (.not. file_exists) then 
                    print *, 'error,file does not exist', &   
                            trim(filename) , ' exiting'
                    call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
                endif           

                ERROR=NF90_OPEN(TRIM(filename), NF90_WRITE, NCID)
                CALL NETCDF_ERR(ERROR, 'OPENING FILE: '//TRIM(filename) )

                ERROR=NF90_INQ_DIMID(NCID, 'xaxis_1', dim_x)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,dim_x,LEN=IDIM_In)
                CALL NETCDF_ERR(ERROR, 'READING xaxis_1' )

                ERROR=NF90_INQ_DIMID(NCID, 'yaxis_1', dim_y)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )
                ERROR=NF90_INQUIRE_DIMENSION(NCID,dim_y,LEN=JDIM_In)
                CALL NETCDF_ERR(ERROR, 'READING yaxis_1' )

                ERROR=NF90_INQ_DIMID(NCID, 'Time', dim_time)
                CALL NETCDF_ERR(ERROR, 'READING Time' )
                ! ERROR=NF90_INQUIRE_DIMENSION(NCID,dim_time,LEN=TDIM_In)
                ! CALL NETCDF_ERR(ERROR, 'READING Time' )

                IF ((IDIM*JDIM) /= (IDIM_In*JDIM_In)) THEN
                    PRINT*,'FATAL ERROR: DIMENSIONS WRONG.'
                    CALL MPI_ABORT(MPI_COMM_WORLD, 88)
                ENDIF

                    !--- define dimensions                
                error = nf90_redef(ncid)  !CALL NCREDF(NCID, RCODE)
                call netcdf_err(error, 'entering define mode' )

                error = nf90_def_dim(ncid, 'snow_levels', 3, dim_sn)
                call netcdf_err(error, 'DEFINING snow_levels DIMENSION' )
                error = nf90_def_dim(ncid, 'snso_levels', 7, dim_snso)
                call netcdf_err(error, 'DEFINING snso_levels DIMENSION' )
                !--- define fields
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
                
                dims_3d(1) = dim_x
                dims_3d(2) = dim_y
                dims_3d(3) = dim_time            
                error = nf90_def_var(ncid, 'weasd', NF90_DOUBLE, dims_3d, id_weasd)
                call netcdf_err(error, 'DEFINING weasd' )
                error = nf90_put_att(ncid, id_weasd, "long_name", "weasd")
                call netcdf_err(error, 'DEFINING weasd LONG NAME' )
                error = nf90_put_att(ncid, id_weasd, "units", "mm")
                call netcdf_err(error, 'DEFINING weasd UNITS' )

                error = nf90_def_var(ncid, 'snwdph_m', NF90_DOUBLE, dims_3d, id_sndpth_m)
                call netcdf_err(error, 'DEFINING sndpth_m' )
                error = nf90_put_att(ncid, id_sndpth_m, "long_name", "sndpth_m")
                call netcdf_err(error, 'DEFINING sndpth_m LONG NAME' )
                error = nf90_put_att(ncid, id_sndpth_m, "units", "mm")
                call netcdf_err(error, 'DEFINING sndpth_m UNITS' )
                
                error = nf90_def_var(ncid, 'snowxy', NF90_DOUBLE, dims_3d, id_snowxy)
                call netcdf_err(error, 'DEFINING snowxy' )
                error = nf90_put_att(ncid, id_snowxy, "long_name", "snowxy")
                call netcdf_err(error, 'DEFINING snowxy LONG NAME' )
                error = nf90_put_att(ncid, id_snowxy, "units", "mm")
                call netcdf_err(error, 'DEFINING snowxy UNITS' )

                error = nf90_def_var(ncid, 'sneqvoxy', NF90_DOUBLE, dims_3d, id_sneqvoxy)
                call netcdf_err(error, 'DEFINING sneqvoxy' )
                error = nf90_put_att(ncid, id_sneqvoxy, "long_name", "sneqvoxy")
                call netcdf_err(error, 'DEFINING sneqvoxy LONG NAME' )
                error = nf90_put_att(ncid, id_sneqvoxy, "units", "mm")
                call netcdf_err(error, 'DEFINING sneqvoxy UNITS' )

                dims_4d(1) = dim_x
                dims_4d(2) = dim_y
                dims_4d(3) = dim_sn
                dims_4d(4) = dim_time
                error = nf90_def_var(ncid, 'tsnoxy', NF90_DOUBLE, dims_4d, id_tsnoxy)
                call netcdf_err(error, 'DEFINING tsnoxy' )
                error = nf90_put_att(ncid, id_tsnoxy, "long_name", "tsnoxy")
                call netcdf_err(error, 'DEFINING tsnoxy LONG NAME' )
                error = nf90_put_att(ncid, id_tsnoxy, "units", "mm")
                call netcdf_err(error, 'DEFINING tsnoxy UNITS' )
                
                error = nf90_def_var(ncid, 'snicexy', NF90_DOUBLE, dims_4d, id_snicexy)
                call netcdf_err(error, 'DEFINING snicexy' )
                error = nf90_put_att(ncid, id_snicexy, "long_name", "snicexy")
                call netcdf_err(error, 'DEFINING snicexy LONG NAME' )
                error = nf90_put_att(ncid, id_snicexy, "units", "mm")
                call netcdf_err(error, 'DEFINING snicexy UNITS' )
                
                error = nf90_def_var(ncid, 'snliqxy', NF90_DOUBLE, dims_4d, id_snliqxy)
                call netcdf_err(error, 'DEFINING snliqxy' )
                error = nf90_put_att(ncid, id_snliqxy, "long_name", "snliqxy")
                call netcdf_err(error, 'DEFINING snliqxy LONG NAME' )
                error = nf90_put_att(ncid, id_snliqxy, "units", "mm")
                call netcdf_err(error, 'DEFINING snliqxy UNITS' )

                dims_4d(1) = dim_x
                dims_4d(2) = dim_y
                dims_4d(3) = dim_snso
                dims_4d(4) = dim_time
                error = nf90_def_var(ncid, 'zsnsoxy', NF90_DOUBLE, dims_4d, id_zsnsoxy)
                call netcdf_err(error, 'DEFINING zsnsoxy' )
                error = nf90_put_att(ncid, id_zsnsoxy, "long_name", "zsnsoxy")
                call netcdf_err(error, 'DEFINING zsnsoxy LONG NAME' )
                error = nf90_put_att(ncid, id_zsnsoxy, "units", "mm")
                call netcdf_err(error, 'DEFINING zsnsoxy UNITS' )
! End def
                error = nf90_enddef(ncid, header_buffer_val,4,0,4)
                call netcdf_err(error, 'End DEFINING HEADER' )
                
                do i = 1, 3
                x_data(i) = float(i)
                enddo
                do i = 1, 7
                y_data(i) = float(i)
                enddo
                error = nf90_put_var( ncid, id_sn, x_data)
                call netcdf_err(error, 'WRITING snow_levels RECORD' )
                error = nf90_put_var( ncid, id_snso, y_data)
                call netcdf_err(error, 'WRITING snow_soil_levels RECORD' )

                dims_strt(1:3) = 1
                dims_end(1) = idim
                dims_end(2) = jdim
                dims_end(3) = 1            
                ERROR=NF90_INQ_VARID(NCID, "sheleg", ID_VAR)
                CALL NETCDF_ERR(ERROR, 'VarID sheleg ID' )
                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)
                DUMMY_out = reshape(SWEANL(myindx,:), (/IDIM, JDIM/))            
                ERROR=NF90_PUT_VAR(NCID, ID_VAR, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing sheleg')
                
                ERROR=NF90_INQ_VARID(NCID, "snwdph", ID_VAR)
                CALL NETCDF_ERR(ERROR, 'VarID snwdph ID')
                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)
                DUMMY_out = reshape(SNDANL(myindx,:), (/IDIM, JDIM/))            
                ERROR=NF90_PUT_VAR(NCID, ID_VAR, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing snwdph')

                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)
                DUMMY_out = reshape(SWEANL(myindx,:), (/IDIM, JDIM/))            
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out(Jdim_xy(ixy), Idim_xy(ixy)) = DUMMY1(ixy)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_weasd, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing sheleg multi-layers')
                
                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)
                DUMMY_out = reshape(SNDANL(myindx,:), (/IDIM, JDIM/))            
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out(Jdim_xy(ixy), Idim_xy(ixy)) = DUMMY2(ixy)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_sndpth_m, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing snwdph multi-layers')

                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out(Jdim_xy(ixy), Idim_xy(ixy)) = DUMMY3(ixy)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_snowxy, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing snowxy multi-layers')

                DUMMY_out = IEEE_VALUE(DUMMY_out, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out(Jdim_xy(ixy), Idim_xy(ixy)) = DUMMY4(ixy)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_sneqvoxy, DUMMY_out(:,:), dims_strt, dims_end)
                CALL NETCDF_ERR(ERROR, 'Writing sneqvoxy multi-layers')

                DUMMY_out2 = IEEE_VALUE(DUMMY_out2, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out2(Jdim_xy(ixy), Idim_xy(ixy), :) = DUMMY5(ixy, :)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_zsnsoxy, DUMMY_out2(:,:,:), &
                start = (/1, 1, 1, 1/), count =(/idim, jdim, 7, 1/))
                CALL NETCDF_ERR(ERROR, 'Writing zsnsoxy multi-layers')

                DUMMY_out3 = IEEE_VALUE(DUMMY_out3, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out3(Jdim_xy(ixy), Idim_xy(ixy),:) = DUMMY6(ixy,:)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_tsnoxy, DUMMY_out3(:,:,:), &
                start = (/1, 1, 1, 1/), count = (/idim, jdim, 3, 1/))
                CALL NETCDF_ERR(ERROR, 'Writing tsnoxy multi-layers')

                DUMMY_out3 = IEEE_VALUE(DUMMY_out3, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out3(Jdim_xy(ixy), Idim_xy(ixy),:) = DUMMY7(ixy,:)
                    endif
                end do  
                ERROR=NF90_PUT_VAR(NCID, id_snicexy, DUMMY_out3(:,:,:), &
                start = (/1, 1, 1, 1/), count = (/idim, jdim, 3, 1/))
                CALL NETCDF_ERR(ERROR, 'Writing snicexy multi-layers')

                DUMMY_out3 = IEEE_VALUE(DUMMY_out3, IEEE_QUIET_NAN)           
                Do ixy=1, vector_length 
                    if (myindx == tile_xy(ixy)) then
                        DUMMY_out3(Jdim_xy(ixy), Idim_xy(ixy),:) = DUMMY8(ixy,:)
                    endif
                end do
                ERROR=NF90_PUT_VAR(NCID, id_snliqxy, DUMMY_out3(:,:,:), &
                start = (/1, 1, 1, 1/), count = (/idim, jdim, 3, 1/))
                CALL NETCDF_ERR(ERROR, 'Writing snliqxy multi-layers')

                
                ERROR = NF90_CLOSE(NCID)    
                CALL NETCDF_ERR(ERROR, 'closing file' )
            End do 

        endif


    end subroutine restart_noah_to_fv3

    subroutine ReadTileInfo(filename, LENSFC, veg_type_landice, mp_start, mp_end, &
                            tile_xy, Idim_xy, Jdim_xy, RLA, RLO, OROG, VETFCS, LANDMASK) 

        use netcdf

        ! type(noahmp_type) :: noahmp
        CHARACTER(LEN=*)  :: filename
        integer           :: LENSFC, vector_length, veg_type_landice, mp_start, mp_end
        integer           :: ncid, dimid, varid, status, IDIM
        real              :: RLA(LENSFC), RLO(LENSFC), OROG(LENSFC), VETFCS(LENSFC)
        Integer           :: tile_xy(LENSFC), Idim_xy(LENSFC), Jdim_xy(LENSFC), &
                             LANDMASK(LENSFC)
        LOGICAL             :: file_exists
        ! CHARACTER(LEN=500)  :: filename
        ! Character(LEN=4)    :: grid_str
        Integer             :: Dummy(LENSFC), i
        
        vector_length = mp_end - mp_start + 1 
        if (vector_length /= LENSFC) then
            print*, "error, vector dimentions vector_length and LENSFC_land don't match"
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif
    
        INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if (.not. file_exists) then 
            print *, 'error,file does not exist', &   
                    trim(filename) , ' exiting'
            call MPI_ABORT(MPI_COMM_WORLD, 10) ! CSD - add proper error trapping?
        endif           

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        CALL NETCDF_ERR(status, 'opening file '//trim(filename))

        status = nf90_inq_varid(ncid, "cube_tile", varid)
        ! if(status /= nf90_noerr) call handle_err(status)
        CALL NETCDF_ERR(status, 'reading cube_tile var id' )
        status = nf90_get_var(ncid, varid, tile_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube tile value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_i", varid)
        CALL NETCDF_ERR(status, 'reading cube i vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Idim_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube i value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "cube_j", varid)
        CALL NETCDF_ERR(status, 'reading cube j vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Jdim_xy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading cube j value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "vegetation_category", varid)
        CALL NETCDF_ERR(status, 'reading vegetation_category vari id' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, Dummy, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading vegetation_category value' )
        ! if(status /= nf90_noerr) call handle_err(status)
        VETFCS = float(Dummy)

        status = nf90_inq_varid(ncid, "latitude", varid)
        CALL NETCDF_ERR(status, 'reading latitude var ID' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)
        status = nf90_get_var(ncid, varid, RLA, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading latitude value' )
        ! if(status /= nf90_noerr) call NETCDF_ERR(status)

        status = nf90_inq_varid(ncid, "longitude", varid)
        CALL NETCDF_ERR(status, 'reading longitude var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, RLO, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading longitude value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "elevation", varid)
        CALL NETCDF_ERR(status, 'reading elevation var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, OROG, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading elevation value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        status = nf90_inq_varid(ncid, "land_mask", varid)
        CALL NETCDF_ERR(status, 'reading land_mask var ID' )
        ! if(status /= nf90_noerr) call handle_err(status)
        status = nf90_get_var(ncid, varid, LANDMASK, &
            start = (/mp_start/), count = (/vector_length/))
        CALL NETCDF_ERR(status, 'reading land_mask value' )
        ! if(status /= nf90_noerr) call handle_err(status)

        do i = 1, LENSFC 
            ! if land, but not land ice, set mask to 1.
            if (NINT(VETFCS(i)) /=  veg_type_landice) then 
                 LANDMASK(i) = 1 
            else 
                 LANDMASK(i) = 0
            endif 
         enddo

        status = nf90_close(ncid)
        CALL NETCDF_ERR(status, 'closing file' )

    end subroutine ReadTileInfo

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

End MODULE M_UFSLAND_SNOW_UPDATE
