module mizuroute_openwq

  USE nrtype
  USE openWQ, only:CLASSWQ_openwq
  !USE data_types, only:gru_hru_doubleVec
  implicit none
  private
  ! Subroutines
  public :: openwq_init
  public :: openwq_run_time_start
  !public :: openwq_run_time_start_go
  !public :: openwq_run_space_step
  !public :: openwq_run_time_end

  ! Global Data for prognostic Variables of HRUs
  !type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes

  contains

  ! Subroutine to initalize the openWQ object
  ! putting it here to keep the SUMMA_Driver clean
subroutine openwq_init(err, message)

  USE globalData, ONLY:openwq_obj
  USE globalData, ONLY: nRch             ! number of reaches in the whoel river network
  !USE globalData,only:gru_struc                               ! gru-hru mapping structures
  !USE globalData,only:prog_meta
  !USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures

  implicit none

  integer(i4b),intent(inout)                      :: err
  character(*),intent(inout)                      :: message         ! error messgage
  !integer(i4b)                                    :: hruCount
  !integer(i4b)                                    :: nCanopy_2openwq
  !integer(i4b)                                    :: nSnow_2openwq
  !integer(i4b)                                    :: nSoil_2openwq
  !integer(i4b)                                    :: nRunoff_2openwq
  !integer(i4b)                                    :: nAquifer_2openwq
  !integer(i4b)                                    :: nYdirec_2openwq     ! number of layers in the y-dir (not used in summa)
  !integer(i4b)                                    :: iGRU, iHRU          ! indices of GRUs and HRUs

  openwq_obj = CLASSWQ_openwq() ! initalize openWQ object

  ! nx -> num of HRUs)
  !hruCount = sum( gru_struc(:)%hruCount )

  ! ny -> this seems to be fixes because SUMMA is based on the HRU concept, so grids are serialized)
  !nYdirec_2openwq = 1

  ! Openwq nz (number of layers)
  !nCanopy_2openwq = 1       ! Canopy has only 1 layer
  !nRunoff_2openwq = 1       ! Runoff has only 1 layer (not a summa variable - openWQ keeps track of this)
  !nAquifer_2openwq = 1      ! GW has only 1 layer
  !nSoil_2openwq = 0         ! Soil may have multiple layers, and gru-hrus may have different values
  !nSnow_2openwq = 0         ! Snow has multiple layers, and gru-hrus may have different values (up to 5 layers)
  !do iGRU = 1, size(gru_struc(:))
  !  do iHRU = 1, gru_struc(iGRU)%hruCount
  !    nSoil_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, nSoil_2openwq )
  !    nSnow_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSnow, nSnow_2openwq )
  !  enddo
  !enddo

  ! intialize openWQ
  err=openwq_obj%decl(    &
    nRch)
  !err=openwq_obj%decl(    &
  !  hruCount,             & ! num HRU
  !  nCanopy_2openwq,      & ! num layers of canopy (fixed to 1)
  !  nSnow_2openwq,        & ! num layers of snow (fixed to max of 5 because it varies)
  !  nSoil_2openwq,        & ! num layers of snoil (variable)
  !  nRunoff_2openwq,      & ! num layers of runoff (fixed to 1)
  !  nAquifer_2openwq,     & ! num layers of aquifer (fixed to 1)
  !  nYdirec_2openwq)             ! num of layers in y-dir (set to 1 because not used in summa)
  
  ! Create copy of state information, needed for passing to openWQ with fluxes that require
  ! the previous time_steps volume
  !call allocGlobal(prog_meta, progStruct_timestep_start, err, message) 

end subroutine openwq_init

subroutine openwq_run_time_start(  &
  openwq_obj)

  !USE summa_type, only: summa1_type_dec            ! master summa data type
  !USE globalData, only: gru_struc
  !USE var_lookup, only: iLookPROG  ! named variables for state variables
  !USE var_lookup, only: iLookTIME  ! named variables for time data structure
  !USE var_lookup, only: iLookATTR  ! named variables for real valued attribute data structure
  !USE multiconst,only:&
  !                      iden_ice,       & ! intrinsic density of ice             (kg m-3)
  !                      iden_water        ! intrinsic density of liquid water    (kg m-3)

  implicit none

  ! Dummy Varialbes
  class(CLASSWQ_openwq), intent(in)   :: openwq_obj
  !type(summa1_type_dec), intent(in)   :: summa1_struc
  ! local variables
  !integer(i4b), intent(in)            :: nSnow_2openwq
  !integer(i4b), intent(in)            :: nSoil_2openwq
  !integer(i4b)                        :: iGRU
  !integer(i4b)                        :: iHRU
  !integer(i4b)                        :: ilay
  !integer(i4b)                        :: iVar
  !integer(i4b)                        :: iDat
  !integer(i4b)                        :: openWQArrayIndex
  !integer(i4b)                        :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  !real(rkind)                         :: airTemp_depVar_summa_K
  !real(rkind)                         :: canopyWatVol_stateVar_summa_m3
  !real(rkind)                         :: aquiferWatVol_stateVar_summa_m3
  !real(rkind)                         :: sweWatVol_stateVar_summa_m3(nSnow_2openwq)
  !real(rkind)                         :: soilWatVol_stateVar_summa_m3(nSoil_2openwq)
  !real(rkind)                         :: soilTemp_depVar_summa_K(nSoil_2openwq)
  !real(rkind)                         :: soilMoist_depVar_summa_frac(nSoil_2openwq)
  !integer(i4b)                        :: soil_start_index ! starting value of the soil in the mLayerVolFracWat(:) array
  integer(i4b)                        :: err
  !real(rkind),parameter               :: valueMissing=-9999._rkind   ! seems to be SUMMA's default value for missing data
  !logical(1)                          :: last_hru_flag

  !summaVars: associate(&
  !    progStruct     => summa1_struc%progStruct             , &
  !    timeStruct     => summa1_struc%timeStruct             , &
  !    attrStruct     => summa1_struc%attrStruct               &
  !)

  !last_hru_flag = .false.
  ! Update dependencies and storage volumes
  ! Assemble the data to send to openWQ

  !openWQArrayIndex = 0 ! index into the arrays that are being passed to openWQ

  !do iGRU = 1, size(gru_struc(:))
  !    do iHRU = 1, gru_struc(iGRU)%hruCount

  !      if (iGRU == size(gru_struc) .and. iHRU == gru_struc(iGRU)%hruCount)then
  !        last_hru_flag = .true.
  !      end if

  !      GeneralVars: associate(&
  !        hru_area_m2                 => attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)                      ,&
  !        Tair_summa_K                => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanairTemp)%dat(1)      ,&
  !        scalarCanopyWat_summa_kg_m2 => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)        ,&
  !        AquiferStorWat_summa_m      => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1)   &
  !      )


        ! ############################
        ! Update unlayered variables and dependencies 
        ! (1 layer only)
        ! ############################

        ! Tair 
        ! (Summa in K)
  !      if(Tair_summa_K /= valueMissing) then
  !        airTemp_depVar_summa_K =  Tair_summa_K 
  !      endif
          
        ! Vegetation
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarCanopyWat [kg m-2], so needs to  to multiply by hru area [m2] and divide by water density
  !      if(scalarCanopyWat_summa_kg_m2 /= valueMissing) then
  !        canopyWatVol_stateVar_summa_m3 = scalarCanopyWat_summa_kg_m2 * hru_area_m2 / iden_water
  !      endif

        ! Aquifer
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarAquiferStorage [m], so needs to  to multiply by hru area [m2] only
  !      if(AquiferStorWat_summa_m /= valueMissing) then
  !        aquiferWatVol_stateVar_summa_m3 = AquiferStorWat_summa_m * hru_area_m2
  !      endif 
        
        ! ############################
        ! Update layered variables and dependenecies
        ! ############################


        ! Snow
  !      if (nSnow_2openwq .gt. 0)then
  !        do ilay = 1, nSnow_2openwq
  !          SnowVars: associate(&
  !            mLayerDepth_summa_m         => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(ilay)         ,&    ! depth of each layer (m)
  !            mLayerVolFracWat_summa_frac => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(ilay)     &    ! volumetric fraction of total water in each layer  (-)
  !          )
            ! Snow
            ! unit for volume = m3 (summa-to-openwq unit conversions needed)
            ! mLayerVolFracIce and mLayerVolFracLiq [-], so needs to  to multiply by hru area [m2] and divide by water density
            ! But needs to account for both ice and liquid, and convert to liquid volumes
  !          if(mLayerVolFracWat_summa_frac /= valueMissing) then

  !            sweWatVol_stateVar_summa_m3(ilay) =                             &
  !              mLayerVolFracWat_summa_frac  * mLayerDepth_summa_m * hru_area_m2
  !          else
  !            sweWatVol_stateVar_summa_m3(ilay) = 0._rkind
  !          endif

  !          end associate SnowVars
  !        enddo
  !      end if

  !      if (nSnow_2openwq == 0)then
  !        soil_start_index = 1
  !      else
  !        soil_start_index = nSnow_2openwq
  !      end if

        ! Soil - needs to start after the snow
  !      do ilay = soil_start_index, nSoil_2openwq + nSnow_2openwq
          
  !        SoilVars: associate(&
  !          mLayerDepth_summa_m         => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(ilay)       ,&    ! depth of each layer (m)
  !          Tsoil_summa_K               => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat(ilay)        ,&
  !          mLayerVolFracWat_summa_frac => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(ilay)   &
  !        )
          ! Tsoil
          ! (Summa in K)
  !        if(Tsoil_summa_K /= valueMissing) then
  !          soilTemp_depVar_summa_K(ilay) = Tsoil_summa_K
  !        endif

  !        soilMoist_depVar_summa_frac(ilay) = 0     ! TODO: Find the value for this varaibles

          ! Soil
          ! unit for volume = m3 (summa-to-openwq unit conversions needed)
          ! mLayerMatricHead [m], so needs to  to multiply by hru area [m2]
  !        if(mLayerVolFracWat_summa_frac /= valueMissing) then
  !          soilWatVol_stateVar_summa_m3(ilay) = mLayerVolFracWat_summa_frac * hru_area_m2 * mLayerDepth_summa_m
  !        endif

  !        end associate SoilVars

  !      enddo

        ! Copy the prog structure
  !      do iVar = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var)
  !        do iDat = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat)
  !          progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:) = progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat(:)
  !        end do
  !      end do

        
          ! add the time values to the array
  !      simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  !      simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  !      simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  !      simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  !      simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute
  
        err=openwq_obj%openwq_run_time_start()
  !      err=openwq_obj%openwq_run_time_start(&
  !            last_hru_flag,                          & 
  !            openWQArrayIndex,                       & ! total HRUs
  !            nSnow_2openwq,                          &
  !            nSoil_2openwq,                          &
  !            simtime,                                &
  !            soilMoist_depVar_summa_frac,            &                    
  !            soilTemp_depVar_summa_K,                &
  !            airTemp_depVar_summa_K,                 &
  !            sweWatVol_stateVar_summa_m3,            &
  !            canopyWatVol_stateVar_summa_m3,         &
  !            soilWatVol_stateVar_summa_m3,           &
  !            aquiferWatVol_stateVar_summa_m3)
        
  !      openWQArrayIndex = openWQArrayIndex + 1 
  !      end associate GeneralVars

  !    end do ! end HRU
  !end do ! end GRU

  !end associate summaVars

end subroutine openwq_run_time_start

end module mizuroute_openwq