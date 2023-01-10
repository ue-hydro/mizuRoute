module mizuroute_openwq

  USE nrtype
  USE openWQ, only:ClassWQ_OpenWQ
  !USE data_types, only:gru_hru_doubleVec
  implicit none
  private
  ! Subroutines
  public :: init_openwq
  !public :: run_time_start
  !public :: run_time_start_go
  !public :: run_space_step
  !public :: run_time_end

  ! Global Data for prognostic Variables of HRUs
  !type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes

  contains

  ! Subroutine to initalize the openWQ object
  ! putting it here to keep the SUMMA_Driver clean
subroutine init_openwq(err, message)

  USE globalData,only:openWQ_obj
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

  openwq_obj = ClassWQ_OpenWQ() ! initalize openWQ object

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

end subroutine init_openwq

end module mizuroute_openwq