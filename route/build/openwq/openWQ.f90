module openwq
   
 USE, intrinsic :: iso_c_binding
 USE nrtype
 private
 public :: ClassWQ_OpenWQ

 include "openWQInterface.f90"

 type ClassWQ_OpenWQ
    private
    type(c_ptr) :: ptr ! pointer to openWQ class

 contains
   !  procedure :: get_num => openWQ_get_num
    procedure :: decl => openWQ_init
    !procedure :: run_time_start => openWQ_run_time_start
    !procedure :: run_space => openWQ_run_space
    !procedure :: run_space_in => openWQ_run_space_in
    !procedure :: run_time_end => openWQ_run_time_end

 end type

 interface ClassWQ_OpenWQ
    procedure create_openwq
 end interface
 contains
    function create_openwq()
        implicit none
        type(ClassWQ_OpenWQ) :: create_openwq
        create_openwq%ptr = create_openwq_c()
    end function

    ! supposed to be decl but needed to openWQ_decl in the interface file
    ! returns integer of either a failure(-1) or success(0)
    integer function openWQ_init(this)                 ! num of layers in y-dir (set to 1 because not used in summa)
   !integer function openWQ_init( &
   !   this,                      & ! openwq object
   !   num_hru,                   & ! num HRU
   !   nCanopy_2openwq,           & ! num layers of canopy (fixed to 1)
   !   nSnow_2openwq,             & ! num layers of snow (fixed to max of 5 because it varies)
   !   nSoil_2openwq,             & ! num layers of snoil (variable)
   !   nRunoff_2openwq,           & ! num layers of runoff (fixed to 1)
   !   nAquifer_2openwq,          & ! num layers of aquifer (fixed to 1)
   !   nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)
      
      implicit none
      class(ClassWQ_OpenWQ) :: this
      !integer(i4b), intent(in) :: num_hru
      !integer(i4b), intent(in) :: nCanopy_2openwq
      !integer(i4b), intent(in) :: nSnow_2openwq
      !integer(i4b), intent(in) :: nSoil_2openwq
      !integer(i4b), intent(in) :: nRunoff_2openwq
      !integer(i4b), intent(in) :: nAquifer_2openwq
      
      !integer(i4b), intent(in) :: nYdirec_2openwq

      openWQ_init = openwq_decl_c(this%ptr)                 ! num of layers in y-dir (set to 1 because not used in summa)
      !openWQ_init = openwq_decl_c(  &
      !   this%ptr,                  & ! openwq object
      !   num_hru,                   & ! num HRU
      !   nCanopy_2openwq,         & ! num layers of canopy (fixed to 1)
      !   nSnow_2openwq,           & ! num layers of snow (fixed to max of 5 because it varies)
      !   nSoil_2openwq,           & ! num layers of snoil (variable)
      !   nRunoff_2openwq,         & ! num layers of runoff (fixed to 1)
      !   nAquifer_2openwq,        & ! num layers of aquifer (fixed to 1)
      !   nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)

    end function

end module openwq