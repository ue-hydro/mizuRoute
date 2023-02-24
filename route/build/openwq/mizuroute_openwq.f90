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
  public :: openwq_run_space_step_basin_in
  public :: openwq_run_space_step
  public :: openwq_run_time_end

  ! Global Data for prognostic Variables of HRUs
  !type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes

  contains

! ####################
! OpenWQ: openwq init
! ####################
subroutine openwq_init(err, message)

  USE globalData, ONLY:openwq_obj
  USE globalData, ONLY: nRch             ! number of reaches in the whoel river network

  implicit none

  integer(i4b),intent(inout)                      :: err
  character(*),intent(inout)                      :: message         ! error messgage

  openwq_obj = CLASSWQ_openwq() ! initalize openWQ object

  ! intialize openWQ
  err=openwq_obj%decl(    &
    nRch)
  
end subroutine openwq_init


! ####################
! OpenWQ: run_time_start
! ####################
subroutine openwq_run_time_start(  &
      openwq_obj)

      USE globalData,        ONLY: modTime        ! previous and current model time
      USE globalData,       ONLY: nRch             ! number of reaches in the whoel river network
      USE globalData,       ONLY : RCHFLX

      implicit none

      ! Dummy Varialbes
      class(CLASSWQ_openwq), intent(in)   :: openwq_obj
      integer(i4b)                         :: iRch
      integer(i4b)                         :: simtime(6) ! 5 time values yy-mm-dd-hh-min
      real(dp)                             :: REACH_VOL_0(nRch)
      integer(i4b)                        :: err

      ! Getting reach volume to update openwq:waterVol_hydromodel
      do iRch = 1, nRch
      REACH_VOL_0(iRch) = RCHFLX(1,iRch)%ROUTE(1)%REACH_VOL(0)
      end do

      ! add the time values to the array
      simtime(1) = modTime(1)%year()
      simtime(2) = modTime(1)%month()
      simtime(3) = modTime(1)%day()
      simtime(4) = modTime(1)%hour()
      simtime(5) = modTime(1)%minute()
      simtime(6) = modTime(1)%sec()

      err=openwq_obj%openwq_run_time_start(     &
      simtime,                                  &
      nRch,                                     &
      REACH_VOL_0)
  
end subroutine openwq_run_time_start

!  OpenWQ space basin/summa
subroutine openwq_run_space_step_basin_in()

      USE globalData,   only: openwq_obj
      USE globalData,       ONLY: modTime        ! previous and current model time
      USE globalData,       ONLY : nRch
      USE globalData,       only : NETOPO
      USE globalData,       only : RCHFLX
      USE globalData,       only : TSEC

      implicit none

      integer(i4b)                           :: iRch      ! variable needed for looping through reaches
      integer(i4b)                           :: simtime(6) ! 5 time values yy-mm-dd-hh-min
      real(dp)                               :: mizuroute_timestep
      integer(i4b)                           :: err
      integer(i4b)                           :: river_network_reaches    = 0
      integer(i4b)                           :: index_s_openwq
      integer(i4b)                           :: index_r_openwq
      integer(i4b)                           :: ix_s_openwq
      integer(i4b)                           :: ix_r_openwq
      integer(i4b)                           :: iy_s_openwq
      integer(i4b)                           :: iy_r_openwq
      integer(i4b)                           :: iz_s_openwq
      integer(i4b)                           :: iz_r_openwq
      real(dp)                               :: wmass_source_openwq
      real(dp)                               :: compt_vol_m3
      real(dp)                               :: flux_m3_sec
      real(dp)                               :: flux_m3_timestep

      ! Getting time
      simtime(1) = modTime(1)%year()
      simtime(2) = modTime(1)%month()
      simtime(3) = modTime(1)%day()
      simtime(4) = modTime(1)%hour()
      simtime(5) = modTime(1)%minute()
      simtime(6) = modTime(1)%sec()

      iy_r_openwq = 1
      iz_r_openwq = 1 

      ! Time step
      mizuroute_timestep = TSEC(1) - TSEC(0)  

      ! ####################################################################
      ! Apply Fluxes
      ! Call RunSpaceStep
      ! ####################################################################

      ! ====================================================
      ! 1 Basin routing: SUMMA to MizuRoute
      ! ====================================================
 
      do iRch = 1, nRch 
            ! *Source*:
            ! PRECIP (external flux, so need call openwq_run_space_in) 
            ! *Recipient*: canopy (only 1 z layer)
            index_r_openwq = river_network_reaches
            ix_r_openwq          = iRch 
            ! *Flux*: the portion of rainfall and snowfall not throughfall
            flux_m3_sec = RCHFLX(1,iRch)%BASIN_QR(0) 
            flux_m3_timestep = flux_m3_sec * mizuroute_timestep
            ! *Call openwq_run_space_in* if wflux_s2r not 0
            err=openwq_obj%openwq_run_space_in(          &
            simtime,                                     &
            'SUMMA_RUNOFF',                              &
            index_r_openwq, ix_r_openwq, iy_r_openwq, iz_r_openwq,  &
            flux_m3_timestep)
      end do

      
end subroutine openwq_run_space_step_basin_in

! OpenWQ space (wihtin mizuroute)
subroutine openwq_run_space_step(segIndex,                            & ! index
                                    REACH_VOL_segIndex,               & ! Volume
                                    Qlocal_in,                        & ! flow in
                                    Qlocal_out)                        ! flow out

      USE globalData,   only: openwq_obj
      USE globalData,       ONLY: modTime        ! previous and current model time
      USE globalData,       ONLY : nRch
      USE globalData,       only : NETOPO
      USE globalData,       only : RCHFLX
      USE globalData,       only : TSEC

      implicit none

      integer(i4b)                           :: iRch      ! variable needed for looping through reaches
      integer(i4b)                           :: simtime(6) ! 5 time values yy-mm-dd-hh-min
      real(dp)                               :: mizuroute_timestep
      integer(i4b)                           :: err

      integer(i4b) :: segIndex                  ! index
      real(dp) :: REACH_VOL_segIndex            ! Volume
      real(dp) :: Qlocal_in                     ! flow in
      real(dp) :: Qlocal_out                    ! flow out

      ! compartment indexes in OpenWQ (defined in the hydrolink)
      integer(i4b)                           :: river_network_reaches    = 0
      integer(i4b)                           :: index_s_openwq
      integer(i4b)                           :: index_r_openwq
      integer(i4b)                           :: ix_s_openwq
      integer(i4b)                           :: ix_r_openwq
      integer(i4b)                           :: iy_s_openwq
      integer(i4b)                           :: iy_r_openwq
      integer(i4b)                           :: iz_s_openwq
      integer(i4b)                           :: iz_r_openwq
      real(dp)                               :: wflux_s2r_openwq
      real(dp)                               :: wmass_source_openwq
      real(dp)                               :: compt_vol_m3
      real(dp)                               :: flux_m3_sec
      real(dp)                               :: flux_m3_timestep


      simtime(1) = modTime(1)%year()
      simtime(2) = modTime(1)%month()
      simtime(3) = modTime(1)%day()
      simtime(4) = modTime(1)%hour()
      simtime(5) = modTime(1)%minute()
      simtime(6) = modTime(1)%sec()

      ! Mizuroute does not have a y-direction, 
      ! so the dimension will always be 1
      iy_s_openwq = 1
      iy_r_openwq = 1
      iz_s_openwq = 1
      iz_r_openwq = 1 

      ! Time step
      mizuroute_timestep = TSEC(1) - TSEC(0)  

      ! ####################################################################
      ! Apply Fluxes
      ! Call RunSpaceStep
      ! ####################################################################

      ! ====================================================
      ! 1 River routing
      ! ====================================================
      ! segIndex to segIndex+1
      index_s_openwq = river_network_reaches
      ix_s_openwq          = segIndex
      compt_vol_m3         = REACH_VOL_segIndex + Qlocal_in ! That's what is received previous iteraction
      wmass_source_openwq  = compt_vol_m3
      ! *Recipient*: 
      index_r_openwq       = river_network_reaches
      ix_r_openwq          = NETOPO(segIndex)%DREACHI
      !ix_r_openwq          = segIndex + 1
      if(ix_r_openwq.eq.-1) return
      ! flux
      flux_m3_timestep = Qlocal_out
      wflux_s2r_openwq = flux_m3_timestep
      ! *Call openwq_run_space* if wflux_s2r_openwq not 0
      err=openwq_obj%openwq_run_space(                          &
      simtime,                                                &
      index_s_openwq, ix_s_openwq, iy_s_openwq, iz_s_openwq,  &
      index_r_openwq, ix_r_openwq, iy_r_openwq, iz_r_openwq,  &
      wflux_s2r_openwq,                                       &
      wmass_source_openwq)

end subroutine openwq_run_space_step


! ####################
! OpenWQ: run_time_end
! ####################
subroutine openwq_run_time_end( &
  openwq_obj)

      USE globalData,        ONLY: modTime        ! previous and current model time

      implicit none

      ! Dummy Varialbes
      class(CLASSWQ_openwq), intent(in)  :: openwq_obj

      ! Local Variables
      integer(i4b)                       :: simtime(6) ! 5 time values yy-mm-dd-hh-min
      integer(i4b)                       :: err ! error control

      ! Get time
      simtime(1) = modTime(1)%year()
      simtime(2) = modTime(1)%month()
      simtime(3) = modTime(1)%day()
      simtime(4) = modTime(1)%hour()
      simtime(5) = modTime(1)%minute()
      simtime(6) = modTime(1)%sec()

      err=openwq_obj%openwq_run_time_end(simtime)           ! minute


end subroutine openwq_run_time_end

end module mizuroute_openwq