! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module terrain_core

   use mpas_framework
   use mpas_timekeeping
   use mpas_stream_manager


   type (MPAS_Clock_type), pointer :: clock


   contains


   !***********************************************************************
   !
   !  function terrain_core_init
   !
   !> \brief   MPAS Core Initialization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs the necessary initialization of a core. This can
   !>  involve things like setting up coeffiecients for remapping and advection,
   !>  but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   function terrain_core_init(domain, startTimeStamp) result(iErr)
   
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      type (mpas_pool_type), pointer :: modelPool
      character (len=StrKIND), pointer :: xtime
      integer :: iErr


      iErr = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, iErr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

      call mpas_pool_get_subpool(domain % blocklist % structs, 'model', modelPool)
      call mpas_pool_get_array(modelPool, 'xtime', xtime)
      xtime = startTimeStamp

      call mpas_stream_mgr_read(domain % streamManager, ierr=iErr)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=iErr)

   end function terrain_core_init


   !***********************************************************************
   !
   !  function terrain_core_run
   !
   !> \brief   MPAS Core Run Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   function terrain_core_run(domain) result(iErr)
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_vector_operations
      use mpas_geometry_utils
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: mesh

      real (kind=RKIND), dimension(:), pointer :: ter

      integer, pointer :: nCells
      integer :: iErr
      integer :: my_id
      integer :: nCell

      iErr = 0
   
      !! Terrain Stuff

      block_ptr => domain % blocklist

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)

      call mpas_pool_get_array(mesh, 'ter', ter)

      my_id = domain % dminfo % my_proc_id    
      do nCell = 1, nCells, 1
        ter(nCell) = my_id
      enddo

      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)

   end function terrain_core_run


   !***********************************************************************
   !
   !  function terrain_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   function terrain_core_finalize(domain) result(iErr)
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function terrain_core_finalize

end module terrain_core
