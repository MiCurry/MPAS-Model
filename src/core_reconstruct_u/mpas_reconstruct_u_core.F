! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module reconstruct_u_core

   use mpas_framework
   use mpas_timekeeping
   use mpas_stream_manager


   type (MPAS_Clock_type), pointer :: clock


   contains

    !-----------------------------------------------------------------------
    !  routine cam_mpas_cell_to_edge_winds
    !
    !> \brief  Projects cell-centered winds to the normal component of velocity on edges
    !> \author Michael Duda
    !> \date   16 January 2020
    !> \details
    !>  Given zonal and meridional winds at cell centers, unit vectors in the east
    !>  and north directions at cell centers, and unit vectors in the normal
    !>  direction at edges, this routine projects the cell-centered winds onto
    !>  the normal vectors.
    !>
    !>  Prior to calling this routine, the halos for the zonal and meridional
    !>  components of cell-centered winds should be updated. It is also critical
    !>  that the east, north, uZonal, and uMerid field are all allocated with
    !>  a "garbage" element; this is handled automatically for fields allocated
    !>  by the MPAS infrastructure.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_cell_to_edge_winds(nEdges, uZonal, uMerid, east, north, edgeNormalVectors, &
                                           cellsOnEdge, uNormal)

       use mpas_kind_types, only : RKIND

       implicit none

       integer, intent(in) :: nEdges
       real(kind=RKIND), dimension(:,:), intent(in) :: uZonal, uMerid
       real(kind=RKIND), dimension(:,:), intent(in) :: east, north, edgeNormalVectors
       integer, dimension(:,:), intent(in) :: cellsOnEdge
       real(kind=RKIND), dimension(:,:), intent(out) :: uNormal

       integer :: iEdge, cell1, cell2


       do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          uNormal(:,iEdge) =  uZonal(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                            + uMerid(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   &
                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell1)   &
                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell1))  &
                            + uZonal(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                            + uMerid(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   &
                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell2)   &
                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell2))
       end do

    end subroutine cam_mpas_cell_to_edge_winds

    !-----------------------------------------------------------------------
    !  routine cam_mpas_compute_unit_vectors
    !
    !> \brief  Computes local unit north, east, and edge-normal vectors
    !> \author Michael Duda
    !> \date   15 January 2020
    !> \details
    !>  This routine computes the local unit north and east vectors at all cell
    !>  centers, storing the resulting fields in the mesh pool as 'north' and
    !>  'east'. It also computes the edge-normal unit vectors by calling
    !>  the mpas_initialize_vectors routine. Before this routine is called,
    !>  the mesh pool must contain 'latCell' and 'lonCell' fields that are valid
    !>  for all cells (not just solve cells), plus any fields that are required
    !>  by the mpas_initialize_vectors routine.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_compute_unit_vectors(domain_ptr)

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_derived_types, only : mpas_pool_type
       use mpas_kind_types, only : RKIND
       use mpas_vector_operations, only : mpas_initialize_vectors

       implicit none

       type (domain_type), intent(inout) :: domain_ptr
       type (mpas_pool_type), pointer :: meshPool
       real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real(kind=RKIND), dimension(:,:), pointer :: east, north
       integer, pointer :: nCells
       integer :: iCell

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'east', east)
       call mpas_pool_get_array(meshPool, 'north', north)

       do iCell = 1, nCells

          east(1,iCell) = -sin(lonCell(iCell))
          east(2,iCell) =  cos(lonCell(iCell))
          east(3,iCell) =  0.0

          ! Normalize
          east(1:3,iCell) = east(1:3,iCell) / sqrt(sum(east(1:3,iCell) * east(1:3,iCell)))

          north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
          north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
          north(3,iCell) =  cos(latCell(iCell))

          ! Normalize
          north(1:3,iCell) = north(1:3,iCell) / sqrt(sum(north(1:3,iCell) * north(1:3,iCell)))

       end do

       call mpas_initialize_vectors(meshPool)

    end subroutine cam_mpas_compute_unit_vectors

   !***********************************************************************
   !
   !  function reconstruct_u_core_init
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
   function reconstruct_u_core_init(domain, startTimeStamp) result(iErr)
   
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

   end function reconstruct_u_core_init


   !***********************************************************************
   !
   !  function reconstruct_u_core_run
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
   function reconstruct_u_core_run(domain) result(iErr)
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_vector_operations
      use mpas_geometry_utils
      use mpas_pool_routines
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCells
      integer :: iErr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp

      real (kind=RKIND), dimension(:,:), pointer :: east, north

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'east', east)
      call mpas_pool_get_array(meshPool, 'north', north)

      iErr = 0

      call cam_mpas_compute_unit_vectors(domain)


!      do while (.not. mpas_is_clock_stop_time(clock))
!         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
!         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
!
!         call mpas_log_write('Current time = '//trim(timeStamp))
!
!         call mpas_advance_clock(clock)
!      end do

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
      call mpas_log_write('Final time = '//trim(timeStamp))

      call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true.)

   end function reconstruct_u_core_run


   !***********************************************************************
   !
   !  function reconstruct_u_core_finalize
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
   function reconstruct_u_core_finalize(domain) result(iErr)
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function reconstruct_u_core_finalize

end module reconstruct_u_core
