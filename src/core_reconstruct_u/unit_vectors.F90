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
    subroutine cam_mpas_compute_unit_vectors()

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_derived_types, only : mpas_pool_type
       use mpas_kind_types, only : RKIND
       use mpas_vector_operations, only : mpas_initialize_vectors

       implicit none

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
