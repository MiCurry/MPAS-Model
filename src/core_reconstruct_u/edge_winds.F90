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
