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

   use iso_c_binding, only : c_float, c_int, c_long, c_loc, c_ptr
   
   interface
      subroutine timer_start(timer_id) bind(C)
         use iso_c_binding, only : c_int
         integer (c_int), intent(in), value :: timer_id
      end subroutine timer_start

      subroutine timer_stop(timer_id, sec, nsec) bind(C)
         use iso_c_binding, only : c_int
         integer (c_int), intent(in), value :: timer_id
         integer (c_int), intent(out) :: sec, nsec
      end subroutine timer_stop
   end interface
 
   type (MPAS_Clock_type), pointer :: clock

   character (len=StrKIND) :: TOPO_DIR_PATH

   type tile_t
      character (len=StrKIND) :: fname
      type (tile_t), pointer :: next => null()
      type (tile_t), pointer :: prev => null()
      real (c_float), dimension(:,:), pointer, contiguous :: tile
      integer, dimension(2) :: x, y ! The tile's range
      integer :: val ! The value to sort on
   end type tile_t

   type hash_table
      type(tile_t), pointer :: ptr => null()
   end type hash_table

   type (tile_t), target :: tiles
   type (hash_table), allocatable, dimension(:, :) :: tile_hash

   integer(KIND=c_int), parameter :: nx = 1200
   integer(KIND=c_int), parameter :: ny = 1200
   integer(KIND=c_int), parameter :: x_halo = 6
   integer(KIND=c_int), parameter :: y_halo = 6
   integer(KIND=c_int), parameter :: word_size = 2

   integer :: add_tile_cnt = 0

   ! TODO: This is probably defined somewhere else?
   real (KIND=RKIND), parameter :: PI = 4.0 * ATAN(1.0)

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

      use iso_c_binding, only : c_float, c_int

      use mpas_log, only : mpas_log_write
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
      integer, dimension(:), allocatable  :: nhs ! Count of pixel values per mesh cell

      integer :: i, j
      integer :: x, y
      integer :: iErr
      integer :: nCell
      integer :: iPoint
      integer :: cnt = 0
      integer :: unallocated_tiles = 0

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

      integer, pointer :: nCells, maxEdges 
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      character(len=StrKIND), pointer :: config_geog_data_path
      character(len=StrKIND), pointer :: config_topo_data

      real (kind=RKIND) :: tval
      real (kind=RKIND) :: plat, plon

      type (tile_t), pointer :: tile => null()

      iErr = 0
   
      block_ptr => domain % blocklist

      call mpas_pool_get_config(block_ptr % configs, 'config_geog_data_path', config_geog_data_path)
      call mpas_pool_get_config(block_ptr % configs, 'config_topo_data', config_topo_data)

      TOPO_DIR_PATH = trim(config_geog_data_path)//trim(config_topo_data)//"/"


      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'maxEdges', maxEdges)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'ter', ter)

      allocate(tile_hash(0:35, 0:17)) ! TODO: Fix these bounds!
      allocate(nhs(nCells))
      nhs(:) = 0
      
      call mpas_timer_start("static_interp") !!! Timing
   
      call mpas_log_write("Loading in terrain tiles..")
      call mpas_timer_start("tile allocation") !!! Timing

      ! Read/allocate all the tiles that we are responsible for
      do nCell = 1, nCells, 1 
         call get_tile(latCell(nCell), lonCell(nCell), tile)
   
         if (modulo(nCell, 50) == 0 ) then
            call mpas_log_write("We allocated 50 cells: $i", intArgs=(/nCell/))
         endif
      enddo

      call mpas_timer_stop("tile allocation") !!! Timing
      call mpas_log_write("About to interpolate terrain data values to mesh...")

      call mpas_timer_start("tile interp")
      ! Process the tiles 
      do x = 0, size(tile_hash, 1) - 1, 1 ! Forwards
         do y = 0, size(tile_hash, 2) - 1, 1 ! Forwards 

         if ( .NOT. associated(tile_hash(x, y) % ptr)) then
            unallocated_tiles = unallocated_tiles + 1
            cycle
         endif

         cnt = cnt + 1

         if ( modulo(cnt, 10) == 0) then
            call mpas_log_write("We proccessed 10 tiles: $i", intArgs=(/cnt/))
         endif


         iPoint = 1
         ! Loop through each pixel data value of each tile
         do i = 1 + 3, nx
            do j = 1 + 3, ny
               tval = tile_hash(x, y) % ptr % tile (j, i)

               ! Calculate the lat, lon of that pixel
               plat = ( ( j + tile_hash(x, y) % ptr % y(1)) - (1 + y_halo)) / 120. - 90. + 1/120. * 1/2.
               plon = ( ( i + tile_hash(x, y) % ptr % x(1)) - (1 + x_halo)) / 120. - 180. + 1/20. * 1/2.

               plat = plat * PI / 180.0
               plon = plon * PI / 180.0
               
               ! Find which mesh cell the pixel is located within
               iPoint = nearest_cell(plat,plon,iPoint,nCells,maxEdges, &
                                     nEdgesOnCell,cellsOnCell, &
                                     latCell,lonCell)

               ! And assign it to its running average
               ter(iPoint) = ter(iPoint) + tval
               nhs(iPoint) = nhs(iPoint) + 1

            enddo
          enddo
         enddo
      enddo
      call mpas_timer_stop("tile interp")

      do nCell = 1,nCells
         if(nhs(nCell) == 0) then
            write(0,*) "We have a divide by zero! nhs(nCell): ", nhs(nCell)
         endif
         ter(nCell) = ter(nCell) / real(nhs(nCell))
      end do


      ! Deallocate all the tiles
      do i = 0, size(tile_hash, 1) - 1 , 1
         do j = 0, size(tile_hash, 2) - 1, 1
            if(associated(tile_hash(i, j) % ptr)) then
               deallocate(tile_hash(i, j) % ptr % tile)
               deallocate(tile_hash(i, j) % ptr )
            endif
         enddo
      enddo

      call mpas_timer_stop("static_interp") !!! Timing
      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)

   end function terrain_core_run


   !***********************************************************************
   !
   !  function terrain_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015 !> \details 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !***********************************************************************
   !
   !  function get_tile
   !
   !> \brief   Retrive the tile associated with coordinates lat and long
   !> \author  Miles Curry
   !> \date    
   !> \details 
   !>    This function is an abstraction for reading values from a tile file.
   !>  It enables the tiles to be read in from disk once and only once.
   !>
   !>  It does this by using a data structure of the `tile_t` type. Because
   !> the we know all about the geogrphaical datafiles from the index file
   !> we know how many datafiles we know to read to read in, and how large they
   !> are.
   !>
   !> So this function calls search_tile for a tile which contains `lat` 
   !> and `lon` and if it is not found to be allocated it calls add_tile
   !> to add that tile. If the tile is found, a pointer is returned in tile.
   !
   !-----------------------------------------------------------------------
   subroutine get_tile(lat, lon, tile)
      use iso_c_binding, only : c_int, c_char, c_short, c_float
      use mpas_log, only : mpas_log_write

      implicit none

      ! Input variables 
      real (kind=RKIND), value :: lat, lon
      type (tile_t), pointer :: tile

      
      ! Search for the tile in the data structure (ds) we have implemented, if
      ! we find the tile in the ds, then the tile_t type will be returned in 
      ! tile.
      !
      ! If we haven't found the tile in the ds, then we need to add the tile
      ! via add tile, which will also give us the tile_t type that conatins
      ! that tile.

      if (search_tile(lat, lon, tile)) then
         ! pass
      else
         write(0,*) add_tile_cnt
         add_tile_cnt = add_tile_cnt + 1
         call add_tile(lat, lon, tile)
      endif

   end subroutine get_tile

   !***********************************************************************
   !
   !  function search_tile
   !
   !> \brief   Searches for an allocated datafile tile
   !> \author  Miles A. Curry
   !> \date    
   !> \details 
   !>  Given a latitude and longitude, this function searches the implemented
   !> data structure to determine if the tile containing that latitude and
   !> longitude have been allocated or not. A pointer to the tile_t type 
   !> containing that tile is returned in the tile argument.
   !>
   !
   !-----------------------------------------------------------------------
   function search_tile(lat, lon, tile)
      implicit none

      ! Input variables
      real (kind=RKIND), value :: lat, lon
      type (tile_t), intent(inout), pointer :: tile
      
      ! Return value
      logical :: search_tile

      ! Local variables
      character (len=StrKIND) :: fname
      integer :: x, y

      ! Calculate the tiles filename, and its start x and y values. We will
      ! use these values to calculate its location in our hash table
      call gen_file(lat, lon, fname, x, y)

      search_tile = .FALSE.

      if (associated(tile_hash( (x-1)/nx , (y-1)/ny ) % ptr)) then
         tile => tile_hash( (x-1)/nx , (y-1)/ny ) % ptr
         search_tile = .TRUE.
         return
      endif

   end function search_tile

   
   !***********************************************************************
   !
   !  function add_tile
   !
   !> \brief   Reads in a tile from disk and allocates it to main memory
   !> \author  Miles A. Curry
   !> \date    
   !> \details 
   !>    This function reads a datafile from disk and allocates it on the heap
   !> /main memory. Given the latitude and longitude of a point, it calculates
   !> which tilefile contains that points and reads in the datafile that contains
   !> that point.
   !
   !-----------------------------------------------------------------------
   subroutine add_tile(lat, lon, tile)

      use mpas_c_interfacing, only : mpas_f_to_c_string
      use iso_c_binding, only : c_char 

      implicit none

      interface
         integer (c_int) function c_get_tile(file, nx, ny, x_halo, y_halo, word_size, tile) bind(C)
            use iso_c_binding, only : c_int, c_char, c_float, c_ptr
            character (c_char), intent(in) :: file
            integer (c_int), intent(in), value :: nx
            integer (c_int), intent(in), value :: ny
            integer (c_int), intent(in), value :: x_halo
            integer (c_int), intent(in), value :: y_halo
            integer (c_int), intent(in), value :: word_size
            type (c_ptr), value :: tile
         end function c_get_tile
      end interface

      ! Input values
      real (kind=RKIND), intent(in), value :: lat, lon
      type (tile_t), intent(inout), pointer :: tile

      ! Local variables
      integer :: x, y
      character (len=StrKIND) :: fname
      character (len=StrKIND, kind=c_char) :: c_fname
      type(c_ptr) :: tile_ptr

      allocate(tile) ! Allocate the new node 
      allocate(tile%tile(nx + x_halo, ny + y_halo)) ! And its array to hold the tile

      ! Build the filename for the call to the c_get_tile(...) 
      call gen_file(lat, lon, fname, x, y)

      write(c_fname, "(A, A)") trim(TOPO_DIR_PATH), trim(fname)
      c_fname = trim(c_fname)//char(0)
      
      tile_ptr = c_loc(tile%tile) ! Pass in the c_loc of our tile array 
      if(c_get_tile(c_fname, nx, ny, x_halo, y_halo, word_size, tile_ptr) == -1) then
          write(0,*) "We had a problem reading from the geogrid!", trim(fname)
          call mpas_log_write("Error reading the geogrid", messageType=MPAS_LOG_ERR)
          stop
      endif

      tile % fname = fname
      tile % x(1) = x
      tile % x(2) = x + nx - 1
      tile % y(1) = y
      tile % y(2) = y + ny - 1

      tile_hash((tile % x(1) - 1 )/ nx, (tile % y(1) - 1) / ny) % ptr => tile

   end subroutine add_tile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !***********************************************************************
   !  function gen_file
   !
   !> \brief   Generate the filename of the datafile containing lat and lon
   !> \input  
   !>    lat -- a latitude value
   !>    lon -- a longitude value
   !>    fname --  intent(out) - String to contain the filename of the datafile
   !>                            that contains the point latitude and longitude
   !>    x -- intent(out) - The x start of the tile file in tile cords
   !>    y -- intent(out) - The y start of the tile file in tile cords
   !
   !-----------------------------------------------------------------------
   subroutine gen_file(lat, lon, fname, x1, y1)
      implicit none

      real (kind=RKIND), value :: lat, lon
      character (len=StrKIND), intent(out) :: fname
      integer, intent(out) :: x1, y1   ! The start of the tile

      character (len=StrKIND), parameter :: FMT1 = "(I5.5, '-', I5.5, '.', I5.5, '-', I5.5)"
      integer, dimension(2) :: x, y    ! The x, y ranges for the constructing the filename

      lat = lat * (180/PI)
      lon = lon * (180/PI)
      if (lon > 180) then
         lon = lon - 360
      end if

      x(1) = int(120 * ( lon + 180 ))
      y(1) = int(120 * ( lat + 90 ))

      x(2) = x(1) + nx
      x(1) = (x(1) - modulo(x(1), nx))  + 1
      x(2) = x(2) - modulo(x(2), nx)

      if (y(1) == 21600 ) then                        ! Ugly edge cases :(
          y(1) = y(1) - 1
          y(1) = (y(1) - modulo(y(1), ny))  + 1
          y(2) = y(1) + ny
          y(2) = y(2) - modulo(y(2), ny)
      else
          y(2) = y(1) + ny
          y(1) = (y(1) - modulo(y(1), ny))  + 1
          y(2) = y(2) - modulo(y(2), ny)
      end if

      write(fname, FMT1) x(1), x(2), y(1), y(2)

      x1 = x(1)
      y1 = y(1)

   end subroutine gen_file


 integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                               nEdgesOnCell, cellsOnCell, latCell, lonCell)
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: target_lat, target_lon
 integer, intent(in) :: start_cell
 integer, intent(in) :: nCells, maxEdges
 integer, dimension(nCells), intent(in) :: nEdgesOnCell
 integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
 real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

 integer :: i
 integer :: iCell
 integer :: current_cell
 real (kind=RKIND) :: current_distance, d
 real (kind=RKIND) :: nearest_distance

 nearest_cell = start_cell
 current_cell = -1

 do while (nearest_cell /= current_cell)
    current_cell = nearest_cell
    current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                       target_lon, 1.0_RKIND)
    nearest_cell = current_cell
    nearest_distance = current_distance
    do i = 1, nEdgesOnCell(current_cell)
       iCell = cellsOnCell(i,current_cell)
       if (iCell <= nCells) then
          d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
          if (d < nearest_distance) then
             nearest_cell = iCell
             nearest_distance = d
          end if
       end if
    end do
 end do

 end function nearest_cell

 real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

!Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
!sphere with given radius.
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
 real (kind=RKIND) :: arg1
 
 arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )

 sphere_distance = 2.*radius*asin(arg1)

 end function sphere_distance

end module terrain_core
