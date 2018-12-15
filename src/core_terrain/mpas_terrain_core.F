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

   use iso_c_binding, only : c_float, c_int, c_long
   
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

   real (kind=RKIND) :: search_time_total_sec = 0, search_time_total_nsec = 0
   real (kind=RKIND) :: add_time_total_sec = 0, add_time_total_nsec=0

   character (len=StrKIND) :: TOPO_DIR_PATH

   type tile_t
      character (len=StrKIND) :: fname
      type (tile_t), pointer :: next => null()
      type (tile_t), pointer :: prev => null()
      real (c_float), dimension(:,:), allocatable :: tile
      integer, dimension(2) :: x, y ! The tile's range
   end type tile_t

   type hash_table
      type(tile_t), pointer :: tile
   end type hash_table

   type (tile_t), target :: tiles
   type (hash_table), allocatable, dimension(:, :) :: tile_hash

   integer(KIND=c_int), parameter :: dx = 1200
   integer(KIND=c_int), parameter :: dy = 1200
   integer(KIND=c_int), parameter :: x_offset = 6
   integer(KIND=c_int), parameter :: y_offset = 6
   integer(KIND=c_int), parameter :: word_size = 2
   integer, parameter :: num_tiles = 649

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

      ! Timing
      type (MPAS_Time_Type) :: startTime
      integer (c_int) :: allocate_time=0, interp_time=1, sec, nsec


      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: mesh

      real (kind=RKIND), dimension(:), pointer :: ter
      integer, dimension(:), allocatable  :: nhs

      integer :: iErr
      integer :: my_id
      integer :: nCell
      integer :: iPoint
      integer :: x, y
      integer :: i, j
      integer :: cnt

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

      integer, pointer :: nCells, nEdges, maxEdges 
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      character(len=StrKIND), pointer :: config_geog_data_path
      character(len=StrKIND), pointer :: config_topo_data

      real (c_float), dimension(dx + x_offset, dy + y_offset) :: tile
      real (kind=RKIND) :: tval
      real (kind=RKIND) :: plat, plon

      type (tile_t), pointer :: cur

      ! TODO: This is probably defined somewhere else?
      real, parameter :: PI = 4.0 * ATAN(1.0)

      iErr = 0
   
      block_ptr => domain % blocklist

      call mpas_pool_get_config(block_ptr % configs, 'config_geog_data_path', config_geog_data_path)
      call mpas_pool_get_config(block_ptr % configs, 'config_topo_data', config_topo_data)

      TOPO_DIR_PATH = trim(config_geog_data_path)//trim(config_topo_data)//"/"

      allocate(tile_hash(0:35, 0:17)) 

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'maxEdges', maxEdges)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'ter', ter)

      my_id = domain % dminfo % my_proc_id
      allocate(nhs(nCells))
      nhs(:) = 0
      
      ! Read in all the tiles that I am responsible for:
      call mpas_log_write("Loading in terrain tiles..")

      call mpas_timer_start("static_interp")

      call timer_start(allocate_time)
      do nCell = 1, nCells, 1
         tile = get_tile(latCell(nCell), lonCell(nCell)) ! TODO: Overload get_tile to also be allocate tile
      enddo
      call timer_stop(allocate_time, sec, nsec)
      
      call mpas_log_write("Allocation finished time was: $i secs + $i nsec for $i cells", intArgs=(/int(sec), int(nsec), nCells/))

      !TODO: if statement and a warning on associated(cur%next)?
      call mpas_log_write("About to interpolate terrain data values to mesh...")

      call timer_start(interp_time)
      cur => tiles
      cnt = 1

      ! Loop through all the tiles
      do while( associated(cur%next) )
         cnt = cnt + 1

         cur=>cur%next
         tile = cur%tile
   
         iPoint = 1
         do i = 1 + x_offset, dx + x_offset
            do j = 1 + y_offset, dy + y_offset
               tval = tile(j, i) ! Grab the terrain value

               ! Calculate the lat, lon of that tile
               plat = ( ( j + cur % y(1)) - (1 + y_offset)) / 120. - 90. + 1/120. * 1/2.
               plon = ( ( i + cur % x(1)) - (1 + x_offset)) / 120. - 180. + 1/20. * 1/2.

               plat = plat * PI / 180.0
               plon = plon * PI / 180.0
               
               ! Find which mesh cell this pixel is for
               iPoint = nearest_cell(plat,plon,iPoint,nCells,maxEdges, &
                                     nEdgesOnCell,cellsOnCell, &
                                     latCell,lonCell)

               ! And assign it to its running average
               ter(iPoint) = ter(iPoint) + tval
               nhs(iPoint) = nhs(iPoint) + 1

            enddo
         enddo
      enddo
      call timer_stop(interp_time, sec, nsec)
      call mpas_log_write("Interp finished time was: $i secs + $i nsec", intArgs=(/int(sec), int(nsec)/))

      do nCell = 1,nCells
         ter(nCell) = ter(nCell) / real(nhs(nCell))
      end do

      call deallocate_tiles() 

      call mpas_timer_stop("static_interp")
      
      call mpas_log_write("Total searchtime: $r s $r ns", realArgs=(/search_time_total_sec, search_time_total_nsec/))
      call mpas_log_write("Total addtime: $r s $r ns", realArgs=(/add_time_total_sec, add_time_total_nsec/))

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



   subroutine init_topo(topo_dir)

      implicit none

      character (len=StrKIND), intent(in) :: topo_dir
      logical :: res

      inquire(file=trim(topo_dir), exist=res)
      if(.not. res) then
        ! Error
        stop 
      end if

      TOPO_DIR_PATH = topo_dir
   end subroutine init_topo


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Abstraction for getting geogrid files
   !
   ! Files/tiles are either:
   !
   !  a.) Have been read in from the disk and are in memory 
   !  or
   !  b.) Are still on the disk.
   !
   ! Get tile will always return a tile, by either finding a tile that is in
   ! memory or by calling add_tile to read it into memory.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function get_tile(lat, lon)
      use iso_c_binding, only : c_int, c_char, c_short, c_float
      use mpas_log, only : mpas_log_write

      implicit none

      ! Input variables 
      real (kind=RKIND), value :: lat, lon
      real (c_float), dimension(dx + x_offset, dy + y_offset) :: tile, get_tile

      integer (c_int) :: search_time = 2, add_time = 3, sec, nsec
      logical :: search_res 

      call timer_start(search_time) 
      search_res = search_tile(lat, lon, tile)
      call timer_stop(search_time, sec, nsec)
      search_time_total_sec = search_time_total_sec + sec + nsec * 10e-9 

      if ( search_res ) then
         ! pass
      else
         call timer_start(add_time)
         call add_tile(lat, lon, tile)
         call timer_stop(add_time, sec, nsec)
         add_time_total_sec = add_time_total_sec + sec + nsec * 10e-9
      endif

      get_tile = tile
   end function get_tile

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! search_tile - Search the tile list for a tile
   !
   ! return: Returns .FALSE. if that file is not in memory and return 
   !         Returns .TRUE. if that file is in memory and sets `tile` to that 
   !         tile.
   !
   ! Tiles are associated as a dictionary as a linked list with their key being
   ! the name of the file.
   !
   ! Given a latitude and longitute one can generate the tile file that contains
   ! that location by call gen_file(...).
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function search_tile(lat, lon, tile)
      implicit none

      ! Input variables
      real (kind=RKIND), value :: lat, lon
      real (c_float), dimension(dx + x_offset, dy + y_offset), optional :: tile
      
      ! Return value
      logical :: search_tile

      ! Local variables
      character (len=StrKIND) :: fname
      type(tile_t), pointer :: cur
      type(tile_t) :: tile_s
      integer :: x, y

      call gen_file(lat, lon, fname, x, y)

      cur => tiles
      search_tile = .FALSE.
      do while(associated(cur%next))
         if(cur%fname == fname) then
            tile = cur%tile 
            search_tile = .TRUE.
            return
         else
            cur => cur%next
         endif
      enddo

   end function search_tile

   !
   ! add_tile - Reads in a tile from disk and adds it to the list
   !
   subroutine add_tile(lat, lon, tile)
      use iso_c_binding, only : c_char 
      implicit none
      interface
         integer (c_int) function c_get_tile(file, dx, dy, x_offset, y_offset, word_size, tile) bind(C)
            use iso_c_binding, only : c_int, c_char, c_float
            character (c_char), intent(in) :: file
            integer (c_int), intent(in), value :: dx
            integer (c_int), intent(in), value :: dy
            integer (c_int), intent(in), value :: x_offset
            integer (c_int), intent(in), value :: y_offset
            integer (c_int), intent(in), value :: word_size
            real (c_float) :: tile(dx, dy)
         end function c_get_tile
      end interface
      ! Input values
      real (kind=RKIND), value :: lat, lon
      real (c_float), dimension(dx + x_offset, dy + y_offset), intent(inout) :: tile

      ! Local variables
      integer :: x, y
      integer :: H ! Hash Value
      character (len=StrKIND) :: fname
      character (len=StrKIND, kind=c_char) :: c_fname
      type(tile_t), pointer :: cur
      type(tile_t), pointer :: tile_s

      allocate(tile_s) ! Allocate the new node and its dynamic array
      allocate(tile_s%tile(dx + x_offset, dy + y_offset))

      ! TODO: I think there is an module to do this for me?
      ! Build the filename for the call to the C c_get_tile(...) function
      call gen_file(lat, lon, fname, x, y)
      write(c_fname, "(A, A)") trim(TOPO_DIR_PATH), trim(fname)
      c_fname = trim(c_fname)//char(0)

      if( c_get_tile(trim(c_fname), dx, dy, x_offset, y_offset, word_size, tile_s%tile) == -1 ) then
          call mpas_log_write("We had a problem reading the geogrid")
          ! Error
          stop
      endif

      tile_s % fname = fname
      tile_s % x(1) = x
      tile_s % x(2) = x + dx
      tile_s % y(1) = y
      tile_s % y(2) = y + dy

      tile_hash((( tile_s%x(1) - 1 ) / dx), (( tile_s%y(1) - 1 ) / dy))%tile => tile_s

      ! Add the new node to front of the list
      if( .NOT. associated(tiles % next ) ) then
         tiles % next => tile_s
         tile_s % prev => null()
      else
          tiles % next % prev => tile_s
          tile_s % next => tiles % next
          tiles % next => tile_s 
      endif
      
      tile = tile_s % tile

   end subroutine add_tile

   subroutine deallocate_tiles()

      implicit none

      type(tile_t), pointer :: cur

      cur => tiles
      do while(associated(cur%next))
         cur => cur%next
         if(associated(cur%prev)) then
            deallocate(cur%prev)
         endif
      enddo

      deallocate(cur)
      tiles%next=> null()

   end subroutine deallocate_tiles

    function check_file(file_path)
      implicit none
      character (len=StrKIND), intent(in) :: file_path
      integer :: check_file
      logical :: res

      inquire(file=trim(file_path), exist=res)
      if(.not. res) then
         write(0, *) "FAILURE!!: This topo file does not exist!: ", trim(file_path)
         check_file = -1
         return
      end if
      check_file = 0
    end function check_file

   subroutine convert_rads(lat, lon)

      implicit none

      real (kind=RKIND) :: lat, lon
      real, parameter :: PI = 4.0 * ATAN(1.0)

      lat = lat * (180/PI)
      lon = lon * (180/PI)
      if (lon > 180) then
         lon = lon - 360
      end if
   end subroutine convert_rads

   function round_t(val, ds)
      integer :: round_t ! Return Value
      integer, intent(in) :: val
      integer :: ds      ! The the size of each tile in a diretion

      round_t = val - modulo(val, ds) 
   end function round_t

   subroutine gen_file(lat, lon, fname, x1, y1)
      implicit none

      real (kind=RKIND), value :: lat, lon
      character (len=StrKIND), intent(out) :: fname
      integer, intent(out) :: x1, y1   ! The start of the tile

      character (len=StrKIND), parameter :: FMT1 = "(I5.5, '-', I5.5, '.', I5.5, '-', I5.5)"
      integer, dimension(2) :: x, y    ! The x, y ranges for the constructing the filename

      call convert_rads(lat, lon)

      x(1) = int(120 * ( lon + 180 ))
      y(1) = int(120 * ( lat + 90 ))

      x(2) = x(1) + dx
      x(1) = round_t(x(1), dx) + 1
      x(2) = round_t(x(2), dx)

      if (y(1) == 21600 ) then                        ! Ugly edge cases :(
          y(1) = round_t(y(1)-1, dy) + 1
          y(2) = y(1) + dy
          y(2) = round_t(y(2), dy)
      else
          y(2) = y(1) + dy
          y(1) = round_t(y(1), dy) + 1
          y(2) = round_t(y(2), dy)
      end if

      write(fname, FMT1) x(1), x(2), y(1), y(2)

      x1 = x(1)
      y1 = y(1)

   end subroutine gen_file

   subroutine calc_indicies(lat, lon, x, y)

      implicit none
      real(kind=RKIND), value :: lat ,lon
      integer, intent(out) :: x, y

      call convert_rads(lat, lon)

      x = int(120 * ( lon + 180 ))
      y = int(120 * ( lat + 90 ))

      x = x - round_t(x, dx) + 1
      y = y - round_t(y, dy) + 1

   end subroutine calc_indicies

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
