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

   character (len=StrKIND) :: TOPO_DIR_PATH

   contains
   
   !***********************************************************************
   !
   !  subroutine  init_topo
   !
   !> \brief    Initalizes varabiles needed to acces a topo file
   !> \author   Miles A. Curry
   !> \date 
   !> \details  Checks to see that `topo_dir` exists and if so sets the global
   !>           variable to build off of it later.
   !>  
   !
   !-----------------------------------------------------------------------
   subroutine init_topo(topo_dir)

      implicit none

      character (len=StrKIND), intent(in) :: topo_dir
      logical :: res

      inquire(directory=trim(topo_dir), exist=res)
      if(.not. res) then
        write(0, *) "ERROR: Topo directory does not exists!"
        write(0, *) "usage: ./static_test <meshFile> <topo_dir>"
        stop 
      end if

      TOPO_DIR_PATH = topo_dir
   end subroutine init_topo

   !***********************************************************************
   !
   !  subroutine check_file(file_path)
   !
   !> \brief    Checks to see that the file at file_path exists
   !> \author   Miles A. Curry
   !> \date
   !> \details  Inquires if file at `file_path` exists. if it does return 0
   !>           and if it does not return 1
   !>  
   !
   !-----------------------------------------------------------------------
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

   !***********************************************************************
   !
   !  subroutine convert_rads(lat, lon)
   !
   !> \brief    convert radians to degrees 
   !> \author   Miles A. Curry
   !> \date
   !>  
   !
   !-----------------------------------------------------------------------
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

   !***********************************************************************
   !
   !  subroutine round_t(val, ds)
   !
   !> \brief    Rounds down a value to match the cloest multiple of ds
   !> \author   Miles A. Curry
   !> \date
   !> \details  Rounds a value to the cloest multiple of ds that `val` is
   !>           greater then and returns that value in `val`
   !>  
   !-----------------------------------------------------------------------
   function round_t(val, ds)
      integer :: round_t ! Return Value
      integer, intent(in) :: val
      integer :: ds      ! The the size of each tile in a diretion

      round_t = val - modulo(val, ds) 
   end function round_t

   !***********************************************************************
   !
   !  subroutine calc_indeices(x, y, dx, dy)
   !
   !> \brief    Calculate the indices for a specific x, y in a single tile
   !> \author   Miles A. Curry
   !> \date
   !> \details  Given the global x, y and the changes in x and y (dx, dy)
   !>           calculate the indicies that contain that x, y in a single tile
   !>           
   !>  
   !-----------------------------------------------------------------------
   subroutine calc_indicies(x, y, dx, dy)

      implicit none
      integer, intent(inout) :: x, y
      integer, intent(in) :: dx, dy

      x = x - round_t(x, dx) + 1
      y = y - round_t(y, dy) + 1

   end subroutine calc_indicies


   !***********************************************************************
   !
   !  subroutine conver_cords(lat, lon, x, y)
   !
   !> \brief    Convert latitude and longitude cords into their x, y equivalent
   !> \author   Miles A. Curry
   !> \date
   !> \details  Converts latitude and longitude via the following equations:
   !>
   !>           Y = 120 ( phi + 90 ) + 1
   !>           X = 120 ( rho + 180 ) + 1
   !>  
   !
   !-----------------------------------------------------------------------
   subroutine convert_cords(lat, lon, x, y)
      implicit none

      real (kind=RKIND), intent(in) :: lat, lon
      integer :: x, y

      x = int(120 * ( lon + 180 ))
      y = int(120 * ( lat + 90 ))

   end subroutine convert_cords

   !***********************************************************************
   !
   !  function    gen_file(x, y, dx, dy) 
   !
   !> \brief    Generates the name of the file in which the lat and lon is stored
   !> \author   Miles A. Curry
   !> \date
   !> \details  
   !>  
   !
   !-----------------------------------------------------------------------
   function gen_file(x, y, dx, dy)
      implicit none

      integer, intent(in) :: x, y
      integer, intent(in) :: dx, dy
      integer, dimension(2) :: xr, yr
      character (len=StrKIND), parameter :: FMT1 = "(I5.5, '-', I5.5, '.', I5.5, '-', I5.5)"
      character (len=StrKIND) :: fname
      character (len=StrKIND) :: gen_file

      xr(1) = x
      xr(2) = x + dx
      xr(1) = round_t(xr(1), dx) + 1
      xr(2) = round_t(xr(2), dx)

      yr(1) = y

      if (y == 21600 ) then                        ! Ugly edge cases :(
          yr(1) = round_t(yr(1)-1, dy) + 1
          yr(2) = yr(1) + dy
          yr(2) = round_t(yr(2), dy)
      else
          yr(1) = y
          yr(2) = y + dx
          yr(1) = round_t(yr(1), dy) + 1
          yr(2) = round_t(yr(2), dy)
      end if


      write(fname, FMT1) xr(1), xr(2), yr(1), yr(2)
      gen_file = fname

   end function gen_file

    function get_tile(lat, lon, dx, dy, x_offset, y_offset, word_size, x, y)
!=============================================================================
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! get_tile
        !
        ! 1. Convert the lat, lon to X and Y
        !    - (Convert the lat,lon to degrees if needed)
        ! 2. From the raw X and Y, determine which file they are contained in
        !   ie round them.
        !     - getFile(X, Y) - which calls checkFile(fname)
        ! 3. Get the entire tile by calling read_tile.c.get_tile(...)
        ! 4. Return the tile
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use iso_c_binding, only : c_int, c_char, c_short, c_float
      use mpas_log, only : mpas_log_write

      implicit none

      interface
        function c_get_tile(file, dx, dy, x_offset, y_offset, word_size, tile) bind(C)
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

        ! Input variables 
        real (kind=RKIND), value :: lat, lon
        integer (c_int), intent(in) :: dx, dy
        integer (c_int), intent(in) :: x_offset, y_offset
        integer (c_int), intent(in) :: word_size
        integer, intent(out) :: x, y

        ! Local Variables
        real (c_float), dimension(dx + x_offset, dy + y_offset) :: tile, get_tile
        character (len=StrKIND), parameter :: FILE_PATH_FMT = "(A, A)"
        character (len=StrKIND) :: fname
        character (len=StrKIND, KIND=c_char) :: c_fname
        integer :: ierr

        ! Convert the lat, lon into geogrid x, y
        call convert_rads(lat, lon)
        call convert_cords(lat, lon, x, y)

        ! Generate the filename
        fname = gen_file(x, y, dx, dy)
        write(c_fname, FILE_PATH_FMT) trim(TOPO_DIR_PATH), trim(fname)
        c_fname = trim(c_fname)//char(0)    ! Trim and add null a character 

        ! Check to see if the file exists
        ierr = check_file(c_fname)
        if(ierr == -1) then 
            call mpas_log_write("ERROR: File does not exist!")
            write(0,*) "This file does not exists!", c_fname
            stop
        end if

        ! Read in the tile
        ierr = c_get_tile(trim(c_fname), dx, dy, x_offset, y_offset, word_size, tile)
        if(ierr == -1) then
            write(0,*) "ERROR: Error when reading the data from the & 
                       &tile file:", trim(c_fname)
            stop
        end if
      
        x = x - round_t(x, 1200) + 1
        y = y - round_t(y, 1200) + 1

        get_tile = tile
    end function get_tile


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

      integer, pointer :: nCells
      integer :: iErr
      integer :: my_id
      integer :: nCell
      integer :: x, y

      integer(KIND=c_int), parameter :: dx = 1200
      integer(KIND=c_int), parameter :: dy = 1200
      integer(KIND=c_int), parameter :: x_offset = 6
      integer(KIND=c_int), parameter :: y_offset = 6
      integer(KIND=c_int), parameter :: word_size = 2

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

      character(len=StrKIND), pointer :: config_geog_data_path
      character(len=StrKIND), pointer :: config_topo_data

      real (c_float), dimension(dx + x_offset, dy + y_offset) :: tile
      real (kind=RKIND) :: ter_value

      iErr = 0
   
      block_ptr => domain % blocklist

      call mpas_pool_get_config(block_ptr % configs, 'config_geog_data_path', config_geog_data_path)
      call mpas_pool_get_config(block_ptr % configs, 'config_topo_data', config_topo_data)

      TOPO_DIR_PATH = trim(config_geog_data_path)//trim(config_topo_data)//"/"

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)

      call mpas_pool_get_array(mesh, 'ter', ter)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)

      my_id = domain % dminfo % my_proc_id    

      do nCell = 1, nCells, 1
         tile = get_tile(latCell(nCell), lonCell(nCell), dx, dy, x_offset, y_offset, word_size, x, y)
         ter(nCell) = tile(y,x) 
         if(modulo(nCell, 100) == 0) then
             call mpas_log_write('YIPPY SKIPPY! TIPY just interoploated another 100 values! Number: $i had val: $r',    &
                                  intArgs=(/nCell/), &
                                  realArgs=(/tile(y,x)/))
         endif
      enddo

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

end module terrain_core
