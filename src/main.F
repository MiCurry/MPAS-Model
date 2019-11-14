program hello

    use netcdf
    use pnetcdf
    use mpi
    use pio
    use pio_types
    use ESMF
    use ESMF_ClockMod


    integer :: ierr, ncid, dimid
    integer (kind=MPI_OFFSET_KIND) :: timeL = 1
    type (iosystem_desc_t) :: io_system
    type (ESMF_Clock) :: clock
    type (ESMF_Time) :: time

#ifdef CORE_ATMOSPHERE
    write(0,*) 'This is MPAS-Atmosphere'
#endif
#ifdef CORE_OCEAN
    write(0,*) 'This is MPAS-Ocean'
#endif

    call MPI_Init(ierr)
    ! NetCDF
    ierr = nf90_create('blah.nc', NCF90_CLOBBER, ncid)
    ierr = nf90_def_dim(ncid, "Time", 1, dimid)
    ierr = nf90_enddef(ncid)
    ierr = nf90_close(ncid)

    ! PNetCDF
    ierr = nf90mpi_create(MPI_COMM_WORLD, 'pnetcdfBlash.nc', NC_CLOBBER, MPI_INFO_NULL, ncid)
    ierr = nf90mpi_def_dim(ncid, "Time", timeL, dimid)
    ierr = nf90mpi_enddef(ncid)
    ierr = nf90mpi_close(ncid)

    ! PIO
    call PIO_init(0, MPI_COMM_WORLD, 1, 0, 1, PIO_rearr_box, io_system)

    call MPI_Finalize(ierr)

    stop

end program hello
