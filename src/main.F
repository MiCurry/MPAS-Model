program hello

    use netcdf
    use pnetcdf
    use mpi

    integer :: ierr, ncid, dimid
    integer (kind=MPI_OFFSET_KIND) :: timeL = 1
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

    call MPI_Finalize(ierr)

    stop

end program hello
