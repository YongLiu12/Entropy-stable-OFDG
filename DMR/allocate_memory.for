c     def of parameter
      Module allocate_memory
      use cons

      contains

      subroutine assign_memory

      implicit none
      include 'mpif.h'
      integer :: ndims, i, j, err_code
      
      
      nxlocal = int(nx/Nx_processor)
      nylocal = int(ny/Ny_Processor)
      ntmpx = nxlocal*nqua2*nequ
      ntmpy = nylocal*nqua2*nequ
      
      
      allocate(x(0:nxlocal+1), y(0:nylocal+1))
      allocate(uc(1:nqua2,1:nequ,0:nxlocal+1,0:nylocal+1,0:6))
      allocate(rhs(1:nqua2,1:nequ,0:nxlocal+1,0:nylocal+1,0:6))
      allocate(uc_old(1:nqua2,1:nequ,0:nxlocal+1,0:nylocal+1))
      allocate(xsdrv_data(1:ntmpx),ysdrv_data(1:ntmpy))
      
      
      ndims = 2
      reorder = .FALSE.
      
      dims(1) = Nx_Processor
      dims(2) = Ny_Processor
      
      periods(1) = .TRUE.
      periods(2) = .TRUE.

!*    Create a grid structure in WORLD group and inquire about own coordinates
      call MPI_CART_CREATE(comm,ndims,dims,periods,reorder,
     * comm_cart,err_code)
      call MPI_COMM_RANK(comm_cart,myid,err_code)
      call MPI_COMM_SIZE(comm_cart,N_Processor,err_code)
      call MPI_CART_COORDS(comm_cart,myid,ndims,coords,err_code)
    
      xcoord = coords(1) + 1
      ycoord = coords(2) + 1
      
!     Loop up the ranks for the neighbors. My process coordinates are (i,j).
!     Neighbors are left(i-1,j), right(i+1,j), bottom(i,j-1), top(i,j+1)

      i = xcoord
      j = ycoord
      
      write(60+myid,*)"myid=",myid,"ix=",xcoord,"jy=",ycoord


!* i:
      coords(2) = j - 1
      if(i .le. 0) then
         write(*,*) ' --- Error: i <= 0!(In Sub.xassign_memory) ---'
         err_code = -1
!        call x_exitmpi
        call MPI_Finalize (IERR)
         stop 
      else if(Nx_Processor .eq. 1) then
         left = MPI_PROC_NULL
      else if(i .eq. 1) then
         coords(1) = Nx_Processor - 1
         call MPI_CART_RANK(comm_cart, coords, left, err_code)
      else
         coords(1) = i - 1 - 1
         call MPI_CART_RANK(comm_cart, coords, left, err_code)
      end if

      if(i .gt. Nx_Processor) then
         write(*, *)' --- Error: i > XBLKNO!(In Sub.xassign_memory) ---'
         err_code = -1 
         !call x_exitmpi
         call MPI_Finalize (IERR)
         stop ' Error: i > XBLKNO!'
      else if(Nx_Processor .eq. 1) then
         right = MPI_PROC_NULL
      else if(i .eq. Nx_Processor) then
         coords(1) = 0
         call MPI_CART_RANK(comm_cart, coords, right, err_code)
      else
         coords(1) = i - 1 + 1 
         call MPI_CART_RANK(comm_cart, coords, right, err_code)
      end if

!* j:
      coords(1) = i - 1
      if(j .le. 0) then
         write(*,*) ' --- Error: j <= 0!(In Sub.xassign_memory) ---'
         err_code = -1
        ! call x_exitmpi
        call MPI_Finalize (IERR)
         stop ' Error: j <= 0!'
      else if(Ny_Processor .eq. 1) then
         bottom = MPI_PROC_NULL
      else if(j .eq. 1) then
         coords(2) = Ny_Processor -1
         call MPI_CART_RANK(comm_cart, coords, bottom, err_code)
      else
         coords(2) = j - 1 - 1
         call MPI_CART_RANK(comm_cart, coords, bottom, err_code)
      end if
 
      if(j .gt. Ny_Processor) then
         write(*, *)' --- Error: j > YBLKNO!(In Sub.xassign_memory) ---'
         err_code = -1
         !call x_exitmpi
         call MPI_Finalize (IERR)
         stop ' Error: j > YBLKNO!'
      else if(Ny_Processor .eq. 1) then
         top = MPI_PROC_NULL
      else if(j .eq. Ny_Processor) then
         coords(2) = 0
         call MPI_CART_RANK(comm_cart, coords, top, err_code)
      else
         coords(2) = j - 1 + 1
         call MPI_CART_RANK(comm_cart, coords, top, err_code)
      end if

      
      write(60+myid,*) "myid=",myid,"ix=",xcoord,"jy=",ycoord
      write(60+myid,*) left, right, bottom, top
      
      
      return
      end subroutine
      

      subroutine delete_memory
      implicit none
      
      deallocate(x, y)
      deallocate(uc, uc_old)
      deallocate(rhs)
      deallocate(xsdrv_data, ysdrv_data)
      
      
      
      return
      end subroutine
      
      
	endmodule allocate_memory
