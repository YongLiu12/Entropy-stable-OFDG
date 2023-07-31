      module initial
	  use cons
      use para
      use exactsolu
      
      contains

c     initial condition
!***************************************************************
!     Initialization
!***************************************************************
      subroutine init 
      
      implicit none
      include 'mpif.h'
      real :: a(0:md)
      real :: xlen, ylen, s1, s2
      integer :: i, j, k
      
c set up the initial condition

	t = 0.0

      xlen = xright - xleft
      ylen = yright - yleft

      dx = xlen / nx
      dy = ylen / ny
      cdx = 1.0 / dx
      cdy = 1.0 / dy
      
      x = 0.0
      y = 0.0
      
      x(0) = xleft + ((xcoord - 1)*nxlocal - 0.5)*dx
      do  i = 1, nxlocal+1
       x(i) = x(i-1) + dx
      enddo
      
      y(0) = yleft + ((ycoord - 1)*nylocal - 0.5)*dy
      do  j = 1, nylocal+1
       y(j) = y(j-1) + dy
      enddo
      

      uc = 0.0
c      amax = 0.0
c      amay = 0.0
      do j = 1, nylocal
      do i = 1, nxlocal
        do k = 1, nqua2
          s1 = x(i) + 0.5*xg2(1,k)*dx
          s2 = y(j) + 0.5*xg2(2,k)*dy
          uc(k,:,i,j,0) = u0(s1,s2)
c          amax = max( amax,abs(d_fluxfunc1(uc(k,i,j,0))) )
c          amay = max( amay,abs(d_fluxfunc2(uc(k,i,j,0))) )
        enddo
      enddo
      enddo
      
      
      if(ic .eq. 7 .and. myid .eq. nroot) then
          
        do k = 1, nqua2
          uc(k,4,1,1,0) = 0.244816*cdx*cdy
        enddo
        
      endif

c	   dt = cfl * dx*dx
      dt = 1.0e-10
      
      
      rhsl = 0.0
      rhsn = 0.0
      rhs = 0.0
      amaxd = 0.0

      kcount = 0
      
      istop = 1
      
      
      return
      endsubroutine
      
      endmodule initial
