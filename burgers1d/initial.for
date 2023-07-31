      module initial
	use cons
      use para
      use exactsolu
      contains

c     initial condition
!***************************************************************
!     L^2 Projection of the intial solution 2sin(x)+1
!***************************************************************
      subroutine init
      real :: a(0:md)
      
c set up the initial condition
c L2 projection
	t = 0.0

      xlen = xright - xleft

      dxuni = xlen / nx
      cdx = 1.0 / dxuni
      
      do  i = 0, nx+1
       x(i) = xleft + dxuni*(i - 0.5)
       do j = 1, nqua
         xx(j,i) = x(i) + 0.5*xg(j)*dxuni
       enddo
      enddo

      uc = 0.0

      do i = 0, nx
        do k = 1, nqua
          uc(k,i,0) = 2.0*burgersexac(xx(k,i),0.0) + 1.0
c          uc(k,i,0) = burgersexac(xx(k,i),0.0) 
        enddo
      enddo

	   dt = cfl * dxuni**2.0
c      dt = 0.01*dxuni
      
      
      rhsl = 0.0
      rhsn = 0.0
      rhs = 0.0
      amaxd = 0.0

      kcount = 0
      
      istop = 1

      return
      endsubroutine
      
      endmodule initial
