c     data input
      module para
	use cons
      use matrix_inverse
      use exactsolu
      
      contains
      
      subroutine initdata

      dimension a(0:md), aic(0:md,0:md)
      real :: s1
      integer :: k, k1, k2
      
c set up the necessary data before setting the initial condition
      
      if(ic .eq. 1) tend = 0.1
      if(ic .eq. 2) tend = 5.0
      
      pi = 4.0*atan(1.0)
      xleft = 0.0
      xright = 2.0*pi

      cfl = 0.6 / (2.0*mp + 1.0)

c	Gauss-Lobatto Points on [-1,1]
      xg = 0.0
      wg = 0.0
          
      xg(1) = - 1.0
      xg(2) = - 0.2**0.5
      xg(3) = - xg(2)
      xg(4) = - xg(1)

      wg(1) = 1.0/6.0
      wg(2) = 5.0/6.0
      wg(3) = wg(2)
      wg(4) = wg(1)
      
c  the Legendre orthogonal basis 1,xi=2(x-x_i)/dx(i), 0.5*(3*xi**2-1), ...
c  evaluate the basis at - 1.0 (bl) and 1.0 (br)
c  bi(k) are the values of intergral of basis function over [-1.0,1.0]

      do k = 0, 5
        bl(k) = poly(k,-1.0)
        br(k) = poly(k,1.0)
        bm(k) = poly(k,0.0)
        do k1 = 1, 200
c          pv(k1,k) = poly(k,xgl(k1))
          pv(k1,k) = poly(k,1.0*(k1-100.0)/100.0)
        enddo
      enddo
      bi(0) = 1.0
      bi(1:5) = 0.0
      

c  the inverse of the mass matrix:

      aic(0,0) = 1.0

      do  mm = 0, 6
        do i = 0, mm
        aic(i,mm) = 1.0/(2.0*i + 1.0)
        enddo
      enddo

      ai(0:mp) = 1.0 / aic(0:mp,mp)
      
      
c  Matrices on Nodal DG method 
      Mass = 0.0
      Mass_inv = 0.0
      do k = 1, nqua
        Mass(k,k) = wg(k)
        Mass_inv(k,k) = 1.0/Mass(k,k)
      enddo
      
      Vand = 0.0
      do k2 = 1, mp+1
        do k1 = 1, nqua
          Vand(k1,k2) = poly(k2-1,xg(k1))
        enddo
      enddo
      
      Vand_gm = 0.0
      do k = 1, mp+1
        Vand_gm(1,k) = poly(k-1,xg(1))
        Vand_gm(2,k) = poly(k-1,xg(nqua))
      enddo
      
      R_gmk = 0.0
      do k = 1, mp+1
        R_gmk(1,1) = 1.0
        R_gmk(2,nqua) = 1.0
      enddo
      
      Mhat = matmul(transpose(Vand),matmul(Mass,Vand))
      call get_mat_inv(Mhat,Mhat_inv,mp+1)
      
      Pk = 0.0
      do k = 0, mp
        Pk(1:k+1,1:nqua,k) = 
     *    matmul(Mhat_inv(1:k+1,1:k+1)
     *    ,matmul(transpose(Vand(1:nqua,1:k+1)),Mass))
      enddo
      Pk(1,1:nqua,-1) = Pk(1,1:nqua,0)
      
      Idmat = 0.0
      do k = 1, nqua
        Idmat(k,k) = 1.0
      enddo
      
c      Bdy_gm = 0.0
c      Bdy_gm(1,1) = -1.0
c      Bdy_gm(nqua,nqua) = 1.0
      Ext_gm(1,1:nqua,1:nqua) 
     *    = matmul(transpose(R_gmk(1:1,1:nqua)),R_gmk(1:1,1:nqua))
      Ext_gm(2,1:nqua,1:nqua) 
     *    = matmul(transpose(R_gmk(2:2,1:nqua)),R_gmk(2:2,1:nqua))
      
      Diffhat = 0.0
      Diffhat(1,2) = 1.0
      do k = 3, mp+1
        Diffhat(k-1,k) = 2.0*k - 3.0
        Diffhat(:,k) = Diffhat(:,k) + Diffhat(:,k-2)
      enddo
      
      Dhatk = 0.0
      do k = 1, mp+1
        Dhatk(k,k,0) = 1.0
      enddo
      do k = 1, mp
        Dhatk(:,:,k) = matmul(Dhatk(:,:,k-1),Diffhat)
      enddo
      
c  Compute the derivatives of each order on the boundary
      do k = 0, mp
        bdy_dx(:,:,k) = matmul( Vand_gm(:,:)
     *    ,matmul(Dhatk(:,:,k),Pk(:,:,mp)) )
      enddo
      

      Diff = -matmul( transpose( R_gmk(1:1,1:nqua) 
     *    + matmul(Vand_gm(1:1,1:mp+1),Pk(:,:,mp)) )
     *    ,R_gmk(1:1,1:nqua) - matmul(Vand_gm(1:1,1:mp+1),Pk(:,:,mp)) )
      Diff = Diff 
     *    + matmul( transpose( R_gmk(2:2,1:nqua) 
     *    + matmul(Vand_gm(2:2,1:mp+1),Pk(:,:,mp)) )
     *    ,R_gmk(2:2,1:nqua) - matmul(Vand_gm(2:2,1:mp+1),Pk(:,:,mp)) )
      Diff = 0.5*matmul(Mass_inv,Diff)
      Diff = Diff + matmul( Vand,matmul(Diffhat,Pk(:,:,mp)) )

      Stiff = matmul(Mass,Diff)


        er1old = 1.0
	  er2old = 1.0
	  eriold = 1.0

	  kcmax = 1e7

      return
      endsubroutine
      
c  Legendre Polynomial on [-1,1], up to degree 6
c  the return value is an array
      recursive function poly(m,xs)
      integer :: m
      real :: xs
            
      if(m .eq. 0) poly = 1.0
      if(m .eq. 1) poly = xs
      if(m .ge. 2) then
        poly = (2.0*m - 1.0)*xs*poly(m-1,xs) - (m - 1.0)*poly(m-2,xs)
        poly = poly/m
      endif
      
      return
      endfunction
      
      
      
      
      endmodule para
