c     data input
      module para
	use cons
      use matrix_inverse
      use exactsolu
      
      contains
      
      subroutine initdata

      dimension a(0:md), aic(0:md,0:md)
      real :: s1
      integer :: k, k1, k2, kk1, kk2, kxg, kyg, ks1, ks2
      integer :: kax(mp+1,mp+1), kay(mp+1,mp+1)
      real :: tmp(mp2,mp2)
      
c set up the necessary data before setting the initial condition
        

      pi = 4.0*atan(1.0)
      gamma = 1.4
      gamma1 = gamma - 1.0
      
      write(60,*) 'Double Mach Reflection' 
      xleft = 0.0
      xright = 3.25
      yleft = 0.0
      yright = 1.0
      tend = 0.2

      cfl = 0.8 / (2.0*mp + 1.0)

!c	Gauss-Legendre Points on [-1,1]
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
      

c  Sort the quadrature points in order
c  Counterclockwise for the boundary quadrature points 
c  lrnbr, btnbr are the indices for the quadrature points
c  start from left-bottom
      do k = 1, nqua
        wg2(k) = wg(1)*wg(k)
        xg2(1,k) = xg(k)
        xg2(2,k) = xg(1)
        lrnbr(k) = nqua + 1 - k
        btnbr(k) = 3*nqua - 1 - k
      enddo
      do k = 2, nqua-1
        wg2(nqua-1+k) = wg(nqua)*wg(k)
        xg2(1,nqua-1+k) = xg(nqua)
        xg2(2,nqua-1+k) = xg(k)
        lrnbr(nqua-1+k) = 4*nqua - 2 - k
        btnbr(nqua-1+k) = 2*nqua - k
      enddo
      do k = 1, nqua
        wg2(2*nqua-2+k) = wg(nqua)*wg(k)
        xg2(1,2*nqua-2+k) = xg(nqua+1-k)
        xg2(2,2*nqua-2+k) = xg(nqua)
        lrnbr(2*nqua-2+k) = 3*nqua - 1 - k
        btnbr(2*nqua-2+k) = nqua + 1 - k
      enddo
      do k = 2, nqua-1
        wg2(3*nqua-3+k) = wg(1)*wg(k)
        xg2(1,3*nqua-3+k) = xg(1)
        xg2(2,3*nqua-3+k) = xg(nqua+1-k)
        lrnbr(3*nqua-3+k) = 2*nqua - k
        btnbr(3*nqua-3+k) = 4*nqua - 2 - k
      enddo
c  Interior points start from left-bottom to right-top
      do k2 = 2, nqua-1
      do k1 = 2, nqua-1
        k = k1 - 1 + (k2 - 2)*(nqua - 2)
        wg2(4*nqua-4+k) = wg(k1)*wg(k2)
        xg2(1,4*nqua-4+k) = xg(k1)
        xg2(2,4*nqua-4+k) = xg(k2)
        lrnbr(4*nqua-4+k) = 4*nqua - 4 + (k2 - 2)*(nqua - 2) 
     *    + nqua - k1
        btnbr(4*nqua-4+k) = 4*nqua - 4 + k1 - 1 
     *    + (nqua - k2 - 1)*(nqua - 2)
      enddo
      enddo
      
      
c  indices for the four boundary points in quadrature
      do k = 1, nqua
        ibpq(k,1) = k
        ibpq(k,2) = nqua - 1 + k
        ibpq(k,3) = 2*nqua - 2 + k
        ibpq(k,4) = 3*nqua - 3 + k
      enddo
      ibpq(nqua,4) = 1
c  Four boundary quadrature points counterclockwise
      do k1 = 1, 4
        xgB(:,:,k1) = xg2(:,ibpq(:,k1))
      enddo
      
        nbr(:,1) = ibpq(nqua:1:-1,3)
        nbr(:,2) = ibpq(nqua:1:-1,4)
        nbr(:,3) = ibpq(nqua:1:-1,1)
        nbr(:,4) = ibpq(nqua:1:-1,2)

      
c     Gauss Points on [-1/2,1/2]
c      xg(1:nqua) = xg(1:nqua)*0.5
c      wg(1:nqua) = wg(1:nqua)*0.5
      
c  the Legendre orthogonal basis 1, 2(x-x_i)/dx(i), 0.5*(3*xi**2-1), ...

      do k = 1, mp2
          
      do k1 = 1, nqua
        do k2 = 1, 4
        bdy(k1,k,k2) = poly2(k,xgB(1,k1,k2),xgB(2,k1,k2))
        enddo
      enddo
      
      do k1 = 1, nqua2
          pv(k1,k) = poly2(k,xg2(1,k1),xg2(2,k1))
      enddo
        
      enddo
      
c read in the inverse of the mass matrix:
      ai(0) = 1.0
      ai(1) = 3.0
      ai(2) = 3.0
      ai(3) = 5.0
      ai(4) = 9.0
      ai(5) = 5.0
      ai(6) = 7.0
      ai(7) = 15.0
      ai(8) = 15.0
      ai(9) = 7.0
      
      
c  Matrices on Nodal DG method 
      Mass = 0.0
      Mass_inv = 0.0
      do k = 1, nqua2
        Mass(k,k) = wg2(k)
        Mass_inv(k,k) = 1.0/Mass(k,k)
      enddo
      
      Vand = 0.0
      do k2 = 1, mp2
        do k1 = 1, nqua2
          Vand(k1,k2) = poly2(k2,xg2(1,k1),xg2(2,k1))
        enddo
      enddo
      
      Vand_gm = 0.0
      do k2 = 1, 4
      do k = 1, mp2
      do k1 = 1, nqua
        Vand_gm(k1,k,k2) = poly2(k,xgB(1,k1,k2),xgB(2,k1,k2))
      enddo
      enddo
      enddo
      
c  R_gmk is an nqua*nqua2 matrix defined on four boundaries
      R_gmk = 0.0
      do k1 = 1, 4
      do k = 1, nqua
        k2 = k + (k1 - 1)*(nqua - 1)
        R_gmk(k,k2,k1) = 1.0
      enddo
      enddo
c      R_gmk(nqua,13,4) = 0.0
      R_gmk(nqua,4*nqua-3,4) = 0.0
      R_gmk(nqua,1,4) = 1.0
      
      Mhat = matmul(transpose(Vand),matmul(Mass,Vand))
      call get_mat_inv(Mhat,Mhat_inv,mp2)
      
      Pk = 0.0
      do k = 0, mp
        k1 = (k + 1)*(k + 2)/2
        Pk(1:k1,1:nqua2,k) = 
     *    matmul(Mhat_inv(1:k1,1:k1)
     *    ,matmul(transpose(Vand(1:nqua2,1:k1)),Mass))
      enddo
      Pk(1,1:nqua2,-1) = Pk(1,1:nqua2,0)
      
      Idmat = 0.0
      do k = 1, nqua2
        Idmat(k,k) = 1.0
      enddo
      
      Bdy_gm = 0.0
      do k = 1, nqua
      Bdy_gm(k,k) = wg(k)
      enddo
      do k1 = 1, 4
      Ext_gm(:,:,k1) 
     *    = matmul( transpose(R_gmk(:,:,k1)) 
     *    , matmul(Bdy_gm,R_gmk(:,:,k1)) )
      enddo
      
      Diffhat = 0.0
      Diffhat(1,2) = 1.0
      do k = 3, mp+1
        Diffhat(k-1,k) = 2.0*k - 3.0
        Diffhat(:,k) = Diffhat(:,k) + Diffhat(:,k-2)
      enddo
c
      kax = 0
      Diffhat2x = 0.0
      do k = 1, mp+1
        kax(k,k) = k*(k + 1)/2
        do k1 = k, mp
        kax(k1+1,k) = kax(k1,k) + k1
        enddo
        do k2 = k, mp+1
        do k1 = k, mp+1
        Diffhat2x(kax(k1,k),kax(k2,k)) = Diffhat(k1-k+1,k2-k+1)
        enddo
        enddo
      enddo
c
      kay = 0
      Diffhat2y = 0.0
      do k = 1, mp+1
        kay(k,k) = k*(k - 1)/2 + 1
        do k1 = k, mp
        kay(k1+1,k) = kay(k1,k) + k1 + 1
        enddo
        do k2 = k, mp+1
        do k1 = k, mp+1
        Diffhat2y(kay(k1,k),kay(k2,k)) = Diffhat(k1-k+1,k2-k+1)
        enddo
        enddo
      enddo
      
      
      Dhatk = 0.0
      Dhatk(:,:,1) = Diffhat
      do k = 2, mp
        do k1 = 3, mp+1
          Dhatk(:,k1,k) = (2.0*k1 - 3.0)*Dhatk(:,k1-1,k-1)
          Dhatk(:,k1,k) = Dhatk(:,k1,k) + Dhatk(:,k1-2,k)
        enddo
      enddo

c  Dhatk2(:,:,m) is the (m1,m2)-th order differential matrix in 2D
      Dhatk2 = 0.0
      Dhatk2(:,:,2) = Diffhat2x
      Dhatk2(:,:,3) = Diffhat2y

      tmp = 0.0
      do k = 1, mp2
        tmp(k,k) = 1.0
      enddo

      do k = 2, mp
        k1 = k*(k + 1)/2 + 1
        k2 = (k + 1)*(k + 2)/2

        do kk = k1, k2
          kk1 = (k + 1)*(k + 2)/2 - kk
          kk2 = k - kk1
c
          Dhatk2(:,:,kk) = tmp
          do ks1 = 1, kk1
            Dhatk2(:,:,kk) = matmul(Dhatk2(:,:,kk),Dhatk2(:,:,2))
          enddo
          do ks2 = 1, kk2
            Dhatk2(:,:,kk) = matmul(Dhatk2(:,:,kk),Dhatk2(:,:,3))
          enddo
        
        enddo

      enddo
      
      
c      Diff = matmul( Idmat + matmul(Vand,Pk),matmul( Ext_gm
c     *    ,Idmat - matmul(Vand,Pk) ) )
      Diff2x = Diff2x + matmul( transpose( R_gmk(:,:,2)
     *    + matmul(Vand_gm(:,:,2),Pk(:,:,mp)) )
     *    , matmul( Bdy_gm(:,:), R_gmk(:,:,2)
     *    - matmul(Vand_gm(:,:,2),Pk(:,:,mp)) ) )

      Diff2x = Diff2x - matmul( transpose( R_gmk(:,:,4)
     *    + matmul(Vand_gm(:,:,4),Pk(:,:,mp)) )
     *    , matmul( Bdy_gm(:,:), R_gmk(:,:,4)
     *    - matmul(Vand_gm(:,:,4),Pk(:,:,mp)) ) )

      Diff2x = 0.5*matmul(Mass_inv,Diff2x)
      Diff2x = Diff2x + matmul( Vand,matmul(Diffhat2x,Pk(:,:,mp)) )
      Stiff2x = matmul(Mass,Diff2x)
c
      Diff2y = Diff2y + matmul( transpose( R_gmk(:,:,3)
     *    + matmul(Vand_gm(:,:,3),Pk(:,:,mp)) )
     *    , matmul( Bdy_gm(:,:), R_gmk(:,:,3)
     *    - matmul(Vand_gm(:,:,3),Pk(:,:,mp)) ) )

      Diff2y = Diff2y - matmul( transpose( R_gmk(:,:,1)
     *    + matmul(Vand_gm(:,:,1),Pk(:,:,mp)) )
     *    , matmul( Bdy_gm(:,:), R_gmk(:,:,1)
     *    - matmul(Vand_gm(:,:,1),Pk(:,:,mp)) ) )

      Diff2y = 0.5*matmul(Mass_inv,Diff2y)
      Diff2y = Diff2y + matmul( Vand,matmul(Diffhat2y,Pk(:,:,mp)) )
      Stiff2y = matmul(Mass,Diff2y)

      
c  bdy_dx(nqua,nqua2,mp2,4) is the differential terms maping from the interior
c  quadrature points to the four boundary quadrature points
      do k1 = 1, 4
        bdy_dx(:,:,1,k1) = matmul(bdy(:,:,k1),Pk(:,:,mp))
      enddo

      do k1 = 1, 4
      do k = 2, mp2
        bdy_dx(:,:,k,k1) = matmul( matmul(Vand_gm(:,:,k1)
     *  ,Dhatk2(:,:,k)),Pk(:,:,mp) )
      enddo
      enddo

      er1old = 1.0
      er2old = 1.0
      eriold = 1.0

	  kcmax = 1e7

      return
      endsubroutine
      
c  Basis functions in 2D
      function poly2(m,xs,ys)
      integer :: m, a(1:2,1:66), k, kk, k1
      real :: xs, ys, poly2
      
      kk = 0
      do k = 0, m
        do k1 = 0, k
        kk = kk + 1
        a(1,kk) = k - k1
        a(2,kk) = k1
        enddo
      enddo
      
      poly2 = poly(a(1,m),xs)*poly(a(2,m),ys)
      
      
      return
      endfunction
      
c  Legendre Polynomial on [-1,1], up to arbitrary degree
c  the return value is an array
      recursive function poly(m,xs)
      integer :: m
      real :: xs, poly
          
      if(m .eq. 0) poly = 1.0
      if(m .eq. 1) poly = xs
      if(m .ge. 2) then
        poly = (2.0*m - 1.0)*xs*poly(m-1,xs) 
     *    - (m - 1.0)*poly(m-2,xs)
        poly = poly/m
      endif
      
      return
      endfunction
      
      
      
      endmodule para
