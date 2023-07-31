      module getrhs
	  use cons
      use boundary
      use para
      use exactsolu
      
      contains

!**********************************************************************
!     the numerical flux for U_t + F(U)_x + G(U)_y = 0 
!**********************************************************************
      subroutine res(m)
      
      implicit none
      include 'mpif.h'
      integer :: i, j, m, ierror, istop0, istop1
      integer :: k, kk, k1, k2
      
      real :: us(1:nqua2,1:nequ)
      real :: f_ent(1:nqua2,1:nqua2,1:nequ,1:2)
      real, dimension(0:md) :: a, b, c, scm_jmp, fx
      real, dimension(1:nqua,1:nequ,1:4) :: flx, flxhat
      real, dimension(1:nequ) :: upos, uneg
      real :: s1, s2, alpha
      real :: den, xmt, ymt, eng, vex, vey, pre, vxm, vym, hm
     *    , qm, cvel, cm, rcm, t0, t1, t2, b1, b2
      real :: abs, max, min
      real, dimension(1:nequ,1:nequ) :: evlx, evrx, evly, evry
      
      real,dimension(0:md,1:nequ) :: sig, sig0
      real :: fjump(0:md,1:nequ,1:4)
      
c      real :: sigma(0:md)
c      dimension :: dampc(0:md,0:nn)
      
      real, allocatable, dimension(:,:) :: amx, amy 
      real, allocatable, dimension(:,:,:,:) :: hh 
      real, allocatable, dimension(:,:) :: w, h, vx, vy
      real, allocatable, dimension(:,:,:,:,:) :: omj, omjx, omjy

      
      allocate(amx(0:nxlocal+1,0:nylocal+1))
      allocate(amy(0:nxlocal+1,0:nylocal+1))
      allocate( w(0:nxlocal+1,0:nylocal+1)
     *    , h(0:nxlocal+1,0:nylocal+1)
     *    , vx(0:nxlocal+1,0:nylocal+1)
     *    , vy(0:nxlocal+1,0:nylocal+1) )

        
c  Periodic B.C.
        
      call bc(m)
      
      
c  Prepare for alpha in LLF flux and Roe mean 
      
      w = 0.0
      h = 0.0
      vx = 0.0
      vy = 0.0
      amx = 0.0
      amy = 0.0
	  amax = 0.0
        amay = 0.0
        istop0 = 1
      
      do j = 0, nylocal+1
      do i = 0, nxlocal+1
          
      if ( i**2 + j**2 .ne. 0 
     *    .and. (i-nxlocal-1)**2 + j**2 .ne. 0 
     *    .and. i**2 + (j-nylocal-1)**2 .ne. 0 
     *    .and. (i-nxlocal-1)**2 
     *    + (j-nylocal-1)**2 .ne. 0 ) then
          
          den = sum(Pk(1,:,mp)*uc(:,1,i,j,m))
          xmt = sum(Pk(1,:,mp)*uc(:,2,i,j,m))
          ymt = sum(Pk(1,:,mp)*uc(:,3,i,j,m))
          eng = sum(Pk(1,:,mp)*uc(:,4,i,j,m))
          t0 = 1.0/den
          vex = xmt*t0
          vey = ymt*t0
          pre = gamma1*( eng - 0.5*den*(vex*vex+vey*vey) )
          
	  if(den .lt. 0.0) then
            write(10000+myid,*) 'density becomes negative' 
            write(10000+myid,*) kcount, m, t, dt
            write(10000+myid,*) i, j, x(i), y(j)
            write(10000+myid,*) den, xmt, ymt, eng
            istop0 = 0
        endif
c          
	  if(pre .lt. 0.0) then
c            write(20000+myid,*) 'pressure becomes negative'
c            write(20000+myid,*) kcount, m, t, dt
c            write(20000+myid,*) myid, i, j, x(i), y(j)
c            write(20000+myid,*) den, xmt, ymt, eng
            pre = epweno
        endif
          
        w(i,j) = abs(den)**0.5
        h(i,j) = (pre + eng)*t0
        vx(i,j) = vex
        vy(i,j) = vey
        

c  compute the alpha in LF flux on cell boundary points 
      do k = 1, nqua2
          
          den = uc(k,1,i,j,m)
          xmt = uc(k,2,i,j,m)
          ymt = uc(k,3,i,j,m)
          eng = uc(k,4,i,j,m)
          t0 = 1.0/den
          vex = xmt*t0
          vey = ymt*t0
          pre = press(uc(k,1:nequ,i,j,m))

          cvel = sqrt( gamma*abs(pre*t0) ) 
          amx(i,j) = max( amx(i,j),abs(vex) + cvel )
          amax = max( amax,amx(i,j) )
          amy(i,j) = max( amy(i,j),abs(vey) + cvel )
          amay = max( amay,amy(i,j) )
          
        enddo

      
      endif
      
      enddo
      enddo
c
      call MPI_Barrier(MPI_Comm_World, IERR)
      
c  negative density occur and stop the simulation

      call MPI_Reduce(istop0, istop1, 1, MPI_INTEGER, MPI_MIN
     * , nroot, comm_cart, IERR)
      if(myid .eq. nroot .and. istop1 .eq. 0) then
        write(60,*) 'density becomes negative'
        stop
      endif
      
c     Find the maximum amax,amay from all the processes:
      call MPI_Reduce(amax, amax0, 1, MPI_REAL8, MPI_MAX, nroot, 
     * comm_cart, IERR)
      call MPI_Reduce(amay, amay0, 1, MPI_REAL8, MPI_MAX, nroot, 
     * comm_cart, IERR)   
      
      call MPI_Bcast(amax0, 1, MPI_REAL8, nroot, comm_cart, IERR)   
      call MPI_Bcast(amay0, 1, MPI_REAL8, nroot, comm_cart, IERR)
      
      
c compute the residue
      
      do j = 1, nylocal
      do i = 1, nxlocal
        
      do k = 1, nqua2
        us(k,:) = uc(k,:,i,j,m)
      enddo
      
      do k2 = 1, nqua2
      do k1 = 1, nqua2
        f_ent(k1,k2,:,1) = fmS(us(k1,:),us(k2,:),1)
        f_ent(k1,k2,:,2) = fmS(us(k1,:),us(k2,:),2)
      enddo
      enddo
      
      do kk = 1, nequ
      do k = 1, nqua2
      
c compute the contribution from the Difference term
      
        rhs(k,kk,i,j,m) = - 2.0*sum( Diff2x(k,1:nqua2)
     *    *f_ent(k,1:nqua2,kk,1) )*2.0*cdx
        rhs(k,kk,i,j,m) = rhs(k,kk,i,j,m)
     *  - 2.0*sum( Diff2y(k,1:nqua2)
     *    *f_ent(k,1:nqua2,kk,2) )*2.0*cdy
      
      enddo
      enddo
      
c compute the contribution from the simultaneous approximation term

      do k = 1, nqua
        flx(k,:,1) = -fluxfunc2(uc(ibpq(k,1),:,i,j,m))
        flx(k,:,2) = fluxfunc1(uc(ibpq(k,2),:,i,j,m))
        flx(k,:,3) = fluxfunc2(uc(ibpq(k,3),:,i,j,m))
        flx(k,:,4) = -fluxfunc1(uc(ibpq(k,4),:,i,j,m))
      enddo
      
      do k = 1, nqua
          
        uneg = uc(nbr(k,1),:,i,j-1,m)
        upos = uc(ibpq(k,1),:,i,j,m)
        alpha = max( amy(i,j-1),amy(i,j) )*1.1
        flxhat(k,:,1) = - 0.5*( fluxfunc2(upos) + fluxfunc2(uneg)
     *    - alpha*(upos - uneg) )
        
        uneg = uc(ibpq(k,2),:,i,j,m)
        upos = uc(nbr(k,2),:,i+1,j,m)
        alpha = max( amx(i,j),amx(i+1,j) )*1.1
        flxhat(k,:,2) = 0.5*( fluxfunc1(upos) + fluxfunc1(uneg)
     *    - alpha*(upos - uneg) )
        
        uneg = uc(ibpq(k,3),:,i,j,m)
        upos = uc(nbr(k,3),:,i,j+1,m)
        alpha = max( amy(i,j),amy(i,j+1) )*1.1
        flxhat(k,:,3) = 0.5*( fluxfunc2(upos) + fluxfunc2(uneg)
     *    - alpha*(upos - uneg) )
        
        uneg = uc(nbr(k,4),:,i-1,j,m)
        upos = uc(ibpq(k,4),:,i,j,m)
        alpha = max( amx(i-1,j),amx(i,j) )*1.1
        flxhat(k,:,4) = - 0.5*( fluxfunc1(upos) + fluxfunc1(uneg)
     *    - alpha*(upos - uneg) )
        
      enddo
        
      do k = 1, nqua2 
      do k1 = 1, nqua
      
        rhs(k,:,i,j,m) = rhs(k,:,i,j,m) + R_gmk(k1,k,1)
     *    *wg(k1)/wg2(k)*(flx(k1,:,1) - flxhat(k1,:,1))*2.0*cdy
        rhs(k,:,i,j,m) = rhs(k,:,i,j,m) + R_gmk(k1,k,2)
     *    *wg(k1)/wg2(k)*(flx(k1,:,2) - flxhat(k1,:,2))*2.0*cdx
        rhs(k,:,i,j,m) = rhs(k,:,i,j,m) + R_gmk(k1,k,3)
     *    *wg(k1)/wg2(k)*(flx(k1,:,3) - flxhat(k1,:,3))*2.0*cdy
        rhs(k,:,i,j,m) = rhs(k,:,i,j,m) + R_gmk(k1,k,4)
     *    *wg(k1)/wg2(k)*(flx(k1,:,4) - flxhat(k1,:,4))*2.0*cdx

      enddo
      enddo
      
      enddo
      enddo
      
      
      deallocate(amx,amy)
      
      allocate( omj(1:nqua2,1:nequ,0:nxlocal,0:nylocal,1:4)
     *    , omjx(1:nqua2,1:nequ,0:nxlocal,0:nylocal,1:4)
     *    , omjy(1:nqua2,1:nequ,0:nxlocal,0:nylocal,1:4)  )
      
      
c  Add the numerical damping

c  compute the jumps for U & V
      
      do j = 0, nylocal
      do i = 0, nxlocal
          
      if(i*i + j*j .ne. 0) then
          
c  characteristic decomposition of F'(U)
          
          t0 = w(i,j) / ( w(i,j) + w(i+1,j) )
          t1 = 1.0 - t0
          vxm = t0 * vx(i,j) + t1 * vx(i+1,j)
          vym = t0 * vy(i,j) + t1 * vy(i+1,j)
          hm = t0 *  h(i,j) + t1 *  h(i+1,j)
          
         call getmtx(vxm,vym,hm,evlx,evrx) 
          
c  characteristic decomposition of G'(U)
          
          t0 = w(i,j) / ( w(i,j) + w(i,j+1) )
          t1 = 1.0 - t0
          vxm = t0 * vx(i,j) + t1 * vx(i,j+1)
          vym = t0 * vy(i,j) + t1 * vy(i,j+1)
          hm = t0 *  h(i,j) + t1 *  h(i,j+1)
          
         call getmty(vxm,vym,hm,evly,evry) 
          
c  Compute the jumps of all derivatives 
        do k = 1, 3!mp2
            
c  jump on (x_{i+1/2},y_{j-1/2}) between K_{i,j}&K_{i+1,j}
        do kk = 1, nequ
          uneg(kk) = sum(uc(:,kk,i,j,m)*bdy_dx(1,:,k,2))
          upos(kk) = sum(uc(:,kk,i+1,j,m)*bdy_dx(nqua,:,k,4))
        enddo
        omj(k,1:nequ,i,j,1) = upos - uneg 
        omjx(k,1:nequ,i,j,1) = matmul(evlx,upos - uneg) 
        
c  jump on (x_{i+1/2},y_{j+1/2}) between K_{i,j}&K_{i+1,j}
        do kk = 1, nequ
          uneg(kk) = sum(uc(:,kk,i,j,m)*bdy_dx(nqua,:,k,2))
          upos(kk) = sum(uc(:,kk,i+1,j,m)*bdy_dx(1,:,k,4))
        enddo
        omj(k,1:nequ,i,j,2) = upos - uneg 
        omjx(k,1:nequ,i,j,2) = matmul(evlx,upos - uneg) 
        
c  jump on (x_{i-1/2},y_{j+1/2}) between K_{i,j}&K_{i,j+1}
        do kk = 1, nequ
          uneg(kk) = sum(uc(:,kk,i,j,m)*bdy_dx(nqua,:,k,3))
          upos(kk) = sum(uc(:,kk,i,j+1,m)*bdy_dx(1,:,k,1))
        enddo
        omj(k,1:nequ,i,j,3) = upos - uneg 
        omjy(k,1:nequ,i,j,3) = matmul(evly,upos - uneg) 
        
c  jump on (x_{i+1/2},y_{j+1/2}) between K_{i,j}&K_{i,j+1}
        do kk = 1, nequ
          uneg(kk) = sum(uc(:,kk,i,j,m)*bdy_dx(1,:,k,3))
          upos(kk) = sum(uc(:,kk,i,j+1,m)*bdy_dx(nqua,:,k,1))
        enddo
        omj(k,1:nequ,i,j,4) = upos - uneg 
        omjy(k,1:nequ,i,j,4) = matmul(evly,upos - uneg) 
      
      enddo
c  end of the computation of the jumps 
        
        endif
          
        
      enddo
        enddo
        
        
      deallocate(w,h,vx,vy)

c compute the jump as the coefficients

      do j = 1, nylocal
      do i = 1, nxlocal

        sig = 0.0
        a = 0.0
        do k1 = 0, 1!mp

          k2 = k1*(k1 + 1)/2 + 1
          do k = k2, k2 + k1

          fjump(k,:,1:4) = 0.0
c     
          fjump(k,:,1) = fjump(k,:,1) 
     *    + omjy(k,1:nequ,i,j-1,3)**2.0
     *    + omjy(k,1:nequ,i,j-1,4)**2.0
          fjump(k,:,2) = fjump(k,:,2) 
     *    + omjx(k,1:nequ,i,j,1)**2.0
     *    + omjx(k,1:nequ,i,j,2)**2.0
          fjump(k,:,3) = fjump(k,:,3) 
     *    + omjy(k,1:nequ,i,j,3)**2.0
     *    + omjy(k,1:nequ,i,j,4)**2.0 
          fjump(k,:,4) = fjump(k,:,4) 
     *    + omjx(k,1:nequ,i-1,j,1)**2.0 
     *    + omjx(k,1:nequ,i-1,j,2)**2.0
          
c         
          do kk = 1, nequ
          sig(k1,kk) = sig(k1,kk) + sum(fjump(k,kk,1:4))/4.0
          enddo
          
          enddo
      
        
        enddo
      
      b = 0.0
      sig0 = 0.0
      do k = 0, 1
        do k1 = 0, k
        sig0(k,:) = sig0(k,:) + sig(k1,:)
        enddo
      enddo
        
      
      a(1) = max( sig0(1,1),sig0(1,2),sig0(1,3),sig0(1,4) )
      a(1) = a(1) ** 0.5
      
      
        amaxd(m) = max( amaxd(m),a(1) )
      
      us = uc(:,:,i,j,m)
c  Pk.b are the coefficients, V.Pk.b are the quadrature point values
      do kk = 1, nequ
      do k = 1, 1  
        rhs(:,kk,i,j,m) = rhs(:,kk,i,j,m) - a(1)*cdx 
     *    *( us(1:nqua2,kk) - matmul( Vand(1:nqua2,1:k) 
     *    ,matmul(Pk(1:k,1:nqua2,k-2),us(1:nqua2,kk)) ) ) 
      enddo
      enddo
      
      
      enddo
      enddo
      
      call MPI_Barrier (MPI_Comm_World, IERR)
      
c     Find the maximum amaxd0 from all the processes:
      call MPI_Reduce(amaxd(3), amaxd0, 1, MPI_REAL8, MPI_MAX, nroot, 
     * comm_cart, IERR)
      call MPI_Bcast(amaxd0, 1, MPI_REAL8, nroot, comm_cart, IERR)
      
      
      
      deallocate(omj,omjx,omjy)
      
      
        
      return
      endsubroutine
      
      


      endmodule getrhs
