      module getrhs
	  use cons
      use boundary
      use exactsolu
      
      contains
      
!**********************************************************************
!     the numerical flux for u_t + uu_x = 0 
!**********************************************************************
      subroutine res(m)

      integer :: i, k, k1, k2
      real :: us(1:nqua), f_ent(1:nqua,1:nqua)
      real, dimension(0:md) :: a, b, scm_jmp, flx
      
      real, allocatable, dimension(:,:) :: omj
      real :: vv(0:md,-1:1), sigma(0:md)
      real :: upos, uneg
      real :: s1, s2, alf
      

c compute the residue
        
c  Periodic B.C.
        
      call bc(m)

      amax = 0.0
      do i = 1, nx
        do k = 1, nqua
          amax = max( amax, abs(uc(k,i,m)) )
        enddo
      enddo
      amax = amax * 1.1

      do  i = 1, nx
        
      do k = 1, nqua
        us(k) = uc(k,i,m)
      enddo
      
      do k2 = 1, nqua
      do k1 = 1, nqua
        f_ent(k1,k2) = fmS(us(k1),us(k2)) 
      enddo
      enddo


      do k = 1, nqua

c compute the contribution from the Difference term

        rhs(k,i,m) = - 2.0*sum(Diff(k,1:nqua)*f_ent(k,1:nqua))

      enddo

c compute the contribution from the simultaneous approximation term

      upos = uc(1,i,m)
      uneg = uc(nqua,i-1,m)
      alf = 0.0
      s1 = sum(Pk(1,1:nqua,0)*uc(1:nqua,i,m))
      s2 = sum(Pk(1,1:nqua,0)*uc(1:nqua,i-1,m))
      alf = max(abs(s1), abs(s2)) * 1.1
      flx(0) = 0.5*(0.5*upos*upos + 0.5*uneg*uneg - alf*(upos - uneg))
      upos = uc(1,i+1,m)
      uneg = uc(nqua,i,m)
      alf = 0.0
      s1 = sum(Pk(1,1:nqua,0)*uc(1:nqua,i,m))
      s2 = sum(Pk(1,1:nqua,0)*uc(1:nqua,i+1,m))
      alf = max(abs(s1), abs(s2)) * 1.1
      flx(1) = 0.5*(0.5*upos*upos + 0.5*uneg*uneg - alf*(upos - uneg))
        rhs(1,i,m) = rhs(1,i,m) - (0.5*us(1)**2.0 - flx(0))/wg(1)
        rhs(nqua,i,m) = rhs(nqua,i,m) 
     *        + (0.5*us(nqua)**2.0 - flx(1))/wg(nqua)

      enddo
        
c take care of the mass matrix

        rhs(1:nqua,1:nx,m) = 2.0*rhs(1:nqua,1:nx,m)*cdx
      

        
c  Add the numerical damping
        
c compute the jump as the coefficients 
        
        allocate(omj(0:mp,0:nx+1))
        vv = 0.0
        
      do i = 1, nx+1
        
c  compute the jumps of each order
        do k = 0, 1
        upos = sum(uc(:,i,m)*bdy_dx(1,:,k))
        uneg = sum(uc(:,i-1,m)*bdy_dx(2,:,k))
        omj(k,i) = (upos - uneg) 
        enddo
      
	enddo
      
c  add the numerical damping to the scheme  

      do i = 1, nx
          
      sigma(0:mp) = omj(0:mp,i)**2.0 + omj(0:mp,i+1)**2.0 
      
      a = 0.0
      do k = 0, 1
        do k1 = 0, k
        a(k) = a(k) + sigma(k1) / (1.0 + k1)
        enddo
      enddo

      a(1) = a(1)**0.5
      
      us(1:nqua) = uc(:,i,m)
c  Pk.b are the coefficients, V.Pk.b are the quadrature point values
      do k = 1, 1
        rhs(1:nqua,i,m) = rhs(1:nqua,i,m) - a(1)*cdx 
     *    *( us(1:nqua) - matmul( Vand(1:nqua,1:k)
     *    ,matmul(Pk(1:k,1:nqua,k-2),us(1:nqua)) ) )
      enddo
        
      
      enddo

      
      deallocate(omj)
        
      return
      endsubroutine
      
      
      


      endmodule getrhs
