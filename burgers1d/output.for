      module output
	use cons
      use exactsolu
      use para
      use boundary
      
      contains

c     L1&L2&L_inf error

      subroutine outp(m,mm)
      integer :: m, mm
	real :: err, uex, ucom, xr, a(0:md), b(0:md)
      real :: perc, vv(0:md,-1:1)

      write(6,*) 'Current time is', t
	write(*,"(a4,1I4)") 'nx=',nx
	write(11+m,*) 'Variables="x""u"'
      write(11+m,"(a16,I5,a3)") 'zone T="nx=',nx,'"'
      write(11+m,*) 'I=',nx!*nqua
	do i = 1, nx
!c
!        do j = 1, nqua
!          ucom = dot_product(uc(0:mp,i,0),pv(j,0:mp))
!          write(11+m,11) xx(j,i), ucom
!        enddo
          write(11+m,11) x(i), sum(Pk(1,:,mp)*uc(:,i,0))
          
      enddo
      
	er1 = 0.0
	er2 = 0.0
	eri = 0.0

c  error for u - u_h
	do i = 1, nx
        do j = 1, 200
            
          xr = x(i) + 1.0*(j - 100.0)/200.0*dxuni
          uex = 2.0*burgersexac(xr - t,2.0*t) + 1.0
c          uex = burgersexac(xx(j,i),t)
          ucom = sum( pv(j,0:mp) * 
     *    matmul(Pk(1:mp+1,1:nqua,mp),uc(1:nqua,i,0)) )
	    err = abs( uex - ucom )
	    eri = max(eri, err)
          er1 = er1 + err/200.0
          er2 = er2 + err**2.0/200.0
        enddo
      
      enddo

	er1 = er1/nx
	er2 = sqrt(er2/nx)

	o1=log(er1old/er1)/log(2.0)
	o2=log(er2old/er2)/log(2.0)
	oi=log(eriold/eri)/log(2.0)
c	write(*,"(a12,es11.3,a21,f8.3)") 'L1_error= ',
c     *er1,'L1_order=',o1
	write(*,"(a12,es11.3,a21,f8.3)") 'L2_error= ',
     *er2,'L2_order=',o2
c	write(*,"(a12,es11.3,a21,f8.3)") 'Linf_error= ',
c     *eri,'Linf_order=',oi
      
	er1old = er1
	er2old = er2
	eriold = eri
      
11    format(6e13.5)   

      
	return  
      endsubroutine      

      subroutine out_ent(p,ent)
      integer :: i, k, p
      real :: ent
      
       ent = 0.0
       do i = 1, nx
         do k = 1, nqua
           ent = ent + wg(k)*entropy(uc(k,i,0))
         enddo
       enddo
       ent = ent*dxuni*0.5
       
       return
      endsubroutine

      subroutine out_exac_ent(p,ent)
      integer :: i, k, p
      real :: ent, uex
      
       ent = 0.0
       do i = 1, nx
         do k = 1, nqua
           uex = 2.0*burgersexac(xx(k,i) - t,2.0*t) + 1.0
           ent = ent + wg(k)*entropy(uex)
         enddo
       enddo
       ent = ent*dxuni*0.5
       
       return
      endsubroutine
      

      subroutine opens
	open(10,file='exac.plt')
c      open(11,file='nx=10.plt')
	open(12,file='nx=20.plt')
	open(13,file='nx=40.plt')
	open(14,file='nx=80.plt')
	open(15,file='nx=160.plt')
	open(16,file='nx=320.plt')
	open(17,file='nx=640.plt')
      return
	endsubroutine

      subroutine closes
	close(10)
c	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	close(17)
      return
      endsubroutine
      
      subroutine plotexac(m)
      
	write(10,*) 'Variables="x""u"'
      write(10,"(a16,I5,a5)") 'zone T="exac"'
      write(10,*) 'I=',801
	do i = 0, 800
          xe = xleft+(xright-xleft)*i/800.0
          s2 = 2.0*burgersexac(xe-tend,2.0*tend) + 1.0        
         write(10,12) xe, s2
      enddo 

12    format(6e13.5)
      
      return
      endsubroutine
      
      endmodule output
      
