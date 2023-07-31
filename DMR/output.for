      module output
	use cons
      use exactsolu
      use para
      use boundary
      
      contains

c     L1&L2&L_inf error

      subroutine outp(m)
      
      implicit none
      include 'mpif.h'
      integer :: m, i, j, k, kk, ierror, iunit, coordstmp(2)
      integer :: ix, jx, irank
      
	real, dimension(1:nequ) :: err, er1local, er2local, erilocal
      real, dimension(1:nequ) :: o1, o2, oi
      real, dimension(1:nequ) :: uex, ucom
      real :: xr, yr, ts, a(0:md), b(0:md)
      real :: perc, vv(0:md,-1:1)
      real :: s1, s2
      
	er1 = 0.0
	er2 = 0.0
	eri = 0.0
      
      er1local = 0.0
      er2local = 0.0
      erilocal = 0.0
      
      if(ic .le. 2) then

c  error for u - u_h
      do j = 1, nylocal
	do i = 1, nxlocal
          
        do k = 1, nqua2

          s1 = x(i) + 0.5*xg2(1,k)*dx
          s2 = y(j) + 0.5*xg2(2,k)*dy
          
          if(ic .eq. 1) then
          uex = u0(s1 - 0.7*t,s2 - 0.3*t)
          elseif(ic .eq. 2) then
          ts = t - 1.0*floor(t)
          uex = u0(s1 - ts,s2 - ts)
          endif
          ucom = uc(k,1:nequ,i,j,0)
	    err = abs( uex - ucom )
	    erilocal = max(erilocal, err)
	    er1local = er1local + err*wg2(k)
	    er2local = er2local + err**2.0*wg2(k)
          
        enddo

      enddo
      enddo

      call MPI_Barrier (MPI_Comm_World, IERROR)
      
      call MPI_Reduce(er1local, er1, nequ, 
     * MPI_REAL8, MPI_SUM, nroot, comm_cart, IERROR)
      
 	call MPI_Reduce(er2local, er2, nequ, 
     * MPI_REAL8, MPI_SUM, nroot, comm_cart, IERROR)

      call MPI_Reduce(erilocal, eri, nequ, 
     * MPI_REAL8, MPI_MAX, nroot, comm_cart, IERROR)

      if (myid .eq. nroot) then
      
	er1 = er1/nx/ny
	er2 = sqrt(er2/nx/ny)
      
	o1 = log(er1old/er1)/log(2.0)
	o2 = log(er2old/er2)/log(2.0)
	oi = log(eriold/eri)/log(2.0)
        write(3000,*) 'nx=',nx,'ny=',ny,'tend=',t
        write(3000,*) 'density'
	  write(3000,"(a12,es12.3,a21,f9.3)") 'L1_error= '
     *    , er1(1),'L1_order=', o1(1)
	  write(3000,"(a12,es12.3,a21,f9.3)") 'L2_error= '
     *    , er2(1),'L2_order=', o2(1)
	  write(3000,"(a12,es12.3,a21,f9.3)") 'Linf_error= '
     *    , eri(1),'Linf_order=', oi(1)
      
      
	er1old = er1
	er2old = er2
	eriold = eri
      
      endif
      
      call outp_data(m)
      
      else
          
        call outp_data(m)
        
      
      endif
      
      
	return  
      endsubroutine  
      
      
      subroutine outp_data(m)
      
      implicit none
      include 'mpif.h'
      integer :: m, i, j, k, kk, ierror, iunit, coordstmp(2)
      integer :: ix, jx, irank
	real, dimension(1:nequ) :: err, uex, ucom
      real :: den, pre, vex, vey
	real :: xr, yr, a(0:6)
      real, allocatable, dimension (:,:,:) :: tmpucom
         
      Iunit = 4000 + myid
         open(Iunit,status='replace')
         
      do j = 1, nylocal
	  do i = 1, nxlocal
          do k = 1, nequ
          ucom(k) =  sum(Pk(1,:,mp)*uc(:,k,i,j,0))
          enddo
          pre = press(ucom)
          vex = ucom(2) / ucom(1)
          vey = ucom(3) / ucom(1)
          write(Iunit,"(4es13.3)") ucom(1), vex, vey, pre
        enddo
      enddo
      
        close(Iunit) 
        
         call MPI_Barrier(MPI_Comm_World, IERROR)
        
      if (myid .eq. nroot) then
      allocate(tmpucom(1:nx,1:ny,1:nequ))

        do ix = 1, Nx_Processor
        do jx = 1, Ny_Processor
          !write(*,*) 'Nx_Processor',Nx_Processor
          !write(*,*) 'Ny_Processor',Ny_Processor
          coordstmp(1) = ix-1
          coordstmp(2) = jx-1
          !write(*,*) coordstmp
          call MPI_Cart_Rank(comm_cart,coordstmp,irank,ierror)

          Iunit = 4000 + irank
          open(Iunit)
             
          do j = 1, nylocal
            do i = 1, nxlocal
             read(Iunit,"(4es13.3)") 
     *       tmpucom((ix-1)*nxlocal+i,(jx-1)*nylocal+j,1:nequ)
            enddo
          enddo
             
          close(Iunit)
             
        enddo
        enddo
                   
          
	write(11+m,*) 'Variables= "x" "y" "rho" "ux" "uy" "pre"'
c      write(11+m,"(a16,I5,a3,a8,I5,a3)") 'zone T="nx=',nx, 'ny=', ny,'"'
      write(11+m,'(a10,a5,I5,a5,I5)') 'zone  ', 'I=',nx, 'J=', ny
      
      do j = 1, ny
	  do i = 1, nx
            
          xr = xleft + dx*(i - 0.5)
          yr = yleft + dy*(j - 0.5)
          write(11+m,"(6es13.3)") xr, yr, tmpucom(i,j,1:nequ)
          
        enddo
      enddo 
!c
!      write(7000+m,*) 'Variables="x""u"'
!      write(7000+m,*) 'zone I=', nx
!      do i = 1, nx
!        xr = xleft + dx*(i - 0.5)
!        write(7000+m,*) xr, tmpucom(i,1,1)
!      enddo
      if( abs(t - tend) .le. 1.0e-12) then 
          
	write(21+m,*) 'Variables="x" "density" "ux" "uy" "pre"'
      write(21+m,'(a10,a5,I5)') 'zone  ', 'T= "256*256 " '
      do i = 1, nx
          xr = xleft + dx*(i - 0.5)
          write(21+m,"(5es13.3)") xr, tmpucom(i,1,1:nequ)
      enddo
      
      endif
      
      deallocate(tmpucom)  
      
      endif
      
      
      return
      endsubroutine


      subroutine out_ent(p,ent)
      implicit none
      integer :: i, j, k, p
      real :: ent
      
       ent = 0.0
      do j = 1, nylocal
      do i = 1, nxlocal
        do k = 1, nqua2
        ent = ent 
     *    + wg2(k)*sum(var_ent(uc(k,:,i,j,0))*uc(k,:,i,j,0))
        enddo
      enddo
      enddo
       
       
       return
      endsubroutine


      subroutine out_res(p,residue)
      implicit none
      integer :: i, j, k, p
      real :: residue
      
        residue = 0.0
        
        do j = 1, nylocal
        do i = 1, nxlocal
          do k = 1, nequ
        residue = residue 
     *    + abs( sum(wg2*(uc(:,k,i,j,0) - uc_old(:,k,i,j))) )/dt
          enddo
        enddo
        enddo
        

        if( isNaN(residue) ) then
        
          write(60+myid,*) 'Residue: NaN occurs. Myid=', myid 
      
          stop
        
        endif
c        
       
       
       return
      endsubroutine
      

      subroutine opens
      implicit none
      
      if (myid .eq. nroot) then

	open(10,file='exac.plt')
	open(11,file='16_16.plt')
	open(12,file='32_32.plt')
	open(13,file='64_64.plt')
	open(14,file='128_128.plt')
	open(15,file='256_256.plt')
      open(16,file='512_512.plt')
      
      end if
      
      
      return
	endsubroutine

      subroutine closes
      
      implicit none
      
      if (myid .eq. nroot) then
	close(10)
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
      close(16)
      
      end if
      
      
      return
      endsubroutine

c
      subroutine plotexac
      implicit none
      integer :: i
      real :: xe, s1(1:nequ)
      
	write(10,*) 'Variables="x""u"'
      write(10,"(a16,I5,a5)") 'zone T="exac"'
      write(10,*) 'I=',801
	do i = 0, 800
          xe = xleft + (xright-xleft)*i/800.0
          s1 = u0(xe-0.7*tend,0.0-0.3*tend)        
         write(10,12) xe, s1(1)
      enddo 

12    format(6e13.5)
      
      return
      endsubroutine
      

      endmodule output
      
