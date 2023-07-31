c  Set dt
      module setdt
	  use cons
      use para
      use exactsolu
      
      contains
      
      subroutine set_dt
      
      implicit none
      include 'mpif.h'
      integer :: i, j, k, ierror
      real :: s1, s2, amaxlocal, amaylocal
      real :: den, xmt, ymt, eng, t0, vex, vey, pre, cvel
      real :: abs, max
	
      amaxlocal = 0.0
      amaylocal = 0.0
      
      do j = 1, nylocal
      do i = 1, nxlocal
        do k = 1, nqua2
          den = uc(k,1,i,j,0)
          xmt = uc(k,2,i,j,0)
          ymt = uc(k,3,i,j,0)
          eng = uc(k,4,i,j,0)
          t0 = 1.0/den
          vex = xmt*t0
          vey = ymt*t0
          pre = press(uc(k,1:nequ,i,j,0))
          cvel = sqrt( gamma*abs(pre*t0) )
          amaxlocal = max( amaxlocal,abs(vex)+cvel )
          amaylocal = max( amaylocal,abs(vey)+cvel )
        enddo
      enddo
      enddo
      
      amaxlocal = amaxlocal * 1.1
      amaylocal = amaylocal * 1.1
      
      call MPI_Barrier (MPI_Comm_World, IERR)

!     Find the maximum amax,amay,amaxd(0) from all the processes:
      call MPI_Reduce(amaxlocal, amax, 1, MPI_REAL8, MPI_MAX,nroot, 
     * comm_cart, IERR)
      call MPI_Reduce(amaylocal, amay, 1, MPI_REAL8, MPI_MAX,nroot,
     * comm_cart, IERR)

!     Send the maximum amax0,amay0,amaxd0 found above to all the processes for using:
      call MPI_Bcast(amax,1,MPI_REAL8,nroot,comm_cart,IERR)
      call MPI_Bcast(amay,1,MPI_REAL8,nroot,comm_cart,IERR)
      
      
      dt = cfl/(amax * cdx + amay*cdy + 0.1*amaxd0*cdx)
      
      
c  dt is problematic and stop the simulation
      if( isNaN(dt) ) then
        
          write(60+myid,*) 'dt is NaN' 
      
          stop
        
      endif
          
          
      if(tend <= t+dt .and. tend > t) then 
          
          dt = tend - t
          
      endif   
        
      if(dt .le. eps .and. t+dt .lt. tend - eps) then 
          
          write(60+myid,*) dt, 'dt is too small' 
          istop = 0
          
      endif

      rhs = 0.0
      rhsl = 0.0
      rhsn = 0.0
      amaxd = 0.0
      
      
          
      return
      endsubroutine
      
      
      endmodule setdt

