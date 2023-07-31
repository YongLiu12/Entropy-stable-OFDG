c  main program: U_t + F(U)_x + G(U)_y = 0 
      
      program main
	use cons
      use para
      use allocate_memory
      use initial
      use setdt
      use getrhs
      use rungekutta
      use output
      
      implicit none
      include 'mpif.h'
      
      integer Status(MPI_STATUS_SIZE)
      
c use discontinuous Galerkin, polynomial truncation for the nonlinear
c terms, to compute 2D conservation laws.


c First test the compressible Euler equations U_t + F(U)_x + G(U)_y = 0

c set up the problem and the initial condition

      integer :: p, ir, i, j, kk, ierror, istop0
      real :: t_start, t_end, tmp
      real :: abs, min
      real :: entlocal, reslocal, tot_ent, tot_res
      
      
      call MPI_init      (                             IERR)
      call MPI_comm_rank (MPI_Comm_World, myid       , IERR)
      call MPI_comm_size (MPI_Comm_World, N_Processor, IERR)
    

      comm = MPI_Comm_World
      
      if(N_Processor .ne. nproc) then
      write(60,*) 'Need ',nproc,' nodes!!!'
      stop 922
      end if
      
      Nx_Processor = 32
      Ny_Processor = 16
      
      call cpu_time(t_start)
      
      call opens
      call initdata
      call plotexac  
      
      if (myid .eq. nroot) then
      write(60,"(a12,f12.6)") 'final time=', tend
      write(60,"(a8,f5.2,a16,f5.2)") 'xleft = ', xleft,
     *'xright = ', xright
      write(60,"(a8,f5.2,a16,f5.2)") 'yleft = ', yleft,
     *'yright = ', yright
      endif 
      
      
	do p = 5, 5
      
	    nx = 26*2**p
          ny = 8*2**p
          
        if (myid .eq. nroot) then
        write(60,"(a8,1I6,a8,1I6)") 'nx= ',nx,'ny= ',ny
        end if
        
        call assign_memory
        
        call init
        
        ir = 0
        
100   continue
        
        uc_old = uc(:,:,:,:,0)
	
      if(kcount .ne. 0) call set_dt
      if(t .ge. tend-1.0e-14) goto 600
      if(kcount .ge. kcmax .or. istop .eq. 0) goto 600

	if(t < tend) then

101   continue
      
c      trk = t + (7.0 - 3.0*ir)*ir/4.0*dt
      
	  call res(ir)
        
        call rk4(ir)
        
	    ir = mod(ir+1,4)
c
        
	  if(ir .ne. 0) goto 101

	endif

	  t = t + dt
        kcount = kcount + 1
        
c  Output the residue and the entropy
        call out_ent(p,entlocal)
        call out_res(p,reslocal)
        
      call MPI_Barrier(MPI_Comm_World, IERR)
        
      call MPI_Reduce(entlocal, tot_ent, 1, MPI_REAL8, MPI_SUM
     *  , nroot, comm_cart, IERR)
      call MPI_Reduce(reslocal, tot_res, 1, MPI_REAL8, MPI_SUM
     *  , nroot, comm_cart, IERR)
       
c       
       if(Myid .eq. nroot .and. kcount .ge. 2) then
           
c  compute the residue
        tot_res = tot_res*dx*dy*0.25
        write(5000+p,*) t, tot_res
           
c  compute the entropy
        tot_ent = tot_ent*dx*dy*0.25
        write(6000+p,*) t, tot_ent
       
       endif
c     
      if(mod(kcount,100*(mp+1)) .eq. 0) then

         call MPI_Barrier (MPI_Comm_World, IERR)
         
          kk = kcount/(100*(mp+1))
          call outp_data(989+kk)
     
      endif
        
        
      if(myid .eq. nroot) then
          
        if(mod(kcount,100*(mp+1)) .eq. 0) then
            
          write(60,*) 'iteration = ', kcount, amax0, amay0, amaxd0
          write(60,*) 'current time = ', t, 'dt = ', dt 
          call cpu_time(tmp)
          write(60,*) 'CPU time = ', tmp - t_start, 's'
          
        endif
        
      endif
        
        
        goto 100

600   continue

c print and compute errors

          if(Myid .eq. nroot) then
              write(60,*) "Completed ESOFDG!"
          endif


        call outp(p)
        
        call delete_memory
      
      enddo      

      call closes
      
      stop
      endprogram

