c     main program
c	u_t + uu_x=0 with periodic boundary condition
      
      program main
	use cons
      use para
      use allocate_memory
      use initial
      use setdt
      use getrhs
      use rungekutta
      use output
      
c use discontinuous Galerkin, polynomial truncation for the nonlinear
c terms, to compute 1D conservation laws.


c First test the nonlinear problem u_t+uu_x=0
c the initial solution is 2sin(x) + 1


c set up the problem and the initial condition

	integer :: p, ir
      real :: tot_ent
      
      call opens
      call initdata
      call plotexac(ic)
      
	do p = 4, 4
	  nx = 8*2**p
        
        call assign_memory
        
        call init
      
        ir = 0
	
100   if(kcount .ne. 0) call set_dt
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
        
        if(mod(kcount,5) .eq. 0) then
c  compute the entropy of the numerical solution
            call out_ent(p,tot_ent)
            write(5000+p,'(2es10.3)') t, tot_ent
c  compute the entropy of the exact solution
            call out_exac_ent(p,tot_ent)
            write(6000+p,'(2es10.3)') t, tot_ent
        endif
        
        if(mod(kcount,10000) .eq. 0) then 
            write(*,*) t, dt, kcount
        endif
        
        
        goto 100

600   continue

c print and compute errors


        call outp(p,ic)
        
        call delete_memory
      
      enddo      

      call closes
      
      stop
      endprogram

