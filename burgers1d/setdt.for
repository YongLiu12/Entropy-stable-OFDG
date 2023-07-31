c  Set dt
      module setdt
	use cons
      use para
      
      contains
      
      subroutine set_dt
      real :: s1
	
c      amax = 1.1
      

      amax = 0.0
      do i = 1, nx
        do j = 1, nqua
          amax = max( amax, abs(uc(j,i,0)) )
        enddo
      enddo
      amax = amax * 1.1
      
      dt = cfl/amax * dxuni
c      dt = 0.01*dxuni
          
      if(tend <= t+dt .and. tend > t) then 
          
          dt = tend - t
          
      endif   
        
      if(dt .le. eps .and. t+dt .lt. tend - eps) then 
          
          write(*,*) dt, 'dt is too small' 
          istop = 0
          
      endif

      rhs = 0.0
      rhsl = 0.0
      rhsn = 0.0
      amaxd = 0.0
      
      
          
      return
      endsubroutine
      
      
      
      
      
      endmodule setdt

