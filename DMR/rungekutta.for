      module rungekutta
      use cons
      use exactsolu
      
      contains
      
      
c     4th order Runge-Kutta
      subroutine rk4(m)
          
        if(m .eq. 0) then
            
        uc(:,:,1:nxlocal,1:nylocal,1) = 
     *    uc(:,:,1:nxlocal,1:nylocal,0)
     *    + 0.5*dt*rhs(:,:,1:nxlocal,1:nylocal,0) 
          

        elseif(m .eq. 1) then
            
        uc(:,:,1:nxlocal,1:nylocal,2) = 
     *    uc(:,:,1:nxlocal,1:nylocal,0) 
     *    + 0.5*dt*rhs(:,:,1:nxlocal,1:nylocal,1)
        
        elseif(m .eq. 2) then
            
        uc(:,:,1:nxlocal,1:nylocal,3) = 
     *    uc(:,:,1:nxlocal,1:nylocal,0) 
     *    + dt*rhs(:,:,1:nxlocal,1:nylocal,2)
        
        elseif(m .eq. 3) then
            
        uc(:,:,1:nxlocal,1:nylocal,0) = 
     *    - uc(:,:,1:nxlocal,1:nylocal,0)/3.0
     *    + uc(:,:,1:nxlocal,1:nylocal,1)/3.0
     *    + 2.0/3.0*uc(:,:,1:nxlocal,1:nylocal,2) 
     *    + uc(:,:,1:nxlocal,1:nylocal,3)/3.0
     *    + dt*rhs(:,:,1:nxlocal,1:nylocal,3)/6.0
     
      endif
      
      
	  return
      endsubroutine
      
      
      
      endmodule rungekutta
      
      
      
      
      
      
      
      
      
      
      
      
