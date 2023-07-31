      module boundary
	use cons
      use exactsolu
      
      contains
      
c	boundary condition

      subroutine bc(m) 

c set up the periodic boundary condition

      uc(:,0,m) = uc(:,nx,m)
      uc(:,nx+1,m) = uc(:,1,m)

      return
      endsubroutine

      endmodule
      
