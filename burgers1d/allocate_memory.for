c     def of parameter
      Module allocate_memory
      use cons

      contains

      subroutine assign_memory
      
      
      allocate(x(0:nx+1), xx(1:nqua,0:nx+1))
      allocate(uc(1:nqua,0:nx+1,0:md))
      allocate(rhs(1:nqua,0:nx+1,0:md))
      allocate(rhsl(1:nqua,0:nx+1,0:md), rhsn(1:nqua,0:nx+1,0:md) )
      

      return
      end subroutine
      
      

      subroutine delete_memory
      
      deallocate(x, xx)
      deallocate(uc)
      deallocate(rhs, rhsl, rhsn)
      
      return
      end subroutine
      
      
	endmodule allocate_memory
