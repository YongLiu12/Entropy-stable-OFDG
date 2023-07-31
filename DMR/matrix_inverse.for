      module matrix_inverse
      use cons
      
      contains
      
      
      subroutine get_mat_inv(mat,mat_inv,nd)
      integer :: nd, k
      real :: mat(1:nd,1:nd), mat_inv(1:nd,1:nd)
      real :: fs(1:nd)
      
      do k = 1, nd
          
        fs = 0.0
        fs(k) = 1.0
        call linearsolver(mat,fs,mat_inv(k,1:nd),nd)
        
      enddo
      
      return
      endsubroutine
      
      
      
*****************   linear solver     *******************
      subroutine linearsolver(mat,b,f,nd)
      integer :: nd, info
      real, dimension(1:nd) :: ffs, b, f
      real :: mat(1:nd,1:nd)
      integer, allocatable, dimension(:) :: ipiv
      real, allocatable, dimension(:,:) :: amtx
      
      allocate(amtx(1:nd,1:nd), ipiv(1:nd))
      
      amtx = mat
      
      ffs = b

      call dgetrf(nd,nd,amtx,nd,ipiv,info)

      call dgetrs('N',nd,1,amtx(1:nd,1:nd),nd,ipiv,ffs(1:nd),nd,info)
      
      f = ffs
      
      deallocate(amtx,ipiv)
      
	return
      endsubroutine

      
      
      endmodule matrix_inverse
      
      
      
      
      
      
      
      
      
      
