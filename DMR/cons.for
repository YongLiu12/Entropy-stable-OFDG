c     def of parameter
      Module cons
      
c  mp is the degree of the polynomial and mp2=(mp+1)*(mp+2)/2
      parameter(md = 30, mp = 2, mp2 = 6)
c  nqua is number of the boundary quadrature points and 
c  nqua2 = nqua*nqua
      parameter(nqua = 4, nqua2 = 16)
      parameter(ic = 3, ic_ent = 2)
      parameter(epweno = 1.0e-4, eps = 1.0e-12)
      parameter(nequ = 4)
c  Nproc: the number of cores; nroot is the root processor
      parameter(nproc = 512, itag1 = 96, itag2 = 97, itag3 = 98
     *        , itag4 = 99, nroot = 0)
      
      real :: ai(0:md)
      real, dimension(nqua,mp2,4) :: bdy
      real, dimension(nqua,nqua2,0:md,4) :: bdy_dx
      real :: xg(1:md), wg(1:md)
      real :: xg2(1:2,1:nqua2), wg2(1:nqua2)
      real :: xgB(1:2,1:nqua,1:4)
      real :: pv(nqua2,0:md)
      
      real, allocatable, dimension(:) :: x, y
      real, allocatable, dimension(:,:,:,:) :: uc_old
      real, allocatable, dimension(:,:,:,:,:) :: uc, rhs
     *    , rhsl, rhsn
      
      
      real :: pi, cfl, dx, dy, cdx, cdy, tend, t, dt, trk
      real :: xleft, xright, yleft, yright 
      real :: amax, amay, amax0, amay0, amaxd(0:md), amaxd0
      integer :: mp2, nqua2, kcmax, kcount, kflux, nx, ny, istop
      real :: gamma, gamma1
      
      real, dimension(1:nequ) :: er1, er2, eri
      real, dimension(1:nequ) :: er1old, er2old, eriold

      real :: acoef(0:10,0:10)
      real :: gmrk, bt1rk, bt2rk, alf1rk, alf2rk
      
c  Matrices and Vectors in Nodal DGM
      real, dimension(nqua2,nqua2) :: Mass, Mass_inv
      real, dimension(mp2,mp2) :: Mhat, Mhat_inv
      real :: Vand(nqua2,mp2)
      real, dimension(nqua2,nqua2) :: Stiff2x, Stiff2y
      real, dimension(nqua2,nqua2) :: Diff2x, Diff2y
      real :: Diffhat(mp+1,mp+1)
      real, dimension(mp2,mp2) :: Diffhat2x, Diffhat2y
      real :: Idmat(nqua2,nqua2)

c  Pk(:,:,m) is the L^2 projection into V^m
      real :: Pk(mp2,nqua2,-1:mp)

c  Bdy_gm is the diagonal matrix with quadrature weights in the diagonal line
c  Vand_gm(:,:,4) is the quadrature point values on four boundaries 
c  of basis functions, R_gmk is the extrapolation matrix projecting 
c  the volume quadrature points to boundary quadrature points, Ext_gm 
c  is to approximate the boundary integral
      real :: Bdy_gm(nqua,nqua), Vand_gm(nqua,mp2,4)
     *    , R_gmk(nqua,nqua2,4), Ext_gm(nqua2,nqua2,4)
      
c  Dhatk(:,:,m) is the differential matrix of m-th derivative in 1D
c  Dhatk2(:,:,m) is the differential matrix of (m1,m2)-th derivative in 2D
      real :: Dhatk(mp+1,mp+1,mp), Dhatk2(mp2,mp2,mp2)
      
c  ibpq(:,4) is to connect the boundary quadrature point and
c  the volume quadrature point
      integer :: ibpq(nqua,4)
c  nbr(:,4) is to connect the boundary quadrature point and 
c  its neighborhood
      integer :: nbr(nqua,4)
c  lrnbr(:), btnbr(:) are to connect the quadrature point between 
c  horizontal and vertical neighboring cells 
      integer :: lrnbr(nqua2), btnbr(nqua2)
      
      
      
c  MPI data   
      integer :: Myid, N_Processor
      integer :: Nx_Processor, Ny_Processor
      integer :: N_First, N_Last
      integer :: left, right, top, bottom
      integer :: nxlocal, nylocal, ntmpx, ntmpy
      integer :: xcoord, ycoord
      integer :: comm_cart
      integer :: comm, comm_world, IERR
      
      integer, dimension (2) :: dims
      integer, dimension (2) :: coords
      logical, dimension (2) :: periods
      logical reorder
      
      real, allocatable, dimension(:) :: xsdrv_data,ysdrv_data 
      
      
      
	endmodule cons
