c     def of parameter
      Module cons
c      implicit none
      
      parameter(md = 10, mp = 2, nqua = 4)
      parameter(ic = 2, ic_ent = 2)
      parameter(epweno = 1.0e-4, eps = 1.0e-12)
      
      real, dimension(0:md) :: bl, br, bi, ai, bm
      real, dimension(1:2,1:nqua,0:md) :: bdy_dx
      real :: xg(1:md), wg(1:md)
      real :: pv(1:200,0:md)
      
      real, allocatable, dimension(:) :: x
      real, allocatable, dimension(:,:) :: xx
      real, allocatable, dimension(:,:,:) :: uc, rhs, rhsl, rhsn
      
      
      real :: pi, cfl, dxuni, cdx, tend, t, dt, trk, xleft, xright
      real :: er1old, er2old, eriold, amax, amaxd(0:md)
      integer :: kcmax, kcount, kflux, nx, istop
      
      
      real :: Mass(nqua,nqua), Stiff(nqua,nqua), Diff(nqua,nqua)
      real :: Bdy_gm(nqua,nqua), Vand(nqua,mp+1), Diffhat(mp+1,mp+1)
      real :: Ext_gm(2,nqua,nqua), Mhat(mp+1,mp+1)
      real :: Idmat(nqua,nqua),Mass_inv(nqua,nqua),Mhat_inv(mp+1,mp+1)

c  Pk(:,:,m) is the L^2 projection into V^m, Vand_gm(2,:) is the two-point 
c  boundary values of basis functions, R_gmk is the two-point values of projecting 
c  the volume quadrature points to boundary quadrature points
      real :: Pk(mp+1,nqua,-1:mp), Vand_gm(2,mp+1), R_gmk(2,nqua)
      
c  Dhatk(:,:,m) is the differential matrix of m-th derivative
      real :: Dhatk(mp+1,mp+1,0:mp)
      
      
	endmodule cons
