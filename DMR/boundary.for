      module boundary
	use cons
      use exactsolu
      
      contains
      
c	boundary condition

      subroutine bc(m) 
      
      implicit none
      include 'mpif.h'
      integer :: kk, err_code
      integer :: topright, topleft, botright, botleft
      integer :: Status(MPI_STATUS_SIZE)
      integer :: k, i, j, m, ierror, m1, m2
      integer :: k1, k2, kyg
      integer :: ns0, ns1
      real :: xr, yr
      real :: den, vex, vey, pre, phy(1:nequ), alf
      
      

c exchanging data 

c      write(700+myid,*) 'exchanging data started'

      call MPI_Barrier(MPI_Comm_World, IERR)
      
c      write(700+myid,*) myid, xcoord, ycoord

c  odd and even processor in x-direction
      if( mod(xcoord,2) .eq. 1 ) then
          
c  left local boundary 
c  save boundary data to ysdrv and send 
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        ysdrv_data(m1:m2) = uc(kk,k,1,1:nylocal,m)
      enddo
      enddo
          
      call MPI_Send(ysdrv_data, ntmpy, MPI_REAL8, left, itag1,
     * comm_cart, IERR)
      
c      write(500+myid,*) myid, xcoord, ycoord, left
      
c  from ysdrv to u
      call MPI_Recv(ysdrv_data, ntmpy, MPI_REAL8, left, itag2,
     * comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        uc(kk,k,0,1:nylocal,m) = ysdrv_data(m1:m2)
      enddo
      enddo
         
c  right local boundary
c  save boundary data to ysdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        ysdrv_data(m1:m2) = uc(kk,k,nxlocal,1:nylocal,m)
      enddo
      enddo
         
      call MPI_Send(ysdrv_data, ntmpy, MPI_REAL8, right, itag2,
     *    comm_cart,IERROR)
         
c  from ysdrv to u
      call MPI_Recv(ysdrv_data, ntmpy, MPI_REAL8, right, itag1,
     *    comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        uc(kk,k,nxlocal+1,1:nylocal,m) = ysdrv_data(m1:m2)
      enddo
      enddo
         
      
      else
      
          
c  right local boundary
c  from ysdrv to u
      call MPI_Recv(ysdrv_data, ntmpy, MPI_REAL8, right, itag1,
     * comm_cart, Status, IERR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        uc(kk,k,nxlocal+1,1:nylocal,m) = ysdrv_data(m1:m2)
      enddo
      enddo
      
c      write(500+myid,*) myid, xcoord, ycoord, right
      
c  save boundary data to ysdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        ysdrv_data(m1:m2) = uc(kk,k,nxlocal,1:nylocal,m)
      enddo
      enddo
          
      call MPI_Send(ysdrv_data, ntmpy, MPI_REAL8, right, itag2,
     * comm_cart, IERROR)
      
      
c  left local boundary
c  from ysdrv to u
      call MPI_Recv(ysdrv_data, ntmpy, MPI_REAL8, left, itag2,
     * comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        uc(kk,k,0,1:nylocal,m) = ysdrv_data(m1:m2)
      enddo
      enddo
      
c  save boundary data to ysdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nylocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nylocal
        ysdrv_data(m1:m2) = uc(kk,k,1,1:nylocal,m)
      enddo
      enddo
      
      call MPI_Send(ysdrv_data, ntmpy, MPI_REAL8, left, itag1, 
     *    comm_cart, IERROR)
    
    
      end if
      
      !call MPI_Barrier (MPI_Comm_World, IERR)
      !write(60+myid,*) myid,'finish x-loop'
      

c      write(800+myid,*) 'y-direction'
c      write(800+myid,*) myid, xcoord, ycoord 
      

c  odd and even processor in y-direction
      if (mod(ycoord,2) .eq. 1 ) then
          
c  bottom local boundary 
c  save boundary data to xsdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        xsdrv_data(m1:m2) = uc(kk,k,1:nxlocal,1,m)
      enddo
      enddo
          
      call MPI_Send(xsdrv_data, ntmpx, MPI_REAL8, bottom, itag3,
     * comm_cart, IERROR)
      
c      write(600+myid,*) xcoord, ycoord, bottom, myid
      
c  from xsdrv to u
      call MPI_Recv(xsdrv_data, ntmpx, MPI_REAL8, bottom, itag4,
     * comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        uc(kk,k,1:nxlocal,0,m) = xsdrv_data(m1:m2)
      enddo
      enddo
         
c  top local boundary
c  save boundary data to xsdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        xsdrv_data(m1:m2) = uc(kk,k,1:nxlocal,nylocal,m)
      enddo
      enddo
         
      call MPI_Send(xsdrv_data, ntmpx, MPI_REAL8, top, itag4,
     *    comm_cart,IERROR)
         
c  from xsdrv to u
      call MPI_Recv(xsdrv_data, ntmpx, MPI_REAL8, top, itag3,
     *    comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        uc(kk,k,1:nxlocal,nylocal+1,m) = xsdrv_data(m1:m2)
      enddo
      enddo
         
      
      else
      
          
c  top local boundary
c  from xsdrv to u
      call MPI_Recv(xsdrv_data, ntmpx, MPI_REAL8, top, itag3,
     * comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        uc(kk,k,1:nxlocal,nylocal+1,m) = xsdrv_data(m1:m2)
      enddo
      enddo
      
c      write(600+myid,*) xcoord, ycoord, top, myid
      
c  save boundary data to xsdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        xsdrv_data(m1:m2) = uc(kk,k,1:nxlocal,nylocal,m)
      enddo
      enddo
          
      call MPI_Send(xsdrv_data, ntmpx, MPI_REAL8, top, itag4,
     * comm_cart, IERROR)
      
      
c  bottom local boundary
c  from xsdrv to u
      call MPI_Recv(xsdrv_data, ntmpx, MPI_REAL8, bottom, itag4,
     * comm_cart, Status, IERROR)
      
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        uc(kk,k,1:nxlocal,0,m) = xsdrv_data(m1:m2)
      enddo
      enddo
      
c  save boundary data to xsdrv and send
      do k = 1, nequ
      do kk = 1, nqua2
        m1 = ( kk - 1 + (k - 1)*nqua2 )*nxlocal + 1
        m2 = ( kk + (k - 1)*nqua2 )*nxlocal
        xsdrv_data(m1:m2) = uc(kk,k,1:nxlocal,1,m)
      enddo
      enddo
      
      call MPI_Send(xsdrv_data, ntmpx, MPI_REAL8, bottom, itag3, 
     *    comm_cart, IERROR)
    
    
      end if

      
      
c set up the boundary condition
      
      if (xcoord .eq. Nx_Processor) then
******************   right boundary : outflow boundary   ********************
      uc(:,:,nxlocal+1,1:nylocal,m) 
     *    = uc(:,:,nxlocal,1:nylocal,m)
      endif

******  (1/6,0)-(0,0)-(0,1)-(x_0(t),1) : post-shock condition   ******
      den = 8.0
      vex = 8.25*cos(pi/6.0)
      vey = - 8.25*sin(pi/6.0)
      pre = 116.5
      phy(1) = den
      phy(2) = den*vex
      phy(3) = den*vey
      phy(4) = pre/gamma1 + 0.5*den*(vex*vex+vey*vey)
      
c  inflow B.C. (0,0) - (0,1)
      if (xcoord .eq. 1) then
      do k = 1, nequ
        uc(:,k,0,0:nylocal+1,m) = phy(k)
      enddo
      endif
      
      
c  inflow B.C. (0,0) - (1/6,0)
      if (ycoord .eq. 1) then
          
      ns0 = floor( 1.0/6.0*cdx )
      if (ns0 .ge. xcoord*nxlocal) then
          ns0 = nxlocal
      do k = 1, nequ
        uc(:,k,1:ns0,0,m) = phy(k)
      enddo
      
      else if(ns0 .gt. (xcoord-1)*nxlocal) then
          ns0 = ns0 - (xcoord - 1)*nxlocal
      
      do k = 1, nequ
        uc(:,k,1:ns0,0,m) = phy(k)
      enddo
        
      end if
      
      end if
      
      
c  inflow B.C. (0,1) - (x_0(t),1)
      if (ycoord .eq. Ny_Processor) then
          
      ns1 = floor( (1.0/6.0 + (1.0 + 20.0*t)/3.0**0.5)*cdx )
      if (ns1 .ge. xcoord*nxlocal) then
          ns1 = nxlocal
      do k = 1, nequ
        uc(:,k,1:ns1,nylocal+1,m) = phy(k)
      enddo
      
      else if (ns1 .gt. (xcoord - 1)*nxlocal) then
          ns1 = ns1 - (xcoord - 1)*nxlocal
      do k = 1, nequ
        uc(:,k,1:ns1,nylocal+1,m) = phy(k)
      enddo
      
      end if
      
      end if
      
**********   (1/6,0)-(4,0) : reflecting boundary condition   **********
      if (ycoord .eq. 1) then
          
          ns0 = floor( 1.0/6.0*cdx )
      if (ns0 .ge. (xcoord-1)*nxlocal 
     *    .and. ns0 .lt. xcoord*nxlocal) then
          ns0 = ns0 - (xcoord - 1)*nxlocal
      do k = 1, nqua2
      uc(k,1,ns0+1:nxlocal+1,0,m) =  uc(btnbr(k),1,ns0+1:nxlocal+1,1,m)
      uc(k,2,ns0+1:nxlocal+1,0,m) =  uc(btnbr(k),2,ns0+1:nxlocal+1,1,m)
      uc(k,3,ns0+1:nxlocal+1,0,m) = -uc(btnbr(k),3,ns0+1:nxlocal+1,1,m)
      uc(k,4,ns0+1:nxlocal+1,0,m) =  uc(btnbr(k),4,ns0+1:nxlocal+1,1,m)
      enddo
      elseif (ns0 .lt. (xcoord-1)*nxlocal) then
          ns0 = 0
      do k = 1, nqua2
      uc(k,1,0:nxlocal+1,0,m) =  uc(btnbr(k),1,0:nxlocal+1,1,m)
      uc(k,2,0:nxlocal+1,0,m) =  uc(btnbr(k),2,0:nxlocal+1,1,m)
      uc(k,3,0:nxlocal+1,0,m) = -uc(btnbr(k),3,0:nxlocal+1,1,m)
      uc(k,4,0:nxlocal+1,0,m) =  uc(btnbr(k),4,0:nxlocal+1,1,m)
      enddo
      
      end if
          
      end if
**************   (x0(t),1)-(4,1) : pre-shock condition   ***************
        den = 1.4
        vex = 0.0
        vey = 0.0
        pre = 1.0
      phy(1) = den
      phy(2) = den*vex
      phy(3) = den*vey
      phy(4) = pre/gamma1 + 0.5*den*(vex*vex + vey*vey)
      
      if (ycoord .eq. Ny_Processor) then
          
          ns1 = floor( (1.0/6.0 + (1.0 + 20.0*t)/3.0**0.5)*cdx )
      if (ns1 .ge. (xcoord-1)*nxlocal 
     *    .and. ns1 .lt. xcoord*nxlocal) then
              ns1 = ns1 - (xcoord - 1)*nxlocal
      do k = 1, nequ
        uc(:,k,ns1+1:nxlocal+1,nylocal+1,m) = phy(k)
      enddo  
      
      else if (ns1 .lt. (xcoord - 1)*nxlocal) then
      do k = 1, nequ
        uc(:,k,1:nxlocal+1,nylocal+1,m) = phy(k)
      enddo
      
      end if
      
      end if
      


      return
      endsubroutine
      
      
      

      endmodule
      
