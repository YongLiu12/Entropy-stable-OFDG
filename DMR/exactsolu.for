      module exactsolu
      use cons
      
      contains


c     initial condition
      function u0(xs,ys)
      implicit none
      real :: xs, ys
      real :: u0(1:nequ), den, vex, vey, pre
      real :: as, alf, s1, s2, bs(0:3)
      
          
      if(xs .lt. 1.0/6.0 + ys/3.0**0.5) then
        den = 8.0
        vex = 8.25*cos(pi/6.0)
        vey = - 8.25*sin(pi/6.0)
        pre = 116.5
      else
        den = 1.4
        vex = 0.0
        vey = 0.0
        pre = 1.0
      endif
      
      u0(1) = den
      u0(2) = den*vex
      u0(3) = den*vey
      u0(4) = pre/gamma1 + 0.5*den*(vex*vex+vey*vey)
      
      
      
	return
      endfunction
      

      function press(uu)
      implicit none
      real :: press
      real :: uu(1:nequ)
      
      press = gamma1*( uu(4) - 0.5*(uu(2)*uu(2) 
     *    + uu(3)*uu(3))/uu(1) )
      
      return
      endfunction
      

      function var_ent(uu)
      real, dimension(1:nequ) :: uu, var_ent
      
        var_ent = d_ent(uu)


      return
      endfunction


      function fmS(uL,uR,m)
      implicit none
      integer :: m
      real, dimension(1:nequ) :: uL, uR, fmS 
      real :: rhoL, rhoR, vexL, vexR, veyL, veyR, preL, preR
     *    , betaL, betaR
      real :: rhom, rho_ln, betam, beta_ln
     *    , prem, pre_ln, vexm, veym, u_sqm
      real :: s1
      
      rhoL = abs(uL(1))
      rhoR = abs(uR(1))
      vexL = uL(2)/uL(1)
      vexR = uR(2)/uR(1)
      veyL = uL(3)/uL(1)
      veyR = uR(3)/uR(1)
      preL = abs(press(uL))
      preR = abs(press(uR))
      betaL = 0.5*rhoL/preL
      betaR = 0.5*rhoR/preR
      
      rhom = 0.5*(rhoL + rhoR)
      vexm = 0.5*(vexL + vexR)
      veym = 0.5*(veyL + veyR)
      betam = 0.5*(betaL + betaR)
      prem = 0.5*rhom/betam
      
      s1 = log(rhoR) - log(rhoL)
      if(abs(s1) .lt. epweno) then
        rho_ln = rhoL*(1.0 + s1/2.0 + s1*s1/6.0 
     *    + s1**3.0/24.0)
      else
        rho_ln = ( rhoR - rhoL ) / s1
      endif
      
      s1 = log(betaR) - log(betaL)
      if(abs(s1) .lt. epweno) then
        beta_ln = betaL*(1.0 + s1/2.0 + s1*s1/6.0 
     *    + s1**3.0/24.0)
      else
        beta_ln = ( betaR - betaL ) / s1
      endif
      
      pre_ln = 0.5*rho_ln/beta_ln
      u_sqm = 2.0*( vexm**2.0 + veym**2.0 )
     *    - 0.5*(vexL*vexL + vexR*vexR) 
     *    - 0.5*(veyL*veyL + veyR*veyR) 
      
      if(m .eq. 1) then
          
!        fmS(1) = rho_ln*vexm
!        fmS(2) = rho_ln*vexm**2.0 + prem
!        fmS(3) = rho_ln*vexm*veym
!        fmS(4) = vexm*( pre_ln/gamma1 
!     *    + rho_ln*u_sqm/2.0 + prem )
!c          
        fmS(1) = rho_ln*vexm
        fmS(2) = prem + vexm*fmS(1)
        fmS(3) = veym*fmS(1)
        fmS(4) = ( 0.5/gamma1/beta_ln 
     *    - 0.25*(vexL*vexL + vexR*vexR) 
     *    - 0.25*(veyL*veyL + veyR*veyR) )*fmS(1) 
     *    + vexm*fmS(2) + veym*fmS(3)
        
      elseif(m .eq. 2) then
      
!        fmS(1) = rho_ln*veym
!        fmS(2) = rho_ln*vexm*veym
!        fmS(3) = rho_ln*veym**2.0 + prem
!        fmS(4) = veym*( pre_ln/gamma1 
!     *    + rho_ln*u_sqm/2.0 + prem )
!c          
        fmS(1) = rho_ln*veym
        fmS(2) = vexm*fmS(1)
        fmS(3) = prem + veym*fmS(1)
        fmS(4) = (1.0/2.0/gamma1/beta_ln 
     *    - 0.25*(vexL*vexL + vexR*vexR) 
     *    - 0.25*(veyL*veyL + veyR*veyR) )*fmS(1) 
     *    + vexm*fmS(2) + veym*fmS(3)
        
      endif
        
      
      return
      endfunction
      
      
      function fluxfunc1(uu)
      implicit none
      real, dimension(1:nequ) :: uu, fluxfunc1
      real :: pre
      
      pre = press(uu(1:nequ))
      
      fluxfunc1(1) = uu(2)
      fluxfunc1(2) = uu(2)**2.0/uu(1) + pre
      fluxfunc1(3) = uu(2)*uu(3)/uu(1)
      fluxfunc1(4) = (uu(4) + pre)*uu(2)/uu(1)
      
      
      return
      endfunction
      
      
      function fluxfunc2(uu)
      implicit none
      real, dimension(1:nequ) :: uu, fluxfunc2
      real :: pre
      
      pre = press(uu(1:nequ))
      
      fluxfunc2(1) = uu(3)
      fluxfunc2(2) = uu(2)*uu(3)/uu(1)
      fluxfunc2(3) = uu(3)**2.0/uu(1) + pre
      fluxfunc2(4) = (uu(4) + pre)*uu(3)/uu(1)
      
      
      return
      endfunction
      
      
      function entropy(uu)
      implicit none
      real :: uu(1:nequ), entropy, den, pre, Ss
      
      den = uu(1)
      pre = press(uu(1:nequ))
      Ss = log( pre*den**(-gamma) )
      
      entropy = - den*Ss/gamma1
      
      
      return
      endfunction
      
      
      function d_ent(uu)
      implicit none
      real, dimension(1:nequ) :: uu, d_ent
      real :: den, vex, vey, pre, Ss, beta
      
      den = uu(1)
      vex = uu(2)/uu(1)
      vey = uu(3)/uu(1)
      pre = press(uu(1:nequ))
      Ss = log( pre*den**(-gamma) )
      beta = 0.5*den/pre
      
      d_ent(1) = (gamma - Ss)/gamma1 - beta*(vex*vex + vey*vey)
      d_ent(2) = 2.0*beta*vex
      d_ent(3) = 2.0*beta*vey
      d_ent(4) = - 2.0*beta
      

      return
      endfunction
      
      
c  factorial
      function stimes(m)
      
      real :: stimes
      integer :: k, m
      
      stimes = 1.0
      do k = 1, m
          stimes = stimes*k
      enddo
      
      return
      endfunction
      
      
      subroutine getmtx(vxm,vym,hm,evl,evr)
      implicit none
      real :: vxm, vym, hm, qm, cm, rcm
      real, dimension(1:nequ,1:nequ) :: evl, evr
      real :: t0, t1, t2, t3, b1, b2
      real :: abs, as(1:4), bs(1:4)
      integer :: k
      
      qm = 0.5 * ( vxm * vxm + vym * vym )
      cm = sqrt( gamma1 * abs( hm - qm ) )
      t0 = vxm * cm

          evr(1,1) = 1.0
          evr(1,2) = 0.0
          evr(1,3) = 1.0
          evr(1,4) = 1.0
          evr(2,1) = vxm - cm
          evr(2,2) = 0.0
          evr(2,3) = vxm
          evr(2,4) = vxm + cm
          evr(3,1) = vym
c          evr(3,2) = 1.0
          evr(3,3) = vym
          evr(3,4) = vym
          evr(4,1) = hm - t0
c          evr(4,2) = vym
          evr(4,3) = qm
          evr(4,4) = hm + t0
          evr(3,2) = -cm
          evr(4,2) = - vym * cm
          
          
          rcm = 1.0 / cm
	    evr = evr * rcm

            b1 = gamma1
            b2 = qm * b1
            t0 = vxm * cm
            t1 = b1 * vxm
            t2 = 0.5 * b1
            t3 = b1 * vym

          evl(1,1) = 0.5 * ( b2 + t0 )
          evl(1,2) = -0.5 * ( t1 + cm )
          evl(1,3) = -0.5 * t3
          evl(1,4) = t2 
c          evl(2,1) = - vym * cm*cm
          evl(2,2) = 0.0
c          evl(2,3) = cm*cm
          evl(2,4) = 0.0
          evl(3,1) = cm*cm - b2
          evl(3,2) = t1
          evl(3,3) = t3
          evl(3,4) = -b1
          evl(4,1) =  0.5 * ( b2 - t0 )
          evl(4,2) = -0.5 * ( t1 - cm )
          evl(4,3) = -0.5 * t3
          evl(4,4) = t2
          evl(2,1) = vym * cm
          evl(2,3) = -cm
          
	    evl = evl * rcm 
          
c          t0 = ( b1/2.0 )**0.25
          t0 = 1.0 / b1
          evl = evl / t0
          evr = evr * t0
          
          
          return
      endsubroutine
      
      
      subroutine getmty(vxm,vym,hm,evl,evr)
      implicit none
      real :: vxm, vym, hm, qm, cm, rcm
      real, dimension(1:nequ,1:nequ) :: evl, evr
      real :: t0, t1, t2, t3, b1, b2
      real :: abs, as(1:4), bs(1:4)
      integer :: k
      
           qm = 0.5 * ( vxm * vxm + vym * vym )
           cm = sqrt( gamma1 * abs( hm - qm ) )
           t0 = vym * cm

          evr(1,1) = 1.0
          evr(1,2) = 1.0
          evr(1,3) = 0.0
          evr(1,4) = 1.0
          evr(2,1) = vxm
          evr(2,2) = vxm
c          evr(2,3) = 1.0
          evr(2,4) = vxm
          evr(3,1) = vym - cm
          evr(3,2) = vym
          evr(3,3) = 0.0
          evr(3,4) = vym + cm
          evr(4,1) = hm - t0
          evr(4,2) = qm
c          evr(4,3) = vxm
          evr(4,4) = hm + t0
          evr(2,3) = -cm
          evr(4,3) = - vxm * cm
          
          rcm = 1.0 / cm
	    evr = evr * rcm

            b1 = gamma1
            b2 = qm * b1
            t0 = vym * cm
            t1 = b1 * vym
            t2 = 0.5 * b1
            t3 = b1 * vxm

          evl(1,1) = 0.5 * ( b2 + t0 )
          evl(1,2) = -0.5 * t3
          evl(1,3) = -0.5 * ( t1 + cm )
          evl(1,4) = t2
          evl(2,1) = cm*cm - b2
          evl(2,2) = t3
          evl(2,3) = t1
          evl(2,4) = -b1
c          evl(3,1) = - vxm * cm*cm
c          evl(3,2) = cm*cm
          evl(3,3) = 0.0
          evl(3,4) = 0.0
          evl(4,1) =  0.5 * ( b2 - t0 )
          evl(4,2) = -0.5 *   t3
          evl(4,3) = -0.5 * ( t1 - cm )
          evl(4,4) = t2
          evl(3,1) = vxm * cm
          evl(3,2) = -cm
          
          evl = evl * rcm
          
c          t0 = ( b1/2.0 )**0.25
          t0 = 1.0/b1
          evl = evl / t0
          evr = evr * t0
          
          
          return
      endsubroutine
      
      
      
      
      endmodule exactsolu
