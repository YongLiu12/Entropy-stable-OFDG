      module exactsolu
      use cons
      contains


c     exact solution for initial condition u(x,0) = sin(x)
	function burgersexac(xe,te)
	real :: x0 , pi , ddx , xs(0:200) , xt , xt0 , xxe , y
	integer :: j
	pi = 4.0 * atan(1.0)
	ddx = pi / 200.0
	xt0 = xe
	xt0 = xt0 - 2.0 * pi * floor((xt0+pi) /2.0 / pi)
	
      xxe = abs(xt0)
	
      do j = 0 , 200
	  xs(j) = ddx * j
	enddo
	do j = 0 , 200
	  xt = sin(xs(j)) * te + xs(j) - xxe
        if(abs(xt) .le. 15.0*ddx) then
        x0 = xs(j)
        do i = 1 , 100
		x0 = x0 - (xt)/(te*cos(x0)+1.0)
		xt = te * sin(x0) + x0 - xxe
		if(abs(xt) .le. 1.0e-14) goto 102
        enddo
        endif
	enddo

	
102	y = sign(1.0,xt0) * sin(x0)
	if(abs(xe) .lt. 1.0e-14) y = 0.0
	if(abs(te) .lt. 1.0e-14) y = sin(xe)
	burgersexac = y

      return
      endfunction
      


      function var_ent(uu)
      real :: uu, var_ent

        var_ent = d_ent(uu)


      return
      endfunction


c  entropy conservative flux
      function fmS(uL,uR)
      real :: uL, uR, fmS
      real :: dpsi(1:4), df(0:3), dU(1:4)
      real :: s1

      if(ic_ent .eq. 1) then
          
      fmS = (uL*uL + uL*uR + uR*uR)/6.0
      
      elseif(ic_ent .eq. 2) then
          
        s1 = var_ent(uR) - var_ent(uL)
          
        if(abs(s1) .lt. epweno) then
            
          df(0) = fluxfunc(uL)
          df(1) = d_fluxfunc(uL)
          df(2) = d2_fluxfunc(uL)
          df(3) = d3_fluxfunc(uL)
          dU(1) = d_ent(uL)
          dU(2) = d2_ent(uL)
          dU(3) = d3_ent(uL)
          dU(4) = d4_ent(uL)
          dpsi(1) = df(0)
          dpsi(2) = df(1)/dU(2)
          dpsi(3) = df(2)/dU(2)**2.0 
     *    - df(1)*dU(3)/dU(2)**3.0
          dpsi(4) = df(3)/dU(2)**3.0 
     *    - 3.0*df(2)*dU(3)/dU(2)**4.0 
     *    - df(1)*( -3.0*dU(3)**2.0/dU(2)**5.0 
     *    + dU(4)/dU(2)**4.0)
          
          fmS = 0.0
          do k = 1, 4
          fmS = fmS + dpsi(k)*s1**(k - 1)/stimes(k)
          enddo
          
        else
            
          fmS = ( (0.05*uR*uR - 0.1*uR + 0.1)*exp(uR) + 0.15*uR**3.0 
     *    - (0.05*uL*uL - 0.1*uL + 0.1)*exp(uL) - 0.15*uL**3.0 ) 
     *     / s1
          
        endif
        
      endif
      
      
      return
      endfunction
      
      function fluxfunc(uu)
      real :: uu, fluxfunc
      
      fluxfunc = 0.5*uu*uu
      
      
      return
      endfunction
      
      function d_fluxfunc(uu)
      real :: uu, d_fluxfunc
      
      d_fluxfunc = uu
      
      
      return
      endfunction
      
      function d2_fluxfunc(uu)
      real :: uu, d2_fluxfunc
      
      d2_fluxfunc = 1.0
      
      
      return
      endfunction
      
      function d3_fluxfunc(uu)
      real :: uu, d3_fluxfunc
      
      d3_fluxfunc = 0.0
      
      
      return
      endfunction
      
c  entropy
      function entropy(uu)
      real :: uu, entropy
      
      if(ic_ent .eq. 1) then
        entropy = 0.5*uu*uu
      elseif(ic_ent .eq. 2) then
        entropy = 0.1*exp(uu) + 0.45*uu*uu
      endif
      
      return
      endfunction
      
      
      function d_ent(uu)
      real :: uu, d_ent

      if(ic_ent .eq. 1) then
        d_ent = uu
      elseif(ic_ent .eq. 2) then
        d_ent = 0.1*exp(uu) + 0.9*uu
      endif


      return
      endfunction
      
      function d2_ent(uu)
      real :: uu, d2_ent

      if(ic_ent .eq. 1) then
        d2_ent = 1.0
      elseif(ic_ent .eq. 2) then
        d2_ent = 0.1*exp(uu) + 0.9
      endif


      return
      endfunction
      
      function d3_ent(uu)
      real :: uu, d3_ent

      if(ic_ent .eq. 1) then
        d3_ent = 0.0
      elseif(ic_ent .eq. 2) then
        d3_ent = 0.1*exp(uu)
      endif


      return
      endfunction
      
      function d4_ent(uu)
      real :: uu, d4_ent

      if(ic_ent .eq. 1) then
        d4_ent = 0.0
      elseif(ic_ent .eq. 2) then
        d4_ent = 0.1*exp(uu)
      endif


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
      
      
      endmodule exactsolu
