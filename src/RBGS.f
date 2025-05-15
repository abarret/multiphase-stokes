c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Function to convert Thn to Ths since Thn + Ths = 1 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      module my_subs

      implicit none

      contains
      double precision function toThs(Thn)
c      
        double precision Thn  
        toThs = 1.d0 - Thn
c
        return
      end function toThs

      end module my_subs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes sol to a linear system on each cell in a patch i.e. box relaxation
c
c       where un, us are vector valued side-centered velocities
c             p is a scalar-valued cell-centered pressure
c       and   thn is a scalar-valued cell-centered network volume fraction
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine rbgs(
     &        dx, ilow0, iup0,
     &        ilow1, iup1, 
     &        un_0, un_1, un_gcw,
     &        us_0, us_1, us_gcw,
     &        p, p_gcw, f_p, f_p_gcw,
     &        f_un_0, f_un_1, f_un_gcw,
     &        f_us_0, f_us_1, f_us_gcw,
     &        thn, thn_gcw, eta_n, eta_s, l_n, l_s,
     &        nu_n, nu_s, xi, w, C, D, red_or_black)
c
      use my_subs
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0  
      integer ilow1,  iup1
      integer un_gcw, us_gcw, p_gcw, f_un_gcw, f_us_gcw
      integer f_p_gcw, thn_gcw, red_or_black
      double precision eta_n, eta_s, l_n, l_s, nu_n, nu_s, xi, w, C, D
c      
      double precision thn(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
c
      double precision un_0(ilow0-un_gcw:iup0+un_gcw+1,
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)
c
      double precision us_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
      double precision p(ilow0-p_gcw:iup0+p_gcw,
     &          ilow1-p_gcw:iup1+p_gcw)

      double precision f_un_0(ilow0-f_un_gcw:iup0+f_un_gcw+1,
     &          ilow1-f_un_gcw:iup1+f_un_gcw)
      double precision f_un_1(ilow0-f_un_gcw:iup0+f_un_gcw,
     &          ilow1-f_un_gcw:iup1+f_un_gcw+1)
c
      double precision f_us_0(ilow0-f_us_gcw:iup0+f_us_gcw+1,
     &          ilow1-f_us_gcw:iup1+f_us_gcw)
      double precision f_us_1(ilow0-f_us_gcw:iup0+f_us_gcw,
     &          ilow1-f_us_gcw:iup1+f_us_gcw+1)
c
      double precision f_p(ilow0-f_p_gcw:iup0+f_p_gcw,
     &          ilow1-f_p_gcw:iup1+f_p_gcw)
c
      double precision thn_lower_x, thn_lower_y
      double precision thn_upper_x, thn_upper_y
      double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
      double precision thn_iphalf_jphalf, thn_iphalf_jmhalf
      double precision ths_lower_x, ths_lower_y
      double precision ths_upper_x, ths_upper_y
      double precision ths_imhalf_jphalf, ths_imhalf_jmhalf
      double precision ths_iphalf_jphalf, ths_iphalf_jmhalf
c    
      double precision A_box(9,9)
      double precision b(9)
      double precision dx_dx, dx_dy, dy_dy
c
      integer i0, i1
c    
      dx_dx = 1.d0 / (dx(0) * dx(0))
      dx_dy = 1.d0 / (dx(0) * dx(1))
      dy_dy = 1.d0 / (dx(1) * dx(1))
      do i1 = ilow1, iup1   
        do i0 = ilow0, iup0  ! same loop as the c++ code (currently this is just GS)
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then
          
          thn_lower_x = 0.5d0*(thn(i0,i1)+thn(i0-1,i1))
          thn_upper_x = 0.5d0*(thn(i0,i1)+thn(i0+1,i1))
          thn_lower_y = 0.5d0*(thn(i0,i1)+thn(i0,i1-1))
          thn_upper_y = 0.5d0*(thn(i0,i1)+thn(i0,i1+1))

          ths_lower_x = toThs(thn_lower_x)
          ths_upper_x = toThs(thn_upper_x)
          ths_lower_y = toThs(thn_lower_y)
          ths_upper_y = toThs(thn_upper_y)

      ! calculate thn at corners
          thn_imhalf_jphalf = 0.25d0*(thn(i0-1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0-1,i1+1))
          thn_imhalf_jmhalf = 0.25d0*(thn(i0,i1)+thn(i0-1,i1)
     &                           +thn(i0,i1-1)+thn(i0-1,i1-1))
          thn_iphalf_jphalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0+1,i1+1))
          thn_iphalf_jmhalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1-1)+thn(i0+1,i1-1))

          ths_imhalf_jphalf = toThs(thn_imhalf_jphalf)
          ths_imhalf_jmhalf = toThs(thn_imhalf_jmhalf)
          ths_iphalf_jphalf = toThs(thn_iphalf_jphalf)
          ths_iphalf_jmhalf = toThs(thn_iphalf_jmhalf)

          A_box(1,1) = D*(
     &      (2.d0 * eta_n - l_n) * dx_dx
     &      * (-thn(i0,i1) - thn(i0-1,i1))
     &      - eta_n * dy_dy
     &      * (thn_imhalf_jmhalf + thn_imhalf_jphalf)
     &      - xi / nu_n * thn_lower_x * ths_lower_x)
     &      + C * thn_lower_x
          A_box(1,2) = D * (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0,i1)
          A_box(1,3) = D * (l_n * thn(i0,i1)
     &      - eta_n * thn_imhalf_jmhalf) * dx_dx
          A_box(1, 4) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_imhalf_jphalf) * dx_dy
          A_box(1,5) = D * xi / nu_s * thn_lower_x * ths_lower_x
          A_box(1,6) = 0.0
          A_box(1,7) = 0.0
          A_box(1,8) = 0.0
          A_box(1,9) = D * (-thn_lower_x) / dx(0)

          A_box(5,1) = D * xi / nu_n * thn_lower_x * ths_lower_x
          A_box(5,2) = 0.0
          A_box(5,3) = 0.0
          A_box(5,4) = 0.0
          A_box(5,5) = D*(
     &      (2.d0 * eta_s - l_s) * dx_dx
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0-1,i1)))
     &      - eta_s * dy_dy
     &      * (ths_imhalf_jmhalf + ths_imhalf_jphalf)
     &      - xi / nu_s * thn_lower_x * ths_lower_x)
     &      + C * ths_lower_x
          A_box(5,6) = D * (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0,i1))
          A_box(5,7) = D * (-eta_s * ths_imhalf_jmhalf
     &      + l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(5,8) = D *(eta_s * ths_imhalf_jphalf
     &      - l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(5,9) = D * (-ths_lower_x) / dx(0)

          A_box(2,1) = D * (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0,i1)
          A_box(2,2) = D * (
     &      (2.d0 * eta_n - l_n) * dx_dx
     &      * (-thn(i0+1,i1) - thn(i0,i1))
     &      - eta_n * dy_dy
     &      * (thn_iphalf_jphalf + thn_iphalf_jmhalf)
     &      - xi / nu_n * thn_upper_x * ths_upper_x)
     &      + C * thn_upper_x
          A_box(2,3) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_iphalf_jmhalf) * dx_dy
          A_box(2,4) = D * (l_n * thn(i0,i1)
     &    - eta_n * thn_iphalf_jphalf) * dx_dy
          A_box(2,6) = D * xi / nu_n * thn_upper_x * ths_upper_x
          A_box(2,5) = 0.0
          A_box(2,7) = 0.0
          A_box(2,8) = 0.0
          A_box(2,9) = D * thn_upper_x / dx(0)
c
          A_box(6,1) = 0.0
          A_box(6,2) = D * xi / nu_s * thn_upper_x * ths_upper_x
          A_box(6,3) = 0.0
          A_box(6,4) = 0.0
          A_box(6,5) = D * (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0,i1))
          A_box(6,6) = D * ((2.d0 * eta_s - l_s) * dx_dx
     &      * (-toThs(thn(i0+1,i1)) - toThs(thn(i0,i1)))
     &      - eta_s * dy_dy
     &      * (ths_iphalf_jphalf + ths_iphalf_jmhalf)
     &      - xi / nu_s * thn_upper_x * ths_upper_x)
     &      + C * ths_upper_x
          A_box(6,7) = D * (eta_s * ths_iphalf_jmhalf
     &      - l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(6,8) = D * (-eta_s * ths_iphalf_jphalf
     &      + l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(6,9) = D * ths_upper_x / dx(0)
c
      ! network at south edge
          A_box(3,1) = D * (l_n * thn(i0,i1)
     &      - eta_n * thn_imhalf_jmhalf) * dx_dy
          A_box(3,2) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_iphalf_jmhalf) * dx_dy
          A_box(3,3) = D*((2.d0 * eta_n - l_n) * dy_dy
     &      * (-thn(i0,i1) - thn(i0,i1-1))
     &      - eta_n * dx_dx
     &      * (thn_iphalf_jmhalf + thn_imhalf_jmhalf)
     &      - xi / nu_n * thn_lower_y * ths_lower_y)
     &      + C * thn_lower_y
          A_box(3,4) = D * (2.d0*eta_n-l_n) * dy_dy
     &      * thn(i0,i1)
          A_box(3,5) = 0.0
          A_box(3,6) = 0.0
          A_box(3,7) = D * xi / nu_s * thn_lower_y * ths_lower_y
          A_box(3,8) = 0.0
          A_box(3,9) = D * (-thn_lower_y) / dx(1)
c
          A_box(7,1) = 0.0
          A_box(7,2) = 0.0
          A_box(7,3) = D * xi / nu_n * thn_lower_y * ths_lower_y
          A_box(7,4) = 0.0
          A_box(7,5) = D * (l_s * toThs(thn(i0,i1))
     &      - eta_s * ths_imhalf_jmhalf) * dx_dy
          A_box(7,6) = D * (-l_s * toThs(thn(i0,i1))
     &      + eta_s * ths_iphalf_jmhalf) * dx_dy
          A_box(7,7) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1-1)))
     &      - eta_s * dx_dx
     &      * (ths_iphalf_jmhalf + ths_imhalf_jmhalf)
     &      - xi / nu_s * thn_lower_y * ths_lower_y)
     &      + C * ths_lower_y
          A_box(7,8) = D * (2.d0*eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1))
          A_box(7,9) = D * (-ths_lower_y) / dx(1)
c
      ! network at north edge
          A_box(4,1) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_imhalf_jphalf) * dx_dy
          A_box(4,2) = D * (l_n * thn(i0,i1)
     &      - eta_n *thn_iphalf_jphalf) * dx_dy
          A_box(4,3) = D * (2.d0 * eta_n - l_n) * dy_dy
     &      * (thn(i0,i1))
          A_box(4,4) = D * ((2.d0 * eta_n - l_n) * dy_dy
     &      * (-thn(i0,i1) - thn(i0,i1+1))
     &      - eta_n * dx_dx
     &      * (thn_iphalf_jphalf + thn_imhalf_jphalf)
     &      - xi / nu_n * thn_upper_y * ths_upper_y)
     &      + C * thn_upper_y
          A_box(4,5) = 0.0
          A_box(4,6) = 0.0
          A_box(4,7) = 0.0
          A_box(4,8) = D * xi / nu_s * thn_upper_y * ths_upper_y
          A_box(4,9) = D * thn_upper_y / dx(1)
c
          A_box(8,1) = 0.0
          A_box(8,2) = 0.0
          A_box(8,3) = 0.0
          A_box(8,4) = D * xi / nu_n * thn_upper_y * ths_upper_y
          A_box(8,5) = D * (-l_s * toThs(thn(i0,i1))
     &      + eta_s * ths_imhalf_jphalf) * dx_dy
          A_box(8,6) = D * (l_s * toThs(thn(i0,i1))
     &      - eta_s * ths_iphalf_jphalf) * dx_dy
          A_box(8,7) = D * (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1))
          A_box(8,8) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1+1)))
     &      - eta_s * dx_dx
     &      * (ths_iphalf_jphalf + ths_imhalf_jphalf)
     &      - xi / nu_s * thn_upper_y * ths_upper_y)
     &      + C * ths_upper_y
          A_box(8,9) = D * ths_upper_y / dx(1)
c
          ! incompressible constrain term at center
          A_box(9,1) = -thn_lower_x / dx(0)
          A_box(9,2) = thn_upper_x / dx(0)
          A_box(9,3) = -thn_lower_y / dx(1)
          A_box(9,4) = thn_upper_y / dx(1)
          A_box(9,5) = -ths_lower_x / dx(0)
          A_box(9,6) = ths_upper_x / dx(0)
          A_box(9,7) = -ths_lower_y / dx(1)
          A_box(9,8) = ths_upper_y / dx(1)
          A_box(9,9) = 0.0
c
          ! network at west edge
          b(1) = f_un_0(i0,i1) + D * (-thn_lower_x / dx(0)
     &      * p(i0-1,i1)
     &      - (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0-1,i1) * un_0(i0-1,i1)
     &      - eta_n * dy_dy
     &      * thn_imhalf_jphalf * un_0(i0,i1+1)
     &      - eta_n * dy_dy
     &      * thn_imhalf_jmhalf * un_0(i0,i1-1)
     &      + eta_n * dx_dy
     &      * thn_imhalf_jphalf * un_1(i0-1,i1+1)
     &      - eta_n * dx_dy
     &      * thn_imhalf_jmhalf * un_1(i0-1,i1)
     &      - l_n * dx_dy * thn(i0-1,i1)
     &      * (un_1(i0-1,i1+1) - un_1(i0-1,i1)))
c
          ! solvent at west edge
          b(5) = f_us_0(i0,i1) + D * (-ths_lower_x / dx(0)
     &      * p(i0-1,i1)
     &      - (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0-1,i1)) * us_0(i0-1,i1)
     &      - eta_s * dy_dy
     &      * ths_imhalf_jphalf * us_0(i0,i1+1)
     &      - eta_s * dy_dy
     &      * ths_imhalf_jmhalf * us_0(i0,i1-1)
     &      + eta_s * dx_dy
     &      * ths_imhalf_jphalf * us_1(i0-1,i1+1)
     &      - eta_s * dx_dy
     &      * ths_imhalf_jmhalf * us_1(i0-1,i1)
     &      - l_s * dx_dy
     &      * toThs(thn(i0-1,i1))
     &      * (us_1(i0-1,i1+1) - us_1(i0-1,i1)))
c
          ! network at east edge
          b(2) = f_un_0(i0+1,i1) + D * (thn_upper_x / dx(0)
     &      * p(i0+1,i1)
     &      - (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0+1,i1) * un_0(i0+2,i1)
     &      - eta_n * dy_dy * thn_iphalf_jphalf
     &      * un_0(i0+1,i1+1)
     &      - eta_n * dy_dy * thn_iphalf_jmhalf
     &      * un_0(i0+1,i1-1)
     &      - eta_n * dx_dy * thn_iphalf_jphalf
     &      * un_1(i0+1,i1+1)
     &      + eta_n * dx_dy * thn_iphalf_jmhalf
     &      * un_1(i0+1,i1)
     &      + l_n * dx_dy * thn(i0+1,i1)
     &      * (un_1(i0+1,i1+1) - un_1(i0+1,i1)))
c
          ! solvent at east edge
          b(6) = f_us_0(i0+1,i1) + D * (ths_upper_x / dx(0)
     &      * p(i0+1,i1)
     &      - (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0+1,i1)) * us_0(i0+2,i1)
     &      - eta_s * dy_dy * ths_iphalf_jphalf
     &      * us_0(i0+1,i1+1)
     &      - eta_s * dy_dy * ths_iphalf_jmhalf
     &      * us_0(i0+1,i1-1)
     &      - eta_s * dx_dy * ths_iphalf_jphalf
     &      * us_1(i0+1,i1+1)
     &      + eta_s * dx_dy * ths_iphalf_jmhalf
     &      * us_1(i0+1,i1)
     &      + l_s * dx_dy * toThs(thn(i0+1,i1))
     &      * (us_1(i0+1,i1+1) - us_1(i0+1,i1)))
c
          ! network at south edge
          b(3) = f_un_1(i0,i1) + D * (-thn_lower_y / dx(1)
     &      * p(i0,i1-1)
     &      - (2.d0 * eta_n - l_n) * dy_dy
     &      * thn(i0,i1-1) * un_1(i0,i1-1)
     &      - eta_n * dx_dx * thn_iphalf_jmhalf
     &      * un_1(i0+1,i1)
     &      - eta_n * dx_dx * thn_imhalf_jmhalf
     &      * un_1(i0-1,i1)
     &      + eta_n * dx_dy * thn_iphalf_jmhalf
     &      * un_0(i0+1,i1-1)
     &      - eta_n * dx_dy * thn_imhalf_jmhalf
     &      * un_0(i0,i1-1)
     &      - l_n * dx_dy * thn(i0,i1-1)
     &      * (un_0(i0+1,i1-1) - un_0(i0,i1-1)))
c
          ! solvent at south edge
          b(7) = f_us_1(i0,i1)+D*(-ths_lower_y / dx(1)
     &      * p(i0,i1-1)
     &      - (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1-1)) * us_1(i0,i1-1)
     &      - eta_s * dx_dx * ths_iphalf_jmhalf
     &      * us_1(i0+1,i1)
     &      - eta_s * dx_dx * ths_imhalf_jmhalf
     &      * us_1(i0-1,i1)
     &      + eta_s * dx_dy * ths_iphalf_jmhalf
     &      * us_0(i0+1,i1-1)
     &      - eta_s * dx_dy * ths_imhalf_jmhalf
     &      * us_0(i0,i1-1)
     &      - l_s * dx_dy * toThs(thn(i0,i1-1))
     &      * (us_0(i0+1,i1-1) - us_0(i0,i1-1)))
c
          ! network at north edge
          b(4) = f_un_1(i0,i1+1) + D * (thn_upper_y / dx(1)
     &      * p(i0,i1+1)
     &      - (2.d0 * eta_n - l_n) * dy_dy * thn(i0,i1+1)
     &      * un_1(i0,i1+2)
     &      - eta_n * dx_dx * thn_iphalf_jphalf
     &      * un_1(i0+1,i1+1)
     &      - eta_n * dx_dx * thn_imhalf_jphalf
     &      * un_1(i0-1,i1+1)
     &      - eta_n * dx_dy * thn_iphalf_jphalf
     &      * un_0(i0+1,i1+1)
     &      + eta_n * dx_dy * thn_imhalf_jphalf
     &      * un_0(i0,i1 + 1)
     &      + l_n * dx_dy * thn(i0,i1+1)
     &      * (un_0(i0+1,i1+1) - un_0(i0,i1 + 1)))
c
          ! solvent at north edge
          b(8) = f_us_1(i0,i1+1) + D * (ths_upper_y / dx(1)
     &      * p(i0,i1+1)
     &      - (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1+1)) * us_1(i0,i1+2)
     &      - eta_s * dx_dx * ths_iphalf_jphalf
     &      * us_1(i0+1,i1+1)
     &      - eta_s * dx_dx * ths_imhalf_jphalf
     &      * us_1(i0-1,i1+1)
     &      - eta_s * dx_dy * ths_iphalf_jphalf
     &      * us_0(i0+1,i1+1)
     &      + eta_s * dx_dy * ths_imhalf_jphalf
     &      * us_0(i0,i1 + 1)
     &      + l_s * dx_dy * toThs(thn(i0,i1+1))
     &      * (us_0(i0+1,i1+1) - us_0(i0,i1 + 1)))
c
          ! pressure at cell center
          b(9) = f_p(i0,i1)
c
          ! solve the system Ax = b, overwriting b with x
          call lu_solve(A_box, b)
c 
          un_0(i0,i1) = (1.d0-w)*un_0(i0,i1) + w*b(1);
          un_0(i0+1,i1) = (1.d0-w)*un_0(i0+1,i1) + w*b(2);
          un_1(i0,i1) = (1.d0-w)*un_1(i0,i1) + w*b(3);
          un_1(i0,i1+1) = (1.d0-w)*un_1(i0,i1+1) + w*b(4);
          us_0(i0,i1) = (1.d0-w)*us_0(i0,i1) + w*b(5);
          us_0(i0+1,i1) = (1.d0-w)*us_0(i0+1,i1) + w*b(6);
          us_1(i0,i1) = (1.d0-w)*us_1(i0,i1) + w*b(7);
          us_1(i0,i1+1) = (1.d0-w)*us_1(i0,i1+1) + w*b(8);
          p(i0,i1) = (1.d0-w)*p(i0,i1) + w*b(9);
c
          end if
        enddo
      enddo
c
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Same as rbgs, but allows for a mask index to indicate which side
c       indices should not be changed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rbgs_mask(
     &        dx, ilow0, iup0,
     &        ilow1, iup1,
     &        un_0, un_1, un_gcw,
     &        us_0, us_1, us_gcw,
     &        p, p_gcw, f_p, f_p_gcw,
     &        f_un_0, f_un_1, f_un_gcw,
     &        f_us_0, f_us_1, f_us_gcw,
     &        thn, thn_gcw, eta_n, eta_s, l_n, l_s,
     &        nu_n, nu_s, xi, w, C, D, red_or_black,
     &        mask_0, mask_1, mask_gcw)
c
      use my_subs
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0
      integer ilow1,  iup1
      integer un_gcw, us_gcw, p_gcw, f_un_gcw, f_us_gcw
      integer f_p_gcw, thn_gcw, red_or_black
      double precision eta_n, eta_s, l_n, l_s, nu_n, nu_s, xi, w, C, D
c
      double precision thn(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw)
c
      double precision un_0(ilow0-un_gcw:iup0+un_gcw+1,
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)
c
      double precision us_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
      double precision p(ilow0-p_gcw:iup0+p_gcw,
     &          ilow1-p_gcw:iup1+p_gcw)

      double precision f_un_0(ilow0-f_un_gcw:iup0+f_un_gcw+1,
     &          ilow1-f_un_gcw:iup1+f_un_gcw)
      double precision f_un_1(ilow0-f_un_gcw:iup0+f_un_gcw,
     &          ilow1-f_un_gcw:iup1+f_un_gcw+1)
c
      double precision f_us_0(ilow0-f_us_gcw:iup0+f_us_gcw+1,
     &          ilow1-f_us_gcw:iup1+f_us_gcw)
      double precision f_us_1(ilow0-f_us_gcw:iup0+f_us_gcw,
     &          ilow1-f_us_gcw:iup1+f_us_gcw+1)
c
      double precision f_p(ilow0-f_p_gcw:iup0+f_p_gcw,
     &          ilow1-f_p_gcw:iup1+f_p_gcw)

      integer mask_gcw
      integer mask_0(ilow0-mask_gcw:iup0+mask_gcw+1,
     &               ilow1-mask_gcw:iup1+mask_gcw)
      integer mask_1(ilow0-mask_gcw:iup0+mask_gcw,
     &               ilow1-mask_gcw:iup1+mask_gcw+1)
c
      double precision thn_lower_x, thn_lower_y
      double precision thn_upper_x, thn_upper_y
      double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
      double precision thn_iphalf_jphalf, thn_iphalf_jmhalf
      double precision ths_lower_x, ths_lower_y
      double precision ths_upper_x, ths_upper_y
      double precision ths_imhalf_jphalf, ths_imhalf_jmhalf
      double precision ths_iphalf_jphalf, ths_iphalf_jmhalf
c
      double precision A_box(9,9)
      double precision b(9)
      double precision dx_dx, dx_dy, dy_dy
c
      integer i0, i1
c
      dx_dx = 1.d0 / (dx(0) * dx(0))
      dx_dy = 1.d0 / (dx(0) * dx(1))
      dy_dy = 1.d0 / (dx(1) * dx(1))
      do i1 = ilow1, iup1
        do i0 = ilow0, iup0  ! same loop as the c++ code (currently this is just GS)
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then

          thn_lower_x = 0.5d0*(thn(i0,i1)+thn(i0-1,i1))
          thn_upper_x = 0.5d0*(thn(i0,i1)+thn(i0+1,i1))
          thn_lower_y = 0.5d0*(thn(i0,i1)+thn(i0,i1-1))
          thn_upper_y = 0.5d0*(thn(i0,i1)+thn(i0,i1+1))

          ths_lower_x = toThs(thn_lower_x)
          ths_upper_x = toThs(thn_upper_x)
          ths_lower_y = toThs(thn_lower_y)
          ths_upper_y = toThs(thn_upper_y)

      ! calculate thn at corners
          thn_imhalf_jphalf = 0.25d0*(thn(i0-1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0-1,i1+1))
          thn_imhalf_jmhalf = 0.25d0*(thn(i0,i1)+thn(i0-1,i1)
     &                           +thn(i0,i1-1)+thn(i0-1,i1-1))
          thn_iphalf_jphalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0+1,i1+1))
          thn_iphalf_jmhalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1-1)+thn(i0+1,i1-1))

          ths_imhalf_jphalf = toThs(thn_imhalf_jphalf)
          ths_imhalf_jmhalf = toThs(thn_imhalf_jmhalf)
          ths_iphalf_jphalf = toThs(thn_iphalf_jphalf)
          ths_iphalf_jmhalf = toThs(thn_iphalf_jmhalf)

          ! network at west edge
          if (mask_0(i0,i1) .eq. 1) then
            A_box(1, 1) = 1.d0
            A_box(1, 2) = 0.d0
            A_box(1, 3) = 0.d0
            A_box(1, 4) = 0.d0
            A_box(1, 5) = 0.d0
            A_box(1, 6) = 0.d0
            A_box(1, 7) = 0.d0
            A_box(1, 8) = 0.d0
            A_box(1, 9) = 0.d0
            A_box(5, 5) = 1.d0
            A_box(5, 6) = 0.d0
            A_box(5, 7) = 0.d0
            A_box(5, 8) = 0.d0
            A_box(5, 1) = 0.d0
            A_box(5, 2) = 0.d0
            A_box(5, 3) = 0.d0
            A_box(5, 4) = 0.d0
            A_box(5, 9) = 0.d0
            b(1) = un_0(i0,i1)
            b(5) = us_0(i0,i1)
          else
            A_box(1,1) = D*(
     &        (2.d0 * eta_n - l_n) * dx_dx
     &        * (-thn(i0,i1) - thn(i0-1,i1))
     &        - eta_n * dy_dy
     &        * (thn_imhalf_jmhalf + thn_imhalf_jphalf)
     &        - xi / nu_n * thn_lower_x * ths_lower_x)
     &        + C * thn_lower_x
            A_box(1,2) = D * (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0,i1)
            A_box(1,3) = D * (l_n * thn(i0,i1)
     &        - eta_n * thn_imhalf_jmhalf) * dx_dy
            A_box(1, 4) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_imhalf_jphalf) * dx_dy
            A_box(1,5) = D * xi / nu_s * thn_lower_x * ths_lower_x
            A_box(1,6) = 0.0
            A_box(1,7) = 0.0
            A_box(1,8) = 0.0
            A_box(1,9) = D * (-thn_lower_x) / dx(0)

            A_box(5,1) = D * xi / nu_n * thn_lower_x * ths_lower_x
            A_box(5,2) = 0.0
            A_box(5,3) = 0.0
            A_box(5,4) = 0.0
            A_box(5,5) = D*(
     &        (2.d0 * eta_s - l_s) * dx_dx
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0-1,i1)))
     &        - eta_s * dy_dy
     &        * (ths_imhalf_jmhalf + ths_imhalf_jphalf)
     &        - xi / nu_s * thn_lower_x * ths_lower_x)
     &        + C * ths_lower_x
            A_box(5,6) = D * (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0,i1))
            A_box(5,7) = D * (-eta_s * ths_imhalf_jmhalf
     &        + l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(5,8) = D *(eta_s * ths_imhalf_jphalf
     &        - l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(5,9) = D * (-ths_lower_x) / dx(0)
c
            b(1) = f_un_0(i0,i1) + D * (-thn_lower_x / dx(0)
     &        * p(i0-1,i1)
     &        - (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0-1,i1) * un_0(i0-1,i1)
     &        - eta_n * dy_dy
     &        * thn_imhalf_jphalf * un_0(i0,i1+1)
     &        - eta_n * dy_dy
     &        * thn_imhalf_jmhalf * un_0(i0,i1-1)
     &        + eta_n * dx_dy
     &        * thn_imhalf_jphalf * un_1(i0-1,i1+1)
     &        - eta_n * dx_dy
     &        * thn_imhalf_jmhalf * un_1(i0-1,i1)
     &        - l_n * dx_dy * thn(i0-1,i1)
     &        * (un_1(i0-1,i1+1) - un_1(i0-1,i1)))
c
          ! solvent at west edge
            b(5) = f_us_0(i0,i1) + D * (-ths_lower_x / dx(0)
     &        * p(i0-1,i1)
     &        - (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0-1,i1)) * us_0(i0-1,i1)
     &        - eta_s * dy_dy
     &        * ths_imhalf_jphalf * us_0(i0,i1+1)
     &        - eta_s * dy_dy
     &        * ths_imhalf_jmhalf * us_0(i0,i1-1)
     &        + eta_s * dx_dy
     &          * ths_imhalf_jphalf * us_1(i0-1,i1+1)
     &        - eta_s * dx_dy
     &        * ths_imhalf_jmhalf * us_1(i0-1,i1)
     &        - l_s * dx_dy
     &        * toThs(thn(i0-1,i1))
     &        * (us_1(i0-1,i1+1) - us_1(i0-1,i1)))
          endif
c
          ! network at east edge
          if (mask_0(i0+1,i1) .eq. 1) then
            A_box(2, 1) = 0.d0
            A_box(2, 2) = 1.d0
            A_box(2, 3) = 0.d0
            A_box(2, 4) = 0.d0
            A_box(2, 5) = 0.d0
            A_box(2, 7) = 0.d0
            A_box(2, 8) = 0.d0
            A_box(2, 6) = 0.d0
            A_box(2, 9) = 0.d0
c
            A_box(6, 5) = 0.d0
            A_box(6, 6) = 1.d0
            A_box(6, 7) = 0.d0
            A_box(6, 8) = 0.d0
            A_box(6, 1) = 0.d0
            A_box(6, 3) = 0.d0
            A_box(6, 4) = 0.d0
            A_box(6, 2) = 0.d0
            A_box(6, 9) = 0.d0

            b(2) = un_0(i0+1,i1)
            b(6) = us_0(i0+1,i1)
          else
            A_box(2,1) = D * (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0,i1)
            A_box(2,2) = D * (
     &        (2.d0 * eta_n - l_n) * dx_dx
     &        * (-thn(i0+1,i1) - thn(i0,i1))
     &        - eta_n * dy_dy
     &        * (thn_iphalf_jphalf + thn_iphalf_jmhalf)
     &        - xi / nu_n * thn_upper_x * ths_upper_x)
     &        + C * thn_upper_x
            A_box(2,3) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_iphalf_jmhalf) * dx_dy
            A_box(2,4) = D * (l_n * thn(i0,i1)
     &        - eta_n * thn_iphalf_jphalf) * dx_dy
            A_box(2,6) = D * xi / nu_n * thn_upper_x * ths_upper_x
            A_box(2,5) = 0.0
            A_box(2,7) = 0.0
            A_box(2,8) = 0.0
            A_box(2,9) = D * thn_upper_x / dx(0)
c
            A_box(6,1) = 0.0
            A_box(6,2) = D * xi / nu_s * thn_upper_x * ths_upper_x
            A_box(6,3) = 0.0
            A_box(6,4) = 0.0
            A_box(6,5) = D * (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0,i1))
            A_box(6,6) = D * ((2.d0 * eta_s - l_s) * dx_dx
     &        * (-toThs(thn(i0+1,i1)) - toThs(thn(i0,i1)))
     &        - eta_s * dy_dy
     &        * (ths_iphalf_jphalf + ths_iphalf_jmhalf)
     &        - xi / nu_s * thn_upper_x * ths_upper_x)
     &        + C * ths_upper_x
            A_box(6,7) = D * (eta_s * ths_iphalf_jmhalf
     &        - l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(6,8) = D * (-eta_s * ths_iphalf_jphalf
     &        + l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(6,9) = D * ths_upper_x / dx(0)

                      ! network at east edge
            b(2) = f_un_0(i0+1,i1) + D * (thn_upper_x / dx(0)
     &        * p(i0+1,i1)
     &        - (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0+1,i1) * un_0(i0+2,i1)
     &        - eta_n * dy_dy * thn_iphalf_jphalf
     &        * un_0(i0+1,i1+1)
     &        - eta_n * dy_dy * thn_iphalf_jmhalf
     &        * un_0(i0+1,i1-1)
     &        - eta_n * dx_dy * thn_iphalf_jphalf
     &        * un_1(i0+1,i1+1)
     &        + eta_n * dx_dy * thn_iphalf_jmhalf
     &        * un_1(i0+1,i1)
     &        + l_n * dx_dy * thn(i0+1,i1)
     &        * (un_1(i0+1,i1+1) - un_1(i0+1,i1)))
c
          ! solvent at east edge
            b(6) = f_us_0(i0+1,i1) + D * (ths_upper_x / dx(0)
     &        * p(i0+1,i1)
     &        - (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0+1,i1)) * us_0(i0+2,i1)
     &        - eta_s * dy_dy * ths_iphalf_jphalf
     &        * us_0(i0+1,i1+1)
     &        - eta_s * dy_dy * ths_iphalf_jmhalf
     &        * us_0(i0+1,i1-1)
     &        - eta_s * dx_dy * ths_iphalf_jphalf
     &        * us_1(i0+1,i1+1)
     &        + eta_s * dx_dy * ths_iphalf_jmhalf
     &        * us_1(i0+1,i1)
     &        + l_s * dx_dy * toThs(thn(i0+1,i1))
     &        * (us_1(i0+1,i1+1) - us_1(i0+1,i1)))
          endif
c
          ! network at south edge
          if (mask_1(i0,i1) .eq. 1) then
            A_box(3, 1) = 0.d0
            A_box(3, 2) = 0.d0
            A_box(3, 3) = 1.d0
            A_box(3, 4) = 0.d0
            A_box(3, 5) = 0.d0
            A_box(3, 6) = 0.d0
            A_box(3, 8) = 0.d0
            A_box(3, 7) = 0.d0
            A_box(3, 9) = 0.d0
c
            A_box(7, 5) = 0.d0
            A_box(7, 6) = 0.d0
            A_box(7, 7) = 1.d0
            A_box(7, 8) = 0.d0
            A_box(7, 1) = 0.d0
            A_box(7, 2) = 0.d0
            A_box(7, 4) = 0.d0
            A_box(7, 3) = 0.d0
            A_box(7, 9) = 0.d0

            b(3) = un_1(i0,i1)
            b(7) = us_1(i0,i1)
          else
            A_box(3,1) = D * (l_n * thn(i0,i1)
     &        - eta_n * thn_imhalf_jmhalf) * dx_dy
            A_box(3,2) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_iphalf_jmhalf) * dx_dy
            A_box(3,3) = D*((2.d0 * eta_n - l_n) * dy_dy
     &        * (-thn(i0,i1) - thn(i0,i1-1))
     &        - eta_n * dx_dx
     &        * (thn_iphalf_jmhalf + thn_imhalf_jmhalf)
     &        - xi / nu_n * thn_lower_y * ths_lower_y)
     &        + C * thn_lower_y
            A_box(3,4) = D * (2.d0*eta_n-l_n) * dy_dy
     &        * thn(i0,i1)
            A_box(3,5) = 0.0
            A_box(3,6) = 0.0
            A_box(3,7) = D * xi / nu_s * thn_lower_y * ths_lower_y
            A_box(3,8) = 0.0
            A_box(3,9) = D * (-thn_lower_y) / dx(1)
c
            A_box(7,1) = 0.0
            A_box(7,2) = 0.0
            A_box(7,3) = D * xi / nu_n * thn_lower_y * ths_lower_y
            A_box(7,4) = 0.0
            A_box(7,5) = D * (l_s * toThs(thn(i0,i1))
     &        - eta_s * ths_imhalf_jmhalf) * dx_dy
            A_box(7,6) = D * (-l_s * toThs(thn(i0,i1))
     &        + eta_s * ths_iphalf_jmhalf) * dx_dy
            A_box(7,7) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1-1)))
     &        - eta_s * dx_dx
     &        * (ths_iphalf_jmhalf + ths_imhalf_jmhalf)
     &        - xi / nu_s * thn_lower_y * ths_lower_y)
     &        + C * ths_lower_y
            A_box(7,8) = D * (2.d0*eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1))
            A_box(7,9) = D * (-ths_lower_y) / dx(1)
c
            b(3) = f_un_1(i0,i1) + D * (-thn_lower_y / dx(1)
     &        * p(i0,i1-1)
     &        - (2.d0 * eta_n - l_n) * dy_dy
     &        * thn(i0,i1-1) * un_1(i0,i1-1)
     &        - eta_n * dx_dx * thn_iphalf_jmhalf
     &        * un_1(i0+1,i1)
     &        - eta_n * dx_dx * thn_imhalf_jmhalf
     &        * un_1(i0-1,i1)
     &        + eta_n * dx_dy * thn_iphalf_jmhalf
     &        * un_0(i0+1,i1-1)
     &        - eta_n * dx_dy * thn_imhalf_jmhalf
     &        * un_0(i0,i1-1)
     &        - l_n * dx_dy * thn(i0,i1-1)
     &        * (un_0(i0+1,i1-1) - un_0(i0,i1-1)))
c
            b(7) = f_us_1(i0,i1)+D*(-ths_lower_y / dx(1)
     &        * p(i0,i1-1)
     &        - (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1-1)) * us_1(i0,i1-1)
     &        - eta_s * dx_dx * ths_iphalf_jmhalf
     &        * us_1(i0+1,i1)
     &        - eta_s * dx_dx * ths_imhalf_jmhalf
     &        * us_1(i0-1,i1)
     &        + eta_s * dx_dy * ths_iphalf_jmhalf
     &        * us_0(i0+1,i1-1)
     &        - eta_s * dx_dy * ths_imhalf_jmhalf
     &        * us_0(i0,i1-1)
     &        - l_s * dx_dy * toThs(thn(i0,i1-1))
     &        * (us_0(i0+1,i1-1) - us_0(i0,i1-1)))
          endif
c
          ! network at north edge
          if (mask_1(i0,i1+1) .eq. 1) then
            A_box(4, 1) = 0.d0
            A_box(4, 2) = 0.d0
            A_box(4, 3) = 0.d0
            A_box(4, 4) = 1.d0
            A_box(4, 5) = 0.d0
            A_box(4, 6) = 0.d0
            A_box(4, 7) = 0.d0
            A_box(4, 8) = 0.d0
            A_box(4, 9) = 0.d0
c
            A_box(8, 5) = 0.d0
            A_box(8, 6) = 0.d0
            A_box(8, 7) = 0.d0
            A_box(8, 8) = 1.d0
            A_box(8, 1) = 0.d0
            A_box(8, 2) = 0.d0
            A_box(8, 3) = 0.d0
            A_box(8, 4) = 0.d0
            A_box(8, 9) = 0.d0

            b(4) = un_1(i0,i1+1)
            b(8) = us_1(i0,i1+1)
          else
            A_box(4,1) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_imhalf_jphalf) * dx_dy
            A_box(4,2) = D * (l_n * thn(i0,i1)
     &        - eta_n *thn_iphalf_jphalf) * dx_dy
            A_box(4,3) = D * (2.d0 * eta_n - l_n) * dy_dy
     &        * (thn(i0,i1))
            A_box(4,4) = D * ((2.d0 * eta_n - l_n) * dy_dy
     &        * (-thn(i0,i1) - thn(i0,i1+1))
     &        - eta_n * dx_dx
     &        * (thn_iphalf_jphalf + thn_imhalf_jphalf)
     &        - xi / nu_n * thn_upper_y * ths_upper_y)
     &        + C * thn_upper_y
            A_box(4,5) = 0.0
            A_box(4,6) = 0.0
            A_box(4,7) = 0.0
            A_box(4,8) = D * xi / nu_s * thn_upper_y * ths_upper_y
            A_box(4,9) = D * thn_upper_y / dx(1)
c
            A_box(8,1) = 0.0
            A_box(8,2) = 0.0
            A_box(8,3) = 0.0
            A_box(8,4) = D * xi / nu_n * thn_upper_y * ths_upper_y
            A_box(8,5) = D * (-l_s * toThs(thn(i0,i1))
     &        + eta_s * ths_imhalf_jphalf) * dx_dy
            A_box(8,6) = D * (l_s * toThs(thn(i0,i1))
     &        - eta_s * ths_iphalf_jphalf) * dx_dy
            A_box(8,7) = D * (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1))
            A_box(8,8) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1+1)))
     &        - eta_s * dx_dx
     &        * (ths_iphalf_jphalf + ths_imhalf_jphalf)
     &        - xi / nu_s * thn_upper_y * ths_upper_y)
     &        + C * ths_upper_y
            A_box(8,9) = D * ths_upper_y / dx(1)
                      ! network at north edge
            b(4) = f_un_1(i0,i1+1) + D * (thn_upper_y / dx(1)
     &        * p(i0,i1+1)
     &        - (2.d0 * eta_n - l_n) * dy_dy * thn(i0,i1+1)
     &        * un_1(i0,i1+2)
     &        - eta_n * dx_dx * thn_iphalf_jphalf
     &        * un_1(i0+1,i1+1)
     &        - eta_n * dx_dx * thn_imhalf_jphalf
     &        * un_1(i0-1,i1+1)
     &        - eta_n * dx_dy * thn_iphalf_jphalf
     &        * un_0(i0+1,i1+1)
     &        + eta_n * dx_dy * thn_imhalf_jphalf
     &        * un_0(i0,i1 + 1)
     &        + l_n * dx_dy * thn(i0,i1+1)
     &        * (un_0(i0+1,i1+1) - un_0(i0,i1 + 1)))
c
          ! solvent at north edge
            b(8) = f_us_1(i0,i1+1) + D * (ths_upper_y / dx(1)
     &        * p(i0,i1+1)
     &        - (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1+1)) * us_1(i0,i1+2)
     &        - eta_s * dx_dx * ths_iphalf_jphalf
     &        * us_1(i0+1,i1+1)
     &        - eta_s * dx_dx * ths_imhalf_jphalf
     &        * us_1(i0-1,i1+1)
     &        - eta_s * dx_dy * ths_iphalf_jphalf
     &        * us_0(i0+1,i1+1)
     &        + eta_s * dx_dy * ths_imhalf_jphalf
     &        * us_0(i0,i1 + 1)
     &        + l_s * dx_dy * toThs(thn(i0,i1+1))
     &        * (us_0(i0+1,i1+1) - us_0(i0,i1 + 1)))
          endif
c
          ! incompressible constrain term at center
          A_box(9, 1) = -thn_lower_x / dx(0)
          A_box(9, 2) = thn_upper_x / dx(0)
          A_box(9, 3) = -thn_lower_y / dx(1)
          A_box(9, 4) = thn_upper_y / dx(1)
          A_box(9, 5) = -toThs(thn_lower_x) / dx(0)
          A_box(9, 6) = toThs(thn_upper_x) / dx(0)
          A_box(9, 7) = -toThs(thn_lower_y) / dx(1)
          A_box(9, 8) = toThs(thn_upper_y) / dx(1)
          A_box(9, 9) = 0.0
          ! pressure at cell center
          b(9) = f_p(i0,i1)
c
          ! solve the system Ax = b, overwriting b with x
          call lu_solve(A_box, b)
c
          if (mask_0(i0,i1) .eq. 0) then
            un_0(i0,i1) = (1.d0-w)*un_0(i0,i1) + w*b(1);
            us_0(i0,i1) = (1.d0-w)*us_0(i0,i1) + w*b(5);
          endif
          if (mask_0(i0+1,i1) .eq. 0) then
            un_0(i0+1,i1) = (1.d0-w)*un_0(i0+1,i1) + w*b(2);
            us_0(i0+1,i1) = (1.d0-w)*us_0(i0+1,i1) + w*b(6);
          endif
          if (mask_1(i0,i1) .eq. 0) then
            un_1(i0,i1) = (1.d0-w)*un_1(i0,i1) + w*b(3);
            us_1(i0,i1) = (1.d0-w)*us_1(i0,i1) + w*b(7);
          endif
          if (mask_1(i0,i1+1) .eq. 0) then
            un_1(i0,i1+1) = (1.d0-w)*un_1(i0,i1+1) + w*b(4);
            us_1(i0,i1+1) = (1.d0-w)*us_1(i0,i1+1) + w*b(8);
          endif
          p(i0,i1) = (1.d0-w)*p(i0,i1) + w*b(9);
c
          end if
        enddo
      enddo
c
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes sol to a linear system on each cell in a patch i.e. box relaxation
c
c       where un, us are vector valued side-centered velocities
c             p is a scalar-valued cell-centered pressure
c       and   thn is a scalar-valued cell-centered network volume fraction
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rbgs_var_xi(
     &        dx, ilow0, iup0,
     &        ilow1, iup1,
     &        un_0, un_1, un_gcw,
     &        us_0, us_1, us_gcw,
     &        p, p_gcw, f_p, f_p_gcw,
     &        f_un_0, f_un_1, f_un_gcw,
     &        f_us_0, f_us_1, f_us_gcw,
     &        thn, thn_gcw, eta_n, eta_s, l_n, l_s,
     &        xi_0, xi_1, xi_gcw,
     &        w, C, D, red_or_black)
c
      use my_subs
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0
      integer ilow1,  iup1
      integer un_gcw, us_gcw, p_gcw, f_un_gcw, f_us_gcw
      integer f_p_gcw, thn_gcw, red_or_black
      double precision eta_n, eta_s, l_n, l_s, w, C, D
c
      double precision thn(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw)
c
      double precision un_0(ilow0-un_gcw:iup0+un_gcw+1,
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)
c
      double precision us_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
      double precision p(ilow0-p_gcw:iup0+p_gcw,
     &          ilow1-p_gcw:iup1+p_gcw)

      double precision f_un_0(ilow0-f_un_gcw:iup0+f_un_gcw+1,
     &          ilow1-f_un_gcw:iup1+f_un_gcw)
      double precision f_un_1(ilow0-f_un_gcw:iup0+f_un_gcw,
     &          ilow1-f_un_gcw:iup1+f_un_gcw+1)
c
      double precision f_us_0(ilow0-f_us_gcw:iup0+f_us_gcw+1,
     &          ilow1-f_us_gcw:iup1+f_us_gcw)
      double precision f_us_1(ilow0-f_us_gcw:iup0+f_us_gcw,
     &          ilow1-f_us_gcw:iup1+f_us_gcw+1)
c
      double precision f_p(ilow0-f_p_gcw:iup0+f_p_gcw,
     &          ilow1-f_p_gcw:iup1+f_p_gcw)

      integer xi_gcw
      double precision xi_0(ilow0-xi_gcw:iup0+xi_gcw+1,
     &                      ilow1-xi_gcw:iup1+xi_gcw)
      double precision xi_1(ilow0-xi_gcw:iup0+xi_gcw,
     &                      ilow1-xi_gcw:iup1+xi_gcw+1)
c
      double precision thn_lower_x, thn_lower_y
      double precision thn_upper_x, thn_upper_y
      double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
      double precision thn_iphalf_jphalf, thn_iphalf_jmhalf
      double precision ths_lower_x, ths_lower_y
      double precision ths_upper_x, ths_upper_y
      double precision ths_imhalf_jphalf, ths_imhalf_jmhalf
      double precision ths_iphalf_jphalf, ths_iphalf_jmhalf
c
      double precision A_box(9,9)
      double precision b(9)
      double precision dx_dx, dx_dy, dy_dy
c
      integer i0, i1
c
      dx_dx = 1.d0 / (dx(0) * dx(0))
      dx_dy = 1.d0 / (dx(0) * dx(1))
      dy_dy = 1.d0 / (dx(1) * dx(1))
      do i1 = ilow1, iup1
        do i0 = ilow0, iup0  ! same loop as the c++ code (currently this is just GS)
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then

          thn_lower_x = 0.5d0*(thn(i0,i1)+thn(i0-1,i1))
          thn_upper_x = 0.5d0*(thn(i0,i1)+thn(i0+1,i1))
          thn_lower_y = 0.5d0*(thn(i0,i1)+thn(i0,i1-1))
          thn_upper_y = 0.5d0*(thn(i0,i1)+thn(i0,i1+1))

          ths_lower_x = toThs(thn_lower_x)
          ths_upper_x = toThs(thn_upper_x)
          ths_lower_y = toThs(thn_lower_y)
          ths_upper_y = toThs(thn_upper_y)

      ! calculate thn at corners
          thn_imhalf_jphalf = 0.25d0*(thn(i0-1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0-1,i1+1))
          thn_imhalf_jmhalf = 0.25d0*(thn(i0,i1)+thn(i0-1,i1)
     &                           +thn(i0,i1-1)+thn(i0-1,i1-1))
          thn_iphalf_jphalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0+1,i1+1))
          thn_iphalf_jmhalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1-1)+thn(i0+1,i1-1))

          ths_imhalf_jphalf = toThs(thn_imhalf_jphalf)
          ths_imhalf_jmhalf = toThs(thn_imhalf_jmhalf)
          ths_iphalf_jphalf = toThs(thn_iphalf_jphalf)
          ths_iphalf_jmhalf = toThs(thn_iphalf_jmhalf)

          A_box(1,1) = D*(
     &      (2.d0 * eta_n - l_n) * dx_dx
     &      * (-thn(i0,i1) - thn(i0-1,i1))
     &      - eta_n * dy_dy
     &      * (thn_imhalf_jmhalf + thn_imhalf_jphalf)
     &      - xi_0(i0,i1))
     &      + C * thn_lower_x
          A_box(1,2) = D * (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0,i1)
          A_box(1,3) = D * (l_n * thn(i0,i1)
     &      - eta_n * thn_imhalf_jmhalf) * dx_dy
          A_box(1, 4) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_imhalf_jphalf) * dx_dy
          A_box(1,5) = D * xi_0(i0,i1)
          A_box(1,6) = 0.0
          A_box(1,7) = 0.0
          A_box(1,8) = 0.0
          A_box(1,9) = D * (-thn_lower_x) / dx(0)

          A_box(5,1) = D * xi_0(i0,i1)
          A_box(5,2) = 0.0
          A_box(5,3) = 0.0
          A_box(5,4) = 0.0
          A_box(5,5) = D*(
     &      (2.d0 * eta_s - l_s) * dx_dx
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0-1,i1)))
     &      - eta_s * dy_dy
     &      * (ths_imhalf_jmhalf + ths_imhalf_jphalf)
     &      - xi_0(i0,i1))
     &      + C * ths_lower_x
          A_box(5,6) = D * (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0,i1))
          A_box(5,7) = D * (-eta_s * ths_imhalf_jmhalf
     &      + l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(5,8) = D *(eta_s * ths_imhalf_jphalf
     &      - l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(5,9) = D * (-ths_lower_x) / dx(0)

          A_box(2,1) = D * (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0,i1)
          A_box(2,2) = D * (
     &      (2.d0 * eta_n - l_n) * dx_dx
     &      * (-thn(i0+1,i1) - thn(i0,i1))
     &      - eta_n * dy_dy
     &      * (thn_iphalf_jphalf + thn_iphalf_jmhalf)
     &      - xi_0(i0+1,i1))
     &      + C * thn_upper_x
          A_box(2,3) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_iphalf_jmhalf) * dx_dy
          A_box(2,4) = D * (l_n * thn(i0,i1)
     &    - eta_n * thn_iphalf_jphalf) * dx_dy
          A_box(2,6) = D * xi_0(i0+1,i1)
          A_box(2,5) = 0.0
          A_box(2,7) = 0.0
          A_box(2,8) = 0.0
          A_box(2,9) = D * thn_upper_x / dx(0)
c
          A_box(6,1) = 0.0
          A_box(6,2) = D * xi_0(i0+1,i1)
          A_box(6,3) = 0.0
          A_box(6,4) = 0.0
          A_box(6,5) = D * (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0,i1))
          A_box(6,6) = D * ((2.d0 * eta_s - l_s) * dx_dx
     &      * (-toThs(thn(i0+1,i1)) - toThs(thn(i0,i1)))
     &      - eta_s * dy_dy
     &      * (ths_iphalf_jphalf + ths_iphalf_jmhalf)
     &      - xi_0(i0+1,i1))
     &      + C * ths_upper_x
          A_box(6,7) = D * (eta_s * ths_iphalf_jmhalf
     &      - l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(6,8) = D * (-eta_s * ths_iphalf_jphalf
     &      + l_s * toThs(thn(i0,i1))) * dx_dy
          A_box(6,9) = D * ths_upper_x / dx(0)
c
      ! network at south edge
          A_box(3,1) = D * (l_n * thn(i0,i1)
     &      - eta_n * thn_imhalf_jmhalf) * dx_dy
          A_box(3,2) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_iphalf_jmhalf) * dx_dy
          A_box(3,3) = D*((2.d0 * eta_n - l_n) * dy_dy
     &      * (-thn(i0,i1) - thn(i0,i1-1))
     &      - eta_n * dx_dx
     &      * (thn_iphalf_jmhalf + thn_imhalf_jmhalf)
     &      - xi_1(i0,i1))
     &      + C * thn_lower_y
          A_box(3,4) = D * (2.d0*eta_n-l_n) * dy_dy
     &      * thn(i0,i1)
          A_box(3,5) = 0.0
          A_box(3,6) = 0.0
          A_box(3,7) = D * xi_1(i0,i1)
          A_box(3,8) = 0.0
          A_box(3,9) = D * (-thn_lower_y) / dx(1)
c
          A_box(7,1) = 0.0
          A_box(7,2) = 0.0
          A_box(7,3) = D * xi_1(i0,i1)
          A_box(7,4) = 0.0
          A_box(7,5) = D * (l_s * toThs(thn(i0,i1))
     &      - eta_s * ths_imhalf_jmhalf) * dx_dy
          A_box(7,6) = D * (-l_s * toThs(thn(i0,i1))
     &      + eta_s * ths_iphalf_jmhalf) * dx_dy
          A_box(7,7) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1-1)))
     &      - eta_s * dx_dx
     &      * (ths_iphalf_jmhalf + ths_imhalf_jmhalf)
     &      - xi_1(i0,i1))
     &      + C * ths_lower_y
          A_box(7,8) = D * (2.d0*eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1))
          A_box(7,9) = D * (-ths_lower_y) / dx(1)
c
      ! network at north edge
          A_box(4,1) = D * (-l_n * thn(i0,i1)
     &      + eta_n * thn_imhalf_jphalf) * dx_dy
          A_box(4,2) = D * (l_n * thn(i0,i1)
     &      - eta_n *thn_iphalf_jphalf) * dx_dy
          A_box(4,3) = D * (2.d0 * eta_n - l_n) * dy_dy
     &      * (thn(i0,i1))
          A_box(4,4) = D * ((2.d0 * eta_n - l_n) * dy_dy
     &      * (-thn(i0,i1) - thn(i0,i1+1))
     &      - eta_n * dx_dx
     &      * (thn_iphalf_jphalf + thn_imhalf_jphalf)
     &      - xi_1(i0,i1+1))
     &      + C * thn_upper_y
          A_box(4,5) = 0.0
          A_box(4,6) = 0.0
          A_box(4,7) = 0.0
          A_box(4,8) = D * xi_1(i0,i1+1)
          A_box(4,9) = D * thn_upper_y / dx(1)
c
          A_box(8,1) = 0.0
          A_box(8,2) = 0.0
          A_box(8,3) = 0.0
          A_box(8,4) = D * xi_1(i0,i1+1)
          A_box(8,5) = D * (-l_s * toThs(thn(i0,i1))
     &      + eta_s * ths_imhalf_jphalf) * dx_dy
          A_box(8,6) = D * (l_s * toThs(thn(i0,i1))
     &      - eta_s * ths_iphalf_jphalf) * dx_dy
          A_box(8,7) = D * (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1))
          A_box(8,8) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &      * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1+1)))
     &      - eta_s * dx_dx
     &      * (ths_iphalf_jphalf + ths_imhalf_jphalf)
     &      - xi_1(i0,i1+1))
     &      + C * ths_upper_y
          A_box(8,9) = D * ths_upper_y / dx(1)

          ! incompressible constrain term at center
          A_box(9,1) = -thn_lower_x / dx(0)
          A_box(9,2) = thn_upper_x / dx(0)
          A_box(9,3) = -thn_lower_y / dx(1)
          A_box(9,4) = thn_upper_y / dx(1)
          A_box(9,5) = -ths_lower_x / dx(0)
          A_box(9,6) = ths_upper_x / dx(0)
          A_box(9,7) = -ths_lower_y / dx(1)
          A_box(9,8) = ths_upper_y / dx(1)
          A_box(9,9) = 0.0
c
          ! network at west edge
          b(1) = f_un_0(i0,i1) + D * (-thn_lower_x / dx(0)
     &      * p(i0-1,i1)
     &      - (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0-1,i1) * un_0(i0-1,i1)
     &      - eta_n * dy_dy
     &      * thn_imhalf_jphalf * un_0(i0,i1+1)
     &      - eta_n * dy_dy
     &      * thn_imhalf_jmhalf * un_0(i0,i1-1)
     &      + eta_n * dx_dy
     &      * thn_imhalf_jphalf * un_1(i0-1,i1+1)
     &      - eta_n * dx_dy
     &      * thn_imhalf_jmhalf * un_1(i0-1,i1)
     &      - l_n * dx_dy * thn(i0-1,i1)
     &      * (un_1(i0-1,i1+1) - un_1(i0-1,i1)))
c
          ! solvent at west edge
          b(5) = f_us_0(i0,i1) + D * (-ths_lower_x / dx(0)
     &      * p(i0-1,i1)
     &      - (2.d0 * eta_s - l_s) * dx_dx
     &      * toThs(thn(i0-1,i1)) * us_0(i0-1,i1)
     &      - eta_s * dy_dy
     &      * ths_imhalf_jphalf * us_0(i0,i1+1)
     &      - eta_s * dy_dy
     &      * ths_imhalf_jmhalf * us_0(i0,i1-1)
     &      + eta_s * dx_dy
     &      * ths_imhalf_jphalf * us_1(i0-1,i1+1)
     &      - eta_s * dx_dy
     &      * ths_imhalf_jmhalf * us_1(i0-1,i1)
     &      - l_s * dx_dy
     &      * toThs(thn(i0-1,i1))
     &      * (us_1(i0-1,i1+1) - us_1(i0-1,i1)))
c
          ! network at east edge
          b(2) = f_un_0(i0+1,i1) + D * (thn_upper_x / dx(0)
     &      * p(i0+1,i1)
     &      - (2.d0 * eta_n - l_n) * dx_dx
     &      * thn(i0+1,i1) * un_0(i0+2,i1)
     &      - eta_n * dy_dy * thn_iphalf_jphalf
     &      * un_0(i0+1,i1+1)
     &      - eta_n * dy_dy * thn_iphalf_jmhalf
     &      * un_0(i0+1,i1-1)
     &      - eta_n * dx_dy * thn_iphalf_jphalf
     &      * un_1(i0+1,i1+1)
     &      + eta_n * dx_dy * thn_iphalf_jmhalf
     &      * un_1(i0+1,i1)
     &      + l_n * dx_dy * thn(i0+1,i1)
     &      * (un_1(i0+1,i1+1) - un_1(i0+1,i1)))
c
          ! solvent at east edge
          b(6) = f_us_0(i0+1,i1) + D * (ths_upper_x / dx(0)
     &      * p(i0+1,i1)
     &      - (2.d0 * eta_s - l_s) / (dx(0)* dx(0))
     &      * toThs(thn(i0+1,i1)) * us_0(i0+2,i1)
     &      - eta_s * dy_dy * ths_iphalf_jphalf
     &      * us_0(i0+1,i1+1)
     &      - eta_s * dy_dy * ths_iphalf_jmhalf
     &      * us_0(i0+1,i1-1)
     &      - eta_s * dx_dy * ths_iphalf_jphalf
     &      * us_1(i0+1,i1+1)
     &      + eta_s * dx_dy * ths_iphalf_jmhalf
     &      * us_1(i0+1,i1)
     &      + l_s * dx_dy * toThs(thn(i0+1,i1))
     &      * (us_1(i0+1,i1+1) - us_1(i0+1,i1)))
c
          ! network at south edge
          b(3) = f_un_1(i0,i1) + D * (-thn_lower_y / dx(1)
     &      * p(i0,i1-1)
     &      - (2.d0 * eta_n - l_n) * dy_dy
     &      * thn(i0,i1-1) * un_1(i0,i1-1)
     &      - eta_n * dx_dx * thn_iphalf_jmhalf
     &      * un_1(i0+1,i1)
     &      - eta_n * dx_dx * thn_imhalf_jmhalf
     &      * un_1(i0-1,i1)
     &      + eta_n * dx_dy * thn_iphalf_jmhalf
     &      * un_0(i0+1,i1-1)
     &      - eta_n * dx_dy * thn_imhalf_jmhalf
     &      * un_0(i0,i1-1)
     &      - l_n * dx_dy * thn(i0,i1-1)
     &      * (un_0(i0+1,i1-1) - un_0(i0,i1-1)))
c
          ! solvent at south edge
          b(7) = f_us_1(i0,i1)+D*(-ths_lower_y / dx(1)
     &      * p(i0,i1-1)
     &      - (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1-1)) * us_1(i0,i1-1)
     &      - eta_s * dx_dx * ths_iphalf_jmhalf
     &      * us_1(i0+1,i1)
     &      - eta_s * dx_dx * ths_imhalf_jmhalf
     &      * us_1(i0-1,i1)
     &      + eta_s * dx_dy * ths_iphalf_jmhalf
     &      * us_0(i0+1,i1-1)
     &      - eta_s * dx_dy * ths_imhalf_jmhalf
     &      * us_0(i0,i1-1)
     &      - l_s * dx_dy * toThs(thn(i0,i1-1))
     &      * (us_0(i0+1,i1-1) - us_0(i0,i1-1)))
c
          ! network at north edge
          b(4) = f_un_1(i0,i1+1) + D * (thn_upper_y / dx(1)
     &      * p(i0,i1+1)
     &      - (2.d0 * eta_n - l_n) * dy_dy * thn(i0,i1+1)
     &      * un_1(i0,i1+2)
     &      - eta_n * dx_dx * thn_iphalf_jphalf
     &      * un_1(i0+1,i1+1)
     &      - eta_n * dx_dx * thn_imhalf_jphalf
     &      * un_1(i0-1,i1+1)
     &      - eta_n * dx_dy * thn_iphalf_jphalf
     &      * un_0(i0+1,i1+1)
     &      + eta_n * dx_dy * thn_imhalf_jphalf
     &      * un_0(i0,i1 + 1)
     &      + l_n * dx_dy * thn(i0,i1+1)
     &      * (un_0(i0+1,i1+1) - un_0(i0,i1 + 1)))
c
          ! solvent at north edge
          b(8) = f_us_1(i0,i1+1) + D * (ths_upper_y / dx(1)
     &      * p(i0,i1+1)
     &      - (2.d0 * eta_s - l_s) * dy_dy
     &      * toThs(thn(i0,i1+1)) * us_1(i0,i1+2)
     &      - eta_s * dx_dx * ths_iphalf_jphalf
     &      * us_1(i0+1,i1+1)
     &      - eta_s * dx_dx * ths_imhalf_jphalf
     &      * us_1(i0-1,i1+1)
     &      - eta_s * dx_dy * ths_iphalf_jphalf
     &      * us_0(i0+1,i1+1)
     &      + eta_s * dx_dy * ths_imhalf_jphalf
     &      * us_0(i0,i1 + 1)
     &      + l_s * dx_dy * toThs(thn(i0,i1+1))
     &      * (us_0(i0+1,i1+1) - us_0(i0,i1 + 1)))
c
          ! pressure at cell center
          b(9) = f_p(i0,i1)
c
          ! solve the system Ax = b, overwriting b with x
          call lu_solve(A_box, b)
c
          un_0(i0,i1) = (1.d0-w)*un_0(i0,i1) + w*b(1);
          un_0(i0+1,i1) = (1.d0-w)*un_0(i0+1,i1) + w*b(2);
          un_1(i0,i1) = (1.d0-w)*un_1(i0,i1) + w*b(3);
          un_1(i0,i1+1) = (1.d0-w)*un_1(i0,i1+1) + w*b(4);
          us_0(i0,i1) = (1.d0-w)*us_0(i0,i1) + w*b(5);
          us_0(i0+1,i1) = (1.d0-w)*us_0(i0+1,i1) + w*b(6);
          us_1(i0,i1) = (1.d0-w)*us_1(i0,i1) + w*b(7);
          us_1(i0,i1+1) = (1.d0-w)*us_1(i0,i1+1) + w*b(8);
          p(i0,i1) = (1.d0-w)*p(i0,i1) + w*b(9);
c
          end if
        enddo
      enddo
c
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Same as rbgs, but allows for a mask index to indicate which side
c       indices should not be changed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rbgs_mask_var_xi(
     &        dx, ilow0, iup0,
     &        ilow1, iup1,
     &        un_0, un_1, un_gcw,
     &        us_0, us_1, us_gcw,
     &        p, p_gcw, f_p, f_p_gcw,
     &        f_un_0, f_un_1, f_un_gcw,
     &        f_us_0, f_us_1, f_us_gcw,
     &        thn, thn_gcw, eta_n, eta_s, l_n, l_s,
     &        xi_0, xi_1, xi_gcw,
     &        w, C, D, red_or_black,
     &        mask_0, mask_1, mask_gcw)
c
      use my_subs
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0
      integer ilow1,  iup1
      integer un_gcw, us_gcw, p_gcw, f_un_gcw, f_us_gcw
      integer f_p_gcw, thn_gcw, red_or_black
      double precision eta_n, eta_s, l_n, l_s, w, C, D
c
      double precision thn(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw)
c
      double precision un_0(ilow0-un_gcw:iup0+un_gcw+1,
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)
c
      double precision us_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
      double precision p(ilow0-p_gcw:iup0+p_gcw,
     &          ilow1-p_gcw:iup1+p_gcw)

      double precision f_un_0(ilow0-f_un_gcw:iup0+f_un_gcw+1,
     &          ilow1-f_un_gcw:iup1+f_un_gcw)
      double precision f_un_1(ilow0-f_un_gcw:iup0+f_un_gcw,
     &          ilow1-f_un_gcw:iup1+f_un_gcw+1)
c
      double precision f_us_0(ilow0-f_us_gcw:iup0+f_us_gcw+1,
     &          ilow1-f_us_gcw:iup1+f_us_gcw)
      double precision f_us_1(ilow0-f_us_gcw:iup0+f_us_gcw,
     &          ilow1-f_us_gcw:iup1+f_us_gcw+1)
c
      double precision f_p(ilow0-f_p_gcw:iup0+f_p_gcw,
     &          ilow1-f_p_gcw:iup1+f_p_gcw)

      integer xi_gcw
      double precision xi_0(ilow0-xi_gcw:iup0+xi_gcw+1,
     &                           ilow1-xi_gcw:iup1+xi_gcw)
      double precision xi_1(ilow0-xi_gcw:iup0+xi_gcw,
     &                           ilow1-xi_gcw:iup1+xi_gcw+1)

      integer mask_gcw
      integer mask_0(ilow0-mask_gcw:iup0+mask_gcw+1,
     &               ilow1-mask_gcw:iup1+mask_gcw)
      integer mask_1(ilow0-mask_gcw:iup0+mask_gcw,
     &               ilow1-mask_gcw:iup1+mask_gcw+1)
c
      double precision thn_lower_x, thn_lower_y
      double precision thn_upper_x, thn_upper_y
      double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
      double precision thn_iphalf_jphalf, thn_iphalf_jmhalf
      double precision ths_lower_x, ths_lower_y
      double precision ths_upper_x, ths_upper_y
      double precision ths_imhalf_jphalf, ths_imhalf_jmhalf
      double precision ths_iphalf_jphalf, ths_iphalf_jmhalf
c
      double precision A_box(9,9)
      double precision b(9)
      double precision dx_dx, dx_dy, dy_dy
c
      integer i0, i1
c
      dx_dx = 1.d0 / (dx(0) * dx(0))
      dx_dy = 1.d0 / (dx(0) * dx(1))
      dy_dy = 1.d0 / (dx(1) * dx(1))
      do i1 = ilow1, iup1
        do i0 = ilow0, iup0  ! same loop as the c++ code (currently this is just GS)
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then

          thn_lower_x = 0.5d0*(thn(i0,i1)+thn(i0-1,i1))
          thn_upper_x = 0.5d0*(thn(i0,i1)+thn(i0+1,i1))
          thn_lower_y = 0.5d0*(thn(i0,i1)+thn(i0,i1-1))
          thn_upper_y = 0.5d0*(thn(i0,i1)+thn(i0,i1+1))

          ths_lower_x = toThs(thn_lower_x)
          ths_upper_x = toThs(thn_upper_x)
          ths_lower_y = toThs(thn_lower_y)
          ths_upper_y = toThs(thn_upper_y)

      ! calculate thn at corners
          thn_imhalf_jphalf = 0.25d0*(thn(i0-1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0-1,i1+1))
          thn_imhalf_jmhalf = 0.25d0*(thn(i0,i1)+thn(i0-1,i1)
     &                           +thn(i0,i1-1)+thn(i0-1,i1-1))
          thn_iphalf_jphalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1+1)+thn(i0+1,i1+1))
          thn_iphalf_jmhalf = 0.25d0*(thn(i0+1,i1)+thn(i0,i1)
     &                           +thn(i0,i1-1)+thn(i0+1,i1-1))

          ths_imhalf_jphalf = toThs(thn_imhalf_jphalf)
          ths_imhalf_jmhalf = toThs(thn_imhalf_jmhalf)
          ths_iphalf_jphalf = toThs(thn_iphalf_jphalf)
          ths_iphalf_jmhalf = toThs(thn_iphalf_jmhalf)

          ! network at west edge
          if (mask_0(i0,i1) .eq. 1) then
            A_box(1, 1) = 1.d0
            A_box(1, 2) = 0.d0
            A_box(1, 3) = 0.d0
            A_box(1, 4) = 0.d0
            A_box(1, 5) = 0.d0
            A_box(1, 6) = 0.d0
            A_box(1, 7) = 0.d0
            A_box(1, 8) = 0.d0
            A_box(1, 9) = 0.d0
            A_box(5, 5) = 1.d0
            A_box(5, 6) = 0.d0
            A_box(5, 7) = 0.d0
            A_box(5, 8) = 0.d0
            A_box(5, 1) = 0.d0
            A_box(5, 2) = 0.d0
            A_box(5, 3) = 0.d0
            A_box(5, 4) = 0.d0
            A_box(5, 9) = 0.d0
            b(1) = un_0(i0,i1)
            b(5) = us_0(i0,i1)
          else
            A_box(1,1) = D*(
     &        (2.d0 * eta_n - l_n) * dx_dx
     &        * (-thn(i0,i1) - thn(i0-1,i1))
     &        - eta_n * dy_dy
     &        * (thn_imhalf_jmhalf + thn_imhalf_jphalf)
     &        - xi_0(i0,i1))
     &        + C * thn_lower_x
            A_box(1,2) = D * (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0,i1)
            A_box(1,3) = D * (l_n * thn(i0,i1)
     &        - eta_n * thn_imhalf_jmhalf) * dx_dy
            A_box(1, 4) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_imhalf_jphalf) * dx_dy
            A_box(1,5) = D * xi_0(i0,i1)
            A_box(1,6) = 0.0
            A_box(1,7) = 0.0
            A_box(1,8) = 0.0
            A_box(1,9) = D * (-thn_lower_x) / dx(0)

            A_box(5,1) = D * xi_0(i0,i1)
            A_box(5,2) = 0.0
            A_box(5,3) = 0.0
            A_box(5,4) = 0.0
            A_box(5,5) = D*(
     &        (2.d0 * eta_s - l_s) * dx_dx
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0-1,i1)))
     &        - eta_s * dy_dy
     &        * (ths_imhalf_jmhalf + ths_imhalf_jphalf)
     &        - xi_0(i0,i1))
     &        + C * ths_lower_x
            A_box(5,6) = D * (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0,i1))
            A_box(5,7) = D * (-eta_s * ths_imhalf_jmhalf
     &        + l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(5,8) = D *(eta_s * ths_imhalf_jphalf
     &        - l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(5,9) = D * (-ths_lower_x) / dx(0)
c
            b(1) = f_un_0(i0,i1) + D * (-thn_lower_x / dx(0)
     &        * p(i0-1,i1)
     &        - (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0-1,i1) * un_0(i0-1,i1)
     &        - eta_n * dy_dy
     &        * thn_imhalf_jphalf * un_0(i0,i1+1)
     &        - eta_n * dy_dy
     &        * thn_imhalf_jmhalf * un_0(i0,i1-1)
     &        + eta_n * dx_dy
     &        * thn_imhalf_jphalf * un_1(i0-1,i1+1)
     &        - eta_n * dx_dy
     &        * thn_imhalf_jmhalf * un_1(i0-1,i1)
     &        - l_n * dx_dy * thn(i0-1,i1)
     &        * (un_1(i0-1,i1+1) - un_1(i0-1,i1)))
c
            ! solvent at west edge
            b(5) = f_us_0(i0,i1) + D * (-ths_lower_x / dx(0)
     &        * p(i0-1,i1)
     &        - (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0-1,i1)) * us_0(i0-1,i1)
     &        - eta_s * dy_dy
     &        * ths_imhalf_jphalf * us_0(i0,i1+1)
     &        - eta_s * dy_dy
     &        * ths_imhalf_jmhalf * us_0(i0,i1-1)
     &        + eta_s * dx_dy
     &        * ths_imhalf_jphalf * us_1(i0-1,i1+1)
     &        - eta_s * dx_dy
     &        * ths_imhalf_jmhalf * us_1(i0-1,i1)
     &        - l_s * dx_dy
     &        * toThs(thn(i0-1,i1))
     &        * (us_1(i0-1,i1+1) - us_1(i0-1,i1)))
          endif
c
          ! network at east edge
          if (mask_0(i0+1,i1) .eq. 1) then
            A_box(2, 1) = 0.d0
            A_box(2, 2) = 1.d0
            A_box(2, 3) = 0.d0
            A_box(2, 4) = 0.d0
            A_box(2, 5) = 0.d0
            A_box(2, 7) = 0.d0
            A_box(2, 8) = 0.d0
            A_box(2, 6) = 0.d0
            A_box(2, 9) = 0.d0
c
            A_box(6, 5) = 0.d0
            A_box(6, 6) = 1.d0
            A_box(6, 7) = 0.d0
            A_box(6, 8) = 0.d0
            A_box(6, 1) = 0.d0
            A_box(6, 3) = 0.d0
            A_box(6, 4) = 0.d0
            A_box(6, 2) = 0.d0
            A_box(6, 9) = 0.d0

            b(2) = un_0(i0+1,i1)
            b(6) = us_0(i0+1,i1)
          else
            A_box(2,1) = D * (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0,i1)
            A_box(2,2) = D * (
     &        (2.d0 * eta_n - l_n) * dx_dx
     &        * (-thn(i0+1,i1) - thn(i0,i1))
     &        - eta_n * dy_dy
     &        * (thn_iphalf_jphalf + thn_iphalf_jmhalf)
     &        - xi_0(i0+1,i1))
     &        + C * thn_upper_x
            A_box(2,3) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_iphalf_jmhalf) * dx_dy
            A_box(2,4) = D * (l_n * thn(i0,i1)
     &      - eta_n * thn_iphalf_jphalf) * dx_dy
            A_box(2,6) = D * xi_0(i0+1,i1)
            A_box(2,5) = 0.0
            A_box(2,7) = 0.0
            A_box(2,8) = 0.0
            A_box(2,9) = D * thn_upper_x / dx(0)
c
            A_box(6,1) = 0.0
            A_box(6,2) = D * xi_0(i0+1,i1)
            A_box(6,3) = 0.0
            A_box(6,4) = 0.0
            A_box(6,5) = D * (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0,i1))
            A_box(6,6) = D * ((2.d0 * eta_s - l_s) * dx_dx
     &        * (-toThs(thn(i0+1,i1)) - toThs(thn(i0,i1)))
     &        - eta_s * dy_dy
     &        * (ths_iphalf_jphalf + ths_iphalf_jmhalf)
     &        - xi_0(i0+1,i1))
     &        + C * ths_upper_x
            A_box(6,7) = D * (eta_s * ths_iphalf_jmhalf
     &        - l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(6,8) = D * (-eta_s * ths_iphalf_jphalf
     &        + l_s * toThs(thn(i0,i1))) * dx_dy
            A_box(6,9) = D * ths_upper_x / dx(0)

                      ! network at east edge
            b(2) = f_un_0(i0+1,i1) + D * (thn_upper_x / dx(0)
     &        * p(i0+1,i1)
     &        - (2.d0 * eta_n - l_n) * dx_dx
     &        * thn(i0+1,i1) * un_0(i0+2,i1)
     &        - eta_n * dy_dy * thn_iphalf_jphalf
     &        * un_0(i0+1,i1+1)
     &        - eta_n * dy_dy * thn_iphalf_jmhalf
     &        * un_0(i0+1,i1-1)
     &        - eta_n * dx_dy * thn_iphalf_jphalf
     &        * un_1(i0+1,i1+1)
     &        + eta_n * dx_dy * thn_iphalf_jmhalf
     &        * un_1(i0+1,i1)
     &        + l_n * dx_dy * thn(i0+1,i1)
     &        * (un_1(i0+1,i1+1) - un_1(i0+1,i1)))
c
          ! solvent at east edge
            b(6) = f_us_0(i0+1,i1) + D * (ths_upper_x / dx(0)
     &        * p(i0+1,i1)
     &        - (2.d0 * eta_s - l_s) * dx_dx
     &        * toThs(thn(i0+1,i1)) * us_0(i0+2,i1)
     &        - eta_s * dy_dy * ths_iphalf_jphalf
     &        * us_0(i0+1,i1+1)
     &        - eta_s * dy_dy * ths_iphalf_jmhalf
     &        * us_0(i0+1,i1-1)
     &        - eta_s * dx_dy * ths_iphalf_jphalf
     &        * us_1(i0+1,i1+1)
     &        + eta_s * dx_dy * ths_iphalf_jmhalf
     &        * us_1(i0+1,i1)
     &        + l_s * dx_dy * toThs(thn(i0+1,i1))
     &        * (us_1(i0+1,i1+1) - us_1(i0+1,i1)))
          endif
c
          ! network at south edge
          if (mask_1(i0,i1) .eq. 1) then
            A_box(3, 1) = 0.d0
            A_box(3, 2) = 0.d0
            A_box(3, 3) = 1.d0
            A_box(3, 4) = 0.d0
            A_box(3, 5) = 0.d0
            A_box(3, 6) = 0.d0
            A_box(3, 8) = 0.d0
            A_box(3, 7) = 0.d0
            A_box(3, 9) = 0.d0
c
            A_box(7, 5) = 0.d0
            A_box(7, 6) = 0.d0
            A_box(7, 7) = 1.d0
            A_box(7, 8) = 0.d0
            A_box(7, 1) = 0.d0
            A_box(7, 2) = 0.d0
            A_box(7, 4) = 0.d0
            A_box(7, 3) = 0.d0
            A_box(7, 9) = 0.d0

            b(3) = un_1(i0,i1)
            b(7) = us_1(i0,i1)
          else
            A_box(3,1) = D * (l_n * thn(i0,i1)
     &        - eta_n * thn_imhalf_jmhalf) * dx_dy
            A_box(3,2) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_iphalf_jmhalf) * dx_dy
            A_box(3,3) = D*((2.d0 * eta_n - l_n) * dy_dy
     &        * (-thn(i0,i1) - thn(i0,i1-1))
     &        - eta_n * dx_dx
     &        * (thn_iphalf_jmhalf + thn_imhalf_jmhalf)
     &        - xi_1(i0,i1))
     &        + C * thn_lower_y
            A_box(3,4) = D * (2.d0*eta_n-l_n) * dy_dy
     &        * thn(i0,i1)
            A_box(3,5) = 0.0
            A_box(3,6) = 0.0
            A_box(3,7) = D * xi_1(i0,i1)
            A_box(3,8) = 0.0
            A_box(3,9) = D * (-thn_lower_y) / dx(1)
c
            A_box(7,1) = 0.0
            A_box(7,2) = 0.0
            A_box(7,3) = D * xi_1(i0,i1)
            A_box(7,4) = 0.0
            A_box(7,5) = D * (l_s * toThs(thn(i0,i1))
     &        - eta_s * ths_imhalf_jmhalf) * dx_dy
            A_box(7,6) = D * (-l_s * toThs(thn(i0,i1))
     &        + eta_s * ths_iphalf_jmhalf) * dx_dy
            A_box(7,7) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1-1)))
     &        - eta_s * dx_dx
     &        * (ths_iphalf_jmhalf + ths_imhalf_jmhalf)
     &        - xi_1(i0,i1))
     &        + C * ths_lower_y
            A_box(7,8) = D * (2.d0*eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1))
            A_box(7,9) = D * (-ths_lower_y) / dx(1)
c
            b(3) = f_un_1(i0,i1) + D * (-thn_lower_y / dx(1)
     &        * p(i0,i1-1)
     &        - (2.d0 * eta_n - l_n) * dy_dy
     &        * thn(i0,i1-1) * un_1(i0,i1-1)
     &        - eta_n * dx_dx * thn_iphalf_jmhalf
     &        * un_1(i0+1,i1)
     &        - eta_n * dx_dx * thn_imhalf_jmhalf
     &        * un_1(i0-1,i1)
     &        + eta_n * dx_dy * thn_iphalf_jmhalf
     &        * un_0(i0+1,i1-1)
     &        - eta_n * dx_dy * thn_imhalf_jmhalf
     &        * un_0(i0,i1-1)
     &        - l_n * dx_dy * thn(i0,i1-1)
     &        * (un_0(i0+1,i1-1) - un_0(i0,i1-1)))
c
            b(7) = f_us_1(i0,i1)+D*(-ths_lower_y / dx(1)
     &        * p(i0,i1-1)
     &        - (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1-1)) * us_1(i0,i1-1)
     &        - eta_s * dx_dx * ths_iphalf_jmhalf
     &        * us_1(i0+1,i1)
     &        - eta_s * dx_dx * ths_imhalf_jmhalf
     &        * us_1(i0-1,i1)
     &        + eta_s * dx_dy * ths_iphalf_jmhalf
     &        * us_0(i0+1,i1-1)
     &        - eta_s * dx_dy * ths_imhalf_jmhalf
     &        * us_0(i0,i1-1)
     &        - l_s * dx_dy * toThs(thn(i0,i1-1))
     &        * (us_0(i0+1,i1-1) - us_0(i0,i1-1)))
          endif
c
          ! network at north edge
          if (mask_1(i0,i1+1) .eq. 1) then
            A_box(4, 1) = 0.d0
            A_box(4, 2) = 0.d0
            A_box(4, 3) = 0.d0
            A_box(4, 4) = 1.d0
            A_box(4, 5) = 0.d0
            A_box(4, 6) = 0.d0
            A_box(4, 7) = 0.d0
            A_box(4, 8) = 0.d0
            A_box(4, 9) = 0.d0
c
            A_box(8, 5) = 0.d0
            A_box(8, 6) = 0.d0
            A_box(8, 7) = 0.d0
            A_box(8, 8) = 1.d0
            A_box(8, 1) = 0.d0
            A_box(8, 2) = 0.d0
            A_box(8, 3) = 0.d0
            A_box(8, 4) = 0.d0
            A_box(8, 9) = 0.d0

            b(4) = un_1(i0,i1+1)
            b(8) = us_1(i0,i1+1)
          else
            A_box(4,1) = D * (-l_n * thn(i0,i1)
     &        + eta_n * thn_imhalf_jphalf) * dx_dy
            A_box(4,2) = D * (l_n * thn(i0,i1)
     &        - eta_n *thn_iphalf_jphalf) * dx_dy
            A_box(4,3) = D * (2.d0 * eta_n - l_n) * dy_dy
     &        * (thn(i0,i1))
            A_box(4,4) = D * ((2.d0 * eta_n - l_n) * dy_dy
     &        * (-thn(i0,i1) - thn(i0,i1+1))
     &        - eta_n * dx_dx
     &        * (thn_iphalf_jphalf + thn_imhalf_jphalf)
     &        - xi_1(i0,i1+1))
     &        + C * thn_upper_y
            A_box(4,5) = 0.0
            A_box(4,6) = 0.0
            A_box(4,7) = 0.0
            A_box(4,8) = D * xi_1(i0,i1+1)
            A_box(4,9) = D * thn_upper_y / dx(1)
c
            A_box(8,1) = 0.0
            A_box(8,2) = 0.0
            A_box(8,3) = 0.0
            A_box(8,4) = D * xi_1(i0,i1+1)
            A_box(8,5) = D * (-l_s * toThs(thn(i0,i1))
     &        + eta_s * ths_imhalf_jphalf) * dx_dy
            A_box(8,6) = D * (l_s * toThs(thn(i0,i1))
     &        - eta_s * ths_iphalf_jphalf) * dx_dy
            A_box(8,7) = D * (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1))
            A_box(8,8) = D * ((2.d0 * eta_s - l_s) * dy_dy
     &        * (-toThs(thn(i0,i1)) - toThs(thn(i0,i1+1)))
     &        - eta_s * dx_dx
     &        * (ths_iphalf_jphalf + ths_imhalf_jphalf)
     &        - xi_1(i0,i1+1))
     &        + C * ths_upper_y
            A_box(8,9) = D * ths_upper_y / dx(1)

                      ! network at north edge
            b(4) = f_un_1(i0,i1+1) + D * (thn_upper_y / dx(1)
     &        * p(i0,i1+1)
     &        - (2.d0 * eta_n - l_n) * dy_dy * thn(i0,i1+1)
     &        * un_1(i0,i1+2)
     &        - eta_n * dx_dx * thn_iphalf_jphalf
     &        * un_1(i0+1,i1+1)
     &        - eta_n * dx_dx * thn_imhalf_jphalf
     &        * un_1(i0-1,i1+1)
     &        - eta_n * dx_dy * thn_iphalf_jphalf
     &        * un_0(i0+1,i1+1)
     &        + eta_n * dx_dy * thn_imhalf_jphalf
     &        * un_0(i0,i1 + 1)
     &        + l_n * dx_dy * thn(i0,i1+1)
     &        * (un_0(i0+1,i1+1) - un_0(i0,i1 + 1)))
c
          ! solvent at north edge
            b(8) = f_us_1(i0,i1+1) + D * (ths_upper_y / dx(1)
     &        * p(i0,i1+1)
     &        - (2.d0 * eta_s - l_s) * dy_dy
     &        * toThs(thn(i0,i1+1)) * us_1(i0,i1+2)
     &        - eta_s * dx_dx * ths_iphalf_jphalf
     &        * us_1(i0+1,i1+1)
     &        - eta_s * dx_dx * ths_imhalf_jphalf
     &        * us_1(i0-1,i1+1)
     &        - eta_s * dx_dy * ths_iphalf_jphalf
     &        * us_0(i0+1,i1+1)
     &        + eta_s * dx_dy * ths_imhalf_jphalf
     &        * us_0(i0,i1 + 1)
     &        + l_s * dx_dy * toThs(thn(i0,i1+1))
     &        * (us_0(i0+1,i1+1) - us_0(i0,i1 + 1)))
          endif
c
          ! incompressible constrain term at center
          A_box(9, 1) = -thn_lower_x / dx(0)
          A_box(9, 2) = thn_upper_x / dx(0)
          A_box(9, 3) = -thn_lower_y / dx(1)
          A_box(9, 4) = thn_upper_y / dx(1)
          A_box(9, 5) = -toThs(thn_lower_x) / dx(0)
          A_box(9, 6) = toThs(thn_upper_x) / dx(0)
          A_box(9, 7) = -toThs(thn_lower_y) / dx(1)
          A_box(9, 8) = toThs(thn_upper_y) / dx(1)
          A_box(9, 9) = 0.0
          ! pressure at cell center
          b(9) = f_p(i0,i1)
c
          ! solve the system Ax = b, overwriting b with x
          call lu_solve(A_box, b)
c
          if (mask_0(i0,i1) .eq. 0) then
            un_0(i0,i1) = (1.d0-w)*un_0(i0,i1) + w*b(1);
            us_0(i0,i1) = (1.d0-w)*us_0(i0,i1) + w*b(5);
          endif
          if (mask_0(i0+1,i1) .eq. 0) then
            un_0(i0+1,i1) = (1.d0-w)*un_0(i0+1,i1) + w*b(2);
            us_0(i0+1,i1) = (1.d0-w)*us_0(i0+1,i1) + w*b(6);
          endif
          if (mask_1(i0,i1) .eq. 0) then
            un_1(i0,i1) = (1.d0-w)*un_1(i0,i1) + w*b(3);
            us_1(i0,i1) = (1.d0-w)*us_1(i0,i1) + w*b(7);
          endif
          if (mask_1(i0,i1+1) .eq. 0) then
            un_1(i0,i1+1) = (1.d0-w)*un_1(i0,i1+1) + w*b(4);
            us_1(i0,i1+1) = (1.d0-w)*us_1(i0,i1+1) + w*b(8);
          endif
          p(i0,i1) = (1.d0-w)*p(i0,i1) + w*b(9);
c
          end if
        enddo
      enddo
c
      end subroutine

c     Solve the system Ax = b.
c     Inputs: the matrix A and vector b.
c     Outputs: LU decomposition of A and solution in b
c     Does not pivot nor checks for zero diagonals. This should be
c     sufficient for our system, which is diagonally dominant
c     (except for the last row).
      subroutine lu_solve(A, b)
      integer N
      parameter (N=9)
      double precision A(0:(N-1),0:(N-1))
      double precision b(0:(N-1))

      integer i, j, k

      do i = 0,(N-1)
        do j = (i+1),(N-1)
          A(j,i) = A(j,i) / A(i,i)

          do k = (i+1),(N-1)
            A(j,k) = A(j,k) - A(j,i) * A(i,k)
          enddo
        enddo
      enddo

      ! Now solve
      do i = 0,(N-1)

        do k = 0, (i-1)
          b(i) = b(i) - A(i,k) * b(k)
        enddo
      enddo

      do i = (N-1),0,-1
        do k = (i+1),(N-1)
          b(i) = b(i) - A(i,k) * b(k)
        enddo

        b(i) = b(i) / A(i,i)
      enddo
      end
