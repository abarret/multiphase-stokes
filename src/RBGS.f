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
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        p_data, p_gcw, f_p_data, f_p_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_gcw, eta_n, eta_s,
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
      double precision eta_n, eta_s, nu_n, nu_s, xi, w, C, D
c      
      double precision thn_data(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
c
      double precision un_data_0(ilow0-un_gcw:iup0+un_gcw+1,  
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_data_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)
c
      double precision us_data_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_data_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
      double precision p_data(ilow0-p_gcw:iup0+p_gcw,
     &          ilow1-p_gcw:iup1+p_gcw)

      double precision f_un_data_0(ilow0-f_un_gcw:iup0+f_un_gcw+1,
     &          ilow1-f_un_gcw:iup1+f_un_gcw)
      double precision f_un_data_1(ilow0-f_un_gcw:iup0+f_un_gcw,
     &          ilow1-f_un_gcw:iup1+f_un_gcw+1)
c
      double precision f_us_data_0(ilow0-f_us_gcw:iup0+f_us_gcw+1,
     &          ilow1-f_us_gcw:iup1+f_us_gcw)
      double precision f_us_data_1(ilow0-f_us_gcw:iup0+f_us_gcw,
     &          ilow1-f_us_gcw:iup1+f_us_gcw+1)
c
      double precision f_p_data(ilow0-f_p_gcw:iup0+f_p_gcw,
     &          ilow1-f_p_gcw:iup1+f_p_gcw)
c
      double precision thn_lower_x, thn_lower_y
      double precision thn_upper_x, thn_upper_y
      double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
      double precision thn_iphalf_jphalf, thn_iphalf_jmhalf
c    
      double precision A_box(9,9)
      double precision b(9)
      integer ipiv(9), info
c
      integer i0, i1
c    
      do i1 = ilow1, iup1   
        do i0 = ilow0, iup0  ! same loop as the c++ code (currently this is just GS)
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then
          
          ! calculate thn at sides
          thn_lower_x = 0.5d0*(thn_data(i0,i1)+thn_data(i0-1,i1))  ! thn(i-1/2, j)
          thn_upper_x = 0.5d0*(thn_data(i0,i1)+thn_data(i0+1,i1))  ! thn(i+1/2, j)
          thn_lower_y = 0.5d0*(thn_data(i0,i1)+thn_data(i0,i1-1))  ! thn(i, j-1/2)
          thn_upper_y = 0.5d0*(thn_data(i0,i1)+thn_data(i0,i1+1))  ! thn(i, j+1/2)

          ! calculate thn at corners
          thn_imhalf_jphalf = 0.25d0*(thn_data(i0-1, i1)+thn_data(i0,i1)
     &                     +thn_data(i0,i1+1)+thn_data(i0-1,i1+1))   ! thn(i-1/2, j+1/2)
          thn_imhalf_jmhalf = 0.25d0*(thn_data(i0, i1)+thn_data(i0-1,i1)
     &                     +thn_data(i0,i1-1)+thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)
          thn_iphalf_jphalf = 0.25d0*(thn_data(i0+1, i1)+thn_data(i0,i1)
     &                     +thn_data(i0,i1+1)+thn_data(i0+1,i1+1))   ! thn(i+1/2, j+1/2)
          thn_iphalf_jmhalf = 0.25d0*(thn_data(i0+1, i1)+thn_data(i0,i1)
     &                     +thn_data(i0,i1-1)+thn_data(i0+1,i1-1))   ! thn(i+1/2, j-1/2)

          ! network at west edge
          A_box(1, 1) = eta_n / (dx(0)*dx(0)) * 
     &    (-thn_data(i0,i1)-thn_data(i0-1,i1))-eta_n/(dx(1)*dx(1))*    
     &    (thn_imhalf_jmhalf + thn_imhalf_jphalf)
     &    - (xi / nu_n) * thn_lower_x*toThs(thn_lower_x) 
          A_box(1,1) = C * thn_lower_x + D * A_box(1, 1)
          A_box(1, 2) = D * eta_n / (dx(0)*dx(0))*thn_data(i0,i1)
          A_box(1, 3) = D * eta_n / (dx(0) * dx(1)) * 
     &               (thn_data(i0,i1)-thn_imhalf_jmhalf)
          A_box(1, 4) = D * eta_n / (dx(0) * dx(1)) *
     &               (thn_imhalf_jphalf - thn_data(i0,i1))
          A_box(1, 5) = D *(xi / nu_n)*thn_lower_x*toThs(thn_lower_x)
          A_box(1, 6) = 0.0
          A_box(1, 7) = 0.0
          A_box(1, 8) = 0.0
          A_box(1, 9) = D * (-thn_lower_x) / dx(0)
c
          A_box(5, 5) = eta_s / (dx(0) * dx(0)) * 
     &    (-toThs(thn_data(i0,i1))-toThs(thn_data(i0-1,i1))) 
     &    - eta_s/(dx(1)*dx(1)) *
     &    (toThs(thn_imhalf_jmhalf)+toThs(thn_imhalf_jphalf)) 
     &    -(xi / nu_s)*thn_lower_x*toThs(thn_lower_x) 
          A_box(5, 5) = C * toThs(thn_lower_x) + D * A_box(5, 5)
          A_box(5, 6) = D * eta_s/(dx(0)*dx(0)) * toThs(thn_data(i0,i1))
          A_box(5, 7) = D * eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_data(i0,i1))-toThs(thn_imhalf_jmhalf))
          A_box(5, 8) = D *eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_imhalf_jphalf) - toThs(thn_data(i0,i1)))
          A_box(5, 1) = D * xi/nu_s * thn_lower_x * toThs(thn_lower_x)
          A_box(5, 2) = 0.0
          A_box(5, 3) = 0.0
          A_box(5, 4) = 0.0
          A_box(5, 9) = D * (-toThs(thn_lower_x)) / dx(0)
c
          ! network at east edge
          A_box(2, 1) = D * eta_n / (dx(0) * dx(0)) * thn_data(i0,i1)
          A_box(2, 2) = eta_n / (dx(0) * dx(0)) * 
     &    (-thn_data(i0+1,i1) - thn_data(i0,i1)) -
     &    eta_n/(dx(1)*dx(1))*(thn_iphalf_jphalf+thn_iphalf_jmhalf) 
     &    -(xi / nu_n) * thn_upper_x * toThs(thn_upper_x) 
          A_box(2, 2) = C * thn_upper_x + D * A_box(2, 2)
          A_box(2, 3) = D * eta_n / (dx(1) * dx(0)) * 
     &    (thn_iphalf_jmhalf - thn_data(i0,i1))
          A_box(2, 4) = D * eta_n / (dx(1) * dx(0)) * 
     &    (thn_data(i0,i1)-thn_iphalf_jphalf)
          A_box(2, 5) = 0.0
          A_box(2, 7) = 0.0
          A_box(2, 8) = 0.0
          A_box(2, 6) = D * xi/nu_n * thn_upper_x * toThs(thn_upper_x)
          A_box(2, 9) = D * thn_upper_x / dx(0)
c
          A_box(6, 5) = D * eta_s/(dx(0)*dx(0)) * toThs(thn_data(i0,i1))
          A_box(6, 6) = eta_s/(dx(0)*dx(0)) * 
     &    (-toThs(thn_data(i0+1,i1))-toThs(thn_data(i0,i1))) 
     &    - eta_s/(dx(1)*dx(1)) * 
     &    (toThs(thn_iphalf_jphalf)+toThs(thn_iphalf_jmhalf)) 
     &    - xi/nu_s * thn_upper_x * toThs(thn_upper_x) 
          A_box(6, 6) = C * toThs(thn_upper_x) + D * A_box(6, 6)
          A_box(6, 7) = D * eta_s / (dx(1)*dx(0)) * 
     &    (toThs(thn_iphalf_jmhalf)-toThs(thn_data(i0,i1)))
          A_box(6, 8) = D * eta_s/(dx(1)*dx(0)) *
     &    (toThs(thn_data(i0,i1))-toThs(thn_iphalf_jphalf))
          A_box(6, 1) = 0.0
          A_box(6, 3) = 0.0
          A_box(6, 4) = 0.0
          A_box(6, 2) = D * xi/nu_s * thn_upper_x * toThs(thn_upper_x)
          A_box(6, 9) = D * toThs(thn_upper_x) / dx(0)
c
          ! network at south edge
          A_box(3, 1) = D * eta_n / (dx(0) * dx(1)) * 
     &    (thn_data(i0,i1)-thn_imhalf_jmhalf)
          A_box(3, 2) = D * eta_n / (dx(0) * dx(1)) * 
     &    (thn_iphalf_jmhalf - thn_data(i0,i1))
          A_box(3, 3) = eta_n / (dx(1) * dx(1)) * 
     &    (-thn_data(i0,i1) - thn_data(i0,i1-1)) 
     &    -eta_n/(dx(0)*dx(0))*(thn_iphalf_jmhalf + thn_imhalf_jmhalf) 
     &    -xi/nu_n * thn_lower_y * toThs(thn_lower_y) 
          A_box(3, 3) = C * thn_lower_y + D * A_box(3, 3)
          A_box(3, 4) = D * eta_n / (dx(1) * dx(1)) * (thn_data(i0,i1))
          A_box(3, 5) = 0.0
          A_box(3, 6) = 0.0
          A_box(3, 8) = 0.0
          A_box(3, 7) = D * xi/nu_n * thn_lower_y * toThs(thn_lower_y)
          A_box(3, 9) = D * (-thn_lower_y) / dx(1)
c
          A_box(7, 5) = D * eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_data(i0,i1)) - toThs(thn_imhalf_jmhalf))
          A_box(7, 6) = D * eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_iphalf_jmhalf) - toThs(thn_data(i0,i1)))
          A_box(7, 7) = eta_s / (dx(1) * dx(1)) * 
     &    (-toThs(thn_data(i0,i1)) - toThs(thn_data(i0,i1-1))) 
     &    - eta_s / (dx(0) * dx(0)) * 
     &    (toThs(thn_iphalf_jmhalf) + toThs(thn_imhalf_jmhalf)) 
     &    - xi/nu_s * thn_lower_y * toThs(thn_lower_y) 
          A_box(7, 7) = C * toThs(thn_lower_y) + D * A_box(7, 7)
          A_box(7, 8) = D * eta_s /(dx(1)*dx(1))*toThs(thn_data(i0,i1))
          A_box(7, 1) = 0.0
          A_box(7, 2) = 0.0
          A_box(7, 4) = 0.0
          A_box(7, 3) = D * xi/nu_s * thn_lower_y * toThs(thn_lower_y)
          A_box(7, 9) = D * (-toThs(thn_lower_y)) / dx(1)
c
          ! network at north edge
          A_box(4, 1) = D * eta_n / (dx(0) * dx(1)) * 
     &    (thn_imhalf_jphalf - thn_data(i0,i1))
          A_box(4, 2) = D * eta_n / (dx(0) * dx(1)) * 
     &    (thn_data(i0,i1)-thn_iphalf_jphalf)
          A_box(4, 3) = D * eta_n / (dx(1) * dx(1)) * (thn_data(i0,i1))
          A_box(4, 4) = eta_n / (dx(1) * dx(1)) * 
     &    (-thn_data(i0,i1) - thn_data(i0,i1+1)) 
     &    - eta_n/(dx(0)*dx(0))*(thn_iphalf_jphalf+thn_imhalf_jphalf) 
     &    - xi/nu_n * thn_upper_y * toThs(thn_upper_y) 
          A_box(4, 4) = C * thn_upper_y + D * A_box(4, 4)
          A_box(4, 5) = 0.0
          A_box(4, 6) = 0.0
          A_box(4, 7) = 0.0
          A_box(4, 8) = D * xi/nu_n * thn_upper_y * toThs(thn_upper_y)
          A_box(4, 9) = D * thn_upper_y / dx(1)
c
          A_box(8, 5) = D * eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_imhalf_jphalf)-toThs(thn_data(i0,i1)))
          A_box(8, 6) = D * eta_s / (dx(0) * dx(1)) * 
     &    (toThs(thn_data(i0,i1))-toThs(thn_iphalf_jphalf))
          A_box(8, 7) = D * eta_s/(dx(1)*dx(1))*toThs(thn_data(i0,i1))
          A_box(8, 8) = eta_s / (dx(1) * dx(1)) * 
     &    (-toThs(thn_data(i0,i1))-toThs(thn_data(i0,i1+1))) 
     &    - eta_s / (dx(0) * dx(0)) * 
     &    (toThs(thn_iphalf_jphalf)+toThs(thn_imhalf_jphalf)) 
     &    - xi/nu_s * thn_upper_y * toThs(thn_upper_y) 
          A_box(8, 8) = C * toThs(thn_upper_y) + D * A_box(8, 8)
          A_box(8, 1) = 0.0
          A_box(8, 2) = 0.0
          A_box(8, 3) = 0.0
          A_box(8, 4) = D * xi/nu_s * thn_upper_y * toThs(thn_upper_y)
          A_box(8, 9) = D * toThs(thn_upper_y) / dx(1)
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
c
          ! network at west edge
          b(1) = f_un_data_0(i0,i1)+D*(-thn_lower_x/dx(0)*
     &        p_data(i0-1,i1) -
     &        eta_n / (dx(0)* dx(0)) * 
     &        thn_data(i0-1,i1) * un_data_0(i0-1,i1) -
     &        eta_n / (dx(1)* dx(1)) * 
     &        thn_imhalf_jphalf * un_data_0(i0,i1+1) -
     &        eta_n / (dx(1)* dx(1)) * 
     &        thn_imhalf_jmhalf * un_data_0(i0,i1-1) +
     &        eta_n / (dx(0)* dx(1)) * 
     &        thn_imhalf_jphalf * un_data_1(i0-1,i1+1) -
     &        eta_n / (dx(0)* dx(1)) * 
     &        thn_imhalf_jmhalf * un_data_1(i0-1,i1) -
     &        eta_n / (dx(0)* dx(1)) * thn_data(i0-1,i1) *
     &        (un_data_1(i0-1,i1+1) - un_data_1(i0-1,i1)))
c
          ! solvent at west edge
          b(5) = f_us_data_0(i0,i1)+D*(-toThs(thn_lower_x) / dx(0) *
     &      p_data(i0-1,i1) -
     &      eta_s / (dx(0)* dx(0)) * 
     &      toThs(thn_data(i0-1,i1)) * us_data_0(i0-1,i1) -
     &      eta_s / (dx(1)* dx(1)) * 
     &      toThs(thn_imhalf_jphalf) * us_data_0(i0,i1+1) -
     &      eta_s / (dx(1)* dx(1)) * 
     &      toThs(thn_imhalf_jmhalf) * us_data_0(i0,i1-1) +
     &      eta_s / (dx(0)* dx(1)) * 
     &      toThs(thn_imhalf_jphalf) * us_data_1(i0-1,i1+1) -
     &      eta_s / (dx(0)* dx(1)) * 
     &      toThs(thn_imhalf_jmhalf) * us_data_1(i0-1,i1) -
     &      eta_s / (dx(0)* dx(1)) * 
     &      toThs(thn_data(i0-1,i1)) *
     &      (us_data_1(i0-1,i1+1) - us_data_1(i0-1,i1)))
c
          ! network at east edge
          b(2) = f_un_data_0(i0+1,i1) + D*(thn_upper_x / dx(0) *
     &        p_data(i0+1,i1) -
     &        eta_n / (dx(0)* dx(0)) * thn_data(i0+1,i1) * 
     &        un_data_0(i0+2,i1) -
     &        eta_n / (dx(1)* dx(1)) * thn_iphalf_jphalf * 
     &        un_data_0(i0+1,i1+1) -
     &        eta_n / (dx(1)* dx(1)) * thn_iphalf_jmhalf * 
     &        un_data_0(i0+1,i1-1) -
     &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jphalf * 
     &        un_data_1(i0+1,i1+1) +
     &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jmhalf * 
     &        un_data_1(i0+1,i1) +
     &        eta_n / (dx(0)* dx(1)) * thn_data(i0+1,i1) *
     &        (un_data_1(i0+1,i1+1) - un_data_1(i0+1,i1)))
c
          ! solvent at east edge
          b(6) = f_us_data_0(i0+1,i1) + D*(toThs(thn_upper_x) / dx(0) *
     &    p_data(i0+1,i1) -
     &    eta_s / (dx(0)* dx(0)) * toThs(thn_data(i0+1,i1)) * 
     &    us_data_0(i0+2,i1) -
     &    eta_s / (dx(1)* dx(1)) * toThs(thn_iphalf_jphalf) * 
     &    us_data_0(i0+1,i1+1) -
     &    eta_s / (dx(1)* dx(1)) * toThs(thn_iphalf_jmhalf) * 
     &    us_data_0(i0+1,i1-1) -
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jphalf) * 
     &    us_data_1(i0+1,i1+1) +
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jmhalf) * 
     &    us_data_1(i0+1,i1) +
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0+1,i1)) *
     &        (us_data_1(i0+1,i1+1) - us_data_1(i0+1,i1)))
c
          ! network at south edge
          b(3) = f_un_data_1(i0,i1)+D*(-thn_lower_y / dx(1) *
     &        p_data(i0,i1-1) -
     &        eta_n / (dx(1)* dx(1)) * thn_data(i0,i1-1) * 
     &        un_data_1(i0,i1-1) -
     &        eta_n / (dx(0)* dx(0)) * thn_iphalf_jmhalf * 
     &        un_data_1(i0+1,i1) -
     &        eta_n / (dx(0)* dx(0)) * thn_imhalf_jmhalf * 
     &        un_data_1(i0-1,i1) +
     &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jmhalf * 
     &        un_data_0(i0+1,i1-1) -
     &        eta_n / (dx(0)* dx(1)) * thn_imhalf_jmhalf * 
     &        un_data_0(i0,i1-1) -
     &        eta_n / (dx(0)* dx(1)) * thn_data(i0,i1-1) *
     &        (un_data_0(i0+1,i1-1) - un_data_0(i0,i1-1)))
c
          ! solvent at south edge
          b(7) = f_us_data_1(i0,i1)+D*(-toThs(thn_lower_y) / dx(1) *
     &    p_data(i0,i1-1) -
     &    eta_s / (dx(1)* dx(1)) * toThs(thn_data(i0,i1-1)) * 
     &    us_data_1(i0,i1-1) -
     &    eta_s / (dx(0)* dx(0)) * toThs(thn_iphalf_jmhalf) * 
     &    us_data_1(i0+1,i1) -
     &    eta_s / (dx(0)* dx(0)) * toThs(thn_imhalf_jmhalf) * 
     &    us_data_1(i0-1,i1) +
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jmhalf) * 
     &    us_data_0(i0+1,i1-1) -
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_imhalf_jmhalf) * 
     &    us_data_0(i0,i1-1) -
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0,i1-1)) *
     &    (us_data_0(i0+1,i1-1) - us_data_0(i0,i1-1)))
c
          ! network at north edge
          b(4) = f_un_data_1(i0,i1+1) + D*(thn_upper_y / dx(1) *
     &        p_data(i0,i1+1) -
     &        eta_n / (dx(1)* dx(1)) * thn_data(i0,i1+1) * 
     &        un_data_1(i0,i1+2) -
     &        eta_n / (dx(0)* dx(0)) * thn_iphalf_jphalf * 
     &        un_data_1(i0+1,i1+1) -
     &        eta_n / (dx(0)* dx(0)) * thn_imhalf_jphalf * 
     &        un_data_1(i0-1,i1+1) -
     &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jphalf * 
     &        un_data_0(i0+1,i1+1) +
     &        eta_n / (dx(0)* dx(1)) * thn_imhalf_jphalf * 
     &        un_data_0(i0,i1 + 1) +
     &        eta_n / (dx(0)* dx(1)) * thn_data(i0,i1+1) *
     &        (un_data_0(i0+1,i1+1) - un_data_0(i0,i1 + 1)))
c
          ! solvent at north edge
          b(8) =  f_us_data_1(i0,i1+1) + D*(toThs(thn_upper_y) / dx(1)*
     &    p_data(i0,i1+1) -
     &    eta_s / (dx(1)* dx(1)) * toThs(thn_data(i0,i1+1)) * 
     &    us_data_1(i0,i1+2) -
     &    eta_s / (dx(0)* dx(0)) * toThs(thn_iphalf_jphalf) * 
     &    us_data_1(i0+1,i1+1) -
     &    eta_s / (dx(0)* dx(0)) * toThs(thn_imhalf_jphalf) * 
     &    us_data_1(i0-1,i1+1) -
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jphalf) * 
     &    us_data_0(i0+1,i1+1) +
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_imhalf_jphalf) * 
     &    us_data_0(i0,i1 + 1) +
     &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0,i1+1)) *
     &    (us_data_0(i0+1,i1+1) - us_data_0(i0,i1 + 1)))
c
          ! pressure at cell center
          b(9) = f_p_data(i0,i1)
c
          ! CALL dgetrf( 9, 9, A_box, 9, ipiv, info)

          !IF( info.EQ.0 ) THEN
          !     print *, 'info is:'
               ! Solve the system A*X = B, overwriting B with X.
          !     CALL dgetrs( 'No transpose', 9, 1, A_box, 9, ipiv, b, 
c    &               9, info )
          !END IF

          ! solve the system Ax = b, overwriting b with x
          call dgesv(9, 1, A_box, 9, ipiv, b, 9, info)
c
          if (info /= 0) then
            ! print *, A_box
            print *, 'ERROR IN DGESV'
            print *, info
          endif
c 
          un_data_0(i0,i1) = (1.d0-w)*un_data_0(i0,i1) + w*b(1);
          un_data_0(i0+1,i1) = (1.d0-w)*un_data_0(i0+1,i1) + w*b(2);
          un_data_1(i0,i1) = (1.d0-w)*un_data_1(i0,i1) + w*b(3);
          un_data_1(i0,i1+1) = (1.d0-w)*un_data_1(i0,i1+1) + w*b(4);
          us_data_0(i0,i1) = (1.d0-w)*us_data_0(i0,i1) + w*b(5);
          us_data_0(i0+1,i1) = (1.d0-w)*us_data_0(i0+1,i1) + w*b(6);
          us_data_1(i0,i1) = (1.d0-w)*us_data_1(i0,i1) + w*b(7);
          us_data_1(i0,i1+1) = (1.d0-w)*us_data_1(i0,i1+1) + w*b(8);
          p_data(i0,i1) = (1.d0-w)*p_data(i0,i1) + w*b(9);
c
          end if
        enddo
      enddo
c
      end subroutine


