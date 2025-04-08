c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Function to convert Thn to Ths since Thn + Ths = 1 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      module converting_to_ths

      implicit none

      contains
      double precision function toThs(Thn)
c      
        double precision Thn  
        toThs = 1.d0 - Thn
c
        return
      end function toThs

      end module converting_to_ths

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Accumulate the forces into respective patch indices for the 
c       network and solvent on a given patch. Assumes ghost cells
c       have been filled for the velocities and volume fraction.
c
c       Specifically, computes
c
c       [ C*thn + A_n + D*xi   -D*xi                ][un]
c       [ -D*xi                 C*ths + A_s + D*xi  ][us]
c       in which
c       A_i = D_*eta_i*div(thn*((grad+grad^T)-div*I))
c      
c       Accumulate the momentum forces for constant coefficient problems 
c       with the network volume fraction only provided at cell centers.
c       In this case, the volume fraction is linearly interpolated to
c       respective sides and nodes when necessary.
   
c       Note that no synchronization is provided on the volume 
c       fraction when linear interpolation is done. 
c
c       un, us, f_un, f_us are vector-valued side-centered velocities
c       and thn is a scalar-valued cell-centered network volume fraction.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
        subroutine m_w_p_p_c_c(
     &        dx, ilow0, iup0,
     &        ilow1, iup1, 
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_gcw, eta_n, eta_s,
     &        l_n, l_s, nu_n, nu_s, xi, C, D)
c
      use converting_to_ths
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0  
      integer ilow1,  iup1
      integer un_gcw, us_gcw, f_un_gcw, f_us_gcw
      integer thn_gcw
      double precision eta_n, eta_s, nu_n, nu_s, xi, C, D
      double precision l_n, l_s
c      
      double precision thn_data(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
c
      double precision un_data_0(ilow0-un_gcw:iup0+un_gcw+1,  ! Adding 1 because  no of sides =  no of cells + 1
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_data_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)                   ! Adding 1 because  no of sides =  no of cells + 1
c
      double precision us_data_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_data_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
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
      double precision thn_lower_x, thn_lower_y
      double precision thn_imh_jph, thn_imh_jmh
      double precision thn_iph_jph, thn_iph_jmh
c    
      double precision dx_dx, dy_dy, dx_dy
c 
      double precision ddx_Thn_dx_un, ddy_Thn_dy_un
      double precision ddy_Thn_dx_vn, ddx_Thn_dy_vn
      double precision ddx_Ths_dx_us, ddy_Ths_dy_us
      double precision ddy_Ths_dx_vs, ddx_Ths_dy_vs
      double precision drag_n, drag_s
c   
      integer i0, i1
c
      dx_dx = dx(0) * dx(0)
      dy_dy = dx(1) * dx(1)
      dx_dy = dx(0) * dx(1)

c     Loop over side-centers in x-dir
      do i1 = ilow1, iup1 
        do i0 = ilow0, iup0 + 1

          ! calculate thn at sides
          thn_lower_x = 0.5d0*(thn_data(i0,i1)+thn_data(i0-1,i1))    ! thn(i-1/2, j)

          ! calculate thn at corners
          thn_imh_jph = 0.25d0*(thn_data(i0-1, i1) + thn_data(i0,i1)
     &                   +thn_data(i0,i1+1) + thn_data(i0-1,i1+1))   ! thn(i-1/2, j+1/2)
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0-1,i1)
     &                   +thn_data(i0,i1-1) + thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)

          ! components of first row (x-component of velocity) of network equation
          ddx_Thn_dx_un = ((2.d0*eta_n-l_n) / dx_dx) *
     &        (thn_data(i0,i1) *
     &              (un_data_0(i0+1,i1)-un_data_0(i0,i1)) 
     &       -thn_data(i0-1,i1) *
     &              (un_data_0(i0,i1)-un_data_0(i0-1,i1)))
          ddy_Thn_dy_un = (eta_n / dy_dy) *
     &        (thn_imh_jph * 
     &              (un_data_0(i0,i1+1) - un_data_0(i0,i1))
     &       -thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1)))
          ddy_Thn_dx_vn = (eta_n / dx_dy) *
     &        (thn_imh_jph *
     &              (un_data_1(i0,i1+1)-un_data_1(i0-1,i1+1))
     &       -thn_imh_jmh *
     &              (un_data_1(i0,i1)-un_data_1(i0-1,i1)))
          ddx_Thn_dy_vn = -(l_n / dx_dy) *
     &        (thn_data(i0,i1) *
     &              (un_data_1(i0,i1+1)-un_data_1(i0,i1))
     &        -thn_data(i0-1,i1) *
     &              (un_data_1(i0-1,i1+1)-un_data_1(i0-1,i1)))
          drag_n = -(xi / nu_n) * thn_lower_x * toThs(thn_lower_x) *
     &              (un_data_0(i0,i1) - us_data_0(i0,i1))

          f_un_data_0(i0,i1) = D * (ddx_Thn_dx_un + ddy_Thn_dy_un 
     &              + ddy_Thn_dx_vn + ddx_Thn_dy_vn + drag_n) 
     &              + C * thn_lower_x * un_data_0(i0,i1)

          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dx_dx) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_0(i0+1,i1)-us_data_0(i0,i1)) 
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_0(i0,i1)-us_data_0(i0-1,i1)))
          ddy_Ths_dy_us = (eta_s / dy_dy) *
     &        (toThs(thn_imh_jph) * 
     &              (us_data_0(i0,i1+1) - us_data_0(i0,i1))
     &       - toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          ddy_Ths_dx_vs = (eta_s / dx_dy) *
     &        (toThs(thn_imh_jph) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0-1,i1+1))
     &       - toThs(thn_imh_jmh) *
     &              (us_data_1(i0,i1)-us_data_1(i0-1,i1)))
          ddx_Ths_dy_vs = -(l_s / dx_dy) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0,i1))
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_1(i0-1,i1+1)-us_data_1(i0-1,i1)))
          drag_s = -(xi / nu_s) * thn_lower_x * toThs(thn_lower_x) *
     &              (us_data_0(i0,i1) - un_data_0(i0,i1))

          f_us_data_0(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us
     &              + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &              + C * toThs(thn_lower_x) * us_data_0(i0,i1)

        end do
      end do
c
c     Loop over side centers in y-dir
      do i1 = ilow1, iup1 + 1
        do i0 = ilow0, iup0
c
          ! calculate thn at (i,j-1/2)
          thn_lower_y = 0.5d0*(thn_data(i0,i1)+thn_data(i0,i1-1))  ! thn(i,j-1/2)

          ! calculate thn at corners
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0,i1-1)
     &                   +thn_data(i0-1,i1) + thn_data(i0-1,i1-1))    ! thn(i-1/2, j-1/2)
          thn_iph_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0,i1-1)
     &                   +thn_data(i0+1,i1) + thn_data(i0+1,i1-1))    ! thn(i+1/2, j-1/2)

          ! components of second row (y-component of network vel) of network equation
          ddy_Thn_dy_un = ((2.d0*eta_n-l_n) / dy_dy) *
     &         (thn_data(i0,i1) * 
     &              (un_data_1(i0,i1+1) - un_data_1(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_1(i0,i1) - un_data_1(i0,i1-1)))

          ddx_Thn_dx_un = (eta_n / dx_dx) *
     &         (thn_iph_jmh * 
     &              (un_data_1(i0+1,i1) - un_data_1(i0,i1)) -
     &          thn_imh_jmh * 
     &              (un_data_1(i0,i1) - un_data_1(i0-1,i1)))

          ddx_Thn_dy_vn = (eta_n / dx_dy) *
     &         (thn_iph_jmh * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0+1,i1-1)) -
     &         thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1))) 
          
          ddy_Thn_dx_vn = -(l_n / dx_dy) *
     &         (thn_data(i0,i1) * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_0(i0+1,i1-1) - un_data_0(i0,i1-1))) 
          
          drag_n = -(xi / nu_n) * thn_lower_y * 
     &        toThs(thn_lower_y) * (un_data_1(i0,i1)-us_data_1(i0,i1))
     
          f_un_data_1(i0,i1) = D * (ddy_Thn_dy_un + ddx_Thn_dx_un 
     &         + ddx_Thn_dy_vn + ddy_Thn_dx_vn + drag_n) 
     &         + C * thn_lower_y * un_data_1(i0,i1) 
          
          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dy_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_1(i0,i1+1) - us_data_1(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_1(i0,i1) - us_data_1(i0,i1-1)))

          ddy_Ths_dy_us = (eta_s / dx_dx) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_1(i0+1,i1) - us_data_1(i0,i1)) -
     &          toThs(thn_imh_jmh) * 
     &              (us_data_1(i0,i1) - us_data_1(i0-1,i1)))

          ddy_Ths_dx_vs = (eta_s / dx_dy) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0+1,i1-1)) -
     &         toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          
          ddx_Ths_dy_vs = -(l_s / dx_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_0(i0+1,i1-1) - us_data_0(i0,i1-1)))
          
          drag_s = -(xi / nu_s) * thn_lower_y * 
     &        toThs(thn_lower_y) * (us_data_1(i0,i1)-un_data_1(i0,i1))
     
          f_us_data_1(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us 
     &         + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &         + C * toThs(thn_lower_y) * us_data_1(i0,i1)

        end do
      end do
c
      end subroutine



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Accumulate the forces into respective patch indices for the 
c       network and solvent on a given patch. Assumes ghost cells
c       have been filled for the velocities and volume fraction.
c
c       Specifically, computes
c
c       [ C*thn + A_n + D*xi   -D*xi                ][un]
c       [ -D*xi                 C*ths + A_s + D*xi  ][us]
c       in which
c       A_i = D*eta_i*div(thn*((grad+grad^T)-div*I))
c      
c       Accumulate the momentum forces for constant coefficient problems with
c       the network volume fraction interpolated to cell nodes & cell sides.
c
c       un, us, f_un, f_us are vector-valued side-centered velocities
c       and thn is a scalar-valued cell-centered network volume fraction.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
        subroutine m_w_p_p_c_c_thn_nodes(
     &        dx, ilow0, iup0,
     &        ilow1, iup1, 
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_nc_data, thn_sc_data_0, thn_sc_data_1, 
     &        thn_gcw, thn_nc_gcw, thn_sc_gcw, 
     &        eta_n, eta_s, l_n, l_s, nu_n, nu_s, xi, C, D)
c
      use converting_to_ths
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0  
      integer ilow1,  iup1
      integer un_gcw, us_gcw, f_un_gcw, f_us_gcw
      integer thn_gcw, thn_nc_gcw, thn_sc_gcw
      double precision eta_n, eta_s, nu_n, nu_s, xi, C, D
      double precision l_n, l_s
c      
      double precision thn_data(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
      double precision thn_nc_data(
     &          ilow0-thn_nc_gcw:iup0+thn_nc_gcw+1,
     &          ilow1-thn_nc_gcw:iup1+thn_nc_gcw+1)           ! cell(i0,i1) bottom left node
      double precision thn_sc_data_0(
     &          ilow0-thn_sc_gcw:iup0+thn_sc_gcw+1,
     &          ilow1-thn_sc_gcw:iup1+thn_sc_gcw) 
      double precision thn_sc_data_1(
     &         ilow0-thn_sc_gcw:iup0+thn_sc_gcw,
     &         ilow1-thn_sc_gcw:iup1+thn_sc_gcw+1) 
c
      double precision un_data_0(ilow0-un_gcw:iup0+un_gcw+1,  ! Adding 1 because  no of sides =  no of cells + 1
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_data_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)                   ! Adding 1 because  no of sides =  no of cells + 1
c
      double precision us_data_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_data_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
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
      double precision thn_lower_x, thn_lower_y
      double precision thn_imh_jph, thn_imh_jmh
      double precision thn_iph_jph, thn_iph_jmh
c    
      double precision dx_dx, dy_dy, dx_dy
c 
      double precision ddx_Thn_dx_un, ddy_Thn_dy_un
      double precision ddy_Thn_dx_vn, ddx_Thn_dy_vn
      double precision ddx_Ths_dx_us, ddy_Ths_dy_us
      double precision ddy_Ths_dx_vs, ddx_Ths_dy_vs
      double precision drag_n, drag_s

c   
      integer i0, i1
c
      dx_dx = dx(0) * dx(0)
      dy_dy = dx(1) * dx(1)
      dx_dy = dx(0) * dx(1)

c     Loop over side-centers in x-dir
      do i1 = ilow1, iup1 
        do i0 = ilow0, iup0 + 1

          ! calculate thn at sides
          thn_lower_x = thn_sc_data_0(i0,i1)    ! thn(i-1/2, j)

          ! calculate thn at corners
          thn_imh_jph = thn_nc_data(i0,i1+1)  ! Upper left
          thn_imh_jmh = thn_nc_data(i0,i1)    ! Lower Left

          ! components of first row (x-component of velocity) of network equation
          ddx_Thn_dx_un = ((2.d0*eta_n - l_n) / dx_dx) *
     &        (thn_data(i0,i1) *
     &              (un_data_0(i0+1,i1)-un_data_0(i0,i1)) 
     &       -thn_data(i0-1,i1) *
     &              (un_data_0(i0,i1)-un_data_0(i0-1,i1)))
          ddy_Thn_dy_un = (eta_n / dy_dy) *
     &        (thn_imh_jph * 
     &              (un_data_0(i0,i1+1) - un_data_0(i0,i1))
     &       -thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1)))
          ddy_Thn_dx_vn = (eta_n / dx_dy) *
     &        (thn_imh_jph *
     &              (un_data_1(i0,i1+1)-un_data_1(i0-1,i1+1))
     &       -thn_imh_jmh *
     &              (un_data_1(i0,i1)-un_data_1(i0-1,i1)))
          ddx_Thn_dy_vn = -(l_n / dx_dy) *
     &        (thn_data(i0,i1) *
     &              (un_data_1(i0,i1+1)-un_data_1(i0,i1))
     &        -thn_data(i0-1,i1) *
     &              (un_data_1(i0-1,i1+1)-un_data_1(i0-1,i1)))
          drag_n = -(xi / nu_n) * thn_lower_x * toThs(thn_lower_x) *
     &              (un_data_0(i0,i1) - us_data_0(i0,i1))

          f_un_data_0(i0,i1) = D * (ddx_Thn_dx_un + ddy_Thn_dy_un 
     &              + ddy_Thn_dx_vn + ddx_Thn_dy_vn + drag_n) 
     &              + C * thn_lower_x * un_data_0(i0,i1)

          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dx_dx) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_0(i0+1,i1)-us_data_0(i0,i1)) 
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_0(i0,i1)-us_data_0(i0-1,i1)))
          ddy_Ths_dy_us = (eta_s / dy_dy) *
     &        (toThs(thn_imh_jph) * 
     &              (us_data_0(i0,i1+1) - us_data_0(i0,i1))
     &       - toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          ddy_Ths_dx_vs = (eta_s / dx_dy) *
     &        (toThs(thn_imh_jph) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0-1,i1+1))
     &       - toThs(thn_imh_jmh) *
     &              (us_data_1(i0,i1)-us_data_1(i0-1,i1)))
          ddx_Ths_dy_vs = -(l_s / dx_dy) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0,i1))
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_1(i0-1,i1+1)-us_data_1(i0-1,i1)))
          drag_s = -(xi / nu_s) * thn_lower_x * toThs(thn_lower_x) *
     &              (us_data_0(i0,i1) - un_data_0(i0,i1))

          f_us_data_0(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us
     &              + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &              + C * toThs(thn_lower_x) * us_data_0(i0,i1)

        end do
      end do
c
c     Loop over side centers in y-dir
      do i1 = ilow1, iup1 + 1   
        do i0 = ilow0, iup0
c
          ! calculate thn at (i,j-1/2)
          thn_lower_y = thn_sc_data_1(i0,i1)  ! thn(i,j-1/2)

          ! calculate thn at corners
          thn_imh_jmh = thn_nc_data(i0,i1)      ! Lower Left
          thn_iph_jmh = thn_nc_data(i0+1,i1)    ! Lower Right

          ! components of second row (y-component of network vel) of network equation
          ddy_Thn_dy_un = ((2.d0*eta_n-l_n) / dy_dy) *
     &         (thn_data(i0,i1) * 
     &              (un_data_1(i0,i1+1) - un_data_1(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_1(i0,i1) - un_data_1(i0,i1-1)))

          ddx_Thn_dx_un = (eta_n / dx_dx) *
     &         (thn_iph_jmh * 
     &              (un_data_1(i0+1,i1) - un_data_1(i0,i1)) -
     &          thn_imh_jmh * 
     &              (un_data_1(i0,i1) - un_data_1(i0-1,i1)))

          ddx_Thn_dy_vn = (eta_n / dx_dy) *
     &         (thn_iph_jmh * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0+1,i1-1)) -
     &         thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1)))
          
          ddy_Thn_dx_vn = -(l_n / dx_dy) *
     &         (thn_data(i0,i1) * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_0(i0+1,i1-1) - un_data_0(i0,i1-1)))
          
          drag_n = -(xi / nu_n) * thn_lower_y * 
     &        toThs(thn_lower_y) * (un_data_1(i0,i1)-us_data_1(i0,i1))
     
          f_un_data_1(i0,i1) = D * (ddy_Thn_dy_un + ddx_Thn_dx_un 
     &         + ddx_Thn_dy_vn + ddy_Thn_dx_vn + drag_n) 
     &         + C * thn_lower_y * un_data_1(i0,i1)
          
          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dy_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_1(i0,i1+1) - us_data_1(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_1(i0,i1) - us_data_1(i0,i1-1)))

          ddy_Ths_dy_us = (eta_s / dx_dx) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_1(i0+1,i1) - us_data_1(i0,i1)) -
     &          toThs(thn_imh_jmh) * 
     &              (us_data_1(i0,i1) - us_data_1(i0-1,i1)))

          ddy_Ths_dx_vs = (eta_s / dx_dy) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0+1,i1-1)) -
     &         toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          
          ddx_Ths_dy_vs = -(l_s / dx_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_0(i0+1,i1-1) - us_data_0(i0,i1-1)))
          
          drag_s = -(xi / nu_s) * thn_lower_y * 
     &        toThs(thn_lower_y) * (us_data_1(i0,i1)-un_data_1(i0,i1))
     
          f_us_data_1(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us 
     &         + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &         + C * toThs(thn_lower_y) * us_data_1(i0,i1)

        end do
      end do
c
      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Accumulate the forces into respective patch indices for the 
c       network and solvent on a given patch with variable drag. 
c       Assumes ghost cells have been filled for the velocities
c       and volume fraction.
c
c       Specifically, computes
c
c       [ C*thn + A_n + D*xi   -D*xi                ][un]
c       [ -D*xi                 C*ths + A_s + D*xi  ][us]
c       in which
c       A_i = D*eta_i*div(thn*((grad+grad^T)-div*I))
c      
c       Accumulate the momentum forces for constant coefficient problems 
c       with the network volume fraction only provided at cell centers.
c       In this case, the volume fraction is linearly interpolated to
c       respective sides and nodes when necessary.
   
c       Note that no synchronization is provided on the volume 
c       fraction when linear interpolation is done. 
c
c       un, us, f_un, f_us are vector-valued side-centered velocities
c       and thn is a scalar-valued cell-centered network volume fraction.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
        subroutine m_w_p_p_var_drag(
     &        dx, ilow0, iup0,
     &        ilow1, iup1, 
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_gcw, eta_n, eta_s,
     &        l_n, l_s,
     &        xi_data_0, xi_data_1, xi_gcw, C, D)
c
      use converting_to_ths
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0  
      integer ilow1,  iup1
      integer un_gcw, us_gcw, f_un_gcw, f_us_gcw
      integer thn_gcw, xi_gcw
      double precision eta_n, eta_s, C, D
      double precision l_n, l_s
c      
      double precision xi_data_0(ilow0-xi_gcw:iup0+xi_gcw+1,  
     &          ilow1-xi_gcw:iup1+xi_gcw) 
      double precision xi_data_1(ilow0-xi_gcw:iup0+xi_gcw,  
     &          ilow1-xi_gcw:iup1+xi_gcw+1) 
c    
      double precision thn_data(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
c
      double precision un_data_0(ilow0-un_gcw:iup0+un_gcw+1,  ! Adding 1 because no of sides = no of cells + 1
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_data_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)                   ! Adding 1 because no of sides = no of cells + 1
c
      double precision us_data_0(ilow0-us_gcw:iup0+us_gcw+1,
     &          ilow1-us_gcw:iup1+us_gcw)
      double precision us_data_1(ilow0-us_gcw:iup0+us_gcw,
     &          ilow1-us_gcw:iup1+us_gcw+1)
c
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
      double precision thn_lower_x, thn_lower_y
      double precision thn_imh_jph, thn_imh_jmh
      double precision thn_iph_jph, thn_iph_jmh
c    
      double precision dx_dx, dy_dy, dx_dy
c 
      double precision ddx_Thn_dx_un, ddy_Thn_dy_un
      double precision ddy_Thn_dx_vn, ddx_Thn_dy_vn
      double precision ddx_Ths_dx_us, ddy_Ths_dy_us
      double precision ddy_Ths_dx_vs, ddx_Ths_dy_vs
      double precision drag_n, drag_s
c   
      integer i0, i1
c
      dx_dx = dx(0) * dx(0)
      dy_dy = dx(1) * dx(1)
      dx_dy = dx(0) * dx(1)

c     Loop over side-centers in x-dir
      do i1 = ilow1, iup1 
        do i0 = ilow0, iup0 + 1

          ! calculate thn at sides
          thn_lower_x = 0.5d0*(thn_data(i0,i1)+thn_data(i0-1,i1))  ! thn(i-1/2, j)

          ! calculate thn at corners
          thn_imh_jph = 0.25d0*(thn_data(i0-1, i1) + thn_data(i0,i1)
     &                   +thn_data(i0,i1+1) + thn_data(i0-1,i1+1))   ! thn(i-1/2, j+1/2)
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0-1,i1)
     &                   +thn_data(i0,i1-1) + thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)

          ! components of first row (x-component of velocity) of network equation
          ddx_Thn_dx_un = ((2.d0*eta_n-l_n) / dx_dx) *
     &        (thn_data(i0,i1) *
     &              (un_data_0(i0+1,i1)-un_data_0(i0,i1)) 
     &       -thn_data(i0-1,i1) *
     &              (un_data_0(i0,i1)-un_data_0(i0-1,i1)))
          ddy_Thn_dy_un = (eta_n / dy_dy) *
     &        (thn_imh_jph * 
     &              (un_data_0(i0,i1+1) - un_data_0(i0,i1))
     &       -thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1)))
          ddy_Thn_dx_vn = (eta_n / dx_dy) *
     &        (thn_imh_jph *
     &              (un_data_1(i0,i1+1)-un_data_1(i0-1,i1+1))
     &       -thn_imh_jmh *
     &              (un_data_1(i0,i1)-un_data_1(i0-1,i1)))
          ddx_Thn_dy_vn = -(l_n / dx_dy) *
     &        (thn_data(i0,i1) *
     &              (un_data_1(i0,i1+1)-un_data_1(i0,i1))
     &        -thn_data(i0-1,i1) *
     &              (un_data_1(i0-1,i1+1)-un_data_1(i0-1,i1)))
     
          drag_n = -xi_data_0(i0,i1) *
     &         (un_data_0(i0,i1) - us_data_0(i0,i1))        

          f_un_data_0(i0,i1) = D * (ddx_Thn_dx_un + ddy_Thn_dy_un 
     &              + ddy_Thn_dx_vn + ddx_Thn_dy_vn + drag_n) 
     &              + C * thn_lower_x * un_data_0(i0,i1)

          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dx_dx) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_0(i0+1,i1)-us_data_0(i0,i1)) 
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_0(i0,i1)-us_data_0(i0-1,i1)))
          ddy_Ths_dy_us = (eta_s / dy_dy) *
     &        (toThs(thn_imh_jph) * 
     &              (us_data_0(i0,i1+1) - us_data_0(i0,i1))
     &       - toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          ddy_Ths_dx_vs = (eta_s / dx_dy) *
     &        (toThs(thn_imh_jph) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0-1,i1+1))
     &       - toThs(thn_imh_jmh) *
     &              (us_data_1(i0,i1)-us_data_1(i0-1,i1)))
          ddx_Ths_dy_vs = -(l_s / dx_dy) *
     &        (toThs(thn_data(i0,i1)) *
     &              (us_data_1(i0,i1+1)-us_data_1(i0,i1))
     &       - toThs(thn_data(i0-1,i1)) *
     &              (us_data_1(i0-1,i1+1)-us_data_1(i0-1,i1)))
          drag_s = -xi_data_0(i0,i1) *
     &              (us_data_0(i0,i1) - un_data_0(i0,i1))

          f_us_data_0(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us
     &              + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &              + C * toThs(thn_lower_x) * us_data_0(i0,i1)

        end do
      end do
c
c     Loop over side centers in y-dir
      do i1 = ilow1, iup1 + 1  
        do i0 = ilow0, iup0
c
          ! calculate thn at (i,j-1/2)
          thn_lower_y = 0.5d0*(thn_data(i0,i1)+thn_data(i0,i1-1))  ! thn(i,j-1/2)

          ! calculate thn at corners
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0,i1-1)
     &                   +thn_data(i0-1,i1) + thn_data(i0-1,i1-1))    ! thn(i-1/2, j-1/2)
          thn_iph_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0,i1-1)
     &                   +thn_data(i0+1,i1) + thn_data(i0+1,i1-1))    ! thn(i+1/2, j-1/2)

          ! components of second row (y-component of network vel) of network equation
          ddy_Thn_dy_un = (2.d0*eta_n-l_n) / dx_dx *
     &         (thn_data(i0,i1) * 
     &              (un_data_1(i0,i1+1) - un_data_1(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_1(i0,i1) - un_data_1(i0,i1-1)))

          ddx_Thn_dx_un = (eta_n / dy_dy) *
     &         (thn_iph_jmh * 
     &              (un_data_1(i0+1,i1) - un_data_1(i0,i1)) -
     &          thn_imh_jmh * 
     &              (un_data_1(i0,i1) - un_data_1(i0-1,i1)))

          ddx_Thn_dy_vn = (eta_n / dx_dy) *
     &         (thn_iph_jmh * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0+1,i1-1)) -
     &         thn_imh_jmh * 
     &              (un_data_0(i0,i1) - un_data_0(i0,i1-1)))
          
          ddy_Thn_dx_vn = -(l_n / dx_dy) *
     &         (thn_data(i0,i1) * 
     &              (un_data_0(i0+1,i1) - un_data_0(i0,i1)) -
     &          thn_data(i0,i1-1) * 
     &              (un_data_0(i0+1,i1-1) - un_data_0(i0,i1-1)))
          
          drag_n = -xi_data_1(i0,i1) *
     &         (un_data_1(i0,i1)-us_data_1(i0,i1))
     
          f_un_data_1(i0,i1) = D * (ddy_Thn_dy_un + ddx_Thn_dx_un 
     &         + ddx_Thn_dy_vn + ddy_Thn_dx_vn + drag_n) 
     &         + C * thn_lower_y * un_data_1(i0,i1) 
          
          ! Solvent equation
          ddx_Ths_dx_us = ((2.d0*eta_s-l_s) / dy_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_1(i0,i1+1) - us_data_1(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_1(i0,i1) - us_data_1(i0,i1-1)))

          ddy_Ths_dy_us = (eta_s / dx_dx) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_1(i0+1,i1) - us_data_1(i0,i1)) -
     &          toThs(thn_imh_jmh) * 
     &              (us_data_1(i0,i1) - us_data_1(i0-1,i1)))

          ddx_Ths_dy_vs = (eta_s / dx_dy) *
     &         (toThs(thn_iph_jmh) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0+1,i1-1)) -
     &         toThs(thn_imh_jmh) * 
     &              (us_data_0(i0,i1) - us_data_0(i0,i1-1)))
          
          ddy_Ths_dx_vs = -(l_s / dx_dy) *
     &         (toThs(thn_data(i0,i1)) * 
     &              (us_data_0(i0+1,i1) - us_data_0(i0,i1)) -
     &          toThs(thn_data(i0,i1-1)) * 
     &              (us_data_0(i0+1,i1-1) - us_data_0(i0,i1-1)))
          
          drag_s = -xi_data_1(i0,i1) *
     &         (us_data_1(i0,i1)-un_data_1(i0,i1))
     
          f_us_data_1(i0,i1) = D * (ddx_Ths_dx_us + ddy_Ths_dy_us 
     &         + ddy_Ths_dx_vs + ddx_Ths_dy_vs + drag_s) 
     &         + C * toThs(thn_lower_y) * us_data_1(i0,i1)

        end do
      end do
c
      end subroutine