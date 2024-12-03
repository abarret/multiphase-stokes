c234567
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Function to convert Thn to Ths since Thn + Ths = 1 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      module convert_to_ths

      implicit none

      contains
      double precision function toThs(Thn)
c      
        double precision Thn  
        toThs = 1.d0 - Thn
c
        return
      end function toThs

      end module convert_to_ths

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This smoother performs red-black Gauss Seidel iterations with successive
c       under-relaxation for the Stokes Block. 
c
c       un, us, f_un, f_us are vector-valued side-centered velocities
c       and thn is a scalar-valued cell-centered network volume fraction.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine velocity_rbgs(
     &        dx, ilow0, iup0,
     &        ilow1, iup1, 
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_gcw, eta_n, eta_s,
     &        nu_n, nu_s, xi, w, C, D, red_or_black)
c
      use convert_to_ths
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilow0,  iup0  
      integer ilow1,  iup1
      integer un_gcw, us_gcw, f_un_gcw, f_us_gcw
      integer thn_gcw, red_or_black
      double precision eta_n, eta_s, nu_n, nu_s, xi, w, C, D
c      
      double precision thn_data(ilow0-thn_gcw:iup0+thn_gcw,
     &          ilow1-thn_gcw:iup1+thn_gcw) 
c
      double precision un_data_0(ilow0-un_gcw:iup0+un_gcw+1,  ! Adding 1 because # of sides = # of cells + 1
     &          ilow1-un_gcw:iup1+un_gcw)
      double precision un_data_1(ilow0-un_gcw:iup0+un_gcw,
     &          ilow1-un_gcw:iup1+un_gcw+1)                   ! Adding 1 because # of sides = # of cells + 1
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
      double precision un_imhalf_j, us_imhalf_j
      double precision un_i_jmhalf, us_i_jmhalf
c
      integer i0, i1
c
      dx_dx = dx(0) * dx(0)
      dy_dy = dx(1) * dx(1)
      dx_dy = dx(0) * dx(1)
c
c     Loop over side-centers in x-dir
      do i1 = ilow1, iup1 
        do i0 = ilow0, iup0 + 1
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then
c
          ! calculate thn at sides
          thn_lower_x = 0.5d0*(thn_data(i0,i1)+thn_data(i0-1,i1))  ! thn(i-1/2, j)

          ! calculate thn at corners
          thn_imh_jph = 0.25d0*(thn_data(i0-1, i1) + thn_data(i0,i1)
     &                   +thn_data(i0,i1+1) + thn_data(i0-1,i1+1))   ! thn(i-1/2, j+1/2)
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0-1,i1)
     &                   +thn_data(i0,i1-1) + thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)

          ! solve for network velocities
          un_imhalf_j = f_un_data_0(i0,i1) - (D * eta_n * 
     &      ((thn_imh_jph-thn_data(i0,i1))/(dx_dy) * un_data_1(i0,i1+1) 
     &      + (thn_data(i0-1,i1))/(dx_dx) * un_data_0(i0-1,i1) 
     &      + (thn_data(i0,i1)/dx_dx * un_data_0(i0+1,i1)) 
     &      + (thn_imh_jph/dy_dy * un_data_0(i0,i1+1)) 
     &      + (thn_imh_jmh/dy_dy * un_data_0(i0,i1-1))
     &      + (thn_data(i0-1,i1) - thn_imh_jph)/(dx_dy) 
     &             * un_data_1(i0-1,i1+1)
     &      + (thn_imh_jmh - thn_data(i0-1,i1))/(dx_dy) 
     &             * un_data_1(i0-1,i1)
     &      + (thn_data(i0,i1)-thn_imh_jmh)/(dx_dy)
     &             * un_data_1(i0,i1)) 
     &      + D * xi * us_data_0(i0,i1)
     &             * thn_lower_x * toThs(thn_lower_x))

          un_imhalf_j = un_imhalf_j/ (D * eta_n * 
     &      (-(thn_data(i0,i1) + thn_data(i0-1,i1)) / (dx_dx) 
     &      - (thn_imh_jph + thn_imh_jmh) / (dy_dy)) 
     &      - D * xi * thn_lower_x * toThs(thn_lower_x)
     &      + C * thn_lower_x)

          un_data_0(i0,i1) = (1.d0-w)*un_data_0(i0,i1) + w*un_imhalf_j

          ! solve for solvent velocities
          us_imhalf_j = f_us_data_0(i0,i1) 
     &     - (D * eta_s * (us_data_1(i0,i1+1) *
     &      (toThs(thn_imh_jph) - toThs(thn_data(i0,i1))) / (dx_dy)
     &      + (toThs(thn_data(i0-1,i1)))/(dx_dx) * us_data_0(i0-1,i1) 
     &      + (toThs(thn_data(i0,i1))/dx_dx * us_data_0(i0+1,i1)) 
     &      + (toThs(thn_imh_jph)/dy_dy * us_data_0(i0,i1+1)) 
     &      + (toThs(thn_imh_jmh)/dy_dy * us_data_0(i0,i1-1))
     &      + (toThs(thn_data(i0-1,i1))-toThs(thn_imh_jph))/(dx_dy)
     &             * us_data_1(i0-1,i1+1)
     &      + (toThs(thn_imh_jmh)-toThs(thn_data(i0-1,i1)))/(dx_dy) 
     &             * us_data_1(i0-1,i1)
     &      + (toThs(thn_data(i0,i1))-toThs(thn_imh_jmh))/(dx_dy) 
     &             * us_data_1(i0,i1)) 
     &      + D * xi * un_data_0(i0,i1) * thn_lower_x
     &             * toThs(thn_lower_x))

          us_imhalf_j = us_imhalf_j/(D * eta_s * 
     &      (-(toThs(thn_data(i0,i1))+toThs(thn_data(i0-1,i1)))/(dx_dx) 
     &       - (toThs(thn_imh_jph) + toThs(thn_imh_jmh)) / (dy_dy)) 
     &       - D * xi * thn_lower_x * toThs(thn_lower_x)
     &       + C * toThs(thn_lower_x))

          us_data_0(i0,i1) = (1.d0-w)*us_data_0(i0,i1) + w*us_imhalf_j
          end if
        end do
      end do
c
c     Loop over side centers in y-dir
      do i1 = ilow1, iup1 + 1
        do i0 = ilow0, iup0
c
          if( mod(i0+i1,2) .EQ. red_or_black ) then
          ! calculate thn at (i,j-1/2)
          thn_lower_y = 0.5d0*(thn_data(i0,i1)+thn_data(i0,i1-1))  ! thn(i, j-1/2)

          ! calculate thn at corners
          thn_imh_jmh = 0.25d0*(thn_data(i0, i1) + thn_data(i0-1,i1)
     &                   +thn_data(i0,i1-1) + thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)
          thn_iph_jmh = 0.25d0*(thn_data(i0+1, i1) + thn_data(i0,i1)
     &                   +thn_data(i0,i1-1) + thn_data(i0+1,i1-1))    ! thn(i+1/2, j-1/2)

          ! solve for network velocities
          un_i_jmhalf = f_un_data_1(i0,i1) - (D * eta_n 
     &       * ((thn_imh_jmh/dx_dx) * un_data_1(i0-1,i1)
     &       + (thn_iph_jmh/dx_dx) * un_data_1(i0+1,i1)
     &       + (thn_data(i0,i1)/dy_dy) * un_data_1(i0,i1+1)
     &       + (thn_data(i0,i1-1)/dy_dy) * un_data_1(i0,i1-1)
     &       + (thn_imh_jmh - thn_data(i0,i1-1))/(dx_dy) 
     &             * un_data_0(i0,i1-1)
     &       + (thn_data(i0,i1-1) - thn_iph_jmh)/(dx_dy) 
     &             * un_data_0(i0+1,i1-1)
     &       + (thn_data(i0,i1) - thn_imh_jmh)/(dx_dy)
     &             * un_data_0(i0,i1)
     &       + (thn_iph_jmh - thn_data(i0,i1))/(dx_dy) 
     &             * un_data_0(i0+1,i1)) 
     &       + D * xi * us_data_1(i0,i1)
     &            * thn_lower_y * toThs(thn_lower_y))

          un_i_jmhalf = un_i_jmhalf/ (D * eta_n *
     &       (-(thn_data(i0,i1) + thn_data(i0,i1-1)) / (dy_dy) 
     &       - (thn_iph_jmh + thn_imh_jmh) / (dx_dx)) 
     &       - D * xi * thn_lower_y * toThs(thn_lower_y)
     &       + C * thn_lower_y)

          un_data_1(i0,i1) = (1.d0-w)*un_data_1(i0,i1) + w*un_i_jmhalf

          ! solve for solvent velocities
          us_i_jmhalf = f_us_data_1(i0,i1) - (D * eta_s 
     &       * ((toThs(thn_imh_jmh)/dx_dx) * us_data_1(i0-1,i1)
     &       + (toThs(thn_iph_jmh)/dx_dx) * us_data_1(i0+1,i1)
     &       + (toThs(thn_data(i0,i1))/dy_dy) * us_data_1(i0,i1+1)
     &       + (toThs(thn_data(i0,i1-1))/dy_dy) * us_data_1(i0,i1-1)
     &       + (toThs(thn_imh_jmh) - toThs(thn_data(i0,i1-1)))/(dx_dy) 
     &             * us_data_0(i0,i1-1)
     &       + (toThs(thn_data(i0,i1-1)) - toThs(thn_iph_jmh))/(dx_dy) 
     &             * us_data_0(i0+1,i1-1)
     &       + (toThs(thn_data(i0,i1)) - toThs(thn_imh_jmh))/(dx_dy)
     &             * us_data_0(i0,i1)
     &       + (toThs(thn_iph_jmh) - toThs(thn_data(i0,i1)))/(dx_dy) 
     &             * us_data_0(i0+1,i1)) 
     &       + D * xi * un_data_1(i0,i1)
     &            * thn_lower_y * toThs(thn_lower_y))

          us_i_jmhalf = us_i_jmhalf/ (D * eta_s * (
     &       - (toThs(thn_data(i0,i1))+toThs(thn_data(i0,i1-1)))/dy_dy
     &       - (toThs(thn_iph_jmh) + toThs(thn_imh_jmh)) / dx_dx)
     &       - D * xi * thn_lower_y * toThs(thn_lower_y)
     &       + C * toThs(thn_lower_y))

          us_data_1(i0,i1) = (1.d0-w)*us_data_1(i0,i1) + w*us_i_jmhalf
          end if
        end do
      end do
c
      end subroutine
