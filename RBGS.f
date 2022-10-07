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
     &        dx, ilower0, iupper0,
     &        ilower1, iupper1, 
     &        un_data_0, un_data_1, un_gcw,
     &        us_data_0, us_data_1, us_gcw,
     &        p_data, p_gcw, f_p_data, f_p_gcw,
     &        f_un_data_0, f_un_data_1, f_un_gcw,
     &        f_us_data_0, f_us_data_1, f_us_gcw,
     &        thn_data, thn_gcw,  
     &        eta_n, eta_s, nu_n, nu_s, xi)
c
      implicit none
cccccccccccccccccccccccccccccccccc INPUTS ccccccccccccccccccccccccccccccccc
      double precision dx(0:1)
      integer ilower0,  iupper0  
      integer iupper1,  ilower1
      integer un_gcw,  us_gcw, p_gcw, f_un_gcw, f_us_gcw, f_p_gcw, thn_gcw
      double precision eta_n, eta_s, nu_n, nu_s, xi

      double precision thn_data(ilower0-thn_gcw:iupper0+thn_gcw,
     &          ilower1-thn_gcw:iupper1+thn_gcw) 
c
      double precision un_data_0(ilower0-un_gcw:iupper0+un_gcw+1,  
     &          ilower1-un_gcw:iupper1+un_gcw)
      double precision un_data_1(ilower0-un_gcw:iupper0+un_gcw,
     &          ilower1-un_gcw:iupper1+un_gcw+1)
c
     double precision us_data_0(ilower0-us_gcw:iupper0+us_gcw+1,
     &          ilower1-us_gcw:iupper1+us_gcw)
      double precision us_data_1(ilower0-us_gcw:iupper0+us_gcw,
     &          ilower1-us_gcw:iupper1+us_gcw+1)
c
     double precision p_data(ilower0-p_gcw:iupper0+p_gcw,
     &          ilower1-p_gcw:iupper1+p_gcw)

     double precision f_un_data_0(ilower0-f_un_gcw:iupper0+f_un_gcw+1,
     &          ilower1-f_un_gcw:iupper1+f_un_gcw)
      double precision f_un_data_1(ilower0-f_un_gcw:iupper0+f_un_gcw,
     &          ilower1-f_un_gcw:iupper1+f_un_gcw+1)
c
     double precision f_us_data_0(ilower0-f_us_gcw:iupper0+f_us_gcw+1,
     &          ilower1-f_us_gcw:iupper1+f_us_gcw)
      double precision f_us_data_1(ilower0-f_us_gcw:iupper0+f_us_gcw,
     &          ilower1-f_us_gcw:iupper1+f_us_gcw+1)
c
     double precision f_p_data(ilower0-f_p_gcw:iupper0+f_p_gcw,
     &          ilower1-f_p_gcw:iupper1+f_p_gcw)
c
    double precision thn_lower_x, thn_lower_y, thn_upper_x, thn_upper_y
    double precision thn_imhalf_jphalf, thn_imhalf_jmhalf
    double precisionthn_iphalf_jphalf, thn_iphalf_jmhalf
    integer ipiv(3), info
c
    integer i0, i1
c
      do i1 = ilower1, iupper1   
        do i0 = ilower0, iupper0  ! same loop as the c++ code
c
          ! calculate thn at sides
          thn_lower_x = 0.5*(thn_data(i0,i1)+thn_data(i0-1,i1))  ! thn(i-1/2, j)
          thn_upper_x = 0.5*(thn_data(i0,i1)+thn_data(i0+1,i1))  ! thn(i+1/2, j)
          thn_lower_y = 0.5*(thn_data(i0,i1)+thn_data(i0,i1-1))  ! thn(i, j-1/2)
          thn_upper_y = 0.5*(thn_data(i0,i1)+thn_data(i0,i1+1))  ! thn(i, j+1/2)
c
          ! calculate thn at corners
          thn_imhalf_jphalf = 0.25*(thn_data(i0-1, i1)+thn_data(i0,i1)
        &                     +thn_data(i0,i1+1)+thn_data(i0-1,i1+1))   ! thn(i-1/2, j+1/2)
          thn_imhalf_jmhalf = 0.25*(thn_data(i0, i1)+thn_data(i0-1,i1)
        &                     +thn_data(i0,i1-1)+thn_data(i0-1,i1-1))   ! thn(i-1/2, j-1/2)
          thn_iphalf_jphalf = 0.25*(thn_data(i0+1, i1)+thn_data(i0,i1)
        &                     +thn_data(i0,i1+1)+thn_data(i0+1,i1+1))   ! thn(i+1/2, j+1/2)
          thn_iphalf_jmhalf = 0.25*(thn_data(i0+1, i1)+thn_data(i0,i1)
        &                     +thn_data(i0,i1-1)+thn_data(i0+1,i1-1))   ! thn(i+1/2, j-1/2)
c
          ! set-up the box matrix A_box
          double precision A_box(9,9)
          double precision b(9,9)
          integer ipiv(3), info
C
          ! network at west edge
          A_box(1, 1) = eta_n / (dx(0)*dx(0)) * 
        &      (-(thn_data(i0,i1) - thn_data(i0-1,i1))       
        &      - eta_n/(dx(1)*dx(1))*(thn_imhalf_jmhalf + thn_imhalf_jphalf)
        &      - (xi / nu_n) * thn_lower_x*toThs(thn_lower_x)
          A_box(1, 2) = eta_n / (dx(0)*dx(0))*thn_data(i0,i1)
          A_box(1, 3) = eta_n / (dx(0) * dx(1)) * 
        &               (thn_data(i0,i1)-thn_imhalf_jmhalf)
          A_box(1, 4) = eta_n / (dx(0) * dx(1)) *
        &               (thn_imhalf_jphalf - thn_data(i0,i1))
          A_box(1, 5) = (xi / nu_n)*thn_lower_x*toThs(thn_lower_x)
          A_box(1, 6) = 0.0
          A_box(1, 7) = 0.0
          A_box(1, 8) = 0.0
          A_box(1, 9) = -thn_lower_x / dx(0)
c
          A_box(5, 5) = eta_s / (dx(0) * dx(0)) * 
        & (-toThs(thn_data(i0,i1))-toThs(thn_data(i0-1,i1))) 
          - eta_s/(dx(1)*dx(1)) *
        & (toThs(thn_imhalf_jmhalf)+toThs(thn_imhalf_jphalf)) 
        &  -(xi / nu_s)*thn_lower_x*toThs(thn_lower_x)
          A_box(5, 6) = eta_s/(dx(0)*dx(0)) * toThs(thn_data(i0,i1))
          A_box(5, 7) =  eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_data(i0,i1))-toThs(thn_imhalf_jmhalf))
          A_box(5, 8) = eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_imhalf_jphalf) - toThs(thn_data(i0,i1)))
          A_box(5, 1) = xi/nu_s * thn_lower_x * toThs(thn_lower_x)
          A_box(5, 2) = 0.0
          A_box(5, 3) = 0.0
          A_box(5, 4) = 0.0
          A_box(5, 9) = -toThs(thn_lower_x) / dx(0)
c
          ! network at east edge
          A_box(2, 1) = eta_n / (dx(0) * dx(0)) * thn_data(i0,i1)
          A_box(2, 2) = eta_n / (dx(0) * dx(0)) * 
        & (-thn_data(i0+1,i1) - thn_data(i0,i1)) -
        & eta_n / (dx(1) * dx(1)) * (thn_iphalf_jphalf + thn_iphalf_jmhalf) 
        & -(xi / nu_n) * thn_upper_x * toThs(thn_upper_x)
          A_box(2, 3) = eta_n / (dx(1) * dx(0)) * 
        & (thn_iphalf_jmhalf - thn_data(i0,i1))
          A_box(2, 4) = eta_n / (dx(1) * dx(0)) * 
        & (thn_data(i0,i1)-thn_iphalf_jphalf)
          A_box(2, 5) = 0.0
          A_box(2, 7) = 0.0
          A_box(2, 8) = 0.0
          A_box(2, 6) = xi/nu_n * thn_upper_x * toThs(thn_upper_x)
          A_box(2, 9) = thn_upper_x / dx(0)
c
          A_box(6, 5) = eta_s/(dx(0)*dx(0)) * toThs(thn_data(i0,i1))
          A_box(6, 6) = eta_s/(dx(0)*dx(0)) * 
        & (-toThs(thn_data)(i0+1,i1))-toThs(thn_data(i0,i1))) 
        & - eta_s/(dx(1)*dx(1)) * 
        & (toThs(thn_iphalf_jphalf)+toThs(thn_iphalf_jmhalf)) 
        & - xi/nu_s * thn_upper_x * toThs(thn_upper_x)
          A_box(6, 7) = eta_s / (dx(1)*dx(0)) * 
        & (toThs(thn_iphalf_jmhalf)-toThs(thn_data(i0,i1)))
          A_box(6, 8) = eta_s/(dx(1)*dx(0)) *
        & (toThs(thn_data(i0,i1))-toThs(thn_iphalf_jphalf))
          A_box(6, 1) = 0.0
          A_box(6, 3) = 0.0
          A_box(6, 4) = 0.0
          A_box(6, 2) = xi/nu_s * thn_upper_x * toThs(thn_upper_x)
          A_box(6, 9) = toThs(thn_upper_x) / dx(0)
c
          ! network at south edge
          A_box(2, 0) = eta_n / (dx(0) * dx(1)) * 
        & (thn_data(i0,i1)-thn_imhalf_jmhalf)
          A_box(2, 1) = eta_n / (dx(0) * dx(1)) * 
        & (thn_iphalf_jmhalf - thn_data(i0,i1))
          A_box(2, 2) = eta_n / (dx(1) * dx(1)) * 
        & (-thn_data(i0,i1) - (*thn_data)(i0,i1-1)) 
        & - eta_n/(dx(0)*dx(0))*(thn_iphalf_jmhalf + thn_imhalf_jmhalf) 
        & - xi/nu_n * thn_lower_y * toThs(thn_lower_y)
          A_box(2, 3) = eta_n / (dx(1) * dx(1)) * (thn_data(i0,i1))
          A_box(2, 4) = 0.0
          A_box(2, 5) = 0.0
          A_box(2, 7) = 0.0
          A_box(2, 6) = xi/nu_n * thn_lower_y * toThs(thn_lower_y)
          A_box(2, 8) = -thn_lower_y / dx(1)
c
          A_box(6, 4) = eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_data(i0,i1)) - toThs(thn_imhalf_jmhalf))
          A_box(6, 5) = eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_iphalf_jmhalf) - toThs(thn_data(i0,i1)))
          A_box(6, 6) = eta_s / (dx(1) * dx(1)) * 
        & (-toThs(thn_data(i0,i1)) - toThs(thn_data(i0,i1-1))) 
        & - eta_s / (dx(0) * dx(0)) * 
        & (toThs(thn_iphalf_jmhalf) + toThs(thn_imhalf_jmhalf)) 
        & - xi/nu_s * thn_lower_y * toThs(thn_lower_y)
          A_box(6, 7) = eta_s /(dx(1)*dx(1))*toThs(thn_data(i0,i1))
          A_box(6, 0) = 0.0
          A_box(6, 1) = 0.0
          A_box(6, 3) = 0.0
          A_box(6, 2) = xi/nu_s * thn_lower_y * toThs(thn_lower_y)
          A_box(6, 8) = -toThs(thn_lower_y) / dx(1)
c
          ! network at north edge
          A_box(4, 1) = eta_n / (dx(0) * dx(1)) * 
        & (thn_imhalf_jphalf - thn_data(i0,i1))
          A_box(4, 2) = eta_n / (dx(0) * dx(1)) * 
        & (thn_data(i0,i1)-thn_iphalf_jphalf)
          A_box(4, 3) = eta_n / (dx(1) * dx(1)) * (thn_data(i0,i1))
          A_box(4, 4) = eta_n / (dx(1) * dx(1)) * 
        & (-thn_data(i0,i1) - thn_data(i0,i1+1)) 
        & - eta_n/(dx(0)*dx(0))*(thn_iphalf_jphalf + thn_imhalf_jphalf) 
        & - xi/nu_n * thn_upper_y * toThs(thn_upper_y)
          A_box(4, 5) = 0.0
          A_box(4, 6) = 0.0
          A_box(4, 7) = 0.0
          A_box(4, 8) = xi/nu_n * thn_upper_y * toThs(thn_upper_y)
          A_box(4, 9) = thn_upper_y / dx(1)
c
          A_box(8, 5) = eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_imhalf_jphalf)-toThs(thn_data(i0,i1)))
          A_box(8, 6) = eta_s / (dx(0) * dx(1)) * 
        & (toThs(thn_data(i0,i1))-toThs(thn_iphalf_jphalf))
          A_box(8, 7) = eta_s/(dx(1)*dx(1))*toThs(thn_data(i0,i1))
          A_box(8, 8) = eta_s / (dx(1) * dx(1)) * 
        & (-toThs(thn_data(i0,i1))-toThs(thn_data(i0,i1+1))) 
        & - eta_s / (dx(0) * dx(0)) * 
        & (toThs(thn_iphalf_jphalf)+toThs(thn_imhalf_jphalf)) 
        & - xi/nu_s * thn_upper_y * toThs(thn_upper_y)
          A_box(8, 1) = 0.0
          A_box(8, 2) = 0.0
          A_box(8, 3) = 0.0
          A_box(8, 4) = xi/nu_s * thn_upper_y * toThs(thn_upper_y)
          A_box(8, 9) = toThs(thn_upper_y) / dx(1)
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
          b(1) = f_un_data_0(i0,i1)-thn_lower_x / dx(0)* p_data(i0-1,i1) -
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
          &        (un_data_1(i0-1,i1+1) - un_data_1(i0-1,i1))
c
          ! solvent at west edge
          b(5) = f_us_data_0(i0,i1)-toThs(thn_lower_x) / dx(0) * 
          &      p_data(i0-1,i1) -
          &      eta_s / (dx(0)* dx(0)) * 
          &      toThs(thn_data(i0-1,i1)) * us_data_0(i0-1,i1) -
          &      eta_s / (dx(1)* dx(1)) * 
          &      toThs(thn_imhalf_jphalf) * us_data_0(i0,i1 + 1) -
          &      eta_s / (dx(1)* dx(1)) * 
          &      toThs(thn_imhalf_jmhalf) * us_data_0(i0,i1 - 1) +
          &      eta_s / (dx(0)* dx(1)) * 
          &      toThs(thn_imhalf_jphalf) * us_data_1(i0-1,i1+1) -
          &      eta_s / (dx(0)* dx(1)) * 
          &      toThs(thn_imhalf_jmhalf) * us_data_1(i0-1,i1) -
          &      eta_s / (dx(0)* dx(1)) * 
          &      toThs(thn_data(i0-1,i1)) *
          &      (us_data_1(i0-1,i1+1) - us_data_1(i0-1,i1))
c
          ! network at east edge
          b(2) = f_un_data_0(i0+1,i1) + thn_upper_x / dx(0)* p_data(i0+1,i1) -
          &        eta_n / (dx(0)* dx(0)) * thn_data(i0+1,i1) * un_data_0(i0+2,i1) -
          &        eta_n / (dx(1)* dx(1)) * thn_iphalf_jphalf * un_data_0(i0+1,i1+1) -
          &        eta_n / (dx(1)* dx(1)) * thn_iphalf_jmhalf * un_data_0(i0+1,i1-1) -
          &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jphalf * un_data_1(i0+1,i1+1) +
          &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jmhalf * un_data_1(i0+1,i1) +
          &        eta_n / (dx(0)* dx(1)) * thn_data(i0+1,i1) *
          &        (un_data_1(i0+1,i1+1) - un_data_1(i0+1,i1))
c
          ! solvent at east edge
          b(6) = f_us_data_0(i0+1,i1) + toThs(thn_upper_x) / dx(0)* p_data(i0+1,i1) -
          &    eta_s / (dx(0)* dx(0)) * toThs(thn_data(i0+1,i1)) * us_data_0(i0+2,i1) -
          &    eta_s / (dx(1)* dx(1)) * toThs(thn_iphalf_jphalf) * us_data_0(i0+1,i1+1) -
          &    eta_s / (dx(1)* dx(1)) * toThs(thn_iphalf_jmhalf) * us_data_0(i0+1,i1-1) -
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jphalf) * us_data_1(i0+1,i1+1) +
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jmhalf) * us_data_1(i0+1,i1) +
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0+1,i1)) *
          &        (us_data_1(i0+1,i1+1) - us_data_1(i0+1,i1))
c
          ! network at south edge
          b(3) = f_un_data_1(i0,i1)-thn_lower_y / dx(1)* p_data(i0,i1-1) -
          &        eta_n / (dx(1)* dx(1)) * thn_data(i0,i1-1) * un_data_1(i0,i1-1) -
          &        eta_n / (dx(0)* dx(0)) * thn_iphalf_jmhalf * un_data_1(i0+1,i1) -
          &        eta_n / (dx(0)* dx(0)) * thn_imhalf_jmhalf * un_data_1(i0-1,i1) +
          &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jmhalf * un_data_0(i0+1,i1-1) -
          &        eta_n / (dx(0)* dx(1)) * thn_imhalf_jmhalf * un_data_0(i0,i1-1) -
          &        eta_n / (dx(0)* dx(1)) * thn_data(i0,i1-1) *
          &            (un_data_0(i0+1,i1-1) - un_data_0(i0,i1-1))
c
          ! solvent at south edge
          b(7) = f_us_data_1(i0,i1)-toThs(thn_lower_y) / dx(1)* p_data(i0,i1-1) -
          &    eta_s / (dx(1)* dx(1)) * toThs(thn_data(i0,i1-1)) * us_data_1(i0,i1-1) -
          &    eta_s / (dx(0)* dx(0)) * toThs(thn_iphalf_jmhalf) * us_data_1(i0+1,i1) -
          &    eta_s / (dx(0)* dx(0)) * toThs(thn_imhalf_jmhalf) * us_data_1(i0-1,i1) +
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jmhalf) * us_data_0(i0+1,i1-1) -
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_imhalf_jmhalf) * us_data_0(i0,i1-1) -
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0,i-1)) *
          &        (us_data_0(i0+1,i1-1) - us_data_0(i0,i1-1))
c
          ! network at north edge
          b(4) = f_un_data_1(i0,i1+1) + thn_upper_y / dx(1)* p_data(i0,i1+1) -
          &        eta_n / (dx(1)* dx(1)) * thn_data(i0,i1+1) * un_data_1(i0,i1+2) -
          &        eta_n / (dx(0)* dx(0)) * thn_iphalf_jphalf * un_data_1(i0+1,i1+1) -
          &        eta_n / (dx(0)* dx(0)) * thn_imhalf_jphalf * un_data_1(i0-1,i1+1) -
          &        eta_n / (dx(0)* dx(1)) * thn_iphalf_jphalf * un_data_0(i0+1,i1+1) +
          &        eta_n / (dx(0)* dx(1)) * thn_imhalf_jphalf * un_data_0(i0,i1 + 1) +
          &        eta_n / (dx(0)* dx(1)) * thn_data(i0,i1+1) *
          &            (un_data_0(i0+1,i1+1) - un_data_0(i0,i1 + 1))
c
          ! solvent at north edge
          b(8) =  f_us_data_1(i0,i1+1) + toThs(thn_upper_y) / dx(1)* p_data(i0,i1+1) -
          &    eta_s / (dx(1)* dx(1)) * toThs((*thn_data)(idx + yp)) * us_data_1(upper_y_idx + yp) -
          &    eta_s / (dx(0)* dx(0)) * toThs(thn_iphalf_jphalf) * us_data_1(upper_y_idx + xp) -
          &    eta_s / (dx(0)* dx(0)) * toThs(thn_imhalf_jphalf) * us_data_1(upper_y_idx - xp) -
          &   eta_s / (dx(0)* dx(1)) * toThs(thn_iphalf_jphalf) * us_data_0(upper_x_idx + yp) +
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_imhalf_jphalf) * us_data_0(i0,i1 + 1) +
          &    eta_s / (dx(0)* dx(1)) * toThs(thn_data(i0,i1+1)) *
          &        (us_data_0(i0+1,i1+1) - us_data_0(i0,i1 + 1))
c
          ! pressure at cell center
          b(9) = f_p_data(i0,i1)
c
          ! solve the system
          call dgesv(9, 1, A, 9, ipiv, b, 9, info)
c
          if (info /= 0) then
            print *, 'ERROR IN DGESV'
          endif
c
          un_data_0(i0,i1) = b(1);
          un_data_0(i0+1,i1) = b(2);
          un_data_1(i0,i1) = b(3);
          un_data_1(i0,i1+1) = b(4);
          us_data_0(i0,i1) = b(5);
          us_data_0(i0+1,i1) = b(6);
          us_data_1(i0,i1) = b(7);
          us_data_1(i0,i1+1) = b(8);
          p_data(i0,i1) = b(9);
c
        enddo
      enddo
c     
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Function to convert Thn to Ths since Thn + Ths = 1 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      double precision function toThs(Thn)
c      
        double precision Thn, Ths
        toThs = 1 - Thn
c
        return
      end 