%%Visualizing Data from IBAMR Staggered Stokes Operator Test


%Load the HDF5 data ...

clear; 
close all;  
clc;

%Compute what I'd expect the staggered stokes operator to give me
dx = 1/32;
dy = 1/32;

x_sc = 0:dx:1;
y_sc = 0:dy:1;



x_cc = 0.5*(x_sc(1:end-1) + x_sc(2:end));
y_cc = 0.5*(y_sc(1:end-1) + y_sc(2:end));
[X_s,Y_s] = meshgrid(x_sc,y_cc);
[X_c,Y_c] = meshgrid(x_cc,y_cc);
[X_v,Y_v] = meshgrid(x_cc,y_sc);

%Exact solutions for plane poisuelle flow are ...
u_exact = @(x,y) 0.5*y.*(1-y);
p_exact = @(x,y) 1-x;
v_exact = @(x,y) 0.*y;
u = zeros(32,33);

%initialize u,p,and v
u = u_exact(X_s,Y_s);
p = p_exact(X_c,Y_c);
v = v_exact(X_c,Y_v);



% add in the ghost nodes
u_g = zeros(34,33);
u_g(2:end-1,:) = u;
u_g(1,:) = -u_g(2,:);
u_g(end,:) = -u_g(end-1,:);

v_g = zeros(33,34);
v_g(:,2:end-1) = v;
v_g(:,1) = -v_g(:,2);
v_g(:,end) = -v_g(:,end-1);

%compute the laplacian of u and v now
lap_u = (32)*(32)*(u_g(1:end-2,2:end-1)-2*u_g(2:end-1,2:end-1)+u_g(3:end,2:end-1)) + ...
    (32)*(32)*(u_g(2:end-1,1:end-2) - 2*u_g(2:end-1,2:end-1) + u_g(2:end-1,3:end));



%compute the gradient of the pressure
D_xp = (32)*(p(:,2:end)-p(:,1:end-1));
D_yp = (32)*(p(2:end,:)-p(1:end-1,:));

Stag_Stokes_u = lap_u + D_xp + u(:,2:end-1);
Stag_Stokes_u = [u(:,1), Stag_Stokes_u, u(:,end)];

%Now interpolate Stag_Stokes onto the cell centers so that I can compare to
%the visit data
Stag_Stokes_cc_u = 0.5*(Stag_Stokes_u(:,1:end-1) + Stag_Stokes_u(:,2:end));












