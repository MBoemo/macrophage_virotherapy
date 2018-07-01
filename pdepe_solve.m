function pdepe_solve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Parameters           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finite difference parameters
maxT = 100; % Solves on time interval [0,maxT]

global x_step;
x_step = 0.1;  % Space step
t_step = 0.01; % Time step


% Finite difference mesh
global xi;
xi = linspace(0,1,1/x_step + 1);       % Space mesh
t = linspace(0,maxT,maxT/t_step + 1);  % Time mesh


% Diffusion parameters
global Sigma_m Sigma_l D_c; 

Sigma_m = 10;
Sigma_l = 30;%10;
D_c = 50;

% Chemotaxis force 
global k;
k = 1;%0.01;


% Growth/death parameters
global c_p p_m_max c_c p_al_max n_0 alpha p_am_max d_a d_m_max d_l_max d_cm_max d_cl_max d_cp_max chi_l m_p;
c_p = 0.6;
c_c = 0.2;
n_0 = 0.2;
alpha = 5;
p_m_max = 1;%0.8; 
p_am_max = 1;
p_al_max = 1;
d_a = 0.01;
d_m_max = 5;
d_l_max = 2;
d_cm_max = 0.5;%0.005; 
d_cl_max = 0.5;%0.005;
d_cp_max = 0.1;%0.01;
chi_l = 2000;%500;
m_p = 0.6;

global k_max a_char;
k_max = 6; 
a_char = 1;


%Boundary condition parameters 
global h_l h_n h_a l_inf n_inf c_inf a_inf; 
h_l = 0.05;%0.001; 
h_n = 5;%100; 
h_a = 0.01;
l_inf = 0;%0.2;
n_inf = 0.2; 
c_inf = 1;
a_inf = 0;


% Vector to store tumour radius values
global R r_counter;
R = zeros(maxT/t_step + 1,1);
R(1) = 13; % Initial radius
r_counter = 1;

global DRdt l_final c_final a_final m_final;

% At t = 0, DRdt = 0
DRdt = 0;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Solve PDEs           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocate solution arrays
m_final = zeros(maxT/t_step,1/x_step + 1);
l_final = zeros(maxT/t_step,1/x_step + 1);
c_final = zeros(maxT/t_step,1/x_step + 1);
a_final = zeros(maxT/t_step,1/x_step + 1);

global m_temp c_temp l_temp a_temp;


wait = waitbar(0,'Solving PDE system...');





for i = 1:maxT/t_step
    
    % pdepe requires a minimum of three time steps, so we introduce an
    % "intermediate" time step t + dt/2
    
      
    sol = pdepe(2,@pde,@pdeIC,@pdeBC,xi,[t(i) t(i)+(t_step/2) t(i)+t_step]);
    
    % Extract solution for each variable
    l_temp = sol(end,:,1);
    m_temp = sol(end,:,2);
    c_temp = sol(end,:,3);
    a_temp = sol(end,:,4);
    
    
    % Controls floating point errors when m is near zero 
    m_temp(m_temp < 0) = 0;
    
    
    % Use it to update R via a Forward Euler approximation
    R(r_counter + 1) = R(r_counter) + t_step*(...
        (1/R(r_counter))*( ((m_temp(1,end) - 1)/(k*m_temp(1,end)))*Sigma_m*(m_temp(1,end) - m_temp(1,end-1))/x_step + ...
        (1/k)*Sigma_l*(l_temp(1,end) - l_temp(1,end-1))/x_step - ...
        (1/k)*l_temp(1,end)*(a_temp(1,end) - a_temp(1,end-1))/x_step) );

       
    
    % Update DRdt
    DRdt = (R(r_counter + 1) - R(r_counter))/t_step;
    
    
    % Update total solution
    l_final(i,:) = l_temp;
    m_final(i,:) = m_temp; 
    c_final(i,:) = c_temp; 
    a_final(i,:) = a_temp;

     
    waitbar(i/(maxT/t_step))

    r_counter = r_counter + 1;
end

close(wait);

% Surface plots
figure;
surf(xi,t(1:end-1),l_final,'EdgeColor','none') 
title('Macrophages')
xlabel('xi (transformed space)')
ylabel('tau (transformed time)')
zlabel('volume fraction')

figure;
surf(xi,t(1:end-1),m_final,'EdgeColor','none') 
title('Tumour Cells')
xlabel('xi (transformed space)')
ylabel('tau (transformed time)')
zlabel('volume fraction')

figure;
surf(xi,t(1:end-1),a_final,'EdgeColor','none') 
title('Chemoattractant')
xlabel('xi (transformed space)')
ylabel('tau (transformed time)')
zlabel('volume fraction')

figure;
surf(xi,t(1:end-1),c_final,'EdgeColor','none') 
title('Oxygen')
xlabel('xi (transformed space)')
ylabel('tau (transformed time)')
zlabel('volume fraction')

figure;
plot(t,R)
title('Tumour Radius')
xlabel('tau (transformed time)')
ylabel('Radius')











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Define PDEs          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,f,s] = pde(x,t,u,DuDx)
global Sigma_l Sigma_m k R r_counter d_a DRdt D_c chi_l;

% Notation:
% l = u(1)
% m = u(2)
% c = u(3)
% a = u(4)



% Time derivative coupling
c = [1;  % l
     1;  % m
     1;  % c
     1]; % a

 
% Flux
f =  [(1/(R(r_counter)^2))*( (Sigma_l/k)*DuDx(1) - chi_l*((1 - u(1))/k)*u(1)*DuDx(4) - (u(1)/k)*(Sigma_m*DuDx(2) + Sigma_l*DuDx(1)) ); % l
     (1/(R(r_counter)^2))*( (Sigma_m/k)*DuDx(2) + chi_l*(u(2)/k)*u(1)*DuDx(4) - (u(2)/k)*(Sigma_m*DuDx(2) + Sigma_l*DuDx(1)) );        % m
     (1/(R(r_counter)^2))*(D_c*DuDx(3));                                                                                               % c
     (1/(R(r_counter)^2))*(DuDx(4))];                                                                                                  % a

 
% Source
s = [(x/R(r_counter))*DRdt*DuDx(1) - u(1)*fun_d_l(u(3));                                                           % l 
     (x/R(r_counter))*DRdt*DuDx(2) + u(2)*fun_p_m_fix( u(2), u(3) ) - u(2)*fun_d_m(u(3)) - u(1)*u(2)*fun_k(u(3));  % m
     (x/R(r_counter))*DRdt*DuDx(3) - fun_d_c( u(2), u(1), 1-u(2)-u(1),u(3));                                       % c
     (x/R(r_counter))*DRdt*DuDx(4) + d_a*(fun_p_a(u(1),u(2),u(3)) - u(4))];                                        % a


 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Initial Conditions   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u0 = pdeIC(x)
global xi l_temp m_temp c_temp a_temp r_counter;

if r_counter == 1
    u0 = [0;     % l
          0.8;   % m
          1;     % c
          0];    % a
else
    u0 = [interp1(xi',l_temp', x);
          interp1(xi',m_temp', x);
          interp1(xi',c_temp', x);
          interp1(xi',a_temp', x)];
end











%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Boundary Conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,ql,pr,qr] = pdeBC(xl,ul,xr,ur,t)
global Sigma_m Sigma_l k R r_counter h_l h_n h_a l_inf n_inf a_inf; 




a_flux = (1/R(r_counter))*h_a*(a_inf - ur(4));

m_flux = ( (1/R(r_counter))*k*ur(2)*h_n*( 1 - ur(1) - ur(2) - n_inf ) )/(Sigma_m*(1 - ur(1) - ur(2)));

l_flux = (1/R(r_counter))*( ((k*h_l)/Sigma_l)*(l_inf - ur(1)) + (1/Sigma_l)*ur(1)*a_flux + (ur(1)*Sigma_m)/(ur(2)*Sigma_l)*m_flux);

pl = [0;  % l
      0;  % m
      0;  % c
      0]; % a

ql = [1;  % l
      1;  % m
      1;  % c
      1]; % a

pr = [-(1/R(r_counter)^2)*( -((1 - ur(1))/k)*ur(1)*a_flux - (ur(1)/k)*Sigma_m*m_flux ) - (1/R(r_counter)^2)*( Sigma_l/k - (ur(1)/k)*Sigma_l )*l_flux; % l
      -(1/R(r_counter)^2)*( (ur(2)/k)*ur(1)*a_flux - (ur(2)/k)*Sigma_l*l_flux ) - (1/R(r_counter)^2)*( Sigma_m/k - (ur(2)/k)*Sigma_m )*m_flux;        % m
      ur(3) - 1;                                                                                                                                      % c
      -a_flux];                                                                                                                                       % a

qr = [1; % l
      1; % m
      0; % c
      R(r_counter)^2];% a
