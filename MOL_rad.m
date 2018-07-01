% Description: Fast solver that implements MOL using ODE45

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Parameters           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finite difference parameters
global maxT;
maxT = 1000; % Solves on time interval [0,maxT]

global xi_step;
xi_step = 0.01;  % Space step

% Diffusion parameters
global Sigma_m Sigma_l D_c; 

Sigma_m = 12;%120;
Sigma_l = 10;%100;
D_c = 20;%120;

% Chemotaxis force 
global k;
k = 1;

% Growth/death parameters
global c_p p_m_max c_c p_al_max n_0 alpha p_am_max d_a d_m_max d_l_max d_cm_max d_cl_max d_cp_max chi_l m_p;
c_p = 0.6;
c_c = 0.2;
n_0 = 0.2;
alpha = 5;
p_m_max = 0.1;
p_am_max = 1;
p_al_max = 1;
d_a = 0.01;
d_m_max = 5.5;
d_l_max = 1;
d_cm_max = 0.5;
d_cl_max = 0.5;
d_cp_max = 0.1;
chi_l =  1100;
m_p = 0.65;

global k_max a_char;
k_max = 1; 
a_char = 1;

%Boundary condition parameters 
global h_l h_n h_a l_inf n_inf c_inf a_inf; 
h_l =  0.05;
h_n =  90;
h_a = 0.01;
l_inf = 0.2;
n_inf = 0.2; 
c_inf = 1;
a_inf = 0;

% Radiation parameters
global rad_start rad_dur rad_strength gamma;
rad_start = 150;       % Time to start radiation
rad_dur = 0.05;       % Duration of radiation
rad_strength = 0.2;%0.1;%0.5;   % Radiation strength
gamma = 0.05;

% Initial radius 
R_0 = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Initial Conditions     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y0 = zeros(1,4*(1/xi_step + 1) + 1);
y0(1:1/xi_step + 1) = 0;                         % l 
y0(1/xi_step + 2:2*(1/xi_step + 1)) = 0.8;       % m
y0(2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)) = 1; % c
y0(3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)) = 0; % a
y0(4*(1/xi_step + 1) + 1) = R_0;                 % R

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Solve            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('OutputFcn',@odeprog,'Events',@odeabort);
[x,Y]=ode15s(@MOL_deriv_rad,[0 maxT],y0,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot Results         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,1:1/xi_step + 1),'EdgeColor','none')
title('Macrophages (l)')
ylabel('Time')
xlabel('Distance from Center')
lighting gouraud

subplot(2,2,2)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,1/xi_step + 2:2*(1/xi_step + 1)),'EdgeColor','none')
title('Tumour Cells (m)')
ylabel('Time')
xlabel('Distance from Center')
lighting gouraud

subplot(2,2,3)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)),'EdgeColor','none')
title('Oxygen (c)')
ylabel('Time')
xlabel('Distance from Center')
lighting gouraud

subplot(2,2,4)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)),'EdgeColor','none')
title('Chemoattractant (a)')
ylabel('Time')
xlabel('Distance from Center')
lighting gouraud

figure;
plot(x,Y(:,4*(1/xi_step + 1) + 1))
title('Radius (R)')
xlabel('Time')
ylabel('Tumour Radius')
lighting gouraud
            