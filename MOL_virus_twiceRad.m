% Description: Fast solver that uses MOL
clc;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Parameters       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finite difference parameters
maxT = 1500; % Solves on time interval [0,maxT]

global xi_step;
xi_step = 0.01;  % Space step

% Diffusion parameters
global Sigma_m Sigma_l D_c; 

Sigma_m = 12;
Sigma_l = 10;
D_c = 20;

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

global k_max;
k_max = 1; 

%Boundary condition parameters 
global h_l h_n h_a l_inf n_inf c_inf a_inf phi_inf h_phi; 
h_l =  0.001;
h_n =  90;
h_a = 0.01;
l_inf = 0.2;
n_inf = 0.2; 
c_inf = 1;
a_inf = 0;
phi_inf = 0;
h_phi = 0.01;

% Macrophage delay
global mac_turn_on;
mac_turn_on = 200;

% Radiation parameters
global rad_start rad_dur rad_strength gamma rad_start2;
rad_start2 = 200;
rad_start = 150;      % Time to start radiation
rad_dur = 0.05;       % Duration of radiation
rad_strength = 0.2;   % Radiation strength
gamma = 0.05;

% Virus parameters
global D_phi p_phi_max r_rep r_phi k_phi c_phi;
D_phi = 1;
p_phi_max = 10;
r_rep = 0.125;
r_phi = 1.5;
k_phi = 0.8;
c_phi = 0.6;

% Initial radius 
R_0 = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Initial Conditions   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y0 = zeros(1,4*(1/xi_step + 1) + 1);
y0(1:1/xi_step + 1) = 0;                           % l 
y0(1/xi_step + 2:2*(1/xi_step + 1)) = 0;           % m_i
y0(2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)) = 0.8; % m_h
y0(3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)) = 1;   % c
y0(4*(1/xi_step + 1) + 1:5*(1/xi_step + 1)) = 0;   % a
y0(5*(1/xi_step + 1) + 1:6*(1/xi_step + 1)) = 0;   % phi
y0(6*(1/xi_step + 1) + 1) = R_0;                   % R




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Solve          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
options = odeset('OutputFcn',@odeprog,'Events',@odeabort);
[t,Y]=ode15s(@MOL_deriv_virus_twiceRad,[0 maxT],y0,options);
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Results      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
subplot(2,2,1)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,1:1/xi_step + 1),'EdgeColor','none')
title('Macrophages (l)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,2)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,1/xi_step + 2:2*(1/xi_step + 1)),'EdgeColor','none')
title('Infected Tumour Cells (m_i)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,3)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)),'EdgeColor','none')
title('Uninfected Tumour Cells (m_h)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,4)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)),'EdgeColor','none')
title('Oxygen (c)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud


figure;
subplot(2,2,1)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,4*(1/xi_step + 1) + 1:5*(1/xi_step + 1)),'EdgeColor','none')
title('Chemoattractant (a)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud


subplot(2,2,2)
surf(linspace(0,1,1/xi_step + 1),t,Y(:,5*(1/xi_step + 1) + 1:6*(1/xi_step + 1)),'EdgeColor','none')
axis square;
title('Virus (phi)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,3)
plot(t,Y(:,6*(1/xi_step + 1) + 1))
axis square;
title('Radius (R)','FontSize',18)
xlabel('Time','FontSize',14)
ylabel('Tumour Radius','FontSize',14)
lighting gouraud   
