% PRODUCES FIGURES 1 AND 2 in Boemo,M. and Byrne,H.M. 2016
% Run plotter.m on the result of this script to produce Figure 2

% DESCRIPTION: Solves the mixture model with a set of parameters that reproduces the
% data from Leek1999.  No treatment is included in this version of the
% model: no radiation or oncolytic virus.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Parameters           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finite difference parameters
global maxT;
maxT = 50; % Solves on time interval [0,maxT]

global xi_step;
xi_step = 0.01;  % Space step

% Diffusion parameters
global Sigma_m Sigma_l D_c; 

Sigma_m = 1000; %500
Sigma_l = 400; %300
D_c = 400; %500

% Chemotaxis force 
global k;
k = 1;

% Growth/death parameters
global c_p p_m_max c_c p_al_max n_0 alpha p_am_max d_a d_m_max d_l_max d_cm_max d_cl_max d_cp_max chi_l m_p;
c_p = 0.6;
c_c = 0.2;
n_0 = 0.2;
alpha = 5;
p_m_max = 1;
p_am_max = 1;
p_al_max = 1;
d_a = 0.01;
d_m_max = 0.3;%5
d_l_max = 1;
d_cm_max = 0.5;
d_cl_max = 0.5;
d_cp_max = 0.1;
chi_l =  2000; % 2000 %1100
m_p = 0.65;

global k_max a_char;
k_max = 1.5; 
a_char = 1;

%Boundary condition parameters 
global h_l h_n h_a l_inf n_inf c_inf a_inf; 
h_l =  5;
h_n =  90; %86
h_a = 0.01;
l_inf = 0.2;
n_inf = 0.2;%0.2; 
c_inf = 1;
a_inf = 0;

% Initial radius 
R_0 = 130; %13


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

tic;
options = odeset('OutputFcn',@odeprog,'Events',@odeabort);
[x,Y]=ode15s(@MOL_deriv,[0 maxT],y0,options);
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot Results         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,2,1)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,1:1/xi_step + 1),'EdgeColor','none')
title('Macrophages (l)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,2)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,1/xi_step + 2:2*(1/xi_step + 1)),'EdgeColor','none')
title('Tumour Cells (m)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,3)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)),'EdgeColor','none')
title('Oxygen (c)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

subplot(2,2,4)
surf(linspace(0,1,1/xi_step + 1),x,Y(:,3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)),'EdgeColor','none')
title('Chemoattractant (a)','FontSize',18)
ylabel('Time','FontSize',14)
xlabel('Distance from Center','FontSize',14)
lighting gouraud

figure;
hold on;
%Leek1999 data from Figure 7.3.1.i (p. 129) of thesis
%(http://europepmc.org/abstract/eth/302449)
lt = [0, 2, 5, 8, 10, 11, 12, 14, 15, 16, 17, 19, 21, 24];
ly = [130, 150, 210, 265,250, 300, 275, 310, 285, 275, 295, 300, 295, 280];
plot(lt,ly);

plot(x,Y(:,4*(1/xi_step + 1) + 1))
title('Radius (R)','FontSize',18)
xlabel('Time','FontSize',14)
ylabel('Tumour Radius','FontSize',14)
lighting gouraud
hold off;      