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
maxT = 500; % Solves on time interval [0,maxT]

global xi_step;
xi_step = 0.01;  % Space step

overall = [ 12;%sigmaM
         10;%sigmaL
         20;%D_c
         1;%k
         0.6;%cp
         0.2;%cc
         0.1;%p_m_max
         1;%p_am_max
         1;%p_al_max
         0.01;%d_a
         5.5;%d_m_max
         1;%d_l_max
         0.5;%d_cm_max
         0.5;%d_cl_max
         0.1;%d_cp_max
         1100;%chi_l
         0.65;%m_p
         0.001;%h_l
         90;%h_n
         0.01];%h_a
     
varNames = ['D_m','D_l','D_c','k','c_p','c_c','p_m^{max}','p_{am}^{max}',...
            'p_{al}^{max}','d_a','d_m^{max}','d_l^{max}','d_{cm}^{max}',...
            'd_{cl}^{max}','d_{cp}^{max}','\chi','m_p','h_l','h_n','h_a'];
figure;
hold on;

for i = 1:20
    
    for j = 1:3
        
        if j == 1
            vars = overall;
            vars(i) = 0.5*overall(i);
            lineSpec = ':';
        elseif j == 2
            vars(i) = overall(i);
            lineSpec = '-';
        elseif j == 3
            vars(i) = 1.5*overall(i);
            lineSpec = '--';
        end

        % Diffusion parameters
        global Sigma_m Sigma_l D_c; 

        Sigma_m = vars(1);
        Sigma_l = vars(2);
        D_c = vars(3);

        % Chemotaxis force 
        global k;
        k = vars(4);

        % Growth/death parameters
        global c_p p_m_max c_c p_al_max n_0 alpha p_am_max d_a d_m_max d_l_max d_cm_max d_cl_max d_cp_max chi_l m_p;
        c_p = vars(5);
        c_c = vars(6);
        n_0 = 0.2;
        alpha = 5;
        p_m_max = vars(7);
        p_am_max = vars(8);
        p_al_max = vars(9);
        d_a = vars(10);
        d_m_max = vars(11);
        d_l_max = vars(12);
        d_cm_max = vars(13);
        d_cl_max = vars(14);
        d_cp_max = vars(15);
        chi_l = vars(16);
        m_p = vars(17);

        global k_max a_char;
        k_max = 1; 
        a_char = 1;

        %Boundary condition parameters 
        global h_l h_n h_a l_inf n_inf c_inf a_inf; 
        h_l =  vars(18);
        h_n =  vars(19);
        h_a = vars(20);
        l_inf = 0.2;
        n_inf = 0.2; 
        c_inf = 1;
        a_inf = 0;

        % Initial radius 
        R_0 = 13; %13


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

        subplot(2,10,i)
        hold on;
        plot(x,Y(:,4*(1/xi_step + 1) + 1),lineSpec,'LineWidth',3)
        title(varNames(i))
        lighting gouraud
        ax1 = gca; % current axes
    end
end