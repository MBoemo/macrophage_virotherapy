function dYdx = MOL_deriv_new_from_MOLvirus(t,Y)

% Finite difference paramters
global xi_step
% Switch paramters
global mac_turn_on;
% Radiation parameters
global rad_strength rad_start rad_dur gamma;
% Virus parameters
global D_phi d_phi p_phi_max r_rep r_uptake r_phi k_phi phi_inf h_phi;
% Diffusion paramters
global Sigma_m Sigma_l D_c; 
% Chemotaxis force 
global k;
% Growth/death parameters
global d_a chi_l m_p c_c;
%Boundary condition parameters 
global h_l h_n h_a l_inf n_inf a_inf; 






% Offsets for ODE matrix
% Batting order: l, m, c, a, R
l = 0*(1/xi_step + 1);
m_i = 1*(1/xi_step + 1);
m_h = 2*(1/xi_step + 1);
c = 3*(1/xi_step + 1);
a = 4*(1/xi_step + 1);
phi = 5*(1/xi_step + 1);
R = 6*(1/xi_step + 1) + 1;

xi_value = linspace(0,1,1/xi_step + 1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Ghost Points         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temporarily set i to its value on the boundary
i = 1/xi_step + 1;

m_i_ghostPoint = Y(m_i + i - 1) + Y(R)*2*xi_step*(...
                    (k*Y(m_i + i)*h_n*(1 - Y(m_h + i) - Y(m_i + i) - Y(l + i) - n_inf))/ ...
                    ((1 - Y(m_i+i) - Y(m_h + i) - Y(l+i))*Sigma_m)...
               );
           
m_h_ghostPoint = Y(m_h + i - 1) + Y(R)*2*xi_step*(...
                    (k*Y(m_h + i)*h_n*(1 - Y(m_h + i) - Y(m_i + i) - Y(l + i) - n_inf))/ ...
                    ((1 - Y(m_i + i) - Y(m_h + i) - Y(l + i))*Sigma_m)...
               );           
           
a_ghostPoint = Y(a + i - 1) + Y(R)*2*xi_step*(...
                    h_a*(a_inf - Y(a + i))...
               );

l_ghostPoint = Y(l + i - 1) + Y(R)*2*xi_step*(...
                    ((k*h_l)/Sigma_l)*(l_inf - Y(l + i)) + ...
                    (chi_l/(Sigma_l*Y(R)))*Y(l + i)*((a_ghostPoint - Y(a + i - 1))/(2*xi_step)) + ...
                    ((Y(l + i)*Sigma_m)/(Sigma_l*( Y(m_i + i) + Y(m_h + i) )*Y(R)))*((m_i_ghostPoint + m_h_ghostPoint - Y(m_i + i - 1) - Y(m_h + i - 1))/(2*xi_step))...
               );
           
phi_ghostPoint = Y(phi + i - 1) + Y(R)*2*xi_step*(...
                    h_phi*(phi_inf - Y(phi + i))...
                 );


DRdt = (1/Y(R))*(...
            (Sigma_l/k)*(l_ghostPoint - Y(l + 1/xi_step))/(2*xi_step) + ...
            (Sigma_m/k)*(1 - 1/(Y(m_i + 1/xi_step + 1) + Y(m_h + 1/xi_step + 1)))*( (m_i_ghostPoint - Y(m_i + 1/xi_step))/(2*xi_step) + (m_h_ghostPoint - Y(m_h + 1/xi_step))/(2*xi_step) ) - ...
            (chi_l/k)*Y(l + 1/xi_step + 1)*(a_ghostPoint - Y(a + 1/xi_step))/(2*xi_step)...
       );

           
                   
           
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Xi = 0 BC (Zero Flux)    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;

% l (Macrophages)
dYdx(l + i, 1) = heaviside(t - mac_turn_on)*(...
                   (6/(xi_step^2*Y(R)^2))*(...
                       (Sigma_l/k)*(Y(l + i + 1) - Y(l + i)) - ...
                       chi_l*((1 - Y(l + 1))/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i)) - ...
                       (Y(l + i)/k)*(Sigma_m*(Y(m_i + i + 1) + Y(m_h + i + 1) - Y(m_i + i) - Y(m_h + i)) + Sigma_l*(Y(l + 1 + 1) - Y(l + i)))...
                   ) - ...
                   Y(l + i)*fun_d_l(Y(c + i))...
                 ) - heaviside(mac_turn_on - t)*Y(l + i);    

% m_i (Infected Tumour Cells)
dYdx(i + m_i, 1) = 0;%(6/(xi_step^2*Y(R)^2))*(...
                   % (Sigma_m/k)*(Y(m_i + i + 1) - Y(m_i + i)) + ...
                   % chi_l*(Y(m_i + i)/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i)) - ...
                   % (Y(l + i)/k)*(Sigma_m*(Y(m_i + i + 1) - Y(m_i + i)) + Sigma_l*(Y(l + i + 1) - Y(l + i)))...
                 %)  - Y(m_i + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) - k_phi*Y(m_i + i) + r_phi*Y(phi + i)*Y(m_h + i);
                 
             
% m_h (Heathly Tumour Cells)
dYdx(i + m_h, 1) = (6/(xi_step^2*Y(R)^2))*(...
                    (Sigma_m/k)*(Y(m_h + i + 1) - Y(m_h + i)) + ...
                    chi_l*(Y(m_h + i)/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i)) - ...
                    (Y(l + i)/k)*(Sigma_m*(Y(m_h + i + 1) - Y(m_h + i)) + Sigma_l*(Y(l + i + 1) - Y(l + i)))...
                 ) + Y(m_h + i)*fun_p_m_fix( Y(m_i + i) + Y(m_h + i), Y(c + i) ) - Y(m_h + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_h + i)*fun_k(Y(c+i))...
                 - r_phi*Y(phi + i)*Y(m_h + i);%...
                 %- heaviside(t - rad_start)*exp(-rad_dur*(t - rad_start))*Y(m_h + i)*rad_strength*fun_S( Y(m_h + i) + Y(m_i + i), m_p );
                              
             
             
         
%c (Oxygen)
dYdx(i + c,1) = (6/(xi_step^2*Y(R)^2))*(...
                    D_c*(Y(c + i + 1) - Y(c + i))...
                ) - fun_d_c( Y(m_h + i) + Y(m_i + i), Y(l + i), 1 - Y(m_i + i) - Y(m_h + i) - Y(l + i), Y(c + i));
         

% a (Chemoattractant)
dYdx(i + a, 1) = (6/(xi_step^2*Y(R)^2))*(...
                    (Y(a + i + 1) - Y(a + i))...
                 ) + d_a*(fun_p_a(Y(l + i),Y(m_i + i)+Y(m_h + i),Y(c + i)) - Y(a + i));%...
                 %+ gamma*heaviside(t - rad_start)*exp(-rad_dur*(t - rad_start))*(Y(m_i + i) + Y(m_h + i))*rad_strength*fun_S( Y(m_i + i) + Y(m_h + i), m_p );
    
% Phi (Virus)
dYdx(i + phi,1) = 0;%(6/(xi_step^2*Y(R)^2))*(...
                  %   D_phi*(Y(phi + i + 1) - Y(phi + i))...
                  %) - d_phi*Y(phi + i) - r_uptake*Y(phi+i)*Y(m_h+i) + p_phi_max*(1 - fun_S(Y(c + i),c_c))*Y(l + i) + r_rep*(Y(m_i + i)*fun_d_m(Y(c + i)) + Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) + k_phi*Y(m_i + i));


              
              
              
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Xi in (0,1) [Interior]   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:1/xi_step 
    
    
    % l (Macrophages)
    dYdx(i + l, 1) = heaviside(t - mac_turn_on)*(...
                         (xi_value(i)/Y(R))*DRdt*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)  + ...
                             (1/Y(R)^2)*(...
                                (2/xi_value(i))*(...
                                    (Sigma_l/k)*(Y(l + i + 1) - Y(l + i -1))/(2*xi_step) - chi_l*((1 - Y(l + i))/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i -1))/(2*xi_step) - (Y(l + i)/k)*(Sigma_m*(Y(m_i + i + 1) + Y(m_h + i + 1) - Y(m_i + i - 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step))...
                                ) + ...
                                (Sigma_l/k)*(Y(l + i + 1) - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2) - chi_l*(1/k)*((Y(i + l + 1) - Y(i + l - 1))/(2*xi_step) - 2*Y(l + i)*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step))*((Y(a + i + 1) - Y(a + i - 1))/(2*xi_step)) - ...
                                chi_l*((1 - Y(l + i))/k)*Y(l + i)*(Y(a + i + 1) - 2*Y(a + i) + Y(a + i -1))/(xi_step^2) - (1/k)*((Y(l + i + 1) - Y(l + i - 1))/(2*xi_step))*(Sigma_m*(Y(m_i + i + 1) + Y(m_h + i + 1) - Y(m_i + i - 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)) - ...
                                (Y(l + i)/k)*(Sigma_l*(Y(l + i + 1) - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2) + Sigma_m*(Y(m_i + i + 1) + Y(m_h + i + 1) - 2*Y(m_i + i) - 2*Y(m_h + i) + Y(m_i + i - 1) + Y(m_h + i - 1))/(xi_step^2))...
                             ) - Y(l + i)*fun_d_l(Y(c + i))...
                     );% - heaviside(mac_turn_on - t)*Y(l + i);
                    
    
                 
                 
    % m_i (Infected Tumour Cells)
    dYdx(i + m_i, 1) = 0;%(xi_value(i)/Y(R))*DRdt*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step) + ...
%                      (1/Y(R)^2)*(...
%                         (2/xi_value(i))*(...
%                             (Sigma_m/k)*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step) + chi_l*(Y(m_i + i)/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step) - (Y(m_i + i)/k)*(Sigma_m*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step))...
%                         ) + ...
%                         (Sigma_m/k)*(Y(m_i + i + 1) - 2*Y(m_i + i) + Y(m_i + i - 1))/(xi_step^2) + chi_l*(Y(a + i + 1) - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)*((Y(m_i + i)*Y(l + i))/k) + ...
%                         chi_l*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step)*((1/k)*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step)*Y(l + i) + (1/k)*Y(m_i + i)*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)) - ...
%                         (1/k)*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step)*(Sigma_m*(Y(m_i + i + 1) - Y(m_i + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)) - ...
%                         (Y(m_i + i)/k)*(Sigma_m*(Y(m_i + i + 1) - 2*Y(m_i + i) + Y(m_i + i - 1))/(xi_step^2) + Sigma_l*(Y(l + i + 1) - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2))...
%                      )  - Y(m_i + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) + r_phi*Y(phi + i)*Y(m_h + i)...
%                  - k_phi*Y(m_i + i);
                 
                 
    % m_h (Healthy Tumour Cells)
    dYdx(i + m_h, 1) = (xi_value(i)/Y(R))*DRdt*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step) + ...
                     (1/Y(R)^2)*(...
                        (2/xi_value(i))*(...
                            (Sigma_m/k)*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step) + chi_l*(Y(m_h + i)/k)*Y(l + i)*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step) - (Y(m_h + i)/k)*(Sigma_m*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step))...
                        ) + ...
                        (Sigma_m/k)*(Y(m_h + i + 1) - 2*Y(m_h + i) + Y(m_h + i - 1))/(xi_step^2) + chi_l*(Y(a + i + 1) - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)*((Y(m_h + i)*Y(l + i))/k) + ...
                        chi_l*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step)*((1/k)*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step)*Y(l + i) + (1/k)*Y(m_h + i)*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)) - ...
                        (1/k)*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step)*(Sigma_m*(Y(m_h + i + 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(Y(l + i + 1) - Y(l + i - 1))/(2*xi_step)) - ...
                        (Y(m_h + i)/k)*(Sigma_m*(Y(m_h + i + 1) - 2*Y(m_h + i) + Y(m_h + i - 1))/(xi_step^2) + Sigma_l*(Y(l + i + 1) - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2))...
                     ) + Y(m_h + i)*fun_p_m_fix( Y(m_h + i) + Y(m_i + i), Y(c + i) ) - Y(m_h + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_h + i)*fun_k(Y(c+i))...
                 - r_phi*Y(phi + i)*Y(m_h + i);%...
                 %- heaviside(t - rad_start)*exp(-rad_dur*(t -
                 %rad_start))*Y(m_h + i)*rad_strength*fun_S( Y(m_h + i) + Y(m_i + i), m_p );
                      
    
    
 
    %c (Oxygen)
    dYdx(i + c,1) = (xi_value(i)/Y(R))*DRdt*(Y(c + i + 1) - Y(c + i - 1))/(2*xi_step) + ...
                    (D_c/Y(R)^2)*(...
                        (2/xi_value(i))*(Y(c + i + 1) - Y(c + i - 1))/(2*xi_step) + ...
                        (Y(c + i + 1) - 2*Y(c + i) + Y(c + i - 1))/(xi_step^2)...
                    ) - fun_d_c( Y(m_i + i) + Y(m_h + i), Y(l + i),1 - Y(m_i + i) - Y(m_h + i) - Y(l + i), Y(c + i));
    
    
    
    % a (Chemoattractant)
    dYdx(i + a, 1) = (xi_value(i)/Y(R))*DRdt*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step) + ...
                    (1/Y(R)^2)*(...
                        (2/xi_value(i))*(Y(a + i + 1) - Y(a + i - 1))/(2*xi_step) + ...
                        (Y(a + i + 1) - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)...
                    ) + d_a*(fun_p_a( Y(l + i), Y(m_i + i) + Y(m_h + i) , Y(c + i) ) - Y(a + i));%...
                 %+ gamma*heaviside(t - rad_start)*exp(-rad_dur*(t - rad_start))*(Y(m_i + i) + Y(m_h + i))*rad_strength*fun_S( Y(m_i + i) + Y(m_h + i), m_p );
      
    
    % phi (Virus)
    dYdx(i + phi, 1) = 0;%(xi_value(i)/Y(R))*DRdt*(Y(phi + i + 1) - Y(phi + i - 1))/(2*xi_step) + ...
%                        (D_phi/Y(R)^2)*(...
%                             (2/xi_value(i))*(Y(phi + i + 1) - Y(phi + i - 1))/(2*xi_step) + ...
%                             (Y(phi + i + 1) - 2*Y(phi + i) + Y(phi + i - 1))/(xi_step^2)...
%                        ) - d_phi*Y(phi + i) - r_uptake*Y(phi+i)*Y(m_h+i) + p_phi_max*(1 - fun_S(Y(c + i),c_c))*Y(l + i) + r_rep*(Y(m_i + i)*fun_d_m(Y(c + i)) + Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) + k_phi*Y(m_i + i));
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Xi = 1 BC and R       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = i + 1;


% l
dYdx(i + l, 1) = heaviside(t - mac_turn_on)*(... 
                    (xi_value(i)/Y(R))*DRdt*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)  + ...
                         (1/Y(R)^2)*(...
                            (2/xi_value(i))*(...
                                (Sigma_l/k)*(l_ghostPoint - Y(l + i -1))/(2*xi_step) - chi_l*((1 - Y(l + i))/k)*Y(l + i)*(a_ghostPoint - Y(a + i - 1))/(2*xi_step) - (Y(l + i)/k)*(Sigma_m*(m_i_ghostPoint + m_h_ghostPoint- Y(m_i + i - 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step))...
                            ) + ...
                            (Sigma_l/k)*(l_ghostPoint - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2) - chi_l*(1/k)*((l_ghostPoint - Y(i + l - 1))/(2*xi_step) - 2*Y(l + i)*(l_ghostPoint - Y(l + i - 1))/(2*xi_step))*((a_ghostPoint - Y(a + i - 1))/(2*xi_step)) - ...
                            chi_l*((1 - Y(l + i))/k)*Y(l + i)*(a_ghostPoint - 2*Y(a + i) + Y(a + i -1))/(xi_step^2) - (1/k)*((l_ghostPoint - Y(l + i - 1))/(2*xi_step))*(Sigma_m*(m_i_ghostPoint + m_h_ghostPoint - Y(m_i + i - 1) - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)) - ...
                            (Y(l + i)/k)*(Sigma_l*(l_ghostPoint - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2) + Sigma_m*(m_i_ghostPoint + m_h_ghostPoint - 2*Y(m_i + i) - 2*Y(m_h + i) + Y(m_i + i - 1) + Y(m_h + i - 1))/(xi_step^2))...
                         ) - Y(l + i)*fun_d_l(Y(c + i))...
                 ) - heaviside(mac_turn_on - t)*Y(l + i);
                 
    
                 
% m_i
dYdx(i + m_i, 1) = 0;%(xi_value(i)/Y(R))*DRdt*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step) + ...
%                      (1/Y(R)^2)*(...
%                         (2/xi_value(i))*(...
%                             (Sigma_m/k)*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step) + chi_l*(Y(m_i + i)/k)*Y(l + i)*(a_ghostPoint - Y(a + i - 1))/(2*xi_step) - (Y(m_i + i)/k)*(Sigma_m*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step))...
%                         ) + ...
%                         (Sigma_m/k)*(m_i_ghostPoint - 2*Y(m_i + i) + Y(m_i + i - 1))/(xi_step^2) + chi_l*(a_ghostPoint - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)*((Y(m_i + i)*Y(l + i))/k) + ...
%                         chi_l*(a_ghostPoint - Y(a + i - 1))/(2*xi_step)*((1/k)*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step)*Y(l + i) + (1/k)*Y(m_i + i)*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)) - ...
%                         (1/k)*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step)*(Sigma_m*(m_i_ghostPoint - Y(m_i + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)) - ...
%                         (Y(m_i + i)/k)*(Sigma_m*(m_i_ghostPoint - 2*Y(m_i + i) + Y(m_i + i - 1))/(xi_step^2) + Sigma_l*(l_ghostPoint - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2))...
%                         )  - Y(m_i + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) + r_phi*Y(phi + i)*Y(m_h + i)...
%                  - k_phi*Y(m_i + i);
                    
                    
% m_h
dYdx(i + m_h, 1) = (xi_value(i)/Y(R))*DRdt*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step) + ...
                     (1/Y(R)^2)*(...
                        (2/xi_value(i))*(...
                            (Sigma_m/k)*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step) + chi_l*(Y(m_h + i)/k)*Y(l + i)*(a_ghostPoint - Y(a + i - 1))/(2*xi_step) - (Y(m_h + i)/k)*(Sigma_m*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step))...
                        ) + ...
                        (Sigma_m/k)*(m_h_ghostPoint - 2*Y(m_h + i) + Y(m_h + i - 1))/(xi_step^2) + chi_l*(a_ghostPoint - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)*((Y(m_h + i)*Y(l + i))/k) + ...
                        chi_l*(a_ghostPoint - Y(a + i - 1))/(2*xi_step)*((1/k)*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step)*Y(l + i) + (1/k)*Y(m_h + i)*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)) - ...
                        (1/k)*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step)*(Sigma_m*(m_h_ghostPoint - Y(m_h + i - 1))/(2*xi_step) + Sigma_l*(l_ghostPoint - Y(l + i - 1))/(2*xi_step)) - ...
                        (Y(m_h + i)/k)*(Sigma_m*(m_h_ghostPoint - 2*Y(m_h + i) + Y(m_h + i - 1))/(xi_step^2) + Sigma_l*(l_ghostPoint - 2*Y(l + i) + Y(l + i - 1))/(xi_step^2))...
                        )  + Y(m_h + i)*fun_p_m_fix(  Y(m_h + i) + Y(m_i + i) , Y(c + i) ) - Y(m_h + i)*fun_d_m(Y(c + i)) - Y(l + i)*Y(m_h + i)*fun_k(Y(c+i))...
                    - r_phi*Y(phi + i)*Y(m_h + i);%...
                   % - heaviside(t - rad_start)*exp(-rad_dur*(t -
                   % rad_start))*Y(m_h + i)*rad_strength*fun_S( Y(m_h + i) + Y(m_i + i), m_p );
                  
                    
    
        
 
%c
dYdx(i + c,1) = 1 - Y(c + i);
    
    
    
% a  
dYdx(i + a, 1) = (xi_value(i)/Y(R))*DRdt*(a_ghostPoint - Y(a + i - 1))/(2*xi_step) + ...
                    (1/Y(R)^2)*(...
                        (2/xi_value(i))*(a_ghostPoint - Y(a + i - 1))/(2*xi_step) + ...
                        (a_ghostPoint - 2*Y(a + i) + Y(a + i - 1))/(xi_step^2)...
                    ) + d_a*(fun_p_a(Y(l + i),Y(m_i + i)+Y(m_h + i),Y(c + i)) - Y(a + i));%...
                 %+ gamma*heaviside(t - rad_start)*exp(-rad_dur*(t - rad_start))*(Y(m_i + i) + Y(m_h + i))*rad_strength*fun_S( Y(m_i + i) + Y(m_h + i), m_p );

                
% phi (Virus)
dYdx(i + phi, 1) = 0;%(xi_value(i)/Y(R))*DRdt*(phi_ghostPoint - Y(phi + i - 1))/(2*xi_step) + ...
%                     (1/Y(R)^2)*(...
%                         (2/xi_value(i))*(phi_ghostPoint - Y(phi + i - 1))/(2*xi_step) + ...
%                         (phi_ghostPoint - 2*Y(phi + i) + Y(phi + i - 1))/(xi_step^2)...
%                    ) - d_phi*Y(phi + i) - r_uptake*Y(phi+i)*Y(m_h+i) + p_phi_max*(1 - fun_S(Y(c + i),c_c))*Y(l + i) + r_rep*(Y(m_i + i)*fun_d_m(Y(c + i)) + Y(l + i)*Y(m_i + i)*fun_k(Y(c+i)) + k_phi*Y(m_i + i));
 
                   
% R
dYdx(R,1) = DRdt;



end





