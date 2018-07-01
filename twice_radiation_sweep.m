clc;clear;close;
figure; hold all;
global mac_turn_on rad_start2;

mac_turn_on = 90;

for rad_start2 = mac_turn_on+1
   
       
%         str_legend = strcat('Second Radiation at ',int2str(rad_start2));
     
        
        MOL_virus_twiceRad;
        
%         plot(t/10,Y(:,6*(1/xi_step + 1) + 1),'DisplayName',str_legend)
%         axis square;
%         title('Radius (R)','FontSize',18)
%         xlabel('Time','FontSize',14)
%         ylabel('Tumour Radius','FontSize',14)
%         lighting gouraud
        
    
    
    
end

hold off;

