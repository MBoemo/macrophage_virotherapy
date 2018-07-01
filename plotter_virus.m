
R = Y(:,end);


for i = 1:length(R)
    pcolor_x(i,:) = 0:R(i)/(1/xi_step):R(i);
    pcolor_y(i,1:1/xi_step+1) = x(i);
    
end




%%%%%%%%%%%%%
figure('DefaultAxesFontSize',20)
subplot(3,2,1)
h = pcolor(pcolor_y,pcolor_x,Y(:,1:1/xi_step + 1))
set(h,'EdgeColor','none')
title('Macrophages (l)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar

subplot(3,2,4)
h = pcolor(pcolor_y,pcolor_x,Y(:,1/xi_step + 2:2*(1/xi_step + 1)))
set(h,'EdgeColor','none')
title('Infected Tumour Cells (m_i)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar

subplot(3,2,3)
h = pcolor(pcolor_y,pcolor_x,Y(:,2*(1/xi_step + 1) + 1:3*(1/xi_step + 1)))
set(h,'EdgeColor','none')
title('Uninfected Tumour Cells (m_h)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar

subplot(3,2,2)
h = pcolor(pcolor_y,pcolor_x,Y(:,5*(1/xi_step + 1) + 1:6*(1/xi_step + 1)))
set(h,'EdgeColor','none')
title('Virus (phi)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar

subplot(3,2,5)
h = pcolor(pcolor_y,pcolor_x,Y(:,4*(1/xi_step + 1) + 1:5*(1/xi_step + 1)))
set(h,'EdgeColor','none')
title('Chemoattractant (a)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar

subplot(3,2,6)
h = pcolor(pcolor_y,pcolor_x,Y(:,3*(1/xi_step + 1) + 1:4*(1/xi_step + 1)))
set(h,'EdgeColor','none')
title('Oxygen (c)','FontSize',20)
xlabel('Time','FontSize',20)
ylabel('Distance from Center','FontSize',20)
lighting gouraud
colorbar


colormap jet
