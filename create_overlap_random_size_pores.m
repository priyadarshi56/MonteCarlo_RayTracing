function [C,r,x_pore,y_pore] = create_overlap_random_size_pores(Lx,Ly,n_pore,r_min,r_max)
%     Lx=100;
%     Ly=50;
%     n_pore=30;
%     r_min=2;
%     r_max=5;
    r=zeros(n_pore,1);
    C=zeros(n_pore,2);
    % C(1,:)=[r+rand*Lx, r+rand*Ly]; % First center 
    for k=1:n_pore
        r(k)=r_min+rand*(r_max-r_min);
        C_temp=[rand*(Lx-1*r(k)), rand*(Ly-1*r(k))]; % Set the center inside domain
        while any(sum((C(1:k-1,:)-C_temp).^2,2)<(0.5*r(k))^2)
                C_temp=[rand*(Lx-1*r(k)), rand*(Ly-1*r(k))];
        end
        C(k,:)=C_temp;
    end
    global nc;
    theta=linspace(0,360,nc);
    x_pore=C(:,1)+r*cosd(theta);
    y_pore=C(:,2)+r*sind(theta);
end
%% ----- Pores Plot --------------------------------------------------------
% figure(14)
% clf
% hold on
% for k=1:n_pore
%     plot(C(k,1),C(k,2),'.g');
%     plot(x_pore(k,:),y_pore(k,:),'r');
%     text(C(k,1),C(k,2),num2str(k),'color',[0 0 0],'fontsize',8,'horizontalAlignment', 'center');        
% end
% rectangle('Position',[0 0 Lx Ly],'linewidth',2,'edgecolor','k');
% % axis([0 Lx*1e9 0 Ly*1e9])
% box on


