function [C,x_pore,y_pore] = create_nonoverlap_random_pores(Lx,Ly,n_pore,pore_r)
%     Lx=100;
%     Ly=50;
%     n_pore=30;
%     r=2;
    C=zeros(n_pore,2);
    rad=zeros(1,n_pore);
    C(1,:)=[pore_r+rand*Lx, pore_r+rand*Ly]; % First center 
    for k=1:n_pore
        C_temp=[rand*(Lx-1*pore_r), rand*(Ly-1*pore_r)];
        while any(sum((C(1:k-1,:)-C_temp).^2,2)<(2*pore_r)^2)
                C_temp=[rand*(Lx-1*pore_r), rand*(Ly-1*pore_r)];
        end
        C(k,:)=C_temp;
        rad(k)=pore_r;
    end
    global nc;
    theta=linspace(0,360,nc);
    x_pore=C(:,1)+pore_r*cosd(theta);
    y_pore=C(:,2)+pore_r*sind(theta);
end
%% ----- Pores Plot -------------------------------------------------------
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
