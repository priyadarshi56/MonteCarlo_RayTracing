function [C,x_pore_rotate,y_pore_rotate] = create_nonoverlap_random_oval_pores(Lx,Ly,n_pore,r)
%     Lx=100e-9;
%     Ly=50e-9;
%     n_pore=30;
%     r=5e-9;
    global nc u;
    theta=linspace(0,360,nc);
    C=zeros(n_pore,2);
    x_pore=zeros(n_pore,nc);
    y_pore=zeros(n_pore,nc);
    x_pore_rotate=zeros(n_pore,nc);
    y_pore_rotate=zeros(n_pore,nc);
    C(1,:)=[r+rand*Lx, r+rand*Ly]; % First center 
    for k=1:n_pore
        C_temp=[rand*(Lx-1*r), rand*(Ly-1*r)];
        while any(sum((C(1:k-1,:)-C_temp).^2,2)<(2*r)^2)
                C_temp=[rand*(Lx-r), rand*(Ly-r)];
        end
        C(k,:)=C_temp;       
        r_major=randi([1 5],1,1)./u; % Major axis
%         e=0.9; % Eccentricity
        r_minor=randi([1 2],1,1)./u; % r_major*sqrt(1-e^2); % Minor axis
        x_pore(k,:)=C(k,1)+r_major*cosd(theta);
        y_pore(k,:)=C(k,2)+r_minor*sind(theta);       
%------ Rotation of axis --------------------------------------------------
        axis_rotate=2*pi*rand; % Random angle in radian
        axis_rotate=axis_rotate*180/pi; % Converting theta in degree from radian
        transformMatrix = [cosd(axis_rotate), sind(axis_rotate); -sind(axis_rotate), cosd(axis_rotate)];
        xAligned = (x_pore(k,:) - C(k,1));
        yAligned = (y_pore(k,:) - C(k,2));
        xyAligned = [xAligned; yAligned]';
        xyRotated = xyAligned * transformMatrix;
        x_pore_rotate(k,:) = xyRotated(:, 1) + C(k,1);
        y_pore_rotate(k,:) = xyRotated(:, 2) + C(k,2);
    end   
end
%% ----- Pores Plot -------------------------------------------------------
% figure(14)
% clf
% subplot(1,2,1)
% hold on
% for k=1:n_pore
%     plot(C(k,1),C(k,2),'.g');
%     plot(x_pore(k,:),y_pore(k,:),'b');
%     text(C(k,1),C(k,2),num2str(k),'color',[0 0 0],'fontsize',8,'horizontalAlignment', 'center');        
% end
% box on
% subplot(1,2,2)
% hold on
% for k=1:n_pore
%     plot(C(k,1),C(k,2),'.g');
%     plot(x_pore_rotate(k,:),y_pore_rotate(k,:),'r');
%     text(C(k,1),C(k,2),num2str(k),'color',[0 0 0],'fontsize',8,'horizontalAlignment', 'center');        
% end
% % rectangle('Position',[0 0 Lx Ly],'linewidth',2,'edgecolor','k');
% % axis([0 Lx*1e9 0 Ly*1e9])
% box on
