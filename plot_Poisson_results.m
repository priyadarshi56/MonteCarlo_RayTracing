%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plot Poisson results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear; clc; close all
plot_geometry='Yes';
plot_Poisson_mesh='No';
plot_initial_condition='No';
plot_Ec='Yes';
plot_n='No';
plot_phi='No';
plot_cutlines='Yes';
LW=2;
FS=14;
MS=2.5;
Ly_cut_line=49e-9; % round(rand*Ly*u)/u;
% Ly_cut_line=round(rand*Ly*u)/u; % 35e-9;
Lx_array=0:dLx:Lx; Ly_array=0:dLy:Ly;
cut_line_index=find(Ly_array>Ly_cut_line,1);

% Ec_cut=-0.2;
% Ec_matrix(Ec_matrix < Ec_cut)=Ec_cut; % replace values in Ec_matrix less than Ec_well_max with Ec_well_max
% [Ec_element]=matrix_to_element(element,node,Npx,Ec_matrix);


if strcmpi(plot_geometry,'Yes')
%----- Geometry -----------------------------------------------------------
    figure
    set(gcf,'units','centimeters','position',[1 12 17 10.2]);
    axes('units','centimeters','position',[2 1.75 14.5 7.5]);    
    hold on; box on; u=1e9; uu=1e3;
    plot(polyshape([0 Lx Lx 0]*u,[0 0 Ly Ly]*u),'facecolor',[0.5 0.5 0.5],'FaceAlpha',0.25); % 0.19 0.55 0.91
    b1=plot([0 0]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
    b2=plot([Lx Lx]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
    b3=plot([0 Lx]*u,[0 0]*u,'-k','linewidth',4); % boundary line in Y-direction
    b4=plot([0 Lx]*u,[Ly Ly]*u,'-k','linewidth',4); % boundary line in Y-direction
    for k=1:n_pore
        poly_pore=polyshape(x_pore(k,:)*u,y_pore(k,:)*u);
        plot(poly_pore,'linestyle','none','linewidth',0.01,'facecolor',[0.4 0.6 1],'facealpha',0.4,'edgecolor',[0.4 0.6 1]);
    end
    set(gca,'linewidth',1.5,'fontsize',FS,'xtick',0:40:1e3,'ytick',0:25:1e3);
    ax=gca;
    ax.XAxis.TickLength=[0 0]; ax.YAxis.TickLength=[0 0]; % Set x- & y-axis tick length
    xlabel('length [nm]');
    ylabel('width [nm]');
    gap=2e-9;
    xlim([0-2*gap Lx+2*gap]*u);
    ylim([0-2*gap Ly+2*gap]*u);
    yline(Ly_array(cut_line_index)*u,'linestyle','-.','linewidth',1.5,'color','k');
    sgtitle(['cut at channel width=',num2str(Ly_array(cut_line_index)*u),'nm',''],'fontsize',FS);
end
if strcmpi(plot_Poisson_mesh,'Yes')
%----- Geometry with meshing details --------------------------------------
    figure
    set(gcf,'units','centimeters','position',[1 12 17 10.2]);
    axes('units','centimeters','position',[2 1.75 14.5 7.5]);    
    hold on; box on; u=1e9; uu=1e3;
    plot(x_cord*u,y_cord*u,'.','markersize',MS/2,'color','k');
%     plot(x_cord_in_pore*u,y_cord_in_pore*u,'.','markersize',4,'color','y');
    plot(x_cord_in_dH*u,y_cord_in_dH*u,'.','markersize',MS,'color',[0.1 0.4 0.7]);
    plot(x_cord_in_both*u,y_cord_in_both*u,'.','markersize',MS,'color',[0.4 0.6 1]);
    set(gca,'linewidth',1.5,'fontsize',FS,'xtick',0:40:1e3,'ytick',0:25:1e3);
    ax=gca;
    ax.XAxis.TickLength=[0 0]; ax.YAxis.TickLength=[0 0]; % Set x- & y-axis tick length
    xlabel('length [nm]');
    ylabel('width [nm]');
    gap=2e-9;
    xlim([0-2*gap Lx+2*gap]*u);
    ylim([0-2*gap Ly+2*gap]*u);

    b1=plot([0 0]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
    b2=plot([Lx Lx]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
    b3=plot([0 Lx]*u,[0 0]*u,'-k','linewidth',4); % boundary line in Y-direction
    b4=plot([0 Lx]*u,[Ly Ly]*u,'-k','linewidth',4); % boundary line in Y-direction
    yline(Ly_array(cut_line_index)*u,'linestyle','-.','linewidth',1.5,'color','k');
    sgtitle(['cut at channel width=',num2str(Ly_array(cut_line_index)*u),'nm',''],'fontsize',FS);
end
if strcmpi(plot_initial_condition,'Yes')
    figure
    set(gcf,'units','centimeters','position',[15 12 23 8]);
    axes('units','centimeters','position',[2.5 2 8 5.5]);
    hold on; box on; u=1e9; uu=1e3;
    surf(x_cord*u,y_cord*u,Ef0_matrix*uu,'linestyle','none','facealpha',0.5); % 'facecolor','b',
    colorbar;
    view(-33,30);
    set(gca,'linewidth',0.2,'fontsize',FS/1.5,'xtick',0:40:1e3,'ytick',0:40:1e3);
    xlabel('channel length [nm]'); 
    ylabel('channel width [nm]');
    zlabel('E_f/E_{redox} [meV]');
    xlim([0 Lx]*u);
    ylim([0 Ly]*u);
    zlim([min([min(Ef0_matrix) min(Ef_equ_matrix)])-0.1 max([max(Ef0_matrix) max(Ef_equ_matrix)])+0.05]*uu);
    
    axes('units','centimeters','position',[14 2 8 5.5]);
    hold on; box on; u=1e9;
    surf(x_cord*u,y_cord*u,ND_matrix,'linestyle','none','facealpha',0.5); % 'facecolor','b',
    colorbar;
    view(-33,30);
    set(gca,'linewidth',0.2,'fontsize',FS/1.5,'xtick',0:40:1e3,'ytick',0:40:1e3);
    xlabel('channel length [nm]'); 
    ylabel('channel width [nm]');
    zlabel('charge density [m^{-3}]');
    xlim([0 Lx]*u);
    ylim([0 Ly]*u);
    zlim([min([min(ND_matrix) min(n0_matrix)])/50 max([max(ND_matrix) max(n0_matrix)])*2]);
    legend('N_D','n_{mobile}','location','best'); legend boxoff
end
% ----- Solved Ec profile -------------------------------------------------
Ec_matrix1=Ec_matrix; 
Ef0_matrix1=Ef0_matrix;
for ky=1:Npy
    for kx=1:Npx
        if y_cord(ky,kx)==y_cord_in_pore(ky,kx) && x_cord(ky,kx)==x_cord_in_pore(ky,kx) && kx~=1 % for pore region after dH layer
            Ec_matrix1(ky,kx)=NaN;
            Ef0_matrix1(ky,kx)=NaN;
            Ef0_matrix1(ky,kx)=NaN;
        end
    end
end
if strcmpi(plot_Ec,'Yes')
    figure
    set(gcf,'units','centimeters','position',[1 12 22 15]);
    hold on; box on; u=1e9; uu=1e3;
%     surf(x_cord*u,y_cord*u,Ef0_matrix1*uu,'linestyle','none','facealpha',0.85,'facecolor',[0.5 0.5 0.5]);
    surf(x_cord*u,y_cord*u,Ec_matrix1*uu,'linestyle','none','facealpha',0.75);
    colorbar;
    colormap hsv;
    colormap(flipud(hsv)); % reverse the color
    set(gca,'linewidth',1,'fontsize',16,'xtick',0:40:1e3,'ytick',0:25:1e3);
    xlabel_handle=xlabel('channel length [nm]'); 
    ylabel_handle=ylabel('channel width [nm]');
    zlabel('E_c [meV]');
    xlim([0 Lx]*u);
    ylim([0 Ly]*u);
%     zlim([-300 200]);
%     view(-50,10);
    view(-36.9142827428624, 27.8950522613367);
    [az,el]=view; % view angle (azimuth and elevation)
    rotation_angle=-az;
    set(xlabel_handle,'rotation',20);

%     zlim([min([min(Ec_matrix) min(Ef0_matrix)])-0.1 max([max(Ec_matrix) max(Ef0_matrix)])+0.05]*uu);
%     legend('E_{f,equ}','E_c','location','best'); legend boxoff
%     for kn=1:length(En)
%         surf(x_cord*u,y_cord*u,(En_aboveEc(kn)*ones(Npy,Npx))*uu,'linestyle','none','facealpha',0.85,'facecolor',[0.5 0.5 0.5]);
%     end
end

%----- Solved charge profile ----------------------------------------------
ND_matrix1=ND_matrix;
pore_matrix1=zeros(Npy,Npx);
n_matrix1=n_matrix;
for ky=1:Npy
    for kx=1:Npx
        if y_cord(ky,kx)==y_cord_in_pore(ky,kx) && x_cord(ky,kx)==x_cord_in_pore(ky,kx) && kx~=1 % for pore region after dH layer
            ND_matrix1(ky,kx)=NaN;
            n_matrix1(ky,kx)=NaN;
            pore_matrix1(ky,kx)=ND_SC;
        end
        if y_cord(ky,kx)==y_cord_in_both(ky,kx) && x_cord(ky,kx)==x_cord_in_both(ky,kx) && kx~=1 % for pore region after dH layer
        
        end
    end
end
if strcmpi(plot_n,'Yes')
    figure
    set(gcf,'units','centimeters','position',[23 12 22 15]);
    hold on; box on
    surf(x_cord*u,y_cord*u,ND_matrix1,'linestyle','none','facealpha',0.4,'facecolor','k');
    surf(x_cord_in_pore*u,y_cord_in_pore*u,pore_matrix1,'linestyle','none','facecolor','b','facealpha',0.4);
    surf(x_cord*u,y_cord*u,n_matrix,'linestyle',':','facealpha',0.5);
    colorbar;
    colormap turbo;
    set(gca,'linewidth',1,'fontsize',FS,'xtick',0:40:1e3,'ytick',0:25:1e3,'ztick',10.^(14:30),'zscale','log');
    
    xlabel_handle=xlabel('channel length [nm]'); 
    ylabel_handle=ylabel('channel width [nm]');
    zlabel('charge density [m^{-3}]');
    xlim([0 Lx]*u);
    ylim([0 Ly]*u);
%     view(-50,10);
    view(-36.9142827428624, 27.8950522613367);
    [az,el]=view; % view angle (azimuth and elevation)
    rotation_angle=-az;
    set(xlabel_handle,'rotation',20);
    
    
   
    zlim([min([min(ND_matrix) min(n_matrix)])/100 max([max(ND_matrix) max(n_matrix)])*1.2]);
    legend('n_{init}','','n_{new}','location','best'); legend boxoff
end
% ----- Poisson potential -------------------------------------------------
figure
set(gcf,'units','centimeters','position',[23 1 22 15]);
hold on; box on; u=1e9; uu=1e3;
surf(x_cord*u,y_cord*u,phi_matrix*uu,'linestyle','none','facealpha',0.75);
colorbar;
colormap hsv;
view(-25,30);
set(gca,'linewidth',1,'fontsize',FS,'xtick',0:40:1e3,'ytick',0:20:1e3);
xlabel('channel length [nm]'); 
ylabel('channel width [nm]');
zlabel('Poisson potential [mV]');
xlim([0 Lx]*u);
ylim([0 Ly]*u);
zlim([min(phi_matrix(:))-0.001 max(phi_matrix(:))+0.001]*uu);
legend('volt','orientation','vertical','location','best','fontsize',16); legend boxoff
close

%----- 1D cutline plots ---------------------------------------------------
if strcmpi(plot_cutlines,'Yes')
    Lx_array=0:dLx:Lx;
    Ly_array=0:dLy:Ly;
    figure
    set(gcf,'units','centimeters','position',[16 2 33 8.5]);
%----- Poisson potential profile ------------------------------------------
    axes('units','centimeters','position',[2.5 2 8 5.5]);
    hold on; box on; u=1e9; uu=1e3;
    plot(Lx_array*u,phi_matrix(cut_line_index,:)*uu,'-b','linewidth',2);
    set(gca,'linewidth',LW,'fontsize',FS,'xtick',0:40:1e3);
    xlabel('channel length [nm]'); ylabel('\phi [mV]');
    xlim([0 Lx]*u); 
    ylim([min(phi_matrix(cut_line_index,:))-0.025 max(phi_matrix(cut_line_index,:))+0.05]*uu);
%----- Solved Ec profile --------------------------------------------------   
    ax2=axes('units','centimeters','position',[13 2 8 5.5]);
    hold on; box on; u=1e9; uu=1e3;
    plot(Lx_array*u,Ec_matrix1(cut_line_index,:)*uu,'-b','linewidth',2);  
%     Ec_poisson_cutline=Ec_matrix(cut_line_index,:);
%     Ec_poisson_cutline(isnan(Ec_poisson_cutline))=0;
%     avg_Ec_poisson_cutline=mean(Ec_poisson_cutline);
%     yline(avg_Ec_poisson_cutline*uu,'linestyle','--','linewidth',2,'color','r');
    plot(Lx_array*u,Ef_equ*ones(1,Npx)*uu,'--k','linewidth',LW);
    set(gca,'linewidth',LW,'fontsize',FS,'xtick',0:40:1e3); % ,'ytick',-1e4:100:1e4);
    xlabel('channel length [nm]'); ylabel('E_c [meV]');
    xlim([0 Lx]*u); 
%     ylim([min([min(Ec_matrix(cut_line_index,:)) Ef_equ])-0.025 max([Ec_matrix(cut_line_index,:) Ef_equ])+0.05]*uu);

% % % %     Test_E_Triangular_Well;
% % % %     yline(En_aboveEc*uu);

%     legend('E_c','E_f','location','best','NumColumns',2,'orientation','horizontal'); legend boxoff
%----- Solved charge density profile --------------------------------------    
    axes('units','centimeters','position',[24 2 8 5.5]);
    hold on; box on; u=1e9;
    plot(Lx_array*u,n_matrix(cut_line_index,:),'-b','linewidth',LW);
    plot(Lx_array*u,ND_matrix(cut_line_index,:),'-k','linewidth',LW);
%     yline(avg_n_poisson,'linestyle','--','linewidth',2,'color','r');
    set(gca,'linewidth',LW,'fontsize',FS,'xtick',0:40:1e3,'ytick',10.^(15:30),'yscale','log');
    xlabel('channel length [nm]'); ylabel('charge density [m^{-3}]');
    xlim([0 Lx]*u); 
    ylim([min([ND_SC min(n_matrix1(cut_line_index,:))])/2 max([ND_SC max(n_matrix1(cut_line_index,:))])*50]);
%     legend('n_{new}','n_{init}','location','best','NumColumns',2,'orientation','horizontal'); legend boxoff   

    sgtitle(['cut at channel width=',num2str(Ly_array(cut_line_index)*u),'nm,' ...
        ' p=',num2str(porosity,2),'% E_{redox}=',num2str(E_redox),'eV',''],'fontsize',FS/1.5);
end
%----- Convergance plot ---------------------------------------------------
figure
set(gcf,'units','centimeters','position',[12 2 16 11]);
plot(n_iter,convergance,'-','linewidth',2);
set(gca,'fontsize',16,'xscale','log','yscale','log'); 
xlabel('# iterations'); 
ylabel('Convergance');
close

% saveas(gcf,'name.fig');
% saveas(gcf,'name.svg');
