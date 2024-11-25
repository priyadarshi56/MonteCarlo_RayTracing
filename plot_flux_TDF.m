%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Ploting electron flux and TDF as a function of E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; % close all % clear

LW1=2; % linewidth
LW2=1.5; % linewidth
MS=30; % markersize
FS=16; % fontsize
u1=1e3; 
if exist('Ec_element','var')==0
    Ec_avg=Ec;
end
% flux_zero=find(flux <= 0);
% E_values=E(flux_zero);
% E=E-E_values(end);
%--------------------------------------------------------------------------
figure
set(gcf,'units','centimeters','position',[16 2 33 8.5]);
ax1=axes('units','centimeters','position',[2.5 2 8 5.5]); hold(ax1,'on'); box(ax1,'on');
ax2=axes('units','centimeters','position',[13.5 2 8 5.5]); hold(ax2,'on'); box(ax2,'on');
ax3=axes('units','centimeters','position',[24 2 8 5.5]); hold(ax3,'on'); box(ax3,'on');
%----- Flux ---------------------------------------------------------------
plot(ax1,E*u1,flux,'linestyle','-','linewidth',LW2,'marker','none','markersize',MS);
set(ax1,'linewidth',LW1,'fontsize',FS);
xlabel(ax1,'E [meV]'); 
ylabel(ax1,'Flux [s^{-1}]');
xlim(ax1,[Ec_avg max(E)]*u1);
%----- DOS ----------------------------------------------------------------
DOS_exponent=round(mean(floor(log10(DOS_3D)))); % extracting exponent part
plot(ax2,E*u1,real(DOS_3D),'linestyle','-','linewidth',LW1,'marker','none','markersize',MS,'color','b');
set(ax2,'linewidth',LW1,'fontsize',FS);
xlabel(ax2,'E [meV]'); 
ylabel(ax2,'DOS [ ]');
xlim(ax2,[Ec_avg max(E)]*u1);
%----- TDF ----------------------------------------------------------------
plot(ax3,E*u1,TDF_MC,'linestyle','-','linewidth',LW2,'marker','none','markersize',MS);
set(ax3,'linewidth',LW1,'fontsize',FS);
xlabel(ax3,'E [meV]'); 
ylabel(ax3,'\Xi [a. u.]');
xlim(ax3,[Ec_avg max(E)]*u1);

%----- Shift TDF...? ------------------------------------------------------
% index_nonzero=find(TDF_MC>0);
% index_zero=1:index_nonzero(1)-1;
% zero_TDF_MC=TDF_MC(index_zero);
% nonzero_TDF_MC=TDF_MC(index_nonzero);
% TDF_MC_shifted=[0 nonzero_TDF_MC];
% shifted_E=linspace(0,max(E),length(index_nonzero)+1);
% plot(ax3,shifted_E*u1,TDF_MC_shifted,'linewidth',LW,'linestyle','-','marker','none','markersize',MS,'color','c');
%--------------------------------------------------------------------------

if exist('porosity','var')
    legend(ax3,['p=',num2str(porosity,2),'%, E_{redox}=',num2str(E_redox),'eV',''],'fontsize',10); legend boxoff
    sgtitle(['p=',num2str(porosity,2),'%, E_{redox}=',num2str(E_redox),'eV',''],'fontsize',FS/2);
else
    legend(ax3,['pristine at Ec=',num2str(Ec),'eV',''],'fontsize',10); legend boxoff
    sgtitle(['pristine at Ec=',num2str(Ec),'eV',''],'fontsize',FS/2);
end

