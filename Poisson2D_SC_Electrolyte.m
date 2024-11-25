%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduction band bending at semiconductor-electrolyte interface (2D)
%   assuming no charge variation in electrolyte region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('\nSolving Poisson for Ec in the domain...');
%----- Co-ordinate extraction for boundary value --------------------------
x_cord_in_pore=zeros(Npy,Npx);
y_cord_in_pore=zeros(Npy,Npx);
x_cord_in_both=zeros(Npy,Npx);
y_cord_in_both=zeros(Npy,Npx);
x_cord_in_dH=zeros(Npy,Npx);
y_cord_in_dH=zeros(Npy,Npx);
n_node1=0;
n_node_in_pore=1;
n_node_in_dH=1;
n_node_in_both=1;
for ky=1:Npy
    for kx=1:Npx
        n_node1=n_node1+1;
        if n_node_in_pore~=length(node_in_pore)
            if n_node1==node_in_pore(n_node_in_pore)
                x_cord_in_pore(ky,kx)=x_cord(ky,kx);
            end
            if n_node1==node_in_pore(n_node_in_pore)
                y_cord_in_pore(ky,kx)=y_cord(ky,kx);
                n_node_in_pore=n_node_in_pore+1;
            end
        end
        if n_node_in_dH~=length(node_in_dH)
            if n_node1==node_in_dH(n_node_in_dH)
                x_cord_in_dH(ky,kx)=x_cord(ky,kx);
            end
            if n_node1==node_in_dH(n_node_in_dH)
                y_cord_in_dH(ky,kx)=y_cord(ky,kx);
                n_node_in_dH=n_node_in_dH+1;
            end
        end
        if n_node_in_both~=length(node_in_both)
            if n_node1==node_in_both(n_node_in_both)
                x_cord_in_both(ky,kx)=x_cord(ky,kx);
            end
            if n_node1==node_in_both(n_node_in_both)
                y_cord_in_both(ky,kx)=y_cord(ky,kx);
                n_node_in_both=n_node_in_both+1;
            end
        end
    end
end
%----- variable definition ------------------------------------------------
Vth=kB/q*T;
Ef_equ=Ef;
Ec_local=0; % dummy CB edge in semiconductor region
E_poisson=linspace(Emin,Emax,5000); % Energy grid (eV) for Poisson iterations
dE=E_poisson(2)-E_poisson(1);
del_Ec_offset=calculate_Ec(ND_SC,m1,T,E_poisson,Ef_equ);
phi_bi=(Ef_equ-E_redox); % built in potential
Nc1=2*((m1*kB*T)/(2*pi*hbar^2))^(3/2); % effective density of CB states in semiconductor region (/m3)
n1=ND_SC; % (ND1-NA1)/2+((ND1-NA1)^2/2+ni1^2)^(1/2); % assuming n=ND
Efn1=Ec_local-Vth*log(Nc1/n1); % calculated Fermi level in SC region
c1=-q/(epsilon_r*epsilon0);
c2=-q/(epsilon_r_sol*epsilon0);
Ec0_matrix=ones(Npy,Npx); % initializing Ec
ND_matrix=ones(Npy,Npx); % initializing ND doping
Ef0_matrix=ones(Npy,Npx); % initializing local Ef 
Ef_equ_matrix=Ef_equ*ones(Npy,Npx); % Equilibrium Ef
n0_matrix=zeros(Npy,Npx); % initializing electron concentration (/m3)
n_cation_matrix=zeros(Npy,Npx);
n_anion_matrix=zeros(Npy,Npx);
rho0_matrix=zeros(Npy,Npx); % initializing complete charge density
%----- Initial energy and charge profile ----------------------------------
for ky=1:Npy
    for kx=1:Npx
        if y_cord(ky,kx)==y_cord_in_pore(ky,kx) && x_cord(ky,kx)==x_cord_in_pore(ky,kx) && kx~=1 % for pore region
            Ec0_matrix(ky,kx)=NaN;
            ND_matrix(ky,kx)=0;
            Ef0_matrix(ky,kx)=E_redox;
            n_cation_matrix(ky,kx)=n_cation;
            n_anion_matrix(ky,kx)=n_anion;
            n0_matrix(ky,kx)=n_cation_matrix(ky,kx)-n_anion_matrix(ky,kx); % mobile charge
            rho0_matrix(ky,kx)=q*(n_cation_matrix(ky,kx)-n_anion_matrix(ky,kx));
        else % for SC region
            Ec0_matrix(ky,kx)=Ec_local;
            ND_matrix(ky,kx)=ND_SC;
            Ef0_matrix(ky,kx)=Ef_equ;
            DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_poisson-Ec0_matrix(ky,kx)))/(pi*hbar) *(q^(3/2)); % DOS in eV
            DOS_3D=real(DOS_3D);
            f=1./(1+exp((E_poisson-Efn1)/Vth)); % Fermi function for electron density
            n0_E=DOS_3D.*f; % Electron density as function of energy (/m3)
            n0_matrix(ky,kx)=sum(n0_E*dE); % Electron density (integral) in /m3 at each (dLy, dLx) point
            rho0_matrix(ky,kx)=q*(ND_matrix(ky,kx)-n0_matrix(ky,kx));
        end
    end
end
%
%----- Iterative calculation ----------------------------------------------
phi_matrix=0*ones(Npy,Npx);
a=0; % NBC at x=0 side
b=0; % NBC at x=end side
c=0; % NBC at y=0 side
d=0; % NBC at y=end side

Ec_matrix=zeros(Npy,Npx);
n_matrix=zeros(Npy,Npx);
rho_matrix=zeros(Npy,Npx);
convergance=zeros();
n_iter=zeros();
count_itr=0;
error=1;
error1=1e-6;
fprintf('\nError to get= %g',error1);
fprintf('\nError=           ');
while error>error1
    count_itr=count_itr+1;
    phi_old_matrix=phi_matrix;
    %----- charge density calculation -------------------------------------
    for ky=1:Npy
        for kx=1:Npx
            if y_cord(ky,kx)==y_cord_in_pore(ky,kx) && x_cord(ky,kx)==x_cord_in_pore(ky,kx) && kx~=1 % for pore region after dH layer
                n_matrix(ky,kx)=0; % *exp(-(phi_bi-phi_old(k))/1); % exp(phi_old(k)/1);
                rho_matrix(ky,kx)=c2*n_matrix(ky,kx);
            elseif y_cord(ky,kx)==y_cord_in_dH(ky,kx) && x_cord(ky,kx)==x_cord_in_dH(ky,kx) && kx~=1 % for H-layer
                n_matrix(ky,kx)=0; % mobile charge
                rho_matrix(ky,kx)=0;
            else
                Ec_matrix(ky,kx)=(phi_bi-phi_old_matrix(ky,kx)+del_Ec_offset);
                DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_poisson-Ec_matrix(ky,kx)))/(pi*hbar) *(q^(3/2)); % DOS in eV
                DOS_3D=real(DOS_3D);
                f=1./(1+exp((E_poisson-Ef_equ)/Vth)); % Fermi function for electron density
                n_E=DOS_3D.*f; % Electron density as function of energy (/m3)
                n_matrix(ky,kx)=sum(n_E*dE); % Electron density (integral) in /m3 at each (dLy, dLx) point
                rho_matrix(ky,kx)=c1*(ND_matrix(ky,kx)-n_matrix(ky,kx));           
            end
        end
    end
    %----- Poisson's solver -----------------------------------------------
    for ky=1:Npy
        for kx=1:Npx
            if ky==1 && kx==1
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky+1,kx)-2*c*dLy) + dLy^2*(2*phi_matrix(ky,kx+1)-2*a*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % 1st point at bottom boundary
            elseif ky==1 && kx==Npx
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky+1,kx)-2*c*dLy) + dLy^2*(2*phi_matrix(ky,kx-1)+2*b*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % last point at bottom boundary
            elseif ky==1 && kx>=2 && kx<Npx
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky+1,kx)-2*c*dLy) + dLy^2*(phi_matrix(ky,kx-1)+phi_matrix(ky,kx+1)) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % in between at bottom boundary
            elseif ky==Npy && kx==1
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky-1,kx)+2*d*dLy) + dLy^2*(2*phi_matrix(ky,kx+1)-2*a*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % 1st point at top boundary
            elseif ky==Npy && kx==Npx
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky-1,kx)+2*d*dLy) + dLy^2*(2*phi_matrix(ky,kx-1)+2*b*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % last point at top boundary
            elseif ky==Npy && kx>=2 && kx<Npx
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(2*phi_matrix(ky-1,kx)+2*d*dLy) + dLy^2*(phi_matrix(ky,kx-1)+phi_matrix(ky,kx+1)) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % in between at top boundary
            elseif ky>=2 && ky<Npy && kx==1
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(phi_matrix(ky-1,kx)+phi_matrix(ky+1,kx)) + dLy^2*(2*phi_matrix(ky,kx+1)-2*a*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % left boundary excluding 1st and last point
            elseif ky>=2 && ky<Npy && kx==Npx
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(phi_matrix(ky-1,kx)+phi_matrix(ky+1,kx)) + dLy^2*(2*phi_matrix(ky,kx-1)+2*b*dLx) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % right boundary excluding 1st and last point
            elseif y_cord(ky,kx)==y_cord_in_pore(ky,kx) && x_cord(ky,kx)==x_cord_in_pore(ky,kx) && kx~=1 % for pore region after dH layer
                phi_matrix(ky,kx)=0;
            else
                phi_matrix(ky,kx)=(1/(2*(dLx^2+dLy^2)))*(dLx^2*(phi_matrix(ky-1,kx)+phi_matrix(ky+1,kx)) + dLy^2*(phi_matrix(ky,kx-1)+phi_matrix(ky,kx+1)) - (dLx^2*dLy^2)*rho_matrix(ky,kx)); % rest other points
            end
        end
    end
    error=max(abs(phi_matrix(:)-phi_old_matrix(:))./abs(phi_old_matrix(:)));
    convergance(count_itr)=error;
    n_iter(count_itr)=count_itr;
    % fprintf('\b\b\b\b');
    % fprintf('%g', error); 
%     fprintf('Error= %g \n',error);

    fprintf(repmat('\b', 1, numel(num2str(error, '%.4e'))));
    fprintf('%.4e', error);
    drawnow;
end
[Ec_element]=matrix_to_element(element,node,Npx,Ec_matrix);
poisson_time=toc;
fprintf('   DONE!');
