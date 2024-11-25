function [MF_avg]=MF_calculator(Lx,dLx,Ly,dLy)
%--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Analytical conductivity calculation --------------------------------
    constants;
    input_parameters;
    Vth=kB/q*T;
    Emin_local=0.0001; Emax_local=0.5; NE_local=100;
    Ef_min=-0.2; Ef_max=0.2; NEf=50;
    E_local=linspace(Emin_local,Emax_local,NE_local); % Energy (in eV)
    Ef_local=linspace(Ef_min,Ef_max,NEf); % Fermi energy (in eV)
    Ec=0;
    %----- Scattering rate calculation ------------------------------------
    cl=mass_density*vs^2; % elastic constant
    c_ap=(pi*kB*T*D_adp^2)/(hbar*cl) * q;
    DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_local-Ec))/(pi*hbar) * (q^(3/2)); % DOS in eV
    DOS_3D=real(DOS_3D);
    APS_rate=c_ap*DOS_3D; % 1/sec
    total_rate=APS_rate;
    total_time=1./total_rate;
    %----- TDF calculation-------------------------------------------------
    k=(sqrt(2*m1*q)/hbar)*sqrt(E_local-Ec); % parabolic relation of E-k (1/m)
    v_k=real((hbar*k)/m1); % speed (m/sec)
    v_E=real(sqrt(2*q*E_local/m1)); % speed (m/sec) when energy in eV
    mfp=v_E.*total_time; % mean free path (meter)
    TDF_analyt=total_time.*v_E.^2.*DOS_3D; % 1/(eV.sec.meter)
    %----- Transport coefficients -----------------------------------------
    dE=E_local(2)-E_local(1);
    sigma_analyt=zeros(1,NEf);
    for kEf=1:NEf
        f=1./(1+exp((E_local-Ef_local(kEf))/Vth)); % Unitless
        dfdE=diff(f)./diff(E_local); % 1/Joule
        factor_TDF_analyt=TDF_analyt(1:end-1).*(-dfdE);
        sigma_analyt(kEf)=q^2*sum(factor_TDF_analyt)*dE/q; % Siemen/m    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Ray Tracing in 2D pristine -----------------------------------------
Nele_local=1000;
n_left_exit=zeros(1,NE_local);
n_right_exit=zeros(1,NE_local);
left_exit_ele=zeros(NE_local,Nele_local);
right_exit_ele=zeros(NE_local,Nele_local);
left_exit_time=zeros(NE_local,Nele_local);
right_exit_time=zeros(NE_local,Nele_local);
sim_prog=0;
fprintf('Looking for MF:  %3d%%\n',sim_prog);
for kE=1:NE_local
    % sim_prog=kE/NE*100; % Simulation progress in percentage
    % fprintf('\b\b\b\b%3.0f%%',sim_prog);
    mfp_E=mfp(kE);
    for kNele=1:Nele_local
        sim_prog=((kE-1)*Nele_local+kNele)/(NE_local*Nele_local)*100; % Simulation progress in percentage
%-------Initialization-----------------------------------------------------
        clearvars x_step y_step r_step t_step v_step x_mfp y_mfp
        E_ele=E_local(kE);
        dx=dLx;
        dy=dLy;
        dr=sqrt(dx^2+dy^2);
        x_covered=0;
        y_covered=rand*Ly;
        ele_exit_flag=false;
        ele_yRefl_flag=false;
        
        n_step=1; % number of step encountered
        x_step(n_step)=x_covered;
        y_step(n_step)=y_covered;
        r_step(n_step)=0;
        v_step(n_step)=0;
        t_step(n_step)=0;
        
        n_mfp=1; % number of mfp encountered
        x_mfp(n_mfp)=x_covered;
        y_mfp(n_mfp)=y_covered;
        E_mfp(n_mfp)=E_ele;
        E_ele=E_ele-Ec;
%----- Ray tracing starts -------------------------------------------------        
        while (x_covered<Lx && x_covered>=0)
            k=(sqrt(2*m1)/hbar)*sqrt(E_ele*q);
            theta=2*pi*rand; % random theta in radian
            theta=theta*180/pi; % converting theta in degree from radian
            kx=k*cosd(theta);
            ky=k*sind(theta);
            vx=(hbar*kx)/m1; % speed in x-direction
            vy=(hbar*ky)/m1; % speed in y-direction
            v=sqrt(vx^2+vy^2);
            % dr=sqrt(dx^2+dy^2);
            if (vy<=0)
                v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
            else
                v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
            end
            if (n_step==1 && vx<0)
                dx=-dr*cosd(v_theta);
                dy=dr*sind(v_theta);
                vx=-vx;
            else
                dx=dr*cosd(v_theta);
                dy=dr*sind(v_theta);
            end
            x_flight=0;
            y_flight=0;
            r_flight=0;
            while (r_flight<=mfp_E)
                n_step=n_step+1;
                x_flight=x_flight+dx;
                y_flight=y_flight+dy;
                r_flight=r_flight+sqrt(dx^2+dy^2);
                x_step(n_step)=x_covered+x_flight;
                y_step(n_step)=y_covered+y_flight;
                r_step(n_step)=sqrt((x_step(n_step)-x_step(n_step-1))^2+(y_step(n_step)-y_step(n_step-1))^2);
                v_step(n_step)=v;
                t_step(n_step)=r_step(n_step)/v;
%---------- Boundary reflection in y-direction ----------------------------
                if (y_covered+y_flight<=0)
                    %----- Update values ----------------------------------
                    n_mfp=n_mfp+1;
                    E_mfp(n_mfp)=E_ele;
                    x_mfp(n_mfp)=x_step(n_step);
                    y_mfp(n_mfp)=0;
                    %----- Next condition ---------------------------------
                    vy=-vy;
                    if (vy<=0)
                        v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
                    else
                        v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
                    end
                    dx=dr*cosd(v_theta);
                    dy=dr*sind(v_theta);
                    x_covered=x_step(n_step);
                    y_covered=0;
                    x_flight=0;
                    y_flight=0;   
                    ele_yRefl_flag=true;
                elseif (y_covered+y_flight>=Ly)
                    %----- Update values ----------------------------------
                    n_mfp=n_mfp+1;
                    E_mfp(n_mfp)=E_ele;
                    x_mfp(n_mfp)=x_step(n_step);
                    y_mfp(n_mfp)=Ly;
                    %----- Next condition ---------------------------------
                    vy=-vy;
                    if (vy<=0)
                        v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
                    else
                        v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
                    end
                    dx=dr*cosd(v_theta);
                    dy=dr*sind(v_theta);
                    x_covered=x_step(n_step);
                    y_covered=Ly;
                    x_flight=0;
                    y_flight=0; 
                    ele_yRefl_flag=true;
                end
%---------- FORWARD MOVEMENT COUNT ----------------------------------------     
                if (x_covered+x_flight>=Lx)
                    %----- Update values ----------------------------------
                    n_mfp=n_mfp+1;
                    x_mfp(n_mfp)=Lx;
                    y_mfp(n_mfp)=y_step(n_step);
                    %----- Counting ---------------------------------------
                    n_right_exit(kE)=n_right_exit(kE)+1;
                    right_exit_ele(kE,kNele)=1;
                    right_exit_time(kE,kNele)=sum(t_step);
                    ele_exit_flag=true;
                    break;
%---------- BACK MOVEMENT COUNT -------------------------------------------
                elseif (x_covered+x_flight<=0)
                    %----- Update values ----------------------------------
                    n_mfp=n_mfp+1;
                    x_mfp(n_mfp)=0;
                    y_mfp(n_mfp)=y_step(n_step);
                    %----- Counting ---------------------------------------
                    n_left_exit(kE)=n_left_exit(kE)+1;
                    left_exit_ele(kE,kNele)=1;
                    left_exit_time(kE,kNele)=sum(t_step);
                    ele_exit_flag=true;
                    break;
                end                   
            end % while for r_flight (mfp)
            x_covered=x_covered+x_flight;
            y_covered=y_covered+y_flight;
            %----- Update values at every mfp -----------------------------
            if (ele_exit_flag==true)
                % Do Nothing
                break;
            else
                n_mfp=n_mfp+1;
                x_mfp(n_mfp)=x_step(n_step);
                y_mfp(n_mfp)=y_step(n_step);
            end
        end % while for x_covered
        fprintf('\b\b\b\b\b\b\b');
        fprintf('%6.2f%%', sim_prog); 
    end % Nele end
end % NE end
%----- MF calculation -----------------------------------------------------
    avg_ToF=sum(right_exit_time,2)./sum(right_exit_ele,2);
    flux=1./avg_ToF;
    flux(isnan(flux))=0;
    DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_local-Ec))/(pi*hbar) * (q^(3/2)); % DOS in eV
    TDF_MC=real(flux)'.*real(DOS_3D);
    dE=E_local(2)-E_local(1);
    sigma_MC=zeros(1,NEf);
    for kEf=1:NEf
        f=1./(1+exp((E_local-Ef_local(kEf))/Vth)); % Unitless
        dfdE=-diff(f)./diff(E_local); % 1/Joule
        factor_TDF1_MC=TDF_MC(1:end-1).*dfdE;
        sigma_MC(kEf)=q^2*sum(factor_TDF1_MC)*dE; % Siemen/m
    end
    MF_sigma=sigma_analyt./sigma_MC;
    MF_avg=mean(MF_sigma);
end

