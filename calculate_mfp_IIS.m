function [mfp] = calculate_mfp_IIS(E_ele,ND)
    constants
    input_parameters
    Vth=kB/q*T;
    %----- Non-polar Acoustic Phonon Scattering (ADP) ---------------------
    if strcmpi(include_ADP,'Yes')
        cl=mass_density*vs^2; % elastic constant
        c_ap=(pi*kB*T*D_adp^2)/(hbar*cl) * q;
        DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*E_ele)/(pi*hbar) * (q^(3/2)); % (1/(eV.m3))
        DOS_3D=real(DOS_3D);
        ADP_rate=c_ap*DOS_3D; % 1/sec
    else
        ADP_rate=zeros(1,length(E_ele));
    end
    %----- Non-polar Optical Phonon Scattering (ODP) ----------------------
    if strcmpi(include_ODP,'Yes')
        w0=hbarw0/hbar;
        N0=1/(exp(hbarw0/Vth)-1); % Bose-Einstein factor
        c_op=(pi*D_odp^2)/(mass_density*w0);
        if strcmpi(include_ODP_abs,'Yes')
            DOS_3D_abs=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_ele+hbarw0))/(pi*hbar) * (q^(3/2)); % (1/(eV.m3))
            DOS_3D_abs=real(DOS_3D_abs);
            ODP_abs_rate=c_op*N0*DOS_3D_abs; % absorption rate (1/sec)
        else
            ODP_abs_rate=zeros(1,length(E_ele));
        end
        if strcmpi(include_ODP_ems,'Yes')
            ODP_ems_rate=zeros(1,length(E_ele));
            DOS_3D_ems=zeros(1,length(E_ele));
            for kE=1:length(E_ele)
                if (E_ele(kE)<hbarw0)
                    ODP_ems_rate(kE)=0; % emission rate (1/sec)
                else
                    DOS_3D_ems(kE)=(m1/(pi*hbar^2)) * sqrt(2*m1*(E_ele(kE)-hbarw0))/(pi*hbar) * (q^(3/2)); % (1/(eV.m3))
                    DOS_3D_ems(kE)=real(DOS_3D_ems(kE));
                    ODP_ems_rate(kE)=c_op*(N0+1)*DOS_3D_ems(kE); % emission rate (1/sec)
                end
            end
        else
            ODP_ems_rate=zeros(1,length(E_ele));
        end
    else
        ODP_abs_rate=zeros(1,length(E_ele));
        ODP_ems_rate=zeros(1,length(E_ele));
    end
    %----- Ionized Impurity Scattering (IIS) ------------------------------
    if strcmpi(include_IIS,'Yes')
    %----- Brooks-Herring Model -------------------------------------------
        ND1=ND; % assuming total ionization
        Ej=E_ele*q; % energy in Joule
        c_IIS=(16*pi*sqrt(2*m1)*epsilon_r^2*epsilon0^2)/(q^4*ND1);
        LDebye=sqrt((epsilon_r*epsilon0*kB*T)/(q^2*ND1)); % Debye length in meter
        gamma2_bh=(8*m1*LDebye^2/hbar^2)*Ej; % unitless
        tau_BH=c_IIS*(log(1+gamma2_bh)-gamma2_bh./(1+gamma2_bh)).^(-1).*Ej.^(3/2);
        IIS_rate_BH=1./(tau_BH); % (1/sec)
    %----- Strong screening -----------------------------------------------
        DOS_3D=(m1/(pi*hbar^2)) * sqrt(2*m1*E_ele)/(pi*hbar) * (q^(3/2)); % (1/(eV.m3))
        DOS_3D=real(DOS_3D);
        tau_SS=(hbar/(pi*ND1)) * ((epsilon0*epsilon_r)^2/(q^2*LDebye^2)^2) .* DOS_3D;
        IIS_rate_SS=1./tau_SS; % (1/sec)
        NI_ref=ND1;
        screen_ratio=ND1/(NI_ref*0.2+ND1) * 0.1;
        IIS_rate_total=1*IIS_rate_BH*(1-screen_ratio)+IIS_rate_SS*(screen_ratio);
    else
        IIS_rate_BH=zeros(1,length(E_ele));
        IIS_rate_total=zeros(1,length(E_ele));
    end
    total_rate=ADP_rate+ODP_abs_rate+ODP_ems_rate+IIS_rate_total;
    total_time=1./total_rate; % time (sec)
    %----- mfp calculation ------------------------------------------------
    k=real((sqrt(2*m1*q)/hbar)*sqrt(E_ele)); % parabolic relation of E-k (1/m)
    v_k=(hbar*k)/m1; % speed (m/sec)
    v_E=real(sqrt(2*q*E_ele/m1)); % speed (m/sec) when energy in eV
    mfp=v_E.*total_time; % mean free path (meter)
end
