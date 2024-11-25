%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduction band edge Ec estimation by iterative method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear; % clc; % close all

function [Ec] = calculate_Ec(n,m,T,E,Ef)
    constants
    % m=0.2*m0; % Effective mass (kg)
    % T=300; % Temperature in Kelvin
    Vth=(kB/q)*T; % Thermal voltage at T
    % n=5e22; % electron density
    % E=linspace(-2,2,2000); % Energy grid (eV)
    dE=E(2)-E(1);
    % Ef=0;
    Ec_guess=0; % initial guess for Ec
    error=1e-6; 
    max_iterations=5000;
    for iter=1:max_iterations
        f=1./(1+exp((E-Ef)/Vth));
        DOS_3D=(m/(pi*hbar^2)) * sqrt(2*m*(E-Ec_guess))/(pi*hbar) *(q^(3/2)); % DOS in eV
        DOS_3D=real(DOS_3D);
        n_E=DOS_3D.*f; % electron density a function of energy (/m3)
        n_E_int=sum(n_E*dE); % electron density (integral) in /m3 OR n1=trapz(E,n1_E)
        if abs(n_E_int-n)<error
            Ec=Ec_guess; % final Ec
            break;
        end
        % Update the initial guess for Ec based on the comparison of n_E_int
        if n_E_int>n
            Ec_guess=Ec_guess+0.0001;
            Ec=Ec_guess;
        else
            Ec_guess=Ec_guess-0.0001;
            Ec=Ec_guess;
        end
    end
    % fprintf('Estimated Ec= %g meV for electron density n= %g per m3 \n',Ec*1e3,n);
end
