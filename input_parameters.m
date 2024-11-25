%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input material and simulation parameters for MC ray-tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
material_parameters;
%----- Material selection -------------------------------------------------
selected_material=silicon;
%----- Scattering rates ---------------------------------------------------
include_ADP='Yes';
include_ODP='No';
if strcmpi(include_ODP,'Yes')
    include_ODP_abs='Yes';
    include_ODP_ems='Yes';
else
    include_ODP_abs='No';
    include_ODP_ems='No';
end
include_IIS='Yes';

include_electrolyte='No';
include_ordered_poly='Yes';
include_gb_doping='No';
%----- MC simulation parameters -------------------------------------------
n_itr=1; % # of MC simulation iterations
Emin=0.00001; % Minimum energy point (eV)
Emax=0.5; % Maximum energy point (eV)
NE=100; % # of energy points
Nele=1000; % # of electrons per energy point
Ec=0; % CB minimum (eV)
Ef_min=-0.2; % minimum Fermi energy level (eV)
Ef_max=0.2; % maximum Fermi energy level (eV)
NEf=40; % # of Fermi energy points
%----- Material parameters ------------------------------------------------
T=300; % Temperature (Kelvin)
m1=selected_material.m_eff*m0; % effective mass (kg)
EA=selected_material.EA; % electron affinity (eV)
Eg=selected_material.Eg; % energy band gap (eV)
mass_density=selected_material.mass_density; % mass density (kg/m3)
vs=selected_material.vs; % sound velocity (m/s)
if strcmpi(include_ADP,'Yes')
    D_adp=selected_material.D_adp; % acoustic deformation potential (eV)
end
if strcmpi(include_ODP,'Yes')
    D_odp=selected_material.D_odp; % optical deformation potential (eV)
    hbarw0=selected_material.hbarw0; % optical phonon energy (eV)
end
if strcmpi(include_IIS,'Yes')
    clearvars Ef_min Ef_max NEf
    ND=0; % doping density in semiconductor region (/m3)
    Ef=0; % reference Fermi level (eV)
    epsilon_r=selected_material.epsilon_r; % relative permittivity of SC (unitless)
    E=linspace(Emin,Emax,NE); % Energy grid (eV)
    [Ec] = calculate_Ec(ND,m1,T,E,Ef);
%     Ec=0;
end
%----- Poisson solver input parameters for electrolyte systems ------------
if strcmpi(include_electrolyte,'Yes')
    clearvars Ec Ef_min Ef_max NEf
    ni=selected_material.ni*1e6;
    if strcmpi(include_IIS, 'Yes')
        ND_SC=ND; % doping density in semiconductor region (/m3)
    else
        ND_SC=0;
    end
    epsilon_r=selected_material.epsilon_r; % relative permittivity of SC (unitless)
    Ef=0; % equilibrium Ef as reference level (eV)
    E_redox=1; % REDOX level in electrolyte (eV)
    n_cation=1e26; % +ve ion concentration in electrolyte region (/m3)
    n_anion=1e26; % -ve ion concentration in electrolyte region (/m3)
    epsilon_r_sol=80; % relative permitivity of electrolyte solution (unitless)
end
%----- Poisson solver input parameters for poly+GB systems ----------------
if strcmpi(include_ordered_poly,'Yes') || strcmpi(include_gb_doping,'Yes')
    m1=selected_material.m_eff*m0; % effective mass in SC region 1 (kg)
    m2=selected_material.m_eff*m0; % effective mass in SC region 2 (kg)
    epsilon_r1=selected_material.epsilon_r; % relative permittivity of SC 1 (unitless)
    epsilon_r2=selected_material.epsilon_r; % relative permittivity of SC 2 (unitless)
    Ec1=0; % CBM in SC bulk region
    Ef1=0; % Fermi level in SC bulk region
    ND1=0; % doping density in SC bulk region (/m3)
    Ec2=0; % CBM in SC poly region
    Ef2=0.1; % Fermi level in SC poly region (change it for different doping levels in ploy region)
    ND2=1e26; % doping density in in SC poly region (/m3)
    Ef=0; % equlibrium Ef as reference level (eV)
    clearvars Ec Ef_min Ef_max
    NEf=1;
end
%--------------------------------------------------------------------------
if isequal(selected_material,silicon)
    select_material='Si';
    clear Bi2Te3;
    clear ATO;
    clear default;
elseif isequal(selected_material,Bi2Te3)
    select_material='Bi2Te3';
    clear silicon;
    clear ATO;
    clear default;
elseif isequal(selected_material,ATO)
    select_material='ATO';
    clear silicon;
    clear Bi2Te3;
    clear default;
else
    select_material='default';
    clear silocon;
    clear Bi2Te3;
    clear ATO;
end
%--------------------------------------------------------------------------
if strcmpi(include_ADP,'Yes') && strcmpi(include_ODP,'No') && strcmpi(include_IIS,'No')
    selected_scattering='ADP';
elseif strcmpi(include_ADP,'No') && strcmpi(include_ODP,'Yes') && strcmpi(include_IIS,'No')
    selected_scattering='ODP';
elseif strcmpi(include_ADP,'Yes') && strcmpi(include_ODP,'Yes') && strcmpi(include_IIS,'No')
    selected_scattering='Ph';
elseif strcmpi(include_ADP,'Yes') && strcmpi(include_ODP,'No') && strcmpi(include_IIS,'Yes')
    selected_scattering='ADP+IIS';
elseif strcmpi(include_ADP,'Yes') && strcmpi(include_ODP,'Yes') && strcmpi(include_IIS,'Yes')
    selected_scattering='Ph+IIS';
elseif strcmpi(include_ADP,'No') && strcmpi(include_ODP,'Yes') && strcmpi(include_IIS,'Yes')
    selected_scattering='ODP+IIS';
end 

if strcmpi(include_electrolyte,'Yes')
    selected_poisson='E_{redox}';
end
