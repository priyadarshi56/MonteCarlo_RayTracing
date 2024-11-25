%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Material parameter database
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Silicon ------------------------------------------------------------
silicon.m_eff=0.9; % effective mass (unitless)
% silicon.m_eff=0.2; % effective mass (unitless)
silicon.ni=1e10; % intrinsic carrier density (/cm3)
silicon.EA=-4.05; % electron affinity (eV)
silicon.Eg=1.1; % energy band gap (eV)
silicon.mass_density=2329; % mass density (kg/m3)
silicon.vs=8540; % sound velocity (m/s)
silicon.D_adp=4.85; % acoustic deformation potential (eV)
% silicon.D_adp=10.7; % acoustic deformation potential (eV)      for electron (FCT)
silicon.epsilon_r=11.7; % relative permittivity (unitless)
silicon.D_odp=8; % optical deformation potential (eV)
silicon.hbarw0=0.063; % optical phonon energy (eV)

%----- Bi2Te3 bulk --------------------------------------------------------
Bi2Te3.m_eff=0.2; % effective mass (unitless)
Bi2Te3.ni=1e19; % intrinsic carrier density (/cm3)
Bi2Te3.EA=-4.2; % electron affinity (eV)
Bi2Te3.Eg=0.1; % energy band gap (eV)
Bi2Te3.mass_density=7700; % mass density (kg/m3)
Bi2Te3.vs=3058; % sound velocity (m/s)
Bi2Te3.D_adp=11; % acoustic deformation potential (eV)
Bi2Te3.epsilon_r=85; % relative permittivity (unitless)
Bi2Te3.D_odp=12; % optical deformation potential (eV)
Bi2Te3.hbarw0=0.06; % phonon frequency ( )

%----- ATO ----------------------------------------------------------------
ATO.m_eff=0.2; % effective mass (unitless)
ATO.ni=1e19; % intrinsic carrier density (/cm3)
ATO.EA=-4.2; % electron affinity (eV)
ATO.Eg=3.6; % energy band gap (eV)
ATO.mass_density=7700; % mass density (kg/m3)
ATO.vs=3058; % sound velocity (m/s)
ATO.D_adp=12; % acoustic deformation potential (eV)
ATO.epsilon_r=40; % relative permittivity (unitless)
ATO.D_odp=9.5; % optical deformation potential (eV)
ATO.hbarw0=0.06; % phonon frequency ( )

%----- Default value ------------------------------------------------------
default.m_eff=0.26; % effective mass (unitless)
default.ni=1e10; % intrinsic carrier density (/cm3)
default.EA=-4.2; % electron affinity (eV)
default.mass_density=2330; % mass density (kg/m3)
default.vs=8540; % sound velocity (m/s)
default.D_adp=9.5; % acoustic deformation potential (eV)
default.epsilon_r=12; % relative permittivity (unitless)
default.D_odp=9.5; % optical deformation potential (eV)
default.hbarw0=0.06; % phonon frequency ( )
