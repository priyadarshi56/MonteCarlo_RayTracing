%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     2D input geometry details for MC ray-tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Domain geometry ---------------------------------------------------
Lx=200e-9; % length of domain (meter)
Ly=100e-9; % widthth of domain (meter)
dLx=1e-9; % step size in x-direction (meter)
dLy=1e-9; % step size in y-direction (meter)
buffer_contact=10e-9; % size of the contact buffer region (highly dopped or like metal)
global nc % # of points on a pore periphery
nc=100;
%% ----- Details of porous geometry ---------------------------------------
include_pores='Yes'; % include pores ? TYPE 'yes' or 'no'
if strcmpi(include_pores,'Yes')
%----- Select only one type of pores arrangement --------------------------
    ordered_pores='Yes'; % TYPE 'yes' or 'no'
    staggered_pores='No'; % TYPE 'yes' or 'no'
    nonoverlap_random_pores='No'; % TYPE 'yes' or 'no'
    nonoverlap_random_size_pores='No'; % TYPE 'yes' or 'no'
    overlap_random_size_pores='No'; % TYPE 'yes' or 'no'
    nonoverlap_random_oval_pores='No'; % TYPE 'yes' or 'no'
    if strcmpi(ordered_pores,'Yes') || strcmpi(staggered_pores,'Yes')
        n_pore_x=10; % # of pores in x-direction
        n_pore_y=5; % # of pores in y-direction
        pore_r=5e-9; % radius of pore
    elseif strcmpi(nonoverlap_random_pores,'Yes')
        n_pore=10; % # of pores
        pore_r=10.095e-9; % radius of pore, for equal size
    elseif strcmpi(nonoverlap_random_size_pores,'Yes')
        n_pore=30; % # of pores
        pore_r_min=5e-9; % minimum radius of a pore
        pore_r_max=10e-9; % maximum radius of a pore
    elseif strcmpi(overlap_random_size_pores,'Yes')
        n_pore=70; % # of pores
        pore_r_min=4e-9; % minimum radius of a pore
        pore_r_max=8e-9; % maximum radius of a pore
    elseif strcmpi(nonoverlap_random_oval_pores,'Yes')
        n_pore=12; % # of pores
        pore_r_min=4e-9; % minimum radius of a pore
        pore_r_max=20e-9; % maximum radius of a pore
    end
    
    if strcmpi(include_electrolyte,'Yes')
        dH_thick=1e-9; % Helmholtz layer thickness
    else
        dH_thick=0;
    end
end
%% ----- Details of nanocrystalline grains --------------------------------
include_random_gb='No'; % include random grains ? TYPE 'yes' or 'no'
include_ordered_gb='No'; % include ordered grain boundries lines ? TYPE 'yes' or 'no'
include_ordered_poly='No'; % include ordered poly region ? TYPE 'yes' or 'no'
if strcmpi(include_random_gb,'Yes')
    n_grain=15; % # grain seeds
elseif strcmpi(include_ordered_gb,'Yes')
    n_grain_x=4; % # of ordered grain boundries in x-direction
    n_grain_y=3; % # of ordered grain boundries in y-direction
    grain_length=15e-9; % length of ordered grains
    grain_width=15e-9; % width of ordered grains
end
if strcmpi(include_ordered_poly,'Yes')
    n_poly_x=4; % # of ordered grains in x-direction
    n_poly_y=3; % # of ordered grains in y-direction
    poly_length=22e-9; % length of poly region
    poly_width=22e-9; % width of poly region
    include_buffer='No';
end
%% ------------------------------------------------------------------------
include_gb_pores='No'; % include combination of GB and pores ? TYPE 'yes' or 'no'

