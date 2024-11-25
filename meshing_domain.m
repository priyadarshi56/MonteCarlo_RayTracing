%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   File to call appropriate geometry and their meshings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Rectangular domain with nodes and elements -------------------------
tic

%----- Meshing domain with rectangular elements and nodes -----------------
d=2; % domain dimension
domain=[0 0;Lx 0;Lx Ly;0 Ly]; % domain coordinates
domain_area=Lx*Ly;
% [n_node,node,n_element,element,Npx,Npy]=create_rectangular_mesh_old(Lx,Ly,dLx,dLy);
[n_node,node,n_element,element,x_cord,y_cord,Npx,Npy] = create_rectangular_mesh(Lx,Ly,dLx,dLy);
%----- Creating pores -----------------------------------------------------
global nc
nc=100;
if strcmpi(include_pores,'Yes')
    if strcmpi(ordered_pores,'Yes')
        n_pore=n_pore_x*n_pore_y; % For ordered pores
        % [C,rad,x_pore,y_pore] = create_ordered_pores(Lx,Ly,n_pore_x,n_pore_y,pore_r);
        [C,rad,x_pore,y_pore,x_pore_both,y_pore_both] = create_ordered_pores_dH(Lx,Ly,n_pore_x,n_pore_y,pore_r,dH_thick);
        [porosity] = calculate_porosity(rad,domain_area);
        selected_geometry='Pores_ord';
    elseif strcmpi(staggered_pores,'Yes')
        [C,rad,x_pore,y_pore,n_pore] = create_staggered_pores(Lx,Ly,n_pore_x,n_pore_y,pore_r);
        [porosity] = calculate_porosity(rad,domain_area);
        selected_geometry='Pores_stg';
    elseif strcmpi(nonoverlap_random_pores,'Yes')
        % [C,x_pore,y_pore] = create_nonoverlap_random_pores(Lx,Ly,n_pore,pore_r);
        [C,rad,x_pore,y_pore,x_pore_both,y_pore_both] = create_nonoverlap_random_pores_dH(Lx,Ly,n_pore,pore_r,dH_thick);
        [porosity] = calculate_porosity(rad,domain_area);
        selected_geometry='Pores_NRP';
    elseif strcmpi(nonoverlap_random_size_pores,'Yes')
        % [C,pore_r,x_pore,y_pore] = create_nonoverlap_random_size_pores(Lx,Ly,n_pore,pore_r_min,pore_r_max);
        [C,pore_r,x_pore,y_pore,x_pore_both,y_pore_both] = create_nonoverlap_random_size_pores_dH(Lx,Ly,n_pore,pore_r_min,pore_r_max,dH_thick);
        [porosity] = calculate_porosity(pore_r,domain_area);
        selected_geometry='Pores_NRSP';
    elseif strcmpi(overlap_random_size_pores,'Yes')
        % [C,pore_r,x_pore,y_pore] = create_overlap_random_size_pores(Lx,Ly,n_pore,pore_r_min,pore_r_max);
        [C,pore_r,x_pore,y_pore,x_pore_both,y_pore_both] = create_overlap_random_size_pores_dH(Lx,Ly,n_pore,pore_r_min,pore_r_max,dH_thick);
        [porosity] = calculate_porosity(pore_r,domain_area);
        selected_geometry='Pores_ORSP';
    elseif strcmpi(nonoverlap_random_oval_pores,'Yes')
        [C,x_pore_rotate,y_pore_rotate] = create_nonoverlap_random_oval_pores(Lx,Ly,n_pore,pore_r); 
        [porosity] = calculate_porosity(r,domain_area);
        selected_geometry='Pores_NROP';
    end
    fprintf('Average porosity= %g %%\n',porosity);
    fprintf(['Meshing ',selected_geometry,' domain...']);
    % [element_in_pore,C_element] = find_element_in_pore(node,element,n_pore,C,x_pore,y_pore); % elements in all pore
    [element_in_pore,node_in_pore,element_in_both,node_in_both,node_in_dH,C_element] = find_element_in_pore_dH(node,element,n_pore,C,x_pore,y_pore,x_pore_both,y_pore_both);
% load('filen_p=40.mat')
    % element_in_pore=sort(element_in_pore);
end
%----- Creating random grain boundaries -----------------------------------
if strcmpi(include_random_gb,'Yes')
    [grain_center_x,grain_center_y,grain_line_xsegments,grain_line_ysegments]=create_grain_boundary(domain,dLx,dLy,n_grain);
    [avg_grain_size]=calculate_grain_size(grain_center_x,grain_center_y);
    fprintf('Average grain size= %g nm\n',avg_grain_size*1e9);
    selected_geometry='rGB';
    fprintf(['Meshing ',selected_geometry,' domain...']);
    for k=1:length(grain_line_xsegments)
        grain_line_p1(k,1)=grain_line_xsegments{k}(1);
        grain_line_p1(k,2)=grain_line_ysegments{k}(1);
        grain_line_p2(k,1)=grain_line_xsegments{k}(end);
        grain_line_p2(k,2)=grain_line_ysegments{k}(end);
    end
    [element_in_grain,grain_line_slope] = find_element_at_grain(grain_line_xsegments,grain_line_ysegments,Lx,Ly,node,Npx,Npy); % elements at grain boundaries
end
%----- Creating ordered grain boundaries ----------------------------------
if strcmpi(include_ordered_gb,'Yes')
    selected_geometry='oGB';
    fprintf(['Meshing ',selected_geometry,' domain...']);
    create_ordered_gb;
end
%----- Creating poly region -----------------------------------------------
if strcmpi(include_ordered_poly,'Yes')
    selected_geometry='Poly';
    fprintf(['Meshing ',selected_geometry,' domain...']);
    create_poly_structure;
    % create_poly_structure_buffer;
    avg_grain_area=poly_length*poly_width*1e18;    
end
%--------------------------------------------------------------------------
%----- Selecting domain name ----------------------------------------------
if strcmpi(include_gb_pores,'Yes')
    selected_geometry='GB+Pores';
elseif strcmpi(include_ordered_poly,'Yes') && strcmpi(include_ordered_gb,'Yes')
    selected_geometry='Poly+oGB';
elseif strcmpi(include_pores,'No') && strcmpi(include_random_gb,'No')...
        && strcmpi(include_ordered_gb,'No') && strcmpi(include_ordered_poly,'No') && strcmpi(include_gb_pores,'No')
    selected_geometry='Pristine';
    fprintf('Pristine meshing...');
end

% fprintf(['Meshing ',selected_geometry,' domain...\n']);
mesh_time=toc;
fprintf('DONE!\nMeshing time=%g seconds\n',mesh_time);

