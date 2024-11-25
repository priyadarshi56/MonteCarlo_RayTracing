%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Ray tracing in two dimensional domain from left to right in rectangular 
%   uniform mesh with nanostructured pores and/or grain boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if strcmpi(include_pores,'Yes')
    refl_count=0;
    if strcmpi(include_electrolyte,'Yes')
        fprintf('\nRay-tracing in porous domain with electrolyte...\n');
    else
        fprintf('\nRay-tracing in porous domain...\n');
    end
elseif strcmpi(include_ordered_poly,'Yes')
    if strcmpi(include_ordered_gb,'Yes')
        fprintf('\nRay-tracing in ordered poly with ordered GB domain with Poisson solution...\n');
        n_reflect=0; n_transmit=0;
    else
        fprintf('\nRay-tracing in ordered poly domain (no scattering) with Poisson solution...\n');
    end
elseif strcmpi(include_random_gb,'Yes')
    if strcmpi(include_gb_doping,'Yes')
        fprintf('\nRay-tracing in random GB domain with Poisson solution...\n');
    else
        fprintf('\nRay-tracing in random GB domain...\n');
    end
    n_reflect=0; n_transmit=0;
elseif strcmpi(include_ordered_gb,'Yes')
    fprintf('\nRay-tracing in ordered GB domain...\n');
    n_reflect=0; n_transmit=0;
elseif strcmpi(include_pores,'Yes') &&  strcmpi(include_random_gb,'Yes')
    fprintf('\nRay-tracing in pores+GB domain...\n');
else
    fprintf('\nRay-tracing in pristine domain...\n');
end
%--------------------------------------------------------------------------
for k_itr=1:n_itr
fprintf('>>>for iteration # %d ...',k_itr);

if exist('Ec_element','var')==1
    Emin=min(min(Ec_element)-0.05,Emin);
    Emax=max(max(Ec_element)+0.1,Emax);
end
E=linspace(Emin,Emax,NE);
u=1e9; % visualization in nm
n_left_exit=zeros(1,NE);
n_right_exit=zeros(1,NE);
left_exit_ele=zeros(NE,Nele);
right_exit_ele=zeros(NE,Nele);
left_exit_time=zeros(NE,Nele);
right_exit_time=zeros(NE,Nele);
sim_prog=0;
fprintf('simulation progress: %2d%%\n',sim_prog);
for kE=1:NE
    if isGuiMode % GUI mode - updating progress on the same line
        sim_prog=kE/NE*100; % Simulation progress in percentage only on NE
        fprintf('\b\b\b\b%3.0f%%',sim_prog); % display in one line
    else % Terminal mode - writing progress on a new line
        fprintf('Remaining energy point= %g \n',NE-kE);
    end
    %----------------------------------------------------------------------    
    if exist('Ec','var')==1 && (E(kE)-Ec<0)
        continue;
    end   
    % if kE>1
    %     delete(g1); delete(g2);
    % end
    for kNele=1:Nele
        % if isGuiMode % GUI mode - updating progress on the same line
        %     sim_prog=((kE-1)*Nele+kNele)/(NE*Nele)*100; % Simulation progress in percentage on NE and Nele
        %     fprintf('\b\b\b\b%3.0f%%',sim_prog); % display in one line
        % else % Terminal mode - writing progress on a new line
        %     fprintf('Energy point= %d, Electron= %d\n',kE,kNele); % display electron number and energy point
        % end        
%--------------------------------------------------------------------------
        % if kNele>1 && kE>=1
        %     delete(g1); delete(g2);
        % end
%--------------------------------------------------------------------------
%-------Initialization-----------------------------------------------------
        clearvars x_step y_step r_step t_step v_step x_mfp y_mfp mfp_E
        dx=dLx*(3/4);
        dy=dLy*(3/4);
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
        found_element = find_element_new(Lx,Ly,Npx,Npy,node,x_covered,y_covered);

        if exist('Ec_element','var')==1
            E_ele=E(kE)-Ec_element(found_element);
        else
            E_ele=E(kE)-Ec;
        end
        if E_ele<=0
            continue; % skip kNele for negative electron energy
        end

        n_mfp=1; % number of mfp encountered
        x_mfp(n_mfp)=x_covered;
        y_mfp(n_mfp)=y_covered;
        E_mfp(n_mfp)=E_ele;
        % fprintf('%g electron starts from element %d at point (%g,%g) nm with Ec= %g meV \n',kNele,found_element,x_covered*u,y_covered*u,E_ele*1e3);
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
%             dr=sqrt(dx^2+dy^2);
            if (vy<=0)
                v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
            else
                v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
            end
            
            if (n_step==1 && vx<0)
                dx=-dr*cosd(v_theta);
                dy=dr*sind(v_theta);
                vx=-vx; % pushed into
                % fprintf('Velocity (%g,%g) at an angle %g w.r.t (+x)-axis, pushed into with new velocity (%g,%g). \n',-vx,vy,v_theta-90,vx,vy);
            else
                dx=dr*cosd(v_theta);
                dy=dr*sind(v_theta);
                % fprintf('Velocity (%g,%g) at an angle %g w.r.t (+x)-axis. \n',vx,vy,v_theta);
            end
            x_flight=0;
            y_flight=0;
            r_flight=0;

            if E_ele<=0
                break; % skip for negative electron energy
            end
            % calculate_mfp_E;
            if strcmpi(include_IIS,'Yes')
                mfp_E=real(calculate_mfp_IIS(E_ele,ND));    % <<<<< check for real
            else
                mfp_E=real(calculate_mfp(E_ele));    % <<<<< check for real
            end
            % fprintf('Incremental step (dx,dy)=(%g,%g) nm, will move a distance= %g nm at each step to cover a mfp= %g nm. \n---------->\n',dx*u,dy*u,dr*u,mfp_E*u);
            while (r_flight<=mfp_E)
                found_element=find_element_new(Lx,Ly,Npx,Npy,node,abs(x_covered+x_flight),abs(y_covered+y_flight));
                if exist('Ec_element','var')==1
                    E_ele=E(kE)-Ec_element(found_element);
                else
                    E_ele=E(kE)-Ec;
                end
                if E_ele<=0
                    break; % skip for negative electron energy
                end
                k=(sqrt(2*m1)/hbar)*sqrt(E_ele*q);
                kx=k*cosd(theta);
                ky=k*sind(theta);
                vx_new=(hbar*kx)/m1; % speed in x-direction
                vy_new=(hbar*ky)/m1; % speed in y-direction
                v=sqrt(vx_new^2+vy_new^2);
                
                n_step=n_step+1;
                x_flight=x_flight+dx;
                y_flight=y_flight+dy;
                r_flight=r_flight+sqrt(dx^2+dy^2);
                x_step(n_step)=x_covered+x_flight;
                y_step(n_step)=y_covered+y_flight;
                r_step(n_step)=sqrt((x_step(n_step)-x_step(n_step-1))^2+(y_step(n_step)-y_step(n_step-1))^2);
                v_step(n_step)=v;
                t_step(n_step)=r_step(n_step)/v; 
                found_element=find_element_new(Lx,Ly,Npx,Npy,node,abs(x_covered+x_flight),abs(y_covered+y_flight));
                % fprintf('In %g move electron in element= %d at point (%g,%g) nm with Ec= %g meV and completing %g nm distance \n',n_step-1,found_element,x_step(n_step)*u,y_step(n_step)*u,E_ele*1e3,r_flight*u);
%---------- Boundary reflection in y-direction ----------------------------
                if (y_covered+y_flight<=0)
                    % fprintf('Found Y-reflection with y=0nm at point (%g,%g) nm \n',x_step(n_step)*u,y_step(n_step)*u);
                    %----- Update values (not required) -------------------
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
                    % fprintf('Velocity (%g,%g) at an angle %g w.r.t (+x)-axis, after reflection at Y=0 nm\n',vx,vy,v_theta);
                    dx=dr*cosd(v_theta);
                    dy=dr*sind(v_theta);
                    % fprintf('After Y=0 reflection, incremental step (dx,dy)=(%g,%g) nm\n',dx*u,dy*u);
                    x_covered=x_step(n_step);
                    y_covered=0;
                    x_flight=0;
                    y_flight=0;   
%                     r_flight=0; % Reset mfp checker
                    ele_yRefl_flag=true;
                    % delete(g1); delete(g2);
                    % g1=plot(x_step*u,y_step*u,'-b','linewidth',1,'marker','.','markersize',8);
                    % g2=plot(x_mfp*u,y_mfp*u,'.y','linewidth',1,'marker','*','markersize',6);            
                    % pause(0.05);
                elseif (y_covered+y_flight>=Ly)
                    % fprintf(1,'Found Y-reflection with y=%gnm at point (%g,%g) nm \n',Ly*u,x_step(n_step)*u,y_step(n_step)*u);
                    %----- Update values (not required) -------------------
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
                    % fprintf('Velocity (%g,%g) at an angle %g w.r.t (+x)-axis, after reflection at Y=%g nm\n',vx,vy,v_theta,Ly*u);
                    dx=dr*cosd(v_theta);
                    dy=dr*sind(v_theta);
                    % fprintf(1,'After Y=%gnm reflection, incremental step (dx,dy)=(%g,%g) nm\n',Ly*u,dx*u,dy*u);
                    x_covered=x_step(n_step);
                    y_covered=Ly;
                    x_flight=0;
                    y_flight=0; 
%                     r_flight=0; % Reset mfp checker
                    ele_yRefl_flag=true;
                    % delete(g1); delete(g2);
                    % g1=plot(x_step*u,y_step*u,'-b','linewidth',1,'marker','.','markersize',8);
                    % g2=plot(x_mfp*u,y_mfp*u,'.y','linewidth',1,'marker','*','markersize',6);
                    % pause(0.05);
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
                    % fprintf('Electron moved forward and exited from right end. >>>>>>>>>> \n');
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
                    % fprintf('Electron moved back and exited from left end. <<<<<<<<<< \n');
                    break;
                end
%%------ REFLECTION from PORE ---------------------------------------------
                if strcmpi(include_pores,'Yes')
                        if any(found_element==element_in_pore)
                            % fprintf('Pore found in element %g at point (%g,%g).\n ',found_element,x_step(n_step)*u,y_step(n_step)*u);
                            % delete(g1); delete(g2);
                            % g1=plot(x_step*u,y_step*u,'-b','linewidth',1,'marker','.','markersize',8);
                            % g2=plot(x_mfp*u,y_mfp*u,'.y','linewidth',1,'marker','*','markersize',6);
                            % pause(0.1);
                            [~,C_element_index]=find(found_element==element_in_pore);
                            C_node_vertex=node(element(C_element(C_element_index),:),:);
                            %----- Update values before reflection --------
    %                         n_step=n_step-1;
                            n_mfp=n_mfp+1;
                            x_mfp(n_mfp)=x_step(n_step);
                            y_mfp(n_mfp)=y_step(n_step);
                            %----- Reflection at pore ---------------------
                            dt_refl=t_step(n_step);
                            [xr,yr,vx,vy]=pore_reflection(x_step(n_step-1),y_step(n_step-1),x_step(n_step),y_step(n_step),dt_refl,C_node_vertex);
                            %----- Update values after reflection ---------
                            n_step=n_step+1;
                            x_step(n_step)=xr;
                            y_step(n_step)=yr;
                            t_step(n_step)=dt_refl;
                            v_step(n_step)=v;
                            %----- Next condition -------------------------
                            if (vy<=0)
                                v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
                            else
                                v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
                            end
                            dx=dr*cosd(v_theta);
                            dy=dr*sind(v_theta);
                            x_flight=x_flight+dx;
                            y_flight=y_flight+dy;
                            r_flight=r_flight+sqrt(dx^2+dy^2);

                            refl_count=refl_count+1;
                            n_pore_refl(kE,kNele)=refl_count;
                            % [found_element,found_node]=find_element(node,element,abs(x_covered+x_flight),abs(y_covered+y_flight)); % NOT required
                            % fprintf('After pore reflection, electon found in element %g at point (%g,%g).\n ',found_element,abs(x_covered+x_flight)*u,abs(y_covered+y_flight)*u);
                        end % Checking for element in pores                    
                end
%---------- PORE reflection section end -----------------------------------
%%------ REFLECTION/TRANSMISSION from GB ----------------------------------
                if strcmpi(include_random_gb,'Yes') || strcmpi(include_ordered_gb,'Yes1') 
                    for k_gb=1:length(element_in_grain)
                        if any(found_element==element_in_grain{k_gb})
                        %----- Update values before interaction -----------
                            n_mfp=n_mfp+1;
                            x_mfp(n_mfp)=x_step(n_step);
                            y_mfp(n_mfp)=y_step(n_step);
                            %----- Interaction at grain boundary ----------
                            dt_refl=t_step(n_step);
                            [xr,yr,vx,vy,reflect,transmit]=grain_boundary_interaction(x_step(n_step-1),y_step(n_step-1),x_step(n_step),y_step(n_step),grain_line_p1(k_gb,:),grain_line_p2(k_gb,:),dt_refl);
                            n_reflect=n_reflect+reflect;
                            n_transmit=n_transmit+transmit;
                             %----- Update values after interaction -------
                            n_step=n_step+1;
                            x_step(n_step)=xr;
                            y_step(n_step)=yr;
                            t_step(n_step)=dt_refl;
                            %----- Next condition -------------------------
                            if (vy<=0)
                                v_theta=360-acosd(vx/(sqrt(vx^2+vy^2)));
                            else
                                v_theta=acosd(vx/(sqrt(vx^2+vy^2)));
                            end
                            dx=dr*cosd(v_theta);
                            dy=dr*sind(v_theta);
                            x_flight=x_flight+dx;
                            y_flight=y_flight+dy;
                            r_flight=r_flight+sqrt(dx^2+dy^2);
                            break;
                        end % Checking for element in GB
                    end % Loop for element in GB
                end
%-----------GB interaction section end ------------------------------------                

            %----- Ploting ------------------------------------------------
            % if n_step>2
            %     delete(g1); delete(g2);
            % end
            % g1=plot(x_step*u,y_step*u,'-b','linewidth',1,'marker','.','markersize',8);
            % g2=plot(x_mfp*u,y_mfp*u,'.r','linewidth',1,'marker','*','markersize',6);
            % pause(0.05);
            end % while for r_flight (mfp)
            x_covered=x_covered+x_flight;
            y_covered=y_covered+y_flight;
            %----- Update values at every mfp -----------------------------
            if (ele_exit_flag==true)
                % Do Nothing
                % delete(g1); delete(g2);
                % g1=plot(x_step*u,y_step*u,'-b','linewidth',1,'marker','.','markersize',8);
                % g2=plot(x_mfp*u,y_mfp*u,'.r','linewidth',1,'marker','*','markersize',6);
                break;
            else
                n_mfp=n_mfp+1;
                x_mfp(n_mfp)=x_step(n_step);
                y_mfp(n_mfp)=y_step(n_step);
                % fprintf(1,'%g MFP completed in %g step, while travelled %g nm\n',n_mfp-1,n_step-1,r_flight*u);
                % fprintf(1,'\n---------- %g MFP started ----------\n',n_mfp);
            end
        end % while for x_covered
%----- Storing values -----------------------------------------------------
        % x_step_store{kE,kNele}=x_step;
        % y_step_store{kE,kNele}=y_step;
        % r_step_store{kE,kNele}=r_step;
        % v_step_store{kE,kNele}=v_step;
        % t_step_store{kE,kNele}=t_step;
        % x_mfp_store{kE,kNele}=x_mfp;
        % y_mfp_store{kE,kNele}=y_mfp;
%------ Ploting ray tracing only ------------------------------------------
        % g1=plot(x_mfp_store{kE,kNele}*u,y_mfp_store{kE,kNele}*u,'-*r','linewidth',1,'markersize',10);
        % g2=plot(x_step_store{kE,kNele}*u,y_step_store{kE,kNele}*u,'-..y','linewidth',1,'markersize',8);
        % fk=1;
    end % Nele end
end % NE end
% ----- FLUX -------------------------------------------------------------
avg_ToF(k_itr,:)=sum(right_exit_time,2)./sum(right_exit_ele,2);
flux(k_itr,:)=1./avg_ToF(k_itr,:);
fprintf(' DONE!\n');
MC_time_taken1(k_itr)=toc;
end
MC_time_taken=sum(MC_time_taken1);
fprintf('\nMC simulation time= %g seconds.\n',MC_time_taken);

