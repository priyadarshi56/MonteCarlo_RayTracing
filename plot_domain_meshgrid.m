%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plot domain geometry with meshgrid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'units','centimeters','position',[24 2 28 16]);
% ax1=axes('units','centimeters','position',[4.8 4.6 26 14]);
hold on; box on; u=1e9; % plot visualization in nm
%----- Rectangular domain -------------------------------------------------
plot(polyshape([0 Lx Lx 0]*u,[0 0 Ly Ly]*u),'facecolor',[0.5 0.5 0.5],'FaceAlpha',0.25); % 0.19 0.55 0.91
b1=plot([0 0]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
b2=plot([Lx Lx]*u,[0 Ly]*u,'-','linewidth',4,'color',[1 0.5 0.8]); % boundary line in X-direction
b3=plot([0 Lx]*u,[0 0]*u,'-k','linewidth',6); % boundary line in Y-direction
b4=plot([0 Lx]*u,[Ly Ly]*u,'-k','linewidth',6); % boundary line in Y-direction
%----- Element plot and indicating number ---------------------------------
light_blue=[0.3 0.3 0.83];
element_limit=5000;
if n_element<element_limit
    for k_element=1:n_element
        hold on;
        plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':','color',light_blue);
        plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':','color',light_blue);
        plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':','color',light_blue);
        plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':','color',light_blue);
        x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
        y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
        text(x*u,y*u,num2str(k_element),'color',light_blue,'fontsize',6,'horizontalAlignment','center');        
    end
    % ----- Node number ploting -------------------------------------------
    for k=1:n_node
        hold on
        plot(node(k,1)*u,node(k,2)*u,'o','markersize',15,'markerfacecolor','y','markeredgecolor','y');
        text(node(k,1)*u,node(k,2)*u,num2str(k),'color',light_blue,'fontsize',6,'horizontalAlignment','center');
    end
end
xtick1=[0:Lx/5:Lx]*u;
set(gca,'linewidth',2,'fontsize',25,'xtick',xtick1); % ,'ytick',0:100:500);
xlabel('length [nm]');
ylabel('width [nm]');
axis([0 Lx*u 0 Ly*u]);
axis equal;
%% ----- Porous structure ploting -----------------------------------------
if strcmpi(include_pores,'Yes') || strcmpi(include_gb_pores,'Yes')
    %----- Pore ploting -------------------------------------------------------
    for k=1:n_pore
        poly_pore=polyshape(x_pore(k,:)*u,y_pore(k,:)*u);
        plot(poly_pore,'linestyle','none','linewidth',0.01,'facecolor',[0.4 0.6 1],'facealpha',0.4,'edgecolor',[0.4 0.6 1]);
        plot(C(k,1)*u,C(k,2)*u,'.k','markersize',5);
        plot(x_pore(k,:)*u,y_pore(k,:)*u,'w','linewidth',1,'marker','.','markersize',2);
%         text(C(k,1)*u,C(k,2)*u,num2str(k),'color','k','fontsize',10,'horizontalAlignment','center');
    end
    %----- Element re-ploting found in pore -------------------------------
    if n_element<element_limit
        air_force_blue=[1 0 0];%[0.36 0.54 0.66];
        plot_element_in_pore=unique(element_in_pore);
        for k_element=1:n_element %length(plot_element_in_pore)
            if (find(plot_element_in_pore==k_element))
                plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':', 'color',air_force_blue);
                plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':', 'color',air_force_blue);
                plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':', 'color',air_force_blue);
                plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':', 'color',air_force_blue);
                x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
                y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
                if any(C_element==k_element)
                    text(x*u,y*u,num2str(k_element),'color','k','fontsize',6,'horizontalAlignment','center');
                else
                    text(x*u,y*u,num2str(k_element),'color',[0.9 0.2 0.2],'fontsize',6,'horizontalAlignment','center');
                end
            end
            %----- Re-ploting center element ------------------------------
            % if any(C_element==k_element) 
            %     pt1=node(element(k_element,1),:);
            %     pt2=node(element(k_element,2),:);
            %     pt3=node(element(k_element,3),:);
            %     pt4=node(element(k_element,4),:);
            %     C_element=polyshape([pt1(1) pt2(1) pt3(1) pt4(1)]*u, [pt1(2) pt2(2) pt3(2) pt4(2)]*u);
            %     plot(C_element,'facecolor','y','faceAlpha',0.2);
            % end
        end
    end
end
%% ----- GB ploting -------------------------------------------------------
if strcmpi(include_random_gb,'Yes')
    %----- Ploting grain boundries ----------------------------------------
    gbc=[0.12 0.2 0.73];
    for k=1:length(grain_line_p1)
        plot([grain_line_p1(k,1) grain_line_p2(k,1)]*u, [grain_line_p1(k,2) grain_line_p2(k,2)]*u,'linestyle','-','linewidth',2,'color',gbc);
    end
    %----- Re-ploting elements found at grain boundary --------------------
    if n_element<element_limit
        re_element=[0.36 0.54 0.66];
        plot_element_in_grain=unique([element_in_grain{:}]);
        for k_element=1:n_element
            if (find(plot_element_in_grain==k_element))
                plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':', 'color',re_element);
                x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
                y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
        %         text(x*u,y*u,num2str(k_element),'color','r','fontsize',10,'horizontalAlignment', 'center');        
            end
        end
    end
end
%% ----- GB poly strcture ploting -----------------------------------------
if strcmpi(include_ordered_poly,'Yes')
    %----- Ploting polygon face -------------------------------------------
    for k=1:n_poly
        poly_shape=polyshape(squeeze(polyVertices(k,:,1))*u,squeeze(polyVertices(k,:,2))*u);
        plot(poly_shape,'linestyle','none','facecolor',[0.1 0.1 1],'facealpha',0.6,'edgecolor',[0.4 0.6 1],'linewidth',0.01)
    end
    %----- Re-ploting elements found in polygon ---------------------------
    if n_element<element_limit
        re_element=[0.36 0.54 0.66];
        for k_element=1:n_element 
            if (find(element_in_poly==k_element))
                plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':', 'color',re_element);
                x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
                y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
            end
        end
    end
end
%% ----- Ordered grain lines ploting --------------------------------------
if strcmpi(include_ordered_gb,'Yes')
    %----- Ploting grain lines --------------------------------------------
    for k=1:length(grain_line_p1)
        plot([grain_line_p1(k,1) grain_line_p2(k,1)]*u, [grain_line_p1(k,2) grain_line_p2(k,2)]*u,'--k','linewidth',2);
    end
    %----- Re-ploting elements found at grain boundary --------------------
    if n_element<element_limit
        re_element=[0.36 0.54 0.66];
        plot_element_in_grain=unique([element_in_grain{:}]);
        for k_element=1:n_element
            if (find(plot_element_in_grain==k_element))
                plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':', 'color',re_element);
                plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':', 'color',re_element);
                x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
                y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
            end
        end
    end
end
%% ----- GB strcture ploting ----------------------------------------------
% if strcmpi(include_gb,'Yes') || strcmpi(include_gb_pores,'Yes')
%     %----- Grain Numbering ----------------------------------------------------
%     disp_n_grain='Yes';
%     if strcmpi(disp_n_grain,'Yes')
%         for k=1:n_grain
%             plot(grain_center_x(k)*u,grain_center_y(k)*u,'s','markersize',15,'markerfacecolor','c','markeredgecolor','k');
%             text(grain_center_x(k)*u,grain_center_y(k)*u,num2str(k),'color','k','fontsize',8,'horizontalAlignment','center');
%         end
%     end
%     %----- Grain Line Ploting -------------------------------------------------
%     for k=1:length(grain_line_p1)
%         hold on
%         plot([grain_line_p1(k,1) grain_line_p2(k,1)]*u, [grain_line_p1(k,2) grain_line_p2(k,2)]*u,'-b','linewidth',2);
%     %     plot(grain_line_xsegments{k}*u, grain_line_ysegments{k}*u,'.y','markersize',12);
%     %     line_x=(grain_line_p1(k,1)+grain_line_p2(k,1))/2; line_y=(grain_line_p1(k,2)+grain_line_p2(k,2))/2;
%     %     text(line_x*u,line_y*u,num2str(k),'color','k','fontsize',20,'horizontalAlignment','center');
%     end
%     %----- Element re-ploting found at grain boundary -------------------------
%     if n_grain<50
%         air_force_blue=[1 0 0];
%         plot_element_in_grain=unique([element_in_grain{:}]);
%         for k_element=1:n_element %length(plot_element_in_grain)
%             if (find(plot_element_in_grain==k_element))
%                 plot([node(element(k_element,1),1),node(element(k_element,2),1)]*u, [node(element(k_element,1),2),node(element(k_element,2),2)]*u,':', 'color',air_force_blue);
%                 plot([node(element(k_element,2),1),node(element(k_element,3),1)]*u, [node(element(k_element,2),2),node(element(k_element,3),2)]*u,':', 'color',air_force_blue);
%                 plot([node(element(k_element,3),1),node(element(k_element,4),1)]*u, [node(element(k_element,3),2),node(element(k_element,4),2)]*u,':', 'color',air_force_blue);
%                 plot([node(element(k_element,4),1),node(element(k_element,1),1)]*u, [node(element(k_element,4),2),node(element(k_element,1),2)]*u,':', 'color',air_force_blue);
%                 x=(node(element(k_element,1),1)+node(element(k_element,2),1)+node(element(k_element,3),1)+node(element(k_element,4),1))/4;
%                 y=(node(element(k_element,1),2)+node(element(k_element,2),2)+node(element(k_element,3),2)+node(element(k_element,4),2))/4;
%         %         text(x*u,y*u,num2str(k_element),'color','r','fontsize',10,'horizontalAlignment', 'center');        
%             end
%         end
%     end
% end

%--------------------------------------------------------------------------
