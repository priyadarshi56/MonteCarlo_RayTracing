function [grain_center_x,grain_center_y,grain_line_xsegments,grain_line_ysegments] = create_grain_boundary(domain,dLx,dLy,n_grain)
%--------- Create Grain Boundary ------------------------------------------
%     n_grain=20;
    Lx=domain(3,1);
    Ly=domain(3,2);
    grain_center_x=rand(1,n_grain)*Lx;
    grain_center_y=rand(1,n_grain)*Ly;
    [vx,vy]=voronoi(grain_center_x,grain_center_y);
    for k=1:length(vx)
        p1=[vx(1,k) vy(1,k)];
        p2=[vx(2,k) vy(2,k)];
        domain1=[domain; domain(1,:)];
        p1_in=inpolygon(p1(1,1),p1(1,2),domain(:,1),domain(:,2));
        p2_in=inpolygon(p2(1,1),p2(1,2),domain(:,1),domain(:,2)); 
        if (p1_in==1 && p2_in==1)
            grain_line_p1(k,:)=p1;
            grain_line_p2(k,:)=p2;
    %         fprintf(1,'Line %g completely INSIDE of the domain. \n',k);
        elseif (p1_in==1 && p2_in==0)
            grain_line_p1(k,:)=p1;
            x1=p1(1);y1=p1(2);x2=p2(1);y2=p2(2);
            for k_check=1:length(domain)
                x3=domain1(k_check,1);y3=domain1(k_check,2);
                x4=domain1(k_check+1,1);y4=domain1(k_check+1,2);
                det1=(y4-y3)*(x4-x1)-(x4-x3)*(y4-y1);
                det2=(y4-y3)*(x4-x2)-(x4-x3)*(y4-y2);
                det3=(y2-y1)*(x2-x3)-(x2-x1)*(y2-y3);
                det4=(y2-y1)*(x2-x4)-(x2-x1)*(y2-y4);
                if ((det1*det2<=0) && (det3*det4<=0)) % Condition for intersection
                    m1=(y1-y2)/(x1-x2); % Slope of line 1
                    m2=(y3-y4)/(x3-x4); % Slope of line 2 (domain)
                    intrcpt1=y1-m1*x1; % Intercept of line 1
                    if (abs(m2)==Inf)
                        x0=x3; % Intercept of line 2 (domain); also abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    else                    
                        intrcpt2=y3-m2*x3; % Intercept of line 2 (domain)
                        x0=(intrcpt2-intrcpt1)/(m1-m2); % Abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    end
    %                 fprintf(1,'Line %g intersect with domain line %g at (%g,%g). \n',k,k_check,x0*u,y0*u);
                    grain_line_p2(k,:)=[x0 y0];
                    break;
                elseif (k_check>=length(domain))
    %                 fprintf(1,'Line %g NOT intersect. \n',k);
                end
            end
        elseif (p1_in==0 && p2_in==1)
            grain_line_p2(k,:)=p2;
            x1=p1(1);y1=p1(2);x2=p2(1);y2=p2(2);
            for k_check=1:length(domain)
                x3=domain1(k_check,1);y3=domain1(k_check,2);
                x4=domain1(k_check+1,1);y4=domain1(k_check+1,2);
                det1=(y4-y3)*(x4-x1)-(x4-x3)*(y4-y1);
                det2=(y4-y3)*(x4-x2)-(x4-x3)*(y4-y2);
                det3=(y2-y1)*(x2-x3)-(x2-x1)*(y2-y3);
                det4=(y2-y1)*(x2-x4)-(x2-x1)*(y2-y4);
                if ((det1*det2<=0) && (det3*det4<=0)) % Condition for intersection
                    m1=(y1-y2)/(x1-x2); % Slope of line 1
                    m2=(y3-y4)/(x3-x4); % Slope of line 2 (domain)
                    intrcpt1=y1-m1*x1; % Intercept of line 1
                    if (abs(m2)==Inf)
                        x0=x3; % Intercept of line 2 (domain); also abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    else                    
                        intrcpt2=y3-m2*x3; % Intercept of line 2 (domain)
                        x0=(intrcpt2-intrcpt1)/(m1-m2); % Abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    end
    %                 fprintf(1,'Line %g intersect with domain line %g at (%g,%g). \n',k,k_check,x0*u,y0*u);
                    grain_line_p1(k,:)=[x0 y0];
                    break;
                elseif (k_check>=length(domain))
    %                 fprintf(1,'Line %g NOT intersect. \n',k);
                end
            end
        elseif (p1_in==0 && p2_in==0)
            x1=p1(1);y1=p1(2);x2=p2(1);y2=p2(2);
            k_pt=1;
            for k_check=1:length(domain)
                x3=domain1(k_check,1);y3=domain1(k_check,2);
                x4=domain1(k_check+1,1);y4=domain1(k_check+1,2);
                det1=(y4-y3)*(x4-x1)-(x4-x3)*(y4-y1);
                det2=(y4-y3)*(x4-x2)-(x4-x3)*(y4-y2);
                det3=(y2-y1)*(x2-x3)-(x2-x1)*(y2-y3);
                det4=(y2-y1)*(x2-x4)-(x2-x1)*(y2-y4);
                if ((det1*det2<=0) && (det3*det4<=0))
                    m1=(y1-y2)/(x1-x2); % Slope of line 1
                    m2=(y3-y4)/(x3-x4); % Slope of line 2 (domain)
                    intrcpt1=y1-m1*x1; % Intercept of line 1
                    if (abs(m2)==Inf)
                        x0=x3; % Intercept of line 2 (domain); also abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    else                    
                        intrcpt2=y3-m2*x3; % Intercept of line 2 (domain)
                        x0=(intrcpt2-intrcpt1)/(m1-m2); % Abscissa of point of intersection
                        y0=m1*x0+intrcpt1; % Ordinate of point of intersection
                    end                
                    if (k_pt==1)
                        grain_line_p1(k,:)=[x0 y0];
                        k_pt=k_pt+1;
                    elseif (k_pt==2)
                        grain_line_p2(k,:)=[x0,y0];
                        break;
                    end
                end
            end
    %         fprintf(1,'Line %g OUTSIDE of the domain. \n',k);
        end
    end
    %----- Eleminating zero points --------------------------------------------
    grain_line_p1=grain_line_p1(any(grain_line_p1,2),:);
    grain_line_p2=grain_line_p2(any(grain_line_p2,2),:);
    %----- Divide a line into segments ----------------------------------------
    for k=1:length(grain_line_p1)
        x1=grain_line_p1(k,1);y1=grain_line_p1(k,2);
        x2=grain_line_p2(k,1);y2=grain_line_p2(k,2);
        r=sqrt((x1-x2)^2+(y1-y2)^2);
        nx=abs(x1-x2)/dLx;
        ny=abs(y1-y2)/dLy;
        n=max(nx,ny);
    %     m1=(y1-y2)/(x1-x2);
    %     n=10; % No. of segments
        deltax=(x2-x1)/n;
        deltay=(y2-y1)/n;
        halfdx=deltax/2;
        halfdy=deltay/2;
        xsegments=x1+halfdx+(0:n-1)*deltax;
        ysegments=y1+halfdy+(0:n-1)*deltay;
        grain_line_xsegments{k}=[grain_line_p1(k,1) xsegments grain_line_p2(k,1)];
        grain_line_ysegments{k}=[grain_line_p1(k,2) ysegments grain_line_p2(k,2)];
    end
end

