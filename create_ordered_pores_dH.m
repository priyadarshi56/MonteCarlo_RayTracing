function [C,rad,x_pore,y_pore,x_pore_both,y_pore_both] = create_ordered_pores_dH(Lx,Ly,n_pore_x,n_pore_y,pore_r,dH_thick)
    n_pore=n_pore_x*n_pore_y;
    global nc
%     nc=30; % # of nodes on a pore
    C=zeros(n_pore,2);
%     x_pore=zeros(nc*n_pore,1);
%     y_pore=zeros(nc*n_pore,1);
    k=1; % Counter
    for k_dLy=1:n_pore_y
        for k_dLx=1:n_pore_x
            C(k,:)=[Lx/n_pore_x*(k_dLx-1/2)  Ly/n_pore_y*(k_dLy-1/2)];  % Center of circle
            rad(k)=pore_r;
            k=k+1;
        end
    end
    theta=linspace(0,360,nc); % linspace(0,2*pi);
    x_pore=C(:,1)+pore_r*cosd(theta); 
    y_pore=C(:,2)+pore_r*sind(theta);
    inner_pore_r=pore_r-dH_thick;
    x_pore_both=C(:,1)+inner_pore_r*cosd(theta); 
    y_pore_both=C(:,2)+inner_pore_r*sind(theta);
end
