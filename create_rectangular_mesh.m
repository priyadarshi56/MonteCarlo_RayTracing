function [n_node,node,n_element,element,x_cord,y_cord,Npx,Npy] = create_rectangular_mesh(Lx,Ly,dLx,dLy)
    Lx_array=0:dLx:Lx;
    Ly_array=0:dLy:Ly;
    Npx=length(Lx_array); % # of node points in the x-direction
    Npy=length(Ly_array); % # of node points in the y-direction
    domain=[0 0;Lx 0;Lx Ly;0 Ly]; % Domain coordinates
    dim=sqrt(length(domain)); % dimension of domain (e.g 2D or 3D)
    x_cord=zeros(Npy,Npx); % node's x-coordinate matrix
    y_cord=zeros(Npy,Npx); % node's y-coordinate matrix
    n_node=Npy*Npx;
    node=zeros(n_node,dim);
    k_node=0;
    for ky=1:Npy
        for kx=1:Npx
            x_cord(ky,kx)=Lx_array(kx);
            y_cord(ky,kx)=Ly_array(ky);
            k_node=k_node+1;
            node(k_node,:)=[Lx_array(kx),Ly_array(ky)];
        end
    end
    element=zeros((Npx-1)*(Npy-1),4);
    n_element=0;
    for ky=1:Npy-1
        for kx=1:Npx-1
            node1=(ky-1)*Npx+kx;
            node2=(ky-1)*Npx+kx+1;
            node3=ky*Npx+kx+1;
            node4=ky*Npx+kx;
            n_element=n_element+1;
            element(n_element,:)=[node1, node2, node3, node4];
        end
    end
end

