function [element_wise] = matrix_to_element(element,node,Npx,matrix_2D)
    % bilinear interpolation
    n_element=length(element);
    element_center=zeros(n_element,2);
    element_wise=zeros(n_element,1);
    nn=1;
    nnn=0;
    for k=1:n_element
        node1=element(k,1);
        node2=element(k,2);
        node3=element(k,3);
        node4=element(k,4);
        x1=node(node1,1); y1=node(node1,2);
        x2=node(node2,1); y2=node(node2,2);
        x3=node(node3,1); y3=node(node3,2);
        x4=node(node4,1); y4=node(node4,2);
        x_c=(x1+x2)/2;
        y_c=(y1+y4)/2;
        element_center(k,:)=[x_c,y_c];
        cc=(x2-x1)*(y3-y2);
        if k<=nn*(Npx-1)
            ky1=nn; ky2=ky1; ky3=ky1+1; ky4=ky3;
            nnn=nnn+1;
            kx1=nnn; kx2=kx1+1; kx3=kx2; kx4=kx1;
        else
            nn=nn+1;
            ky1=nn; ky2=ky1; ky3=ky1+1; ky4=ky3;
            nnn=1;
            kx1=nnn; kx2=kx1+1; kx3=kx2; kx4=kx1;
        end
        % hold on
        % p1=plot(x_cord(ky1,kx1)*u,y_cord(ky1,kx1)*u,'or');
        % p2=plot(x_cord(ky2,kx2)*u,y_cord(ky2,kx2)*u,'or');
        % p3=plot(x_cord(ky3,kx3)*u,y_cord(ky3,kx3)*u,'or');
        % p4=plot(x_cord(ky4,kx4)*u,y_cord(ky4,kx4)*u,'or');
        phi1=matrix_2D(ky1,kx1);
        phi2=matrix_2D(ky2,kx2);
        phi3=matrix_2D(ky3,kx3);
        phi4=matrix_2D(ky4,kx4);
        element_wise(k)=(1/cc)*(phi1*(x2-x_c)*(y3-y_c)+phi2*(x_c-x1)*(y3-y_c)+phi3*(x_c-x1)*(y_c-y2)+phi4*(x2-x_c)*(y_c-y1)); % bilinear interpolation
        % element_wise(k)=(1/4)*(phi1+phi2+phi3+phi4);
    % delete(p1); delete(p2); delete(p3); delete(p4);
    end
    % x_cord_element_center=reshape(element_center(:,1),Npy-1,Npx-1);
    % y_cord_element_center=reshape(element_center(:,2),Npy-1,Npx-1);
end
