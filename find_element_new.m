function [found_element] = find_element_new(Lx,Ly,Npx,Npy,node,x,y)
    % x=abs(x_covered+x_flight); % rand*Lx;
    % y=abs(y_covered+y_flight); % rand*Ly;
    for k_dLx=1:Npx-1
        if x<0
            kx_element=1;
            break;
        elseif x>Lx
            kx_element=Npx-1;
            break;
        elseif (x>=node(k_dLx,1) && (x<=node(k_dLx+1,1)))
            kx_element=k_dLx;
            break;
        end
    end
    for k_dLy=1:Npy-1
        if y<0
            ky_element=1;
            break;
        elseif y>Ly
            ky_element=Npy-1;
            break;
        elseif (y>=node((k_dLy-1)*Npx+1,2) && (y<=node(k_dLy*Npx+1,2)))
            ky_element=k_dLy;
            break;
        end
    end
    found_element=kx_element+(Npx-1)*(ky_element-1);
end  
