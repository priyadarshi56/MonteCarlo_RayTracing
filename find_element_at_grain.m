function [element_in_grain,grain_line_slope] = find_element_at_grain(grain_line_xsegments,grain_line_ysegments,Lx,Ly,node,Npx,Npy)
    for k_line=1:length(grain_line_xsegments)
        x=grain_line_xsegments{k_line};
        y=grain_line_ysegments{k_line};
        k_count=1;
        for k_pt=1:length(x)
            [found_element] = find_element_new(Lx,Ly,Npx,Npy,node,x(k_pt),y(k_pt));
            % [found_element,found_node] = find_element_new(node,element,Npx,Npy,x(k_pt),y(k_pt));
            found_element_in_grain(k_count)=found_element;
            m_slope(k_count)=(y(1)-y(end))/(x(1)-x(end)); % Slope of grain boundary
            k_count=k_count+1;
            % node_on_grain(k_count,:)=element(found_element,:); % giving 4 node numbers
        end
%     found_element_in_grain=unique(found_element_in_grain);
    element_in_grain{k_line}=found_element_in_grain;
    grain_line_slope{k_line}=m_slope;
    end
end
