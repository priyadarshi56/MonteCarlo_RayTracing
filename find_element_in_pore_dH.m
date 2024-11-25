function [element_in_pore,node_in_pore,element_in_both,node_in_both,node_in_dH,C_element] = find_element_in_pore_dH(node,element,n_pore,C,x_pore,y_pore,x_pore_both,y_pore_both)
%----- Finding elements having center of pores ----------------------------
    C_pore_in_element=zeros(length(C),1);
    C_pore_in_node=zeros(length(C),4);
    for k_element=1:length(element)
        in_element=inpolygon(C(:,1),C(:,2),node(element(k_element,[1 2 3 4]),1),node(element(k_element,[1 2 3 4]),2));
        C_pore_in_element(in_element)=k_element;
        C_pore_in_node(in_element,:)=element(C_pore_in_element(in_element),:);
    end
%----- Finding elements in all pores --------------------------------------
    k_element_pore=1; % Counter for number of element in a pore
    k_element_in_both=1; % Counter for number of element in both pore
    k_C_element=1; % length(C_pore_in_element);
    % tic;
    for k_n_pore=1:n_pore
        node_in_pore=inpolygon(node(:,1),node(:,2),x_pore(k_n_pore,:),y_pore(k_n_pore,:));
        node_in_pore=find(node_in_pore~=0);   
        for k_node_in_pore=1:length(node_in_pore)            
            for k_element=1:length(element)
                element_in=find(element(k_element,:)==node_in_pore(k_node_in_pore));
                if (element_in~=0)
                    element_in_pore1(k_element_pore)=k_element;
                    C_element1(k_element_pore)=C_pore_in_element(k_n_pore); % Center for corresponding elements
                    k_element_pore=k_element_pore+1;
                end
            end
        end
        k_C_element=k_C_element+1;
%--------------------------------------------------------------------------
        node_in_both=inpolygon(node(:,1),node(:,2),x_pore_both(k_n_pore,:),y_pore_both(k_n_pore,:));
        node_in_both=find(node_in_both~=0);
        for k_node_in_pore=1:length(node_in_both)            
            for k_element=1:length(element)
                element_in_both1=find(element(k_element,:)==node_in_both(k_node_in_pore));
                if (element_in_both1~=0)
                    element_in_both(k_element_in_both)=k_element;
        %         %     C_element(k_element_in_both,1)=C_pore_in_element(k_C_element,1); % Center for corresponding elements
                    k_element_in_both=k_element_in_both+1;
                end
            end
        end
    end
%----- Sorting unique elements and corresponding center of pore -----------
    unique_C_element=unique(C_element1,'stable');
    element_in_pore=[];
    C_element=[];    
    for k=1:length(unique(C_element1,'stable'))
        index=unique_C_element(k)==C_element1;
        element_in_pore=[element_in_pore, unique(element_in_pore1(index))];
        C_element=[C_element, unique_C_element(k)*ones(1,length(unique(element_in_pore1(index))))]; 
    end
%--------------------------------------------------------------------------
%----- Specially for Poisson solver ---------------------------------------
%----- Finding nodes in pores ---------------------------------------------
    sort_element_in_pore=unique(element_in_pore1);
    node_in_pore1=element(sort_element_in_pore(:),:);
    node_in_pore=reshape(node_in_pore1,numel(node_in_pore1),1);
    node_in_pore=unique(sort(node_in_pore));
%----- Finding nodes in both pores ----------------------------------------
    sort_element_in_both=unique(element_in_both);
    node_in_both1=element(sort_element_in_both(:),:);
    node_in_both=reshape(node_in_both1,numel(node_in_both1),1);
    node_in_both=unique(sort(node_in_both));
%----- Finding nodes in the dH ring ---------------------------------------
    node_in_dH=setdiff(node_in_pore,node_in_both);
end



