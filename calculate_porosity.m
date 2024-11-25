function [porosity] = calculate_porosity(r,domain_area)
    pore_area=zeros(1,length(r));
    for k=1:length(r)
        pore_area(k)=pi*r(k)^2;
    end
    pores_area=sum(pore_area);
    porosity=(pores_area/domain_area)*100;
end