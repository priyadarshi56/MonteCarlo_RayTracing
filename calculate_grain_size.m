function [avg_grain_size] = calculate_grain_size(grain_center_x,grain_center_y)
    grain_center=[grain_center_x' grain_center_y'];
    dtt=delaunayTriangulation(grain_center);
    connecting_lines=size(dtt);
    connecting_line_length=zeros(connecting_lines);
    for k=1:connecting_lines(1)
        vertex_number=dtt(k,:);
        vertex_cordinate1=grain_center(vertex_number(1),:);
        vertex_cordinate2=grain_center(vertex_number(2),:);
        vertex_cordinate3=grain_center(vertex_number(3),:);
        connecting_line_length(k,1)=sqrt((vertex_cordinate1(1)-vertex_cordinate2(1))^2+(vertex_cordinate1(2)-vertex_cordinate2(2))^2);
        connecting_line_length(k,2)=sqrt((vertex_cordinate2(1)-vertex_cordinate3(1))^2+(vertex_cordinate2(2)-vertex_cordinate3(2))^2);
        connecting_line_length(k,3)=sqrt((vertex_cordinate3(1)-vertex_cordinate1(1))^2+(vertex_cordinate3(2)-vertex_cordinate2(2))^2);
    end
    avg_grain_size=mean(mean(connecting_line_length));
end