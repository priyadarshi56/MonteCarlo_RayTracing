function [xr,yr,vx,vy,reflect,transmit] = grain_boundary_interaction(xi1,yi1,xi2,yi2,grain_line_p1,grain_line_p2,dt_refl)
%     x1=grain_line_p1(1);y1=grain_line_p1(2);
%     x2=grain_line_p2(1);y2=grain_line_p2(2);
%----- Tangent calculation ------------------------------------------------
    T1=grain_line_p1;
    T2=grain_line_p2;
    mt=(T1(2)-T2(2))/(T1(1)-T2(1)); % Slope of grain boundary line/tangent in radian
    T1T2_angle=atand(mt*180/pi); % Slope of grain boundary line/tangent in degree
    T1T2_len=sqrt((T1(1)-T2(1))^2+(T1(2)-T2(2))^2); % norm(T1T2); % Length of tangent
%----- Incident ray -------------------------------------------------------
    Inc1=[xi1, yi1];
    Inc2=[xi2, yi2];
    mi=(yi2-yi1)/(xi2-xi1); % Slope of incident ray
    %----- Random decision on reflection & transmission -------------------
    sel=rand;
    if (sel<=0.5) % Transmission >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        transmit=1; reflect=0;
        xr=xi2+xi2-xi1;
        yr=yi2+yi2-yi1;
        vx=(xi2-xi1)/dt_refl;
        vy=(yi2-yi1)/dt_refl;
    elseif (sel>0.5) % Reflection >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        reflect=1; transmit=0;
%----- Finding perpendicular to tangent -----------------------------------
        mp=-1/mt; % Slope of normal
        P1=[xi2 yi2]; % First point on perpendicular
        P1P2_incpt=P1(2)-mp*P1(1); % Y-intercept of perpendicular
        P2(1)=P1(1)-(T1T2_len/4)/sqrt(1+mp^2); % Abscissa of point2 of perpendicular
        P2(2)=mp*P2(1)+P1P2_incpt; % Ordinate of point2 of perpendicular
        xp1=P1(1);yp1=P1(2);xp2=P2(1);yp2=P2(2);
        mp1=(yp1-yp2)/(xp1-xp2); %
        P1P2_angle=atand(mp*180/pi); % Slope of grain boundary line/tangent in degree
%----- Reflected ray ------------------------------------------------------
        theta_i=atand(abs((mi-mp)/(1+mi*mp))); % Angle of incidence
        xir_mid=(xi1+mp^2*xp1+mp*yi1-mp*yp1)/(1+mp^2); % Mid point of line joining i & r
        yir_mid=yp1+mp*(xir_mid-xp1); % Mid point of line joining i & r
        xr=2*xir_mid-xi1; % Abscissa of reflected point
        yr=2*yir_mid-yi1; % Ordinate of reflected point
        len_r=sqrt((xr-xi2)^2+(yr-yi2)^2); % Length of reflected ray
        mr=(yr-yi2)/(xr-xi2); % Slope of reflected ray
%----- Velocity -----------------------------------------------------------
       vx=(xr-xi2)/dt_refl;
       vy=(yr-yi2)/dt_refl;
    end
end


%----- Tangent calculation from differential dy/dx ------------------------
%     T1=grain_line_p1;
%     T2=grain_line_p2;
%     m_grain=(T1(2)-T2(2))/(T1(1)-T2(1)); % Slope of grain boundary line/tangent
%     T1T2_angle=atand(m_grain*180/pi); % Slope of grain boundary line/tangent in degree
%     %----- Finding perpendicular to tangent -----------------------------------
%     T1T2=T2-T1; % Tangent vector from T2 to T1
%     T1T2_len=norm(T1T2); % Length of perpendicular = Length of tangent
%     T1T2=T1T2/T1T2_len; % Normalizing T1T2 to have unit length
%     T1T2_perp=T1T2*[0 -1;1 0]; % Perpendicular vector to the tangent vector
%     T1T2_mid=(T1+T2)/2; % Mid point of T1 and T2
%     P1=T1+T1T2_mid/2*T1T2_perp; % First point on perpendicular
%     P2=T1-T1T2_mid/2*T1T2_perp; % Second point on perpendicular
%     xp1=P1(1);yp1=P1(2);xp2=P2(1);yp2=P2(2);
%     mp=(yp2-yp1)/(xp2-xp1);
%     % mp=-1/mt; % Slope of normal 
%     mp*mt;
%     %----- Reflected ray ------------------------------------------------------
%     theta_i=atand(abs((mi-mp)/(1+mi*mp))); % Angle of incidence
%     xir_mid=(xi1+mp^2*xp1+mp*yi1-mp*yp1)/(1+mp^2); % Mid point of line joining i & r
%     yir_mid=yp1+mp*(xir_mid-xp1); % Mid point of line joining i & r
%     xr=2*xir_mid-xi1; % Reflected point
%     yr=2*yir_mid-yi1; % Reflected point
%     len_r=sqrt((xr-xi2)^2+(yr-yi2)^2); % Length of reflected ray
%     mr=(yr-yi2)/(xr-xi2); % Slope of reflected ray
%     %----- Velocity
%     r_angle=mi;
%     vx=(xr-xi2)/dt;
%     vy=(yr-yi2)/dt;
