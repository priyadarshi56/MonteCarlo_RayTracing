function [xr,yr,vx,vy] = pore_reflection(xi1,yi1,xi2,yi2,dt_refl,C_node_vertex)
    % xi1=x_step(n_step-1); yi1=y_step(n_step-1); xi2=x_step(n_step); yi2=y_step(n_step);
    % dt_refl; C_node_vertex;
    global nc;
    C_refl(1)=(C_node_vertex(1,1)+C_node_vertex(3,1))/2; % x-center of element
    C_refl(2)=(C_node_vertex(1,2)+C_node_vertex(3,2))/2; % y-center of element
%----- Tangent calculation from nearest point on the pore -------------------     
    xt1=xi2; yt1=yi2; % first point on tangent
    r=sqrt((C_refl(1)-xt1)^2+(C_refl(2)-yt1)^2); % radius
    theta=linspace(0,360,nc*2);
    x_pore=C_refl(1)+r*cosd(theta);
    y_pore=C_refl(2)+r*sind(theta);    
    distance=sqrt((x_pore-xt1').^2+(y_pore-yt1').^2); % finding nearest point to (xt1,yt1) on the pore
    [~,loc_t2]=min(distance,[],2); % location for the closest point on the pore (i.e. those with min distance)
    xt2=x_pore(loc_t2); % second x-point on tangent
    yt2=y_pore(loc_t2); % second y-point on tangent
    tangent_len=sqrt((xt2-xt1)^2+(yt2-yt1)^2); % length of tangent
    mt=(yt2-yt1)/(xt2-xt1); % slope of tangent (unitless)
    mt_radian=atan(mt); % slope of tangent (radian)
    mt_theta=mt_radian*180/pi; % rad2deg(T1T2_radian); % slope of tangent (degree)
%----- Finding perpendicular to tangent -----------------------------------
    mp=-1/mt; % slope of perpendicular (unitless)
    mp_radian=atan(mp); % slope of perpendicular (radian)
    mp_theta=mp_radian*180/pi; % slope of perpendicular (degree)
    xp1=xi2; yp1=yi2; % first point on perpendicular
    perp_len=tangent_len; % length of perpendicular
    perp_incpt=yp1-mp*xp1; % Y-intercept of perpendicular
    xp2=xp1-perp_len/sqrt(1+mp^2); % absciss of point2 of perpendicular
    yp2=mp*xp2+perp_incpt; % ordinate of point2 of perpendicular 
%----- Incident ray -------------------------------------------------------
    mi=(yi2-yi1)/(xi2-xi1); % slope of incident ray (unitless)
    mi_radian=atan(mi); % slope of incident ray (radian)
    mi_theta=mi_radian*180/pi; % angle of incident ray (degree)
    len_i=sqrt((xi2-xi1)^2+(yi2-yi1)^2); % length of incident ray
%----- Reflected ray ------------------------------------------------------
    % theta_i=atand(abs((mi-mp)/(1+mi*mp))); % angle of incidence
    xir_mid=(xi1+mp^2*xp1+mp*yi1-mp*yp1)/(1+mp^2); % mid point of line joining i & r
    yir_mid=yp1+mp*(xir_mid-xp1); % mid point of line joining i & r
    xr=2*xir_mid-xi1; % x-reflected point
    yr=2*yir_mid-yi1; % y-reflected point
    len_r=sqrt((xr-xi2)^2+(yr-yi2)^2); % length of reflected ray
    % len_i==len_r; % probably be exactly equal
%----- Velocity calculation -----------------------------------------------
    vx=(xr-xi2)/dt_refl;
    vy=(yr-yi2)/dt_refl;
end
