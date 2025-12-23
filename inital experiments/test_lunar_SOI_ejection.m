format compact
clear
clc
clf reset

%----------

error_tolerance = 1e-11;
unit_vec_size = 5e7;

G = 6.6743e-11;

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; %m

mu = G*earth_mass;

lunar_orbital_radius = 	363300e3; %m
lunar_eccentricity = 0.0549;
longitude_of_ascending_node = 0;
argument_of_perihelion = 0;
inclination = 0; %no inc we're in the lunar plane

r_initial = [
cosd(longitude_of_ascending_node), -sind(longitude_of_ascending_node), 0
sind(longitude_of_ascending_node), cosd(longitude_of_ascending_node), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(inclination), -sind(inclination)
0, sind(inclination), cosd(inclination)
]*[
cosd(argument_of_perihelion), -sind(argument_of_perihelion), 0
sind(argument_of_perihelion), cosd(argument_of_perihelion), 0
0, 0, 1
]*[0; lunar_orbital_radius; 0];

semi_major_axis = norm(r_initial)/(1-lunar_eccentricity);
inital_velocity = sqrt((mu/semi_major_axis)*((1+lunar_eccentricity)/(1-lunar_eccentricity)));

r_dot_initial = [
cosd(longitude_of_ascending_node), -sind(longitude_of_ascending_node), 0
sind(longitude_of_ascending_node), cosd(longitude_of_ascending_node), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(inclination), -sind(inclination)
0, sind(inclination), cosd(inclination)
]*[
cosd(argument_of_perihelion), -sind(argument_of_perihelion), 0
sind(argument_of_perihelion), cosd(argument_of_perihelion), 0
0, 0, 1
]*[-inital_velocity; 0; 0];

state_initial = [r_initial(1:3).',r_dot_initial(1:3).'];
orbital_period = sqrt(((semi_major_axis^3)/mu)*(2*pi)^2);

orbit_timespan = linspace(0,orbital_period * (1+2e-3),6e3);

[lunar_timerange, lunar_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),orbit_timespan,state_initial,odeset('Reltol',error_tolerance));
lunar_orbit_statematrix = [lunar_timerange,lunar_state_matrix];

lunar_SOI_radius = semi_major_axis*(1-lunar_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);


%% ejecting out of lunar SOI

SOI_ejection_longitude = 0; %degrees, parameter

r_SOI_eject = [
cosd(SOI_ejection_longitude), -sind(SOI_ejection_longitude), 0
sind(SOI_ejection_longitude), cosd(SOI_ejection_longitude), 0
0, 0, 1    
]*[0;-lunar_SOI_radius;0];

v_SOI_exit = [20, 0, 0]; %test case

SOI_eject_semimajor = (lunar_SOI_radius + lunar_radius) / 2;
SOI_eject_raw_period = sqrt(( ( SOI_eject_semimajor^3 ) / (lunar_mass*G) ) * (2*pi)^2);

state_initial_SOIeject = [(r_SOI_eject).',v_SOI_exit];
orbit_timespan_SOIeject = linspace(SOI_eject_raw_period,0,5e3);

[timerange_SOIeject, state_matrix_SOIeject] = ode45(@(timerange_SOIeject, state_matrix_SOIeject) g_dynamics_twobody(timerange_SOIeject, state_matrix_SOIeject, lunar_mass), orbit_timespan_SOIeject, state_initial_SOIeject, odeset('Reltol',error_tolerance));
SOIeject_orbit_statematrix = [timerange_SOIeject,state_matrix_SOIeject];

eject_good = false;
altitude_current = nan;
altitude_prev = nan;
for n=2:height(state_matrix_SOIeject)
    lunar_altitude_SOIeject = norm(SOIeject_orbit_statematrix(n-1,2:4));
    altitude_current = lunar_altitude_SOIeject;
    lunar_altitude_series(n-1) = altitude_current;
    if lunar_altitude_SOIeject < lunar_radius
        eject_good = true;
        ind_launch = n;
        break
    end
    altitude_prev = altitude_current;
end

if sum(lunar_altitude_series > lunar_SOI_radius) > 10
    %patched conic approximation does not apply if we exit SOI during transfer
    eject_good = false;
end

if eject_good
    SOIeject_orbit_statematrix(ind_launch:end,:) = [];

    surface_interp = interp1([altitude_current,altitude_prev],[0,1],lunar_radius);
    surface_eject_time = interp1([0,1],[SOIeject_orbit_statematrix(end,1),SOIeject_orbit_statematrix(end-1,1)],surface_interp);
    for n=2:7
        surface_eject_state(n) = interp1([0,1],[SOIeject_orbit_statematrix(end,n),SOIeject_orbit_statematrix(end-1,n)],surface_interp);
    end
    
    SOIeject_orbit_statematrix(:,1) = SOIeject_orbit_statematrix(:,1) - surface_eject_time;
    SOIeject_orbit_statematrix(end,:) = surface_eject_state;
    
    SOIeject_time = SOIeject_orbit_statematrix(1,1);

    SOIeject_orbit_statematrix = flipud(SOIeject_orbit_statematrix);
    %add on lunar influence so we're in the earth reference frame
    for n=1:height(SOIeject_orbit_statematrix)
        frametransform_inds = find_closest_inds(lunar_orbit_statematrix(:,1),SOIeject_orbit_statematrix(n,1));
        lunar_pos_spec = interp_statematrix_timebetween(lunar_orbit_statematrix, SOIeject_orbit_statematrix(n,1), frametransform_inds);
        SOIeject_orbit_statematrix(n,:) = SOIeject_orbit_statematrix(n,:) + lunar_pos_spec;
    end

    inds_lunarstate_ejecttime = find_closest_inds(lunar_orbit_statematrix(:,1),SOIeject_time);
    SOIeject_lunar_statematrix = interp_statematrix_timebetween(lunar_orbit_statematrix,SOIeject_time,inds_lunarstate_ejecttime);
else
    SOIeject_orbit_statematrix = repelem(nan,7);
    surface_eject_state = repelem(nan,7);
    SOIeject_lunar_statematrix = repelem(nan,7);
end


%% target MEO
MEO_radius = 5e7 + earth_radius
MEO_target_velocity = sqrt(mu/MEO_radius);
MEO_target_TA = 80; %degrees, parameter

target_r = [
cosd(MEO_target_TA), -sind(MEO_target_TA), 0
sind(MEO_target_TA), cosd(MEO_target_TA), 0
0, 0, 1
]*[0; -MEO_radius; 0];
target_r = target_r.';

SOIearth_inital_r = SOIeject_orbit_statematrix(end,2:4);

target_v = [
cosd(MEO_target_TA), -sind(MEO_target_TA), 0
sind(MEO_target_TA), cosd(MEO_target_TA), 0
0, 0, 1
]*[MEO_target_velocity; 0; 0];
target_v = target_v.';

[v_1, v_2, delta_t, dv_match] = construct_transfer(SOIearth_inital_r, target_r, 0.9, mu, target_v); %dt will be a param

state_initial_earthtransit = [SOIearth_inital_r,v_1];
orbit_timespan_earthtransit = linspace(0,delta_t,5e3);

if eject_good

    [timerange_earthtransit, state_matrix_earthtransit] = ode45(@(timerange_earthtransit, state_matrix_earthtransit) g_dynamics_twobody(timerange_earthtransit, state_matrix_earthtransit, earth_mass), orbit_timespan_earthtransit, state_initial_earthtransit, odeset('Reltol',error_tolerance));
    earthtransit_statematrix = [timerange_earthtransit,state_matrix_earthtransit];
    
    SOI_switchover_v_error = norm(earthtransit_statematrix(1,5:7) - SOIeject_orbit_statematrix(end,5:7))

    %find thrust req
    impulse_assumption_angle = 5;
    r_final = earthtransit_statematrix(end,2:4);
    impulse_assumption_time = nan;
    for n=height(earthtransit_statematrix):-4:1
        r_spec = earthtransit_statematrix(n,2:4);
        theta_r = acosd(dot(r_final,r_spec) / (norm(r_final)*norm(r_spec)));
        if theta_r > impulse_assumption_angle
            impulse_assumption_time = earthtransit_statematrix(end,1) - earthtransit_statematrix(n,1);
            break
        end
    end

    impulse_assumption_time
    dv_match
end

%% plots

hold on
grid on
axis equal padded
ax = gca;
ax.FontSize = 20;
set(gcf, 'Color', [1,1,1])

%plot earth
earth_map = imread('earth_map.jpg');
earth_map = flipud(earth_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(0,0,0,earth_radius,earth_radius,earth_radius,70);
earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(earth_obj,'cdata',earth_map,'facecolor','texturemap')

%plot luna
%luna at release
lunar_map = imread('lunar_map.jpg');
lunar_map = flipud(lunar_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(r_initial(1),r_initial(2),r_initial(3),lunar_radius,lunar_radius,lunar_radius,80);
lunar_obj_1 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj_1,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)
%luna at SOI eject
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(SOIeject_lunar_statematrix(2),SOIeject_lunar_statematrix(3),SOIeject_lunar_statematrix(4),lunar_radius,lunar_radius,lunar_radius,80);
lunar_obj_2 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj_2,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)


%lunar orbit
plot3(lunar_state_matrix(:,1),lunar_state_matrix(:,2),lunar_state_matrix(:,3),"k")

%surface ejection point
scatter3(SOIeject_orbit_statematrix(1,2),SOIeject_orbit_statematrix(1,3),SOIeject_orbit_statematrix(1,4),30,"k","filled")

%SOI ejection path
plot3(SOIeject_orbit_statematrix(:,2),SOIeject_orbit_statematrix(:,3),SOIeject_orbit_statematrix(:,4),"-b")
scatter3(SOIeject_orbit_statematrix(end,2),SOIeject_orbit_statematrix(end,3),SOIeject_orbit_statematrix(end,4),30,"k","filled")

%lunar SOI
plot(nsidedpoly(1000, 'Center', [r_initial(1),r_initial(2)], 'Radius', lunar_SOI_radius), 'FaceColor', [repelem(0.9,3)])
plot(nsidedpoly(1000, 'Center', [SOIeject_lunar_statematrix(2),SOIeject_lunar_statematrix(3)], 'Radius', lunar_SOI_radius), 'FaceColor', [repelem(0.9,3)])
% xlim([r_initial(1)-1e7,r_initial(1)+1e7])
% ylim([r_initial(2)-1e7,r_initial(2)+1e7])

%target orbit
plot(nsidedpoly(1000, 'Center', [0,0], 'Radius', MEO_radius),FaceAlpha=0)
scatter3(target_r(1),target_r(2),target_r(3),"k","filled")
target_v_unit = target_v./norm(target_v);
quiver3(target_r(1),target_r(2),target_r(3), target_v_unit(1)*unit_vec_size,target_v_unit(2)*unit_vec_size,target_v_unit(3)*unit_vec_size,"k")

%earth transit
plot3(earthtransit_statematrix(:,2),earthtransit_statematrix(:,3),earthtransit_statematrix(:,4),"r")

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')


function closest_inds = find_closest_inds(a,b)
    [~,closest_inds] = mink( abs(a-b), 2);
    closest_inds = sort(closest_inds).';
end

function state_out = interp_statematrix_timebetween(statematrix,time,between_inds)
    for n=2:width(statematrix)
        state_out(n) = interp1([statematrix(between_inds(1),1),statematrix(between_inds(2),1)], [statematrix(between_inds(1),n),statematrix(between_inds(2),n)], time,"linear","extrap");
    end
end


function output = g_dynamics_twobody(t,state,primary_mass)
    G = 6.67e-11;
    velocity = state(4:6);
    r = state(1:3); 
    r_vector = norm(r);
    r_unit_vector = r/r_vector;
    force_gravity = r_unit_vector*(-G*primary_mass/(r_vector^2));
    acceleration = force_gravity;
    output = [velocity; acceleration];
end