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

orbit_timespan = linspace(0,orbital_period*1.1,3e3);

[lunar_timerange, lunar_state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),orbit_timespan,state_initial,odeset('Reltol',error_tolerance));
lunar_orbit_statematrix = [lunar_timerange,lunar_state_matrix];

lunar_SOI_radius = semi_major_axis*(1-lunar_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);


%% ejecting out of lunar SOI

SOI_ejection_longitude = 125; %degrees, parameter

r_SOI_enter = [
cosd(SOI_ejection_longitude), -sind(SOI_ejection_longitude), 0
sind(SOI_ejection_longitude), cosd(SOI_ejection_longitude), 0
0, 0, 1    
]*[0;-lunar_SOI_radius;0];

r_SOIluna_enter = r_SOI_enter.' + lunar_orbit_statematrix(1,2:4);

%% target MEO
MEO_radius = 5e7 + earth_radius
MEO_target_velocity = sqrt(mu/MEO_radius);
MEO_target_TA = 90; %degrees, parameter

target_r = [
cosd(MEO_target_TA), -sind(MEO_target_TA), 0
sind(MEO_target_TA), cosd(MEO_target_TA), 0
0, 0, 1
]*[0; -MEO_radius; 0];
target_r = target_r.';

target_v = [
cosd(MEO_target_TA), -sind(MEO_target_TA), 0
sind(MEO_target_TA), cosd(MEO_target_TA), 0
0, 0, 1
]*[MEO_target_velocity; 0; 0];
target_v = target_v.';

[v_1, v_2, delta_t, dv_match] = construct_transfer(r_SOIluna_enter, target_r, 0.5, mu, target_v); %dt will be a param

state_initial_earthtransit = [r_SOIluna_enter,v_1];
orbit_timespan_earthtransit = linspace(0,delta_t,1e3);

[timerange_earthtransit, state_matrix_earthtransit] = ode45(@(timerange_earthtransit, state_matrix_earthtransit) g_dynamics_twobody(timerange_earthtransit, state_matrix_earthtransit, earth_mass), orbit_timespan_earthtransit, state_initial_earthtransit, odeset('Reltol',error_tolerance));
earthtransit_statematrix = [timerange_earthtransit,state_matrix_earthtransit];

%find thrust req
impulse_assumption_angle = 10;
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

%making sure the earth orbit doesn't cross the lunar SOI
no_SOI_reentry = true;
for n=1:height(earthtransit_statematrix)
    orbital_r = norm(earthtransit_statematrix(n,2:4));
    if orbital_r > (lunar_orbital_radius - lunar_SOI_radius) || n==1
        lunar_check_statematrix(n,:) = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), earthtransit_statematrix(n,1), earthtransit_statematrix(n,1));
        if norm(lunar_check_statematrix(2:4) - earthtransit_statematrix(n,2:4)) < lunar_SOI_radius
            error("re-enters lunar SOI")
            no_SOI_reentry = false;
            break
        end
    end

end

%% running back in lunar SOI to surface

SOI_enter_stateinital = earthtransit_statematrix(1,2:7) - lunar_check_statematrix(1,2:7);
SOI_backtrack_max_dt = sqrt(( ( (lunar_SOI_radius/2)^3 ) / (lunar_mass*G) ) * (2*pi)^2); %max dt allowed for backsolve back down to lunar surface
SOI_backtrack_timespan = linspace(SOI_backtrack_max_dt,0,5e3);

[timerange_lunarbacktrack, state_matrix_lunarbacktrack] = ode45(@(timerange_lunarbacktrack, state_matrix_lunarbacktrack) g_dynamics_twobody(timerange_lunarbacktrack, state_matrix_lunarbacktrack, lunar_mass), SOI_backtrack_timespan, SOI_enter_stateinital, odeset('Reltol',error_tolerance));
lunarbacktrack_statematrix = [timerange_lunarbacktrack,state_matrix_lunarbacktrack];

%do we hit the moon?
ind_eject = 0;
reaches_surface = false;
orbital_r_prev = 0;
for n=1:height(lunarbacktrack_statematrix)
    orbital_r_lunar = norm(lunarbacktrack_statematrix(n,2:4));
    lunar_altitude_series(n) = orbital_r_lunar;
    if orbital_r_lunar < lunar_radius
        ind_eject = n;
        reaches_surface = true;
        break
    end
    orbital_r_prev = orbital_r_lunar;
end
if reaches_surface
    closest_altitude = 0;
else
    closest_altitude = min(lunar_altitude_series) %if we miss the surface, how close do we get anyway?
end

stays_in_SOI = true;
if reaches_surface
    lunarbacktrack_statematrix(ind_eject+1:end,:) = [];
    surface_interp = interp1([orbital_r_lunar, orbital_r_prev],[0,1],lunar_radius);
    surface_eject_time = interp1([0,1],[lunarbacktrack_statematrix(end,1),lunarbacktrack_statematrix(end-1,1)],surface_interp);
    lunar_eject_statematrix = interp_statematrix_timebetween(lunarbacktrack_statematrix, lunarbacktrack_statematrix(:,1), surface_eject_time, surface_eject_time);
    lunarbacktrack_statematrix(end,1:7) = lunar_eject_statematrix;
    
    %finding exact revolution time
    [~,ind_rev] = mink(abs(lunar_orbit_statematrix(:,2)),5);
    ind_rev = sort(ind_rev);
    lunar_rev_time = interp1([lunar_orbit_statematrix(ind_rev(end-1),2), lunar_orbit_statematrix(ind_rev(end),2)], [lunar_orbit_statematrix(ind_rev(end-1),1), lunar_orbit_statematrix(ind_rev(end),1)], 0);
    
    lunar_backtrack_time_lunanorm(:,1) = repelem(lunar_rev_time, height(height(lunar_orbit_statematrix))) + (lunarbacktrack_statematrix(:,1)-lunarbacktrack_statematrix(1,1));
    lunarbacktrack_statematrix(:,1) = lunarbacktrack_statematrix(:,1)-lunarbacktrack_statematrix(1,1);
    
    %putting backtrack statematrix in geo coord system
    for n=1:height(lunarbacktrack_statematrix)
        lunar_statematrix_during_backtrack(n,:) = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), lunar_backtrack_time_lunanorm(n,1), lunar_backtrack_time_lunanorm(n,1));
        lunar_backtrack_statematrix_geocoords(n,:) = lunarbacktrack_statematrix(n,2:7) + lunar_statematrix_during_backtrack(n,2:7);
        if norm(lunar_backtrack_statematrix_geocoords(n,2:4) - lunar_statematrix_during_backtrack(n,2:4)) < lunar_SOI_radius && n < height(lunarbacktrack_statematrix)-10
            stays_in_SOI = false;
            error("exits lunar SOI prematurely")
            break
        end
    end
    
    lunar_backtrack_statematrix_geocoords = [lunarbacktrack_statematrix(:,1), lunar_backtrack_statematrix_geocoords];
    
    eject_v = lunarbacktrack_statematrix(end,5:7);
    
    lunar_theta_at_release = acosd(dot([1,0,0],-lunar_statematrix_during_backtrack(end,5:7)) / (norm([1,0,0])*norm(-lunar_statematrix_during_backtrack(end,5:7))))

    eject_r_true = [
    cosd(lunar_theta_at_release), -sind(lunar_theta_at_release), 0
    sind(lunar_theta_at_release), cosd(lunar_theta_at_release), 0
    0, 0, 1
    ]*[lunarbacktrack_statematrix(end,2:4).'];

    eject_r_true = eject_r_true.';

    %finding ejection longitude

    plane_cross = cross([0,-1,0],eject_r_true);
    plane_sign = sign(dot(plane_cross, [0,0,1] )) * norm(plane_cross);
    ejection_longitude = atan2d(plane_sign,dot( [0,-1,0], eject_r_true))

    ejection_elevation = 90-acosd(dot(lunarbacktrack_statematrix(end,2:4),eject_v) / (norm(lunarbacktrack_statematrix(end,2:4))*norm(eject_v)));

    total_transit_statematrix = [flipud(lunar_backtrack_statematrix_geocoords(2:end,:)), repelem(1,height(lunar_backtrack_statematrix_geocoords)-1).'; earthtransit_statematrix, repelem(2,height(earthtransit_statematrix)).'];

end

check_matrix = [ %if any of these are low its not a valid traj
reaches_surface
stays_in_SOI
no_SOI_reentry;
];

return_scalar = [
dv_match
delta_t
ejection_longitude
ejection_elevation
];

return_vector = [
eject_r_true
eject_v
];

total_transit_statematrix;

%% plots

subplot(1,2,1)
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
%luna at SOI entry
lunar_map = imread('lunar_map.jpg');
lunar_map = flipud(lunar_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_statematrix_during_backtrack(end,2),lunar_statematrix_during_backtrack(end,3),lunar_statematrix_during_backtrack(end,4),lunar_radius,lunar_radius,lunar_radius,80);
lunar_obj_1 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj_1,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)

%lunar orbit
plot3(lunar_state_matrix(:,1),lunar_state_matrix(:,2),lunar_state_matrix(:,3),"k")

%lunar SOI
%plot(nsidedpoly(1000, 'Center', [r_initial(1),r_initial(2)], 'Radius', lunar_SOI_radius), 'FaceColor', [repelem(0.9,3)])
% xlim([r_initial(1)-1e8,r_initial(1)+1e8])
% ylim([r_initial(2)-1e8,r_initial(2)+1e8])

%target orbit
plot(nsidedpoly(1000, 'Center', [0,0], 'Radius', MEO_radius),FaceAlpha=0)
scatter3(target_r(1),target_r(2),target_r(3),"k","filled")
target_v_unit = target_v./norm(target_v);
quiver3(target_r(1),target_r(2),target_r(3), target_v_unit(1)*unit_vec_size,target_v_unit(2)*unit_vec_size,target_v_unit(3)*unit_vec_size,"k")

%earth transit
plot3(earthtransit_statematrix(:,2),earthtransit_statematrix(:,3),earthtransit_statematrix(:,4),"r")

%lunar backtrack
plot3(lunar_backtrack_statematrix_geocoords(:,2),lunar_backtrack_statematrix_geocoords(:,3),lunar_backtrack_statematrix_geocoords(:,4),"b",LineWidth=1)
%scatter3(eject_r_true(1),eject_r_true(2),eject_r_true(3),10,"k","filled")

subplot(1,2,2)
hold on
grid on
axis equal
plot3(lunar_state_matrix(:,1),lunar_state_matrix(:,2),lunar_state_matrix(:,3),"k")
xlim([r_initial(1)-1e7,r_initial(1)+1e7])
ylim([r_initial(2)-1e7,r_initial(2)+1e7])

[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(r_initial(1),r_initial(2),r_initial(3),lunar_radius,lunar_radius,lunar_radius,80);
lunar_obj_1 = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj_1,'cdata',lunar_map,'facecolor','texturemap','FaceAlpha',0.5)

plot3(lunarbacktrack_statematrix(:,2) + lunar_orbit_statematrix(1,2),lunarbacktrack_statematrix(:,3)+lunar_orbit_statematrix(1,3),lunarbacktrack_statematrix(:,4)+lunar_orbit_statematrix(1,4),"b")

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')



% for n=1:100:height(total_transit_statematrix)
%     subplot(1,2,1)
%     scatter3(total_transit_statematrix(n,2),total_transit_statematrix(n,3),total_transit_statematrix(n,4),"k")
%     drawnow()
% end

function closest_inds = find_closest_inds(a,b)
    [~,closest_inds] = mink( abs(a-b), 2);
    closest_inds = sort(closest_inds).';
end

function state_out = interp_statematrix_timebetween(statematrix, timeseries, time_close, time_target)
    between_inds = find_closest_inds(timeseries, time_close);
    for n=1:width(statematrix)
        state_out(n) = interp1([statematrix(between_inds(1),1),statematrix(between_inds(2),1)], [statematrix(between_inds(1),n),statematrix(between_inds(2),n)], time_target, "linear","extrap");
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