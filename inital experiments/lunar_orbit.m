format compact
clear
clc
clf reset

%----------

error_tolerance = 1e-13;

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
inclination = 5.15;

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

orbit_timespan = linspace(0,orbital_period * (1+2e-3),3e3);

[timerange, state_matrix] = ode45(@(timerange, state_matrix)g_dynamics_twobody(timerange, state_matrix,earth_mass),orbit_timespan,state_initial,odeset('Reltol',error_tolerance));
lunar_orbit_statematrix = [timerange,state_matrix];

lunar_SOI_radius = semi_major_axis*(1-lunar_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);


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
lunar_map = imread('lunar_map.jpg');
lunar_map = flipud(lunar_map);
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(r_initial(1),r_initial(2),r_initial(3),lunar_radius,lunar_radius,lunar_radius,70);
lunar_obj = surf(surface_map_x,surface_map_y,surface_map_z,'EdgeColor','none'); 
set(lunar_obj,'cdata',lunar_map,'facecolor','texturemap')

%lunar orbit
plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3),"k")

%lunar SOI
plot(nsidedpoly(1000, 'Center', [r_initial(1),r_initial(2)], 'Radius', lunar_SOI_radius), 'FaceColor', [repelem(0.8,3)])

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')



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