format compact
clear
clc
clf reset

%----------

%plotting simple results

lunar_orbital_radius = 	363300e3; %m
libration_angle = 6.87;
lunar_radius = 1.7374e6;

filename = "PSO_result.mat";
load(filename)

pso_results_bestdv = struct();

for n=1:length(target_altitude_range)
    
    best_dv = inf;
    
    for m=1:runs_per_alt
        tmp_struct = getfield(pso_result_struct, "altitude_"+string(n),"run_"+string(m));
        tmp_scalar = tmp_struct.scalar;
        tmp_vector = tmp_struct.vector;

        tmp_statematrix = tmp_struct.statematrix;

        %finding extra dv to account for libration
        approx_traversal_orbitalplane = norm(tmp_statematrix(1,2:4)) - norm(tmp_statematrix(end,2:4));
        
        t_a = approx_traversal_orbitalplane + lunar_radius;
        t_b = lunar_radius;
        t_c = sqrt( t_a^2 + t_b^2 - 2*t_a*t_b*cosd(libration_angle) );
        departure_inc_change = asind( t_a*sind(libration_angle)/t_c );
        arrival_inc_change = 180 - (180-departure_inc_change) - libration_angle;

        v_depart = norm(tmp_vector(2,:));        
        extra_dv_departure = tand(departure_inc_change)*v_depart;
        extra_dv_arrival = tand(arrival_inc_change)*norm(tmp_statematrix(end,5:7));

        dv_margins = [extra_dv_departure, extra_dv_arrival];

        dv_profile(n,m) = tmp_scalar(1);
        eject_long_profile(n,m) = tmp_scalar(3);
        eject_elevation_profile(n,m) = tmp_scalar(4);
        min_acc_profile(n,m) = tmp_scalar(5);
        target_TA_profile(n,m) = tmp_scalar(6);
        eject_v_profile(n,m) = norm(tmp_vector(2,:));

        %writing to best dv matrix
        if tmp_scalar(1) < best_dv
            best_dv = tmp_scalar(1);
            pso_results_bestdv = setfield(pso_results_bestdv,"altitude_"+n,tmp_struct);
            pso_results_bestdv = setfield(pso_results_bestdv,"altitude_"+n,"dv_margins",dv_margins);
        end

    end
end

save("bestdv_results.mat","pso_results_bestdv","target_altitude_range","lunar_orbit_statematrix")

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, width(dv_profile)));


hold on
grid on
axis tight
set(gcf, 'Color', [1,1,1])
for n=1:width(dv_profile)
    scatter(target_altitude_range./1e3, dv_profile(:,n)./1e3, "filled", MarkerFaceColor=[cmap(n,:)], MarkerEdgeColor=[0,0,0])
end

set(gca,"xscale","log")
set(gca,"yscale","log")
ylim("padded")

xlabel("target orbit altitude (km)", Interpreter="latex", FontSize=20)
ylabel("rendezvous $\Delta$v (km/s)", Interpreter="latex", FontSize=20)
%ylabel("ejection longitude (degrees)", Interpreter="latex", FontSize=20)
%ylabel("ejection elevation (degrees)", Interpreter="latex", FontSize=20)
%ylabel("payload T/W ratio", Interpreter="latex", FontSize=20)
%ylabel("mass driver exit velocity (km/s)", Interpreter="latex", FontSize=20)

ax = gca;
ax.FontSize = 20;
set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')