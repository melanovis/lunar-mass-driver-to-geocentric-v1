format compact
clear
clc
clf reset

%----------

load("bestdv_results.mat")

for n=1:length(target_altitude_range)
    tmp_struct = getfield(pso_results_bestdv, "altitude_"+string(n));
    tmp_scalar = tmp_struct.scalar;
    tmp_vector = tmp_struct.vector;

    dv_profile(n) = tmp_scalar(1);
    eject_long_profile(n) = tmp_scalar(3);
    eject_elevation_profile(n) = tmp_scalar(4);
    min_acc_profile(n) = tmp_scalar(5);
    target_TA_profile(n) = tmp_scalar(6);
    eject_v_profile(n) = norm(tmp_vector(2,:));

    dv_bounds(n) = sum(tmp_struct.dv_margins); %to account for libration
end

hold on
grid on
ax = gca;
ax.FontSize = 20;
set(gcf, 'Color', [1,1,1])

xlabel("target orbit altitude (km)", Interpreter="latex", FontSize=20)
set(gca,"XScale","log")

% plot(nan,nan,"r")
% scatter(target_altitude_range./1e3,dv_profile./1e3,"k","filled")
% for n=1:length(target_altitude_range)
%     plot([target_altitude_range(n)/1e3,target_altitude_range(n)/1e3],[dv_profile(n)/1e3, (dv_profile(n)+dv_bounds(n))/1e3 ],"r")
%     scatter(target_altitude_range(n)/1e3, (dv_profile(n)+dv_bounds(n))/1e3,"rx")
% end
% set(gca,"yScale","log")
% xlim padded
% ylabel("payload rocket $\Delta$v (km/s)", Interpreter="latex", FontSize=20)
% legend(" $\space$ additional $\Delta$v to account for lunar libration", Interpreter="latex", FontSize=20)
% legend boxoff

% scatter(target_altitude_range./1e3, eject_v_profile./1e3,"k","filled")
% ylim padded
% xlim padded
% ylabel("mass driver exit velocity (km/s)", Interpreter="latex", FontSize=20)

% scatter(target_altitude_range./1e3, eject_long_profile,"k","filled")
% ylim padded
% xlim padded
% ylabel("launch longitude (degrees)", Interpreter="latex", FontSize=20)

% scatter(target_altitude_range./1e3, eject_elevation_profile,"k","filled")
% ylim padded
% xlim padded
% ylabel("launch elevation (degrees)", Interpreter="latex", FontSize=20)
% set(gca,"yScale","log")

% scatter(target_altitude_range./1e3, min_acc_profile./9.81,"k","filled")
% ylim padded
% xlim padded
% ylabel("min payload rocket T/W ratio", Interpreter="latex", FontSize=20)
% set(gca,"yScale","log")

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
