format compact
clear
clc
clf reset

%----------

record_video = true;

animation_length = 60; %secs
fps = 60;
aspect_ratio = 1920/1080;

focus_group = 3;

unit_vec_size = 1e7;

if record_video
    v = VideoWriter("transits_test", 'MPEG-4');
    v.FrameRate = fps;
    open(v);
end

earth_mass = 5.972e24;
G = 6.6743e-11;
mu = G*earth_mass;
error_tolerance = 1e-12;
lunar_SOI_radius = 5.7922e7;

load("bestdv_results.mat")

%finding sim bounds 
min_bounds = [repelem(inf,4)];
max_bounds = [repelem(0,4)];
for n=1:length(target_altitude_range)
    tmp_struct = getfield(pso_results_bestdv, "altitude_"+string(n));
    tmp_statematrix = tmp_struct.statematrix;
    for b = 1:4
        min_bounds(b) = min([min_bounds(b),tmp_statematrix(:,b).']);
        max_bounds(b) = max([max_bounds(b),tmp_statematrix(:,b).']);
    end
end
for b = 2:4
    min_bounds(b) = min([min_bounds(b),lunar_orbit_statematrix(:,b).']);
    max_bounds(b) = max([max_bounds(b),lunar_orbit_statematrix(:,b).']);
end

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; 

[~,ind_rev] = mink(abs(lunar_orbit_statematrix(:,2)),5);
ind_rev = sort(ind_rev);
lunar_rev_time = interp1([lunar_orbit_statematrix(ind_rev(end-1),2), lunar_orbit_statematrix(ind_rev(end),2)], [lunar_orbit_statematrix(ind_rev(end-1),1), lunar_orbit_statematrix(ind_rev(end),1)], 0);

%making lunar backtracked statematrix
[~,inds_cycle] = mink(abs(lunar_orbit_statematrix(:,1)-lunar_rev_time),2);
lunar_state_finalcycle = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), lunar_rev_time, lunar_rev_time);
lunar_orbit_statematrix(inds_cycle:end,:) = [];
lunar_orbit_statematrix = [lunar_orbit_statematrix; lunar_state_finalcycle];
lunar_orbit_backtracked = [-(lunar_rev_time - lunar_orbit_statematrix(:,1)), lunar_orbit_statematrix(:,2:7)];
lunar_orbit_statematrix = [lunar_orbit_backtracked(1:end-1,:); lunar_orbit_statematrix];
lunar_orbit_statematrix_original = lunar_orbit_statematrix;
for n=1:2
    forward_time = max(lunar_orbit_statematrix(:,1)) + 1e-10 + lunar_orbit_statematrix_original(:,1) + abs(min(lunar_orbit_statematrix_original(:,1)));
    additional_lunar_statematrix = [forward_time, lunar_orbit_statematrix_original(:,2:7)];
    lunar_orbit_statematrix = [lunar_orbit_statematrix; additional_lunar_statematrix];
end

group_1 = [62:100];
group_2 = [1:26];
group_3 = [27:61];
group_inds(1,1:length(group_1)) = group_1;
group_inds(2,1:length(group_2)) = group_2;
group_inds(3,1:length(group_3)) = group_3;
group_colours = [
255, 132, 0
0, 149, 255
81, 255, 0
]./255;

group_spec = group_inds(focus_group,:);
group_spec(group_spec==0)=[];

new_timebounds = [inf,-inf];
for n=1:length(target_altitude_range)
    if ismember(n,group_spec)
        tmp_struct = getfield(pso_results_bestdv, "altitude_"+string(n));
        tmp_statematrix = tmp_struct.statematrix;
        new_timebounds(1) = min([ new_timebounds(1), min([tmp_statematrix(:,1)]) ]);
        new_timebounds(2) = max([ new_timebounds(2), max([tmp_statematrix(:,1)]) ]);
    end
end

timeseries = linspace(new_timebounds(1)-10,new_timebounds(2)+10,round(animation_length*fps));
timeseries = sort([timeseries,0]);

filename = "station_orbits.mat";

if ~exist(filename, 'file')
    %getting simple target orbit statematrixes
    target_station_orbits = struct();
    for n=1:1:length(target_altitude_range)
    
        target_radius = target_altitude_range(n)+earth_radius;
    
        tmp_struct = getfield(pso_results_bestdv, "altitude_"+string(n));
        tmp_statematrix = tmp_struct.statematrix;
    
        station_v_tangent = [
        0, -1, 0
        1, 0, 0
        0, 0, 1
        ]*[tmp_statematrix(end,2:4).'];
        station_v_tangent_norm = station_v_tangent/norm(station_v_tangent);
        station_v_tangent = station_v_tangent_norm .* sqrt(mu/norm(tmp_statematrix(end,2:4)));
    
        state_station_final = [tmp_statematrix(end,2:4), station_v_tangent.'];
        timerange_station = linspace(tmp_statematrix(end,1), min_bounds(1)*4,25e3);
        [timerange_targetstation, state_matrix_targetstation] = ode45(@(timerange_targetstation, state_matrix_targetstation) g_dynamics_twobody(timerange_targetstation, state_matrix_targetstation, earth_mass), timerange_station, state_station_final, odeset('Reltol',error_tolerance));
        statematrix_station = [timerange_targetstation, state_matrix_targetstation]; 
        statematrix_station = flipud(statematrix_station);
    
        target_station_orbits = setfield(target_station_orbits, "altitude_"+string(n),statematrix_station);
        n
    end
    save(filename, "target_station_orbits");
end
load(filename);

magnitude_window = 120;
magnitude_series = zeros(1,magnitude_window);
Kp = 0.2;
Ki = 0.1;
Kd = 0.1;
filter_tau = 15;
intergral_x = 0;
error_prev_x = 0;
intergral_y = 0;
error_prev_y = 0;

for ind_t = 1:length(timeseries)

    time_spec = timeseries(ind_t);

    scatter(nan,nan,"w")
    hold on
    grid on
    axis equal
    ax = gca;
    ax.FontSize = 20;
    set(gcf, 'Color', [1,1,1])

    transform_obj = hgtransform;

    plot3(lunar_orbit_statematrix(:,2),lunar_orbit_statematrix(:,3),lunar_orbit_statematrix(:,4),"k","parent",transform_obj);

    rocket_points_obj = [];
    station_points_obj = [];
    target_orbit_obj = [];
    rocket_fade_tail = [];
    rocket_v_vector = [];
    rocket_label = [];

    %lunar position
    lunar_state_spec = interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix(:,1), time_spec, time_spec);
    lunar_obj = scatter3(lunar_state_spec(2),lunar_state_spec(3),lunar_state_spec(4),18,"k","filled","parent",transform_obj);
    lunar_SOI_obj = plot(nsidedpoly(1000, 'Center', [lunar_state_spec(2),lunar_state_spec(3)], 'Radius', lunar_SOI_radius), 'FaceColor', [repelem(0.95,3)],"parent",transform_obj);
    
    %plot earth
    earth_obj = scatter3(0,0,0,18,"b","filled");

    ind_focus = 1;
    focus_points = [1e-3,2e-3,3e-3];
    plot_good = false;
    for n=1:length(target_altitude_range)

        [ind_row,~] = find(group_inds==n);
        spec_color = group_colours(ind_row,:);

        spec_station_struct = getfield(target_station_orbits, "altitude_"+string(n));
        spec_station_statematrix = spec_station_struct;
        station_plot_good = false;

        spec_rocket_struct = getfield(pso_results_bestdv, "altitude_"+string(n));
        spec_rocket_statematrix = spec_rocket_struct.statematrix;

        if time_spec < spec_rocket_statematrix(end,1)
            target_orbit_color = [repelem(0.8,3)];
        else
            target_orbit_color = [repelem(0.1,3)];
        end

        %target_orbit_obj(n) = plot3(spec_station_statematrix(:,2),spec_station_statematrix(:,3),spec_station_statematrix(:,4),linewidth=0.1, color = target_orbit_color);
        plot3(spec_station_statematrix(:,2),spec_station_statematrix(:,3),spec_station_statematrix(:,4),linewidth=0.1, color = target_orbit_color,Parent=transform_obj);

        if time_spec > spec_station_statematrix(1,1) && time_spec < spec_station_statematrix(end,1)
            station_plot_good = true;
        end
    
        if station_plot_good
            % interp and plot the station
            spec_interp_station = interp_statematrix_timebetween(spec_station_statematrix, spec_station_statematrix(:,1), time_spec, time_spec);
            station_points_obj(n) = scatter3(spec_interp_station(2),spec_interp_station(3),spec_interp_station(4)+1e3,5,"k","filled","parent",transform_obj);
        end

        rocket_plot_good = false;

        if time_spec > spec_rocket_statematrix(1,1) && time_spec < spec_rocket_statematrix(end,1)
            rocket_plot_good = true;
        end
        time_prev = timeseries(max([ind_t-1,1]));
        time_next = timeseries(min([ind_t+1,length(timeseries)]));
        plot_v_vector_start = false; 
        plot_v_vector_end = false;
        if rocket_plot_good && time_prev < spec_rocket_statematrix(1,1)
            plot_v_vector_start = true;
        end
        if rocket_plot_good && time_next > spec_rocket_statematrix(end,1)
            plot_v_vector_end = true;
        end

        if rocket_plot_good
            spec_interp_rocket = interp_statematrix_timebetween(spec_rocket_statematrix, spec_rocket_statematrix(:,1), time_spec, time_spec);
            
            rocket_points_obj(n) = scatter3(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4),4,"filled","parent",transform_obj,MarkerFaceColor=spec_color);
            if ind_row == focus_group
                rocket_label(n) = text(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4)," MD-"+string(n),fontsize=7);
            end

            [~,ind_tail_end] = min(abs(spec_rocket_statematrix(:,1) - (time_spec - 2*24*3600 ) ));
            [~,inds_tail_start] = mink( abs(spec_rocket_statematrix(:,1)-time_spec), 2);
            fade_tail_statematrix = [spec_rocket_statematrix(ind_tail_end:min(inds_tail_start),:); spec_interp_rocket];

            rocket_fade_tail(n) =  patch([fade_tail_statematrix(:,2); nan],[fade_tail_statematrix(:,3); nan],[fade_tail_statematrix(:,4); nan],'black','EdgeColor','black',...
                'FaceVertexAlphaData',linspace(0,0.5,height(fade_tail_statematrix)+1).','AlphaDataMapping','none',...
                'EdgeAlpha','interp',"parent",transform_obj);

            if plot_v_vector_end || plot_v_vector_start
                if plot_v_vector_start
                    spec_v_vector = spec_interp_rocket(5:7)/norm(spec_interp_rocket(5:7));
                end
                if plot_v_vector_end
                    spec_v_tmp = [spec_interp_station(5:7) - spec_interp_rocket(5:7)];
                    spec_v_vector = spec_v_tmp/norm(spec_v_tmp);
                end
                quiver3(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4),spec_v_vector(1)*unit_vec_size,spec_v_vector(2)*unit_vec_size,spec_v_vector(3)*unit_vec_size,"r",MaxHeadSize=0.9,Parent=transform_obj);
            end
        end
        
        %finding bounds of focus group
        if rocket_plot_good && ind_row == focus_group
            focus_points(ind_focus,:) = [spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4)];
            ind_focus = ind_focus+1;
            plot_good = true;
        end
    end

    focus_centroid = [mean(focus_points(:,1)),mean(focus_points(:,2)),mean(focus_points(:,3))];
    max_mag = -inf;
    for n=1:height(focus_points)
        mag_spec = norm(focus_points(n,:) - focus_centroid);
        if mag_spec > max_mag
            max_mag = mag_spec;
        end
    end
    max_mag = max_mag*3;
    
    magnitude_series = [max_mag, magnitude_series];
    if length(magnitude_series) > magnitude_window
        magnitude_series(magnitude_window:end) = [];
    end
    max_mag = median(magnitude_series);
    max_mag = max([max_mag,3e7]);

    if ind_t == 1
        centroid_x_ref = lunar_state_spec(2);
        centroid_y_ref = lunar_state_spec(3)+1e7;
        centroid_x = centroid_x_ref;
        centroid_y = centroid_y_ref;
    end

    %pid for plot magnitude
    centroid_x_ref = centroid_x_ref + (1/filter_tau)*(focus_centroid(1) - centroid_x_ref);
    error_x = centroid_x_ref - centroid_x;
    intergral_x = intergral_x + error_x;
    intergral_x = max(min(intergral_x, 1), -1);
    derivative_x = error_x - error_prev_x;
    mag_ux = Ki*intergral_x + Kp*error_x + Kd*derivative_x;
    centroid_x = centroid_x + mag_ux;
    error_prev_x = error_x;

    centroid_y_ref = centroid_y_ref + (1/filter_tau)*(focus_centroid(2) - centroid_y_ref);
    error_y = centroid_y_ref - centroid_y;
    intergral_y = intergral_y + error_y;
    intergral_y = max(min(intergral_y, 1), -1);
    derivative_y = error_x - error_prev_y;
    mag_uy = Ki*intergral_y + Kp*error_y + Kd*derivative_y;
    centroid_y = centroid_y + mag_uy;
    error_prev_y = error_y;

    if ind_t>1
        xlim([centroid_x-max_mag,centroid_x+max_mag])
        ylim([centroid_y-max_mag/aspect_ratio,centroid_y+max_mag/aspect_ratio])
    end

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    
    plot(nsidedpoly(1000, 'Center', [lunar_state_spec(2),lunar_state_spec(3)], 'Radius', lunar_radius), 'FaceColor', [repelem(0.5,3)],"parent",transform_obj)
    plot(nsidedpoly(1000, 'Center', [0,0], 'Radius', earth_radius), 'FaceColor', [0,0,1],"parent",transform_obj)

    if time_spec > 0
        string_additional = "t+";
    else
        string_additional = "t";
    end
    legend(string_additional + string( round( time_spec/3600 ,2)) + " hours SOI departure", Interpreter="latex", FontSize=17 ,location="northwest")
    legend boxoff

    % xlim([lunar_state_spec(2)-1.778e8,lunar_state_spec(2)+1.778e8])
    % ylim([lunar_state_spec(3)-1e8,lunar_state_spec(3)+1e8])

    lunar_r = lunar_state_spec(2:4);
    lunar_v = lunar_state_spec(5:7);
    lunar_v = lunar_v/norm(lunar_v);
    lunar_v_cross = cross([-1,0,0],lunar_v);
    lunar_v_sign = sign(dot(lunar_v_cross, [0,0,1] )) * norm(lunar_v_cross);
    lunar_theta = atan2d(lunar_v_sign, dot( [-1,0,0], lunar_v) );

    lunar_r_shifted = [
    cosd(-lunar_theta), -sind(-lunar_theta), 0
    sind(-lunar_theta), cosd(-lunar_theta), 0
    0, 0, 1
    ]*lunar_r.';

    %for lunar-following
    % xlim([lunar_r_shifted(1)-1.778e8,lunar_r_shifted(1)+1.778e8])
    % ylim([lunar_r_shifted(2)-1e8,lunar_r_shifted(2)+1e8])
    % rotation_lunarfollow = makehgtform('zrotate', deg2rad(-lunar_theta));
    % set(transform_obj, 'Matrix', rotation_lunarfollow)

    hold off

    if plot_good
        drawnow()
        if record_video
            frame = getframe(gcf);
            %frame.cdata = imresize(frame.cdata, [1080, nan]);
            writeVideo(v,frame)
        end
        fprintf("completed frame %i of %i.\n",ind_t,animation_length*fps)
    else
        fprintf("not plotting! - %i\n",ind_t)
    end
end

if record_video
    close(v);
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);

%--------------

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
