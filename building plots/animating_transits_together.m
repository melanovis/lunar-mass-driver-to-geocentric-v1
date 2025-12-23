format compact
clear
clc
clf reset

%----------

record_video = false;

animation_length = 4; %secs
fps = 60;

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

filename = "bestdv_results.mat";
load(filename)

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

timeseries = linspace(min_bounds(1)-10,max_bounds(1),round(animation_length*fps));
timeseries = sort([timeseries,0]);

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
1,0,0
0,1,0
0,0,1
];

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

    for n=1:length(target_altitude_range)

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
            station_points_obj(n) = scatter3(spec_interp_station(2),spec_interp_station(3),spec_interp_station(4)+1e3,2,"k","filled","parent",transform_obj);
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
            
            rocket_points_obj(n) = scatter3(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4),3,"k","filled","parent",transform_obj);
            rocket_label(n) = text(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4)," "+string(n));

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
                %rocket_v_vector(n) = quiver3(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4),spec_v_vector(1)*unit_vec_size,spec_v_vector(2)*unit_vec_size,spec_v_vector(3)*unit_vec_size,"r",MaxHeadSize=0.9);
                quiver3(spec_interp_rocket(2),spec_interp_rocket(3),spec_interp_rocket(4),spec_v_vector(1)*unit_vec_size,spec_v_vector(2)*unit_vec_size,spec_v_vector(3)*unit_vec_size,"r",MaxHeadSize=0.9,Parent=transform_obj);
            end
        end
        
    end

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    
    plot(nsidedpoly(1000, 'Center', [lunar_state_spec(2),lunar_state_spec(3)], 'Radius', lunar_radius), 'FaceColor', [repelem(0.5,3)],"parent",transform_obj)

    if time_spec > 0
        string_additional = "t+";
    else
        string_additional = "t";
    end
    legend(string_additional + string( round( time_spec/3600 ,2)) + " hours SOI departure", Interpreter="latex", FontSize=17 ,location="northwest")
    legend boxoff

    xlim([lunar_state_spec(2)-1.778e8,lunar_state_spec(2)+1.778e8])
    ylim([lunar_state_spec(3)-1e8,lunar_state_spec(3)+1e8])

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
    drawnow()

    if record_video
        frame = getframe(gcf);
        %frame.cdata = imresize(frame.cdata, [1080, nan]);
        writeVideo(v,frame)
    end

    % delete(rocket_points_obj);
    % delete(station_points_obj);
    % delete(target_orbit_obj);
    % delete(rocket_fade_tail);
    % delete(rocket_v_vector);
    % delete(rocket_label);

    fprintf("completed frame %i of %i.\n",ind_t,animation_length*fps)
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