format compact
clear
clc
clf reset

%----------

load("bestdv_results.mat")
load("station_orbits.mat")

G = 6.6743e-11;
earth_mass = 5.972e24;
mu = G*earth_mass;
earth_radius = 6.371e6; 

[~,ind_rev] = min(abs(lunar_orbit_statematrix(:,2))); %fine for this case
ref_vector = lunar_orbit_statematrix(ind_rev,2:4);
ref_vector = ref_vector/norm(ref_vector); 

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


lunar_TA_series = TA_from_statematrix(lunar_orbit_statematrix,ref_vector);

lunar_orbit_statematrix = [lunar_orbit_statematrix,lunar_TA_series.'];
lunar_orbit_statematrix_prior = lunar_orbit_statematrix;

for n=1:length(target_altitude_range)

    if n>1
        lunar_segment_statematrix_downsample = [];
        spec_statematrix_downsample = [];
        spec_statematrix = [];
        lunar_orbit_statematrix = lunar_orbit_statematrix_prior;
    end

    spec_statematrix = getfield(target_station_orbits,"altitude_"+string(n));
    spec_TA_series = TA_from_statematrix(spec_statematrix,ref_vector).';
    spec_statematrix = [spec_statematrix, spec_TA_series];

    TA_arrival = rem(spec_TA_series(end),360);

    time_arrival = spec_statematrix(end,1);

    [~,inds_a] = mink(abs(lunar_orbit_statematrix(:,1)-time_arrival),2);
    lunar_TA_at_arrival = interp1([lunar_orbit_statematrix(inds_a(1),1), lunar_orbit_statematrix(inds_a(2),1)], [lunar_orbit_statematrix(inds_a(1),8),lunar_orbit_statematrix(inds_a(1),8)] ,time_arrival);

    lunar_TA_at_arrival = rem(lunar_TA_at_arrival,360);
    arrival_TA_difference = TA_arrival - lunar_TA_at_arrival;

    ind_repeat = 2;
    TA_diff_downsample = [0];

    speccut_shiftback = 0;
    while speccut_shiftback < min(spec_statematrix(:,8))
        speccut_shiftback = speccut_shiftback+360;
    end
    spec_statematrix(:,8) = spec_statematrix(:,8) - speccut_shiftback;
    [~,ind_cut] = mink(abs(spec_statematrix(:,8)-360),2);
    spec_statematrix(ind_cut(2):end,:)=[];

    [~,ind_cut] = mink(abs(lunar_orbit_statematrix(:,8)-360),2);
    lunar_orbit_statematrix(ind_cut(2):end,:)=[];
    
    while max(TA_diff_downsample) < 720

        spec_original = spec_statematrix;
        spec_statematrix_additional = [spec_original(end,1) + (spec_original(:,1) - spec_original(1,1)), spec_original(:,2:7), spec_original(end,8) + (spec_original(:,8) - spec_original(1,8)) ];
        spec_statematrix = [spec_statematrix; spec_statematrix_additional]; 

        lunar_original = lunar_orbit_statematrix;
        lunar_statematrix_additional = [lunar_original(end,1) + (lunar_original(:,1) - lunar_original(1,1)), lunar_original(:,2:7), lunar_original(end,8) + (lunar_original(:,8) - lunar_original(1,8)) ];
        lunar_orbit_statematrix = [lunar_orbit_statematrix; lunar_statematrix_additional]; 

        %get segment of lunar statematrix on same time as spec statematrix
        ind_downsample = 1;
        for ind_t = 1:10:height(spec_statematrix)
            lunar_segment_statematrix_downsample(ind_downsample,1:8) =  interp_statematrix_timebetween(lunar_orbit_statematrix, lunar_orbit_statematrix, spec_statematrix(ind_t,1), spec_statematrix(ind_t,1));
            spec_statematrix_downsample(ind_downsample,:) = spec_statematrix(ind_t,:);
            ind_downsample = ind_downsample+1;
        end

        TA_diff_downsample = spec_statematrix_downsample(:,8) - lunar_segment_statematrix_downsample(:,8); 
    end

    if arrival_TA_difference<0
        arrival_TA_difference = 360+arrival_TA_difference;
    end

    back_time = interp1(TA_diff_downsample,spec_statematrix_downsample(:,1),arrival_TA_difference,"linear");
    forward_time = interp1(TA_diff_downsample,spec_statematrix_downsample(:,1),arrival_TA_difference+360,"linear");

    fprintf("----------\n")
    n
    alignment_time = abs(back_time - forward_time)

    alignment_time_series(n) = alignment_time;

    approx = sqrt( (4*(pi^2)*((earth_radius + target_altitude_range(n))^3)) / mu)
end
save("alignment_time.mat","alignment_time_series")

hold on
grid on
axis tight
set(gcf, 'Color', [1,1,1])
set(gca,"xscale","log")
set(gca,"yscale","log")

scatter(target_altitude_range./1e3,alignment_time_series./3600,"k","filled")
ax = gca;
ax.FontSize = 20;
xlim padded
ylim padded
xlabel("target orbit altitude (km)", Interpreter="latex", FontSize=20)
ylabel("launch cycle period (hours)", Interpreter="latex", FontSize=20)

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')



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

function TA_series = TA_from_statematrix(statematrix,ref_vector)

    theta_prev = 1;
    d_theta_prev = 0;
    inc = 0;
    inc_sign = 1;
    for n=1:height(statematrix)
        r_spec = statematrix(n,2:4);
        r_norm = r_spec/norm(r_spec);

        a = cross(r_norm,ref_vector);
        b = sign(dot(a,[0,0,1])) * norm(a);
        theta = atan2d(b,dot(r_norm,ref_vector));

        d_theta = abs(theta - theta_prev);
        if d_theta > 180
            inc = inc - 360;
        end

        TA_series(n) = theta + inc;
        theta_prev = theta;
        d_theta_prev = d_theta;
    end

    TA_series = -TA_series;

end