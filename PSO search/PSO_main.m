format compact
clear
clc
clf reset

%----------

target_altitude_range = logspace(log10(1e6),log10(2.75e8), 100); %m
runs_per_alt = 5;

filename = "PSO_result.mat";

if ~exist(filename, 'file')
    ind_run = 1;
    ind_alt = 1;
    pso_result_struct = struct();
    save(filename,"ind_run","ind_alt","pso_result_struct")
end
load(filename)

if ~gcp().Connected %start up the cores
    delete(gcp('nocreate'));
    parpool('local',8);
end


ind_run = 1;
while ind_alt <= length(target_altitude_range)
    
    while ind_run <= runs_per_alt
        load(filename);

        [global_best_scalar, global_best_vector, global_best_statematrix,global_best_checkgood,lunar_orbit_statematrix] = run_pso( target_altitude_range(ind_alt) );
        
        fprintf("completed run %i of alt index %i.\n",ind_run,ind_alt);

        pso_result_struct = setfield(pso_result_struct,"altitude_"+string(ind_alt),"run_"+string(ind_run),"scalar", global_best_scalar );
        pso_result_struct = setfield(pso_result_struct,"altitude_"+string(ind_alt),"run_"+string(ind_run),"vector", global_best_vector );
        pso_result_struct = setfield(pso_result_struct,"altitude_"+string(ind_alt),"run_"+string(ind_run),"statematrix", global_best_statematrix );
        pso_result_struct = setfield(pso_result_struct,"altitude_"+string(ind_alt),"run_"+string(ind_run),"check_good", global_best_checkgood );
        ind_run = ind_run + 1;
        save(filename,"ind_run","ind_alt","pso_result_struct","target_altitude_range","runs_per_alt","lunar_orbit_statematrix"); %lunar orbit statematrix never changes, not the best way to do this but its fine for this application
    end
    ind_alt = ind_alt + 1;
    ind_run = 1;
    save(filename,"ind_run","ind_alt","pso_result_struct","target_altitude_range","runs_per_alt","lunar_orbit_statematrix");
end


sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);

%-----

function [global_best_scalar, global_best_vector, global_best_statematrix,global_best_checkgood,lunar_orbit_statematrix] = run_pso(target_altitude)

    error_tolerance = 1e-11;

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
    inclination = 0;
    
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
    
    % input_vector = rand(1,3);
    % target_altitude = 1e8;
    % [fitness, check_matrix, return_scalar, return_vector, total_transit_statematrix] = PSO_handler(input_vector, target_altitude, lunar_SOI_radius, lunar_orbit_statematrix);
    
    n_DOF = 3;
    var_size = [1, n_DOF]; %solution matrix size
    var_min = 0;
    var_max = 1;
    
    %clerc kennedy construction coefficient function
    phi_1 = 2.05;
    phi_2 = 2.05;
    phi = phi_1+phi_2;
    kappa = 1;
    chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
    
    population = 8*10;
    w = chi;
    w_damp = 0.7;
    w_original = w_damp;
    c1 = chi*phi_1;
    c2 = chi*phi_2;
    max_velocity = (var_max-var_min)*0.25;
    min_velocity = -max_velocity;
    
    %initalise particle template
    empty_particle.position = [];
    empty_particle.velocity = [];
    empty_particle.fitness = [];
    empty_particle.best.position = [];
    empty_particle.best.fitness = [];
    
    %initalise global best
    global_best.fitness = -inf;
    
    particle = repmat(empty_particle, population, 1);
    
    parfor n=1:population
    
        particle(n).position = unifrnd(var_min,var_max,var_size);
        input_vector = particle(n).position;
        [fitness, check_matrix, return_scalar, return_vector, total_transit_statematrix] = PSO_handler(input_vector, target_altitude, lunar_SOI_radius, lunar_orbit_statematrix);
    
        particle(n).fitness = fitness;
    
        particle(n).velocity = zeros(var_size); 
        
        particle(n).best.position = particle(n).position;
        particle(n).best.fitness = particle(n).fitness;
    end
    
    
    for n=1:population
        if particle(n).best.fitness > global_best.fitness
            global_best = particle(n).best;
        end
        particle(n).velocity = rand(var_size).*max_velocity;
    end
    
    fprintf("-----------\n")
    
    iter = 1;
    best_fit_prev = -inf;
    ag_countdown = 5;
    global_best_fitness = -inf;
    steps_forward = 0;
    
    global_best_scalar = nan;
    global_best_vector = nan;
    global_best_statematrix = nan;
    global_best_checkgood = false;
    
    while iter <= 100 
    
        parfor n=1:population
            
            particle(n).velocity = w*particle(n).velocity ...
                + c1*rand(var_size).*(particle(n).best.position - particle(n).position) ...
                + c2*rand(var_size).*(global_best.position - particle(n).position);
            
            particle(n).velocity = max(particle(n).velocity, min_velocity);
            particle(n).velocity = min(particle(n).velocity, max_velocity);
    
            particle(n).position = particle(n).position + particle(n).velocity;
            
            particle(n).position = max(particle(n).position, var_min);
            particle(n).position = min(particle(n).position, var_max);
    
            input_vector = particle(n).position;
            [fitness, check_matrix, return_scalar, return_vector, total_transit_statematrix] = PSO_handler(input_vector, target_altitude, lunar_SOI_radius, lunar_orbit_statematrix);
    
            particle(n).fitness = fitness;
    
            %update personal best
            if particle(n).fitness > particle(n).best.fitness
                particle(n).best.position = particle(n).position;
                particle(n).best.fitness = particle(n).fitness;
            end
            
            %summarise in matrix
            position_summary(n,:) = particle(n).position;
        end
    
        update_plot = false;
        for n=1:population
            if particle(n).best.fitness > global_best.fitness
                global_best = particle(n).best;
                if global_best.fitness > best_fit_prev
                    
                    input_vector = particle(n).position;
                    [fitness, check_matrix, return_scalar, return_vector, total_transit_statematrix] = PSO_handler(input_vector, target_altitude, lunar_SOI_radius, lunar_orbit_statematrix);
    
                    if all(check_matrix)
                        %update global best
                        global_best_fitness = global_best.fitness;
                        best_fit_prev = global_best.fitness;
                        steps_forward = steps_forward+1;
                        global_best_scalar = return_scalar;
                        global_best_vector = return_vector;
                        global_best_statematrix = total_transit_statematrix;
                        global_best_checkgood = true;

                        update_plot = true;
                    end
                end
            end
        end
    
        bestfitnesss(iter) = global_best.fitness;
    
        if iter > 20
            if std(bestfitnesss(end-3:end)) < 5e-3
                ag_countdown = ag_countdown-1;
            else
                ag_countdown = ceil(rand()*30);
            end
            if ag_countdown <= 0
                fitness_score = [];
                for n=1:population
                    fitness_score(n) = particle(n).fitness;
                end
                [~,ind_sort] = sort(fitness_score);
        
                ind_ag = ind_sort(1:round(length(ind_sort)*0.9));
    
                for n=1:length(ind_ag)
                    particle(ind_ag(n)).position = rand(1,n_DOF)*var_max;
                    particle(ind_ag(n)).velocity = rand(1,n_DOF)*max_velocity;
                end
                
                for n=1:population
                    if ~ismember(n, ind_sort(end-3:end))
                        particle(n).position = clamp(particle(n).position + (rand(1,n_DOF)-0.5).*1e-4, 0,1);
                        particle(n).velocity = clamp(particle(n).velocity + (rand(1,n_DOF)-0.5).*1e-2, min_velocity,max_velocity);
                    end
                end
    
                w = w_original;
                ag_countdown = 10;
                fprintf("agitating.\n")
            end
        end
    
        w = w * w_damp;
    
        fprintf("completed iter %i, best fitness: %3.4f, best dv: %3.4f km/s.\n",iter,bestfitnesss(iter),global_best_scalar(1)/1e3)
    
        iter = iter+1;

        % if update_plot
        %     scatter(nan,nan)
        %     hold on
        %     grid on
        %     axis equal
        %     plot3(lunar_orbit_statematrix(:,2),lunar_orbit_statematrix(:,3),lunar_orbit_statematrix(:,4),"k")
        %     plot3(global_best_statematrix(:,2),global_best_statematrix(:,3),global_best_statematrix(:,4),"r")
        %     hold off
        %     drawnow()
        % end

    end 
end


function b = clamp(a,l,u)
    a(a<l) = l;
    a(a>u) = u;
    b=a;
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

