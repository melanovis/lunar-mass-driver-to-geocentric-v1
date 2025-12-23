function [v_1, v_2, delta_t, dv_total]  = construct_and_optimise_transfer(r1, r2, mu, arrivebody_v)

dv_tolerance = 1e-3;

dt_seq = [0,0.5,1];
dv_prev = inf;
dv_total = 0;

for m=1:1e3

    for n=1:length(dt_seq)
        [v_1, v_2, delta_t, dv_spec] = construct_transfer(r1, r2, dt_seq(n), mu, arrivebody_v);
        dv_seq(n) = dv_spec;
    end
    [~,ind_sort] = sort(dv_seq);
    dt_newbounds = dt_seq(flip(ind_sort(1:2)));
    dt_newmid = mean([dt_newbounds]);

    %fprintf("---\n")
    dv_total = dv_seq(2);

    if length(unique(dv_seq)) ~= length(dv_seq)
        %duplicate elements, cannot converge, need to shuffle
        dt_newbounds = rand(1,2);
    end

    if abs(dv_total-dv_prev) < dv_tolerance
        % m
        % delta_t
        % dv_total
        break
    end
    dv_prev = dv_total;
    dt_seq = sort([dt_newmid, dt_newbounds]);
end


end

