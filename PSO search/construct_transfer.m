function [v_1, v_2, delta_t, dv_total] = construct_transfer(r1, r2, dt_normalised, mu, arrivebody_v)
    
[v_1, v_2, delta_t , a_n] = lambert_solve(r1, r2, dt_normalised, mu);

%need to make sure we're arriving in the same direction as the station
if ~(dot(arrivebody_v,v_2) >= 0)
    v_1 = -v_1;
    v_2 = -v_2;
    transfer_period_full = sqrt(((4*pi^2)/mu)*a_n^3);
    delta_t = transfer_period_full - delta_t;
end

capture_dv = abs(norm(arrivebody_v) - norm(v_2));

dv_total = capture_dv;

end

