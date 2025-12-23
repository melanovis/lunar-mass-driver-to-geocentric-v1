format compact
clear
clc
clf reset

%------

table_raw = string(readcell("craters_raw.txt"));

min_crater_diameter = 80; %km

% crater_labels = [];
% crater_latlong = [];
% crater_diameter = [];

ind_c = 1;
for n=1:height(table_raw)
    if str2double(table_raw(n,4)) >= min_crater_diameter
        
        crater_labels(ind_c,1) = table_raw(n,1);

        if contains(table_raw(n,2),"N")
            crater_lat = str2double( erase(table_raw(n,2),"N") );
        else
            crater_lat = -str2double( erase(table_raw(n,2),"S") );
        end

        if contains(table_raw(n,3),"W")
            crater_long = -str2double( erase(table_raw(n,3),"W") );
        else
            crater_long = str2double( erase(table_raw(n,3),"E") );
        end

        crater_latlong(ind_c,1:2) = single( [crater_lat,crater_long] );
        crater_diameter(ind_c,1) = single( str2double(table_raw(n,4)) );

        ind_c = ind_c+1;
    end
end

save("craters_formatted.mat","crater_labels","crater_latlong","crater_diameter")
