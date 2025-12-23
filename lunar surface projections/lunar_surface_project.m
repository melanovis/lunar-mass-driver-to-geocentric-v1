format compact
clear
clc
clf reset


%------

load("bestdv_results.mat")
load("craters_formatted.mat")

hold on
grid on
axis tight equal

lunar_obliquity = 6.687; 

xlim([-180,180])
ylim([-90,90])

surfacemap = imread('lunar_map.png');
surfacemap = flipud(surfacemap);

wm = image(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

set(gca,'ydir','normal')
uistack(wm,'down')

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, length(target_altitude_range)));

for n=1:length(target_altitude_range)
    tmp_scalar = getfield(pso_results_bestdv,"altitude_"+string(n),"scalar");
    long_series(n) = abs(tmp_scalar(3));
end

[~,ind_longsort] = sort(long_series);
%ind_longsort = flip(ind_longsort);
heightrange = linspace(10,79,length(target_altitude_range));
for n=1:length(target_altitude_range)
    label_height(ind_longsort(n)) = heightrange(n);
end

flip_mask = rem([1:length(label_height)],2)==1;
label_height(flip_mask) = -label_height(flip_mask);

manual_tuning = [
100,-2
49,1
43,-1
59,-1
97,-1
53,-1
47,-1.5
49,-2
5,1
85,1
77,1
65,1
63,1
70,-0.5
23,-0.5
1,-0.5
17,-0.5
19,-0.5
22,0.5
2,1
16,-0.5
14,-0.75
22,0.75
6,0.5
80,1
24,-0.75
3,0.1
46,0.5
98,-0.5
];

for n=1:height(manual_tuning)
    label_height(manual_tuning(n,1)) = label_height(manual_tuning(n,1)) + manual_tuning(n,2);
end

rotation_series(flip_mask) = -30;
rotation_series(~flip_mask) = 30;

% for n=1:height(crater_labels)
%     if crater_latlong(n,2) < 170
%         scatter(crater_latlong(n,2),crater_latlong(n,1),"rx")
%         text(crater_latlong(n,2),crater_latlong(n,1)," "+string(crater_labels(n)),fontsize = 9,Color=[1,1,1])
%     end
% end

for n=1:length(target_altitude_range)
    
    tmp_scalar = getfield(pso_results_bestdv,"altitude_"+string(n),"scalar");
    tmp_vector = getfield(pso_results_bestdv,"altitude_"+string(n),"vector");
    
    spec_dv = tmp_scalar(1);
    spec_long = tmp_scalar(3);
    spec_elev = tmp_scalar(4);
    spec_eject_v = norm(tmp_vector(2,:));
    elev_label = compose("%3.2e", spec_elev);

    if spec_long>0
        rot_modifer = -1;
        align = "right";
    else
        rot_modifer = 1;
        align = "left";
    end

    colour_ease = [1,1,1];
    if spec_long > 150
        label_height(n) = label_height(n)-10;
    end

    plot([spec_long,spec_long],[lunar_obliquity.*cosd(spec_long),label_height(n)],linewidth=0.5, color = colour_ease)
    text(spec_long, label_height(n),"MD-"+string(n)+", "+string(round(spec_eject_v./1e3,2))+" km/s, "+elev_label+char(176), FontSize=7, HorizontalAlignment = align, Rotation = rotation_series(n)*rot_modifer, color = colour_ease)

    scatter(spec_long, lunar_obliquity.*cosd(spec_long), 120, "kx");
    scatter(spec_long, lunar_obliquity.*cosd(spec_long), 30, "r","filled","markeredgecolor","k","markerfacecolor",[cmap(n,:)]);
end

colormap(cmap)
clim([min(target_altitude_range),max(target_altitude_range)]./1e3)
%set(gca,"colorscale","log")
h = colorbar;
set(get(h,'label'),'string','target orbit altitude (km)', Interpreter="latex", FontSize=20);

ax = gca;
ax.FontSize = 20;
set(gcf, 'Color', [1,1,1])
xlabel("longitude (degrees)", Interpreter="latex", FontSize=20)
ylabel("latitude (degrees)", Interpreter="latex", FontSize=20)

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

