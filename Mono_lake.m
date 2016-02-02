%% Tufa growth in Mono Lake
% Computational Modelling course GEOL5700-004 - project 2
% KKochanski, 31-Jan-2016

% For each point on the 2D grid,
%  zV(x,y,t) is the height of the lake floor
%  zT(x,y,t) is the height of tufa growth above the rock
%  zW(t) is the water depth
% All heights are elevations

% %|--------Notes to reviewer---------------
% This code models tufa growth in Mono Lake, CA.
% *What are tufa?* Mono Lake is alkali, but many of the springs that feed
% it are acidic. When the waters mix, they precipitate solid CO3,
% which forms near-vertical towers.
% *Why Mono Lake?* The lake has a fabulous history of water piracy. It was
% drained almost dry in the mid-1900s to water LA. Many once-buried
% tufa are now exposed on the shore.
% *What's this model?* I model the fluctuations of water elevation in the 
% basin, and use this to grow a set of randomly-distributed tufa towers.
% *What can I adjust?* The relative growth and collapse rates of the towers
% (lines 60-85). The behavior is simplest with zero collapse chance.

clear all
set(gca, 'fontsize', 14)
%% Parameters
% Set spatial grid size
Lx = 2500; %m
dx = 5; %m
% Calculate grids of x and y values
xvector = 0:dx:Lx;
Nx = length(xvector);

% Initial topography: 
% imagine lake as a parabola with curvature alpha...
alpha = 10^-4; % 
zV0 = ((xvector-Lx/2).^2).*alpha;
% ...reaching max depth zVmax...
zVmin = 118; %m (depth of lake bottom below max level)
zV0 = zV0 - zVmin;
% ...and maximum height zVmax
zVmax = 2000; %m
zV0 = zV0 + zVmax;
zV0 = min(zV0, zVmax);

% Random locations of acidic (tufa-building) springs
num_springs = 15;                     
tufa_i = round((Nx-1)*rand(num_springs, 1))+1; %number of nodes
tufa_elevations = zV0(tufa_i);

% Ask user to check initial topography
hold on
plot(xvector,zV0, 'k', 'linewidth', 2)
plot(xvector(tufa_i), tufa_elevations, '.r', 'markersize', 15)
title('Initial topography with locations of tufa springs')
disp('Showing initial topography... press any key to accept.')
pause()

%% How do tufa grow and collapse?
% 1. Tufa grow, but only if they are underwater.
%      I am assuming that tufa growth rate underwater oscillates yearly
%      because of changes in stream flow rate (highest in spring)
%      In reality, there will also be long term trends as streams die,
%      get diverted for agriculture, etc.

growth_const = 0.1; %m/yr
tufa_growth_rate = @(zT, zW, t) max(0, zW-zT)/(zW-zT).* ... 1 for underwater tufa, else 0
                              growth_const.*(cos(2*pi*t-pi/2)+1);

% 2. Tufa erosion in water and in air
%      Tufa towers are tall and narrow. I am assuming that, once they 
%      exceed a stable height, they are at risk of catastrophic collapse.
% Towers lower than threshold are wider than tall - these never collapse.
% Above this, collapse rate increases linearly.

total_collapse_const = 0.001; % 1/m/year
collapse_threshold   = 1; %m

total_collapse_chance = @(zT)  (zT-collapse_threshold).*total_collapse_const;

%% Time series for lake surface elevation
% (sorry reviewer, this is long!)
start_time = 1800; % years
end_time   = 2016; % years
dt         = 0.1;  % years
tvector    = start_time:dt:end_time;
% Long-term changes in water elevation
%  simplified from data on monolake.org
zW_managed = 1954;
t = start_time + dt;
while t < 1940
    % natural increase in water level
    zW_managed = [zW_managed, zW_managed(end) + 0.1*dt];
    t = t + dt;
end; 
while t < 1982
    % uncontrolled draining to water LA
    zW_managed = [zW_managed, zW_managed(end) - 0.86*dt];
    t = t + dt;
end; 
while t < 1985
    % temporary restriction ordering minimum stream flow
    zW_managed = [zW_managed, zW_managed(end) + 0.81*dt];
    t = t + dt;
end; 
while t < 1994
    % restriction lifts
    zW_managed = [zW_managed, zW_managed(end) - 0.24*dt];
    t = t + dt;
end; 
while t < 1999
    % Water resources control board gets strict
    zW_managed = [zW_managed, zW_managed(end) + 0.28*dt];
    t = t + dt;
end; 
while t < 2012
    % more or less stable
    zW_managed = [zW_managed, zW_managed(end)];
    t = t + dt;
end; 
while t <= end_time
    % recent California drought
    zW_managed = [zW_managed, zW_managed(end) - 0.28*dt];
    t = t + dt;
end;
% Yearly oscillation (currently regulated as max +/- 6ft or 1.8m)
zW_yearly  = 0.3.*cos(2.*pi.*tvector);
zWvector = zW_yearly + zW_managed;

% Plot the historic water levels
figure(2);
set(gca,'fontsize', 14)
plot(tvector, zWvector, 'b', 'linewidth', 2)
title('Model of historic water elevations in Mono Lake')
xlabel('Year'); ylabel('Water elevation (m)')
disp('Plotting historic water levels. Press any key to accept.')
pause()

%% Finally, the model!
% Initialize a vector of tufa tower heights
tufa_heights = zeros(num_springs, length(tvector)+1); %above bedrock

% Start the model!
for i = 1:length(tvector)
    % get current time
    t  = tvector(i);
    % get current water level
    zW = zWvector(i);
    % grow & maybe collapse the tufa towers
    for j = 1:num_springs
        zT = tufa_heights(j,i) + tufa_elevations(j);
        if rand() < dt*total_collapse_chance(tufa_heights(j,i)) ; 
            zT = zT - tufa_heights(j,i)/2;
        end
        zT = zT + dt*tufa_growth_rate(zT,zW,t);
        tufa_heights(j,i+1) = zT - tufa_elevations(j);
    end
end
tufa_heights = tufa_heights(:,2:end);

%% Now we know topography, water levels, tufa heights. Animate!
figure()
disp('Showing water depth and tufa growth through time...')
set(gca, 'fontsize', 14)
for i = 1000:5:length(tvector)
    plot([0,Lx],[zWvector(i), zWvector(i)], '-b') % water
    hold on
    area(xvector, zV0);                           % land
    axis([0, Lx, min(zV0), max(zV0)])
    for j = 1:num_springs
        tufa_x = xvector(tufa_i(j));
        plot([tufa_x, tufa_x],...
            [tufa_elevations(j), tufa_elevations(j)+tufa_heights(j,i)],...
            'k', 'linewidth', 1.5)
    end
    title(sprintf('Year %d', round(tvector(i))))
    xlabel('Horizontal distance (m)')
    ylabel('Elevation (m)')
    pause(0.1)
    hold off
end

figure()
set(gca, 'fontsize', 14)
plot(tvector,tufa_heights');
xlabel('Time (years)'); ylabel('Height of tufa(m)')
title('Tufa height over time')