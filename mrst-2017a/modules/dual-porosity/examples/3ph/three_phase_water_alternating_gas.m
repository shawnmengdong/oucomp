%% Set up and run a five-spot problem with water-alternating-gas (WAG) drive
% One approach to hydrocarbon recovery is to inject gas or water. In the
% water-alternating-gas approach, wells change between injection of gas and
% water to improve the sweep efficiency.
%
clear;
clc;
close all;
% We begin by loading the required modules
mrstModule clear
mrstModule add ad-core ad-blackoil ad-props mrst-gui dual-porosity

[ state, model, schedule ] = three_phase_water_alternating_gas_setup();

%% Simulate the schedule using a fully implicit solver
% For comparison purposes, we also solve the fully implicit case
[~, statesFIMP] = simulateScheduleAD(state, model, schedule);

for i = 1:numel(statesFIMP)
    figure(1);
    p = plotCellData(model.G,statesFIMP{i}.swm)
    p.EdgeAlpha = 0;
    c = colorbar;
    grid on
    drawnow
    
    figure(2);
    p = plotCellData(model.G,1-statesFIMP{i}.swm-statesFIMP{i}.sgm)
    p.EdgeAlpha = 0;
    c = colorbar;
    grid on
    drawnow
    
    figure(3);
    p = plotCellData(model.G,statesFIMP{i}.sgm)
    p.EdgeAlpha = 0;
    c = colorbar;
    grid on
    drawnow
    
    figure(4);
    p = plotCellData(model.G,statesFIMP{i}.s(:,1))
    p.EdgeAlpha = 0;
    c = colorbar;
    caxis([0,1])
    grid on
    drawnow
    
    figure(5);
    p = plotCellData(model.G,1-statesFIMP{i}.s(:,1)-statesFIMP{i}.s(:,3))
    p.EdgeAlpha = 0;
    c = colorbar;
    caxis([0,1])
    grid on
    drawnow
    
    figure(6);
    p = plotCellData(model.G,statesFIMP{i}.s(:,3))
    p.EdgeAlpha = 0;
    c = colorbar;
    caxis([0,1])
    grid on
    drawnow
   
    pause(0.5)
end

% plotToolbar(model.G, statesFIMP);
