
%{

"Passive Transmission Line Modeling using Vector Fitting with Real Poles,
           Non-Negative Least Squares and Clarke Transformation"
                         Reference: Pending

                              Authors:
                  - José de Jesús Reyes Ramírez
                - Eduardo Salvador Bañuelos Cabral
                  - José Alberto Gutiérrez Robles
                    - José de Jesús Nuño Ayón

           Institution: University of Guadalajara, México. 
                      Last Version: 25/08/26

This repository provides MATLAB implementations of the proposed
passive transmission-line modeling methodology using VF-RP and NNLS.
The code includes the Clarke transformation for modal decomposition,
passive model synthesis, and time-domain simulation. Three-phase and 
double-circuit three-phase transmission-line cases are included.

"Stage 1: Parameter calculation"
"Stage 2: Rational fitting using VF-RP"
"Stage 3: Time-domain simulation"
"Stage 4: ATP file generation"
"Stage 5: Validation with ATP-Draw"

Each one can be executed independently by running its corresponding section.

%}

%% 
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%                     Stage 1: Parameter calculation
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

clc; clear; close all; format long; Figures();

[Fitting_data, Parameters] = Parameter_calculation_data();
return

%% 
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%                  Stage 2: Rational fitting with VFPR                    
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

clc; clear; close all; format long; Figures();

[Values_ab0, Fitting_data] = Rational_fitting();

%%
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%                    Stage 3: Time-domain simulation
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

clc; clear; close all; format long; Figures();

[Phase_voltages, Data, Matrices] = Time_domain_solution();

%%
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%                      Stage 4: ATP file generation
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

Generate_ATPDraw();

%%
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
%                     Stage 5: Validation with ATP-Draw
% ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

%{ 
To compare the ATP-Draw solution with the SS, EMT-type, and NLT solutions, 
first run the generated solution.atp file using the ATP Launcher and save 
the simulation results as ATP_data.MAT. Then, run this section to load the 
ATP-Draw results and plot all four solutions together.
%}

clc; close all;

Phase_voltages.ATP = load("ATP_data.MAT"); % Load the ATP-Draw simulation results

figure; hold on;

for k = 1:3
    h_NLT = plot(Data.Time_data.t, Phase_voltages.V_NLT(k,:), '-k');
	h_SS  = plot(Data.Time_data.t, Phase_voltages.SS(k,:),    '--b');
	h_EMT = plot(Data.Time_data.t, Phase_voltages.EMT(k,:),   '-.c');
end

h_ATP = plot(Phase_voltages.ATP.t, Phase_voltages.ATP.vPhasea, 'r:');
	    plot(Phase_voltages.ATP.t, Phase_voltages.ATP.vPhaseb, 'r:');
	    plot(Phase_voltages.ATP.t, Phase_voltages.ATP.vPhasec, 'r:');

legend([h_NLT, h_SS, h_EMT, h_ATP], 'NLT', 'SS', 'EMT', 'ATP');
xlabel('Time (s)'); ylabel('Voltage (V)');

% ─────────────────────────────────────────────────────────────────────────