clear variables
%% Subsystem 2: Deck
%% Latin Hypercube Sampling
latin_hypercube = lhsdesign(40,4, 'criterion', 'correlation');

deck_width_min = 120;
deck_width_max = 240;
deck_width_range = deck_width_max - deck_width_min;
deck_thickness_min = 2;
deck_thickness_max = 32;
deck_thickness_range = deck_thickness_max - deck_thickness_min;
battery_case_thickness_min = 4;
battery_case_thickness_max = 24;
battery_case_range = battery_case_thickness_max - battery_case_thickness_min;
material_density_choice = [2700; 2810; 1700; 7800];

lhs_deck_width = deck_width_min +((latin_hypercube(:,1))*deck_width_range);
lhs_deck_thickness = deck_thickness_min +((latin_hypercube(:,2))*deck_thickness_range);
lhs_battery_case_thickness = battery_case_thickness_min +((latin_hypercube(:,3))*battery_case_range);
material_density_lhs = ceil(4*latin_hypercube(:,4));
material_density_sim = material_density_choice(material_density_lhs);

%Final sample values used for simulatons:
lhs_deck_width = [237.95; 228.23; 173.12; 220.57; 181.83; 185.04; 193.80; 159.37; 147.14; 197.99; 199.20; 156.64; 232.08; 143.65; 234.87; 176.71; 166.25; 213.46; 164.31; 188.41; 155.89; 126.52; 227.92; 145.90; 191.48; 205.49; 121.51; 150.41; 124.29; 218.49; 210.27; 169.01; 208.49; 129.53; 224.14; 177.36; 136.41; 134.72; 139.55; 203.03];
lhs_deck_thickness = [9.09; 16.28; 27.97; 16.23; 13.45; 28.76; 9.75; 21.13; 18.61; 17.94; 15.20; 6.11; 4.56; 30.87; 17.20; 11.20; 26.60; 24.16; 7.61; 5.64; 23.16; 22.13; 20.38; 8.56; 29.28; 19.66; 31.81; 22.60; 13.18; 3.83; 2.11; 24.70; 12.10; 6.75 ;2.85; 27.09; 29.82; 10.58; 25.36; 14.41]; 
lhs_battery_case_thickness = [8.92; 8.38; 4.19; 4.78; 12.70; 10.38; 11.53; 15.65; 16.66; 16.30; 9.68; 19.28; 18.51; 14.95; 21.27; 22.65; 20.53; 7.68; 17.19; 21.68; 17.92; 23.88; 13.22; 7.38; 18.25; 23.10; 6.17; 11.48; 5.90; 6.80; 10.54; 9.20; 22.39; 12.33; 14.37; 13.81; 20.06; 19.72; 5.34; 15.18];
material_density = [7800; 2810; 2810; 1700; 1700; 1700; 1700; 2700; 2700; 2810; 2810; 2700; 2700; 1700; 7800; 7800; 2700; 2810; 7800; 7800; 7800; 7800; 2700; 1700; 2810; 1700; 2700; 2810; 2700; 2700; 1700; 2810; 2700; 7800; 7800; 7800; 1700; 2810; 1700; 2810];

%Scatter plot to check the spread of random data points for simulation inputs
figure(1)

scatter(latin_hypercube(:,1), latin_hypercube(:,2), 'filled')
title('Latin Hypercube Values for Variable 1 vs 2')
xlabel('Latin Hypercube Deck Width')
ylabel('Latin Hypercube Deck Thickness') 

%% Inputting Raw and Normalised Data of Simulations from CSV

T = csvread('Simulation_Data_2.csv',1,1,[1,1,40,18]);

% Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Material Density 
% Material: Magnesium Alloy 
input_mag_alloy = T(1:10,[1,3,5,7]);
input_mag_alloy_norm = T(1:10,[2,4,6,8]);

% Outputs from Solidworks Simulation - Mass, Yield Strength Safety Factor
% Material: Magnesium Alloy
output_mag_alloy = T(1:10,[9,15]);
output_mag_alloy_norm = T(1:10,[10,16]);

% Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Material Density 
% Material: Aluminium Alloy 6061
input_al_alloy_6061 = T(11:20,[1,3,5,7]);
input_al_alloy_6061_norm = T(11:20,[2,4,6,8]);

% Outputs from Solidworks Simulation - Mass, Yield Strength Safety Factor
% Material: Aluminium Alloy 6061
output_al_alloy_6061 = T(11:20,[9,15]);
output_al_alloy_6061_norm = T(11:20,[10,16]);

% Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Material Density 
% Material: Aluminium Alloy 7075
input_al_alloy_7075 = T(21:30,[1,3,5,7]);
input_al_alloy_7075_norm = T(21:30,[2,4,6,8]);
    
% Outputs from Solidworks Simulation - Mass, Yield Strength Safety Factor
% Material: Aluminium Alloy 7075
output_al_alloy_7075 = T(21:30,[9,15]);
output_al_alloy_7075_norm = T(21:30,[10,16]);

% Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Material Density 
% Material: Carbon Steel
input_carbon_steel = T(31:40,[1,3,5,7]);
input_carbon_steel_norm = T(31:40,[2,4,6,8]);
    
% Outputs from Solidworks Simulation - Mass, Yield Strength Safety Factor
% Material: Carbon Steel
output_carbon_steel = T(31:40,[9,15]);
output_carbon_steel_norm = T(31:40,[10,16]);

%% Plotting Simulation Data

% Scatter plot visually looking at relations between inputs and mass

figure(2)

X = input_al_alloy_6061(:,1);
Y = input_al_alloy_6061(:,2);
Z = input_al_alloy_6061(:,3);
s = output_al_alloy_6061(:,1); %marker Size proportional to mass

scatter3(X,Y,Z,s/60, 'filled')
title('Aluminium Alloy (Marker Size proportional to mass)')
xlabel('Deck Width mm')
ylabel('Deck Thickness mm')
zlabel('Battery Case Thickness mm')


%% Calculating Linear Regression equations

[rsq_mass_mag_alloy, beta_mass_mag_alloy] = regression2(input_mag_alloy_norm(:,1:3), output_mag_alloy_norm(:,1));
[rsq_sf_mag_alloy, beta_sf_mag_alloy] = regression(input_mag_alloy_norm(:,1:3), output_mag_alloy_norm(:,2));
[rsq_mass_mag_alloy_unnorm, beta_mass_mag_alloy_unnorm] = regression2(input_mag_alloy(:,1:3), output_mag_alloy(:,1));

[rsq_mass_al_alloy_6061, beta_mass_al_alloy_6061] = regression2(input_al_alloy_6061_norm(:,1:3), output_al_alloy_6061_norm(:,1));
[rsq_sf_al_alloy_6061, beta_sf_al_alloy_6061] = regression(input_al_alloy_6061_norm(:,1:3), output_al_alloy_6061_norm(:,2));
[rsq_mass_al_alloy_6061_unnorm, beta_mass_al_alloy_6061_unnorm] = regression2(input_al_alloy_6061(:,1:3), output_al_alloy_6061(:,1));

[rsq_mass_al_alloy_7075, beta_mass_al_alloy_7075] = regression2(input_al_alloy_7075_norm(:,1:3), output_al_alloy_7075_norm(:,1));
[rsq_sf_al_alloy_7075, beta_sf_al_alloy_7075] = regression(input_al_alloy_7075_norm(:,1:3), output_al_alloy_7075_norm(:,2));
[rsq_mass_al_alloy_7075_unnorm, beta_mass_al_alloy_7075_unnorm] = regression2(input_mag_alloy(:,1:3), output_al_alloy_7075(:,1));

[rsq_mass_carbon_steel, beta_mass_carbon_steel] = regression2(input_carbon_steel_norm(:,1:3), output_carbon_steel_norm(:,1));
[rsq_sf_carbon_steel, beta_sf_carbon_steel] = regression(input_carbon_steel_norm(:,1:3), output_carbon_steel_norm(:,2));
[rsq_mass_carbon_steel_unnorm, beta_mass_carbon_steel_unnorm] = regression2(input_carbon_steel(:,1:3), output_carbon_steel(:,1));

%% Plotting Regression Data

% plotting surface of mass regression equation with against 2 variables
figure(3)

m1 = @(x1, x2) ((beta_mass_al_alloy_6061(1,1)*x1) + (beta_mass_al_alloy_7075(2,1)*x2));

fsurf(m1)
title('Aluminium Alloy: Normalised Regression Equation of 2 of 4 variables vs Mass')
xlabel('Normalised Deck Width')
ylabel('Normalised Deck Thickness')
zlabel('Normalised Mass')

figure(4)

% plotting surface of mass regression equation with simulation data points
% with color showing whether data passed the static study safety factor

N = 10;
X = rand(N,1);
Y = 10;
% Initialize a blue map
colorMap = [zeros(N, 1), zeros(N, 1), ones(N,1)];
% If y > 0.5, make the markers red.
for k = 1 : N
	if output_al_alloy_6061(k,2) > 1
		colorMap(k, :) = [0,1,0]; % Green
	else
		colorMap(k, :) = [1,0,0]; % Red
	end
end

fsurf(m1)
hold on;
scatter3(input_al_alloy_6061_norm(:,1),input_al_alloy_6061_norm(:,2),output_al_alloy_6061_norm(:,1), 100* ones(length(output_carbon_steel(:,1)), 1), colorMap, 'filled');
title('Aluminium Alloy: Normalised Mass Regression with simulation raw data showing safety factor ')
xlabel('Normalised Deck Width')
ylabel('Normalised Deck Thickness')
zlabel('Normalised Mass')

%% Optimisation: fmincon interior point and sqp & Gentic Algorithm

B = csvread('Simulation_Data_2.csv',58,2,[58,2,65,5]);
N = csvread('Simulation_Data_2.csv',50,2,[50,2,57,5]);

% Upper Bounds of Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Density 
ub_inputs_mag_alloy = B(1,1:4) ;
ub_inputs_al_6061 = B(2,1:4);
ub_inputs_al_7075 = B(3,1:4);
ub_inputs_steel = B(4,1:4);

% Normalised Lower Bounds of Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Density 
lb_inputs_mag_alloy = B(5,1:4);
lb_inputs_al_6061 = B(6,1:4);
lb_inputs_al_7075 = B(7,1:4);
lb_inputs_steel = B(8,1:4);

% Averages of Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Density 
av_inputs_mag_alloy = N(1,1:4);
av_inputs_al_6061 = N(2,1:4);
av_inputs_al_7075 = N(3,1:4);
av_inputs_steel = N(4,1:4);

% Standard Deviations of Inputs: Deck Width, Deck Thickness, Battery Case Thickness, Density 
sd_inputs_mag_alloy = N(5,1:4);
sd_inputs_al_6061 = N(6,1:4);
sd_inputs_al_7075 = N(7,1:4);
sd_inputs_steel = N(8,1:4);

beta_1 = beta_mass_mag_alloy(1,1);
beta_2 = beta_mass_mag_alloy(2,1);
beta_3 = beta_mass_mag_alloy(3,1);
beta_sf1 = beta_sf_mag_alloy(1,1);
beta_sf2 = beta_sf_mag_alloy(2,1);
beta_sf3 = beta_sf_mag_alloy(3,1);

mass_fmincon1 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3); %+ beta_4*x(4);
mass_ga1 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);

x0 = ub_inputs_mag_alloy(1:3);
% A = [-beta_sf1, -beta_sf2, -beta_sf3, -beta_sf4];
A = [-beta_sf1, -beta_sf2, -beta_sf3];
b = 1.1; % normalised safety factor
Aeq = [];
beq = [];
lb = lb_inputs_mag_alloy(1:3);
ub = ub_inputs_mag_alloy(1:3);
LB = lb;
UB = ub;
nonlcon = [];
options = optimoptions('fmincon','Algorithm','sqp');

[opti_mass_mag_alloy_sqp, mass_mag_alloy_sqp] = fmincon(mass_fmincon1,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[opti_mass_mag_alloy_intpoint, mass_mag_alloy_intpoint] = fmincon(mass_fmincon1,x0,A,b,Aeq,beq,lb,ub);
[opti_mass_mag_alloy_ga, mass_mag_alloy_ga] = ga(mass_ga1,3,A,b,Aeq,beq,LB,UB);

% For Robustness

% rng default 
% gs = GlobalSearch;
% mass_gs1 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);
% problem = createOptimProblem('fmincon','x0',x0,'objective',mass_gs1,'lb',lb,'ub',ub,'A', A, 'b', b);
% opti_mass_mag_alloy_gs = run(gs,problem)

beta_1 = beta_mass_al_alloy_6061(1,1);
beta_2 = beta_mass_al_alloy_6061(2,1);
beta_3 = beta_mass_al_alloy_6061(3,1);
beta_sf1 = beta_sf_al_alloy_6061(1,1);
beta_sf2 = beta_sf_al_alloy_6061(2,1);
beta_sf3 = beta_sf_al_alloy_6061(3,1);

mass_fmincon2 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);
mass_ga2 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);

x0 = ub_inputs_al_6061(1:3);
A = [-beta_sf1, -beta_sf2, -beta_sf3];
b = 1; % normalised safety factor
Aeq = [];
beq = [];
lb = lb_inputs_al_6061(1:3);
ub = ub_inputs_al_6061(1:3);
LB = lb;
UB = ub;

[opti_mass_al_alloy_6061_sqp, mass_al_alloy_6061_sqp] = fmincon(mass_fmincon2,x0,A,b,Aeq,beq,lb,ub, nonlcon,options);
[opti_mass_al_alloy_6061_intpoint, mass_al_alloy_6061_intpoint] = fmincon(mass_fmincon2,x0,A,b,Aeq,beq,lb,ub);
[opti_mass_al_alloy_6061_ga, mass_al_alloy_6061_ga] = ga(mass_ga2,3,A,b,Aeq,beq,LB,UB);

beta_1 = beta_mass_al_alloy_7075(1,1);
beta_2 = beta_mass_al_alloy_7075(2,1);
beta_3 = beta_mass_al_alloy_7075(3,1);
beta_sf1 = beta_sf_al_alloy_7075(1,1);
beta_sf2 = beta_sf_al_alloy_7075(2,1);
beta_sf3 = beta_sf_al_alloy_7075(3,1);

mass_fmincon3 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);
mass_ga3 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);

x0 = lb_inputs_al_7075(1:3);
A = [-beta_sf1, -beta_sf2, -beta_sf3];
b = 1.64; % normalised safety factor
Aeq = [];
beq = [];
lb = lb_inputs_al_7075(1:3);
ub = ub_inputs_al_7075(1:3);
LB = lb;
UB = ub;

[opti_mass_al_alloy_7075_sqp, mass_al_alloy_7075_sqp] = fmincon(mass_fmincon3,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[opti_mass_al_alloy_7075_intpoint, mass_al_alloy_7075_intpoint] = fmincon(mass_fmincon3,x0,A,b,Aeq,beq,lb,ub);
[opti_mass_al_alloy_7075_ga, mass_al_alloy_7075_ga] = ga(mass_ga3,3,A,b,Aeq,beq,LB,UB);


beta_1 = beta_mass_carbon_steel(1,1);
beta_2 = beta_mass_carbon_steel(2,1);
beta_3 = beta_mass_carbon_steel(3,1);
beta_sf1 = beta_sf_carbon_steel(1,1);
beta_sf2 = beta_sf_carbon_steel(2,1);
beta_sf3 = beta_sf_carbon_steel(3,1);

mass_fmincon4 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3); %+ beta_4*x(4);
mass_ga4 = @(x) beta_1*x(1) + beta_2*x(2) + beta_3*x(3);

x0 = ub_inputs_steel(1:3);
A = [-beta_sf1, -beta_sf2, -beta_sf3];
b = 1.11; % normalised safety factor
Aeq = [];
beq = [];
lb = lb_inputs_steel(1:3);
ub = ub_inputs_steel(1:3);
LB = lb;
UB = ub;

[opti_mass_carbon_steel_sqp, mass_carbon_steel_sqp] = fmincon(mass_fmincon4,x0,A,b,Aeq,beq,lb,ub, nonlcon,options);
[opti_mass_carbon_steel_intpoint, mass_carbon_steel_intpoint] = fmincon(mass_fmincon4,x0,A,b,Aeq,beq,lb,ub);
[opti_mass_carbon_steel_ga, mass_carbon_steel_ga] = ga(mass_ga4,3,A,b,Aeq,beq,LB,UB);

%% Removing Normalisation & Finding Actual Mass

opti_mass_mag_alloy_x = zeros(1,4);
opti_mass_al_alloy_6061_x = zeros(1,4);
opti_mass_al_alloy_7075_x = zeros(1,4);
opti_mass_carbon_steel_x = zeros(1,4);
m1 = zeros(1,4);
m2 = zeros(1,4);
m3 = zeros(1,4);
m4 = zeros(1,4);
m = [1e+3, 0.33e+3, -0.33e+3, 0];
mass_mag_alloy = 0;
mass_al_alloy_6061 = 0;
mass_al_alloy_7075 = 0;
mass_carbon_steel = 0;

opti_mass_mag_alloy_norm = opti_mass_mag_alloy_intpoint;
opti_mass_al_alloy_6061_norm = opti_mass_al_alloy_6061_intpoint;
opti_mass_al_alloy_7075_norm = opti_mass_al_alloy_7075_intpoint;
opti_mass_carbon_steel_norm = opti_mass_carbon_steel_intpoint;

for i = 1:3
    m1(1,i) = opti_mass_mag_alloy_norm(i) * beta_mass_mag_alloy(i,1);
    opti_mass_mag_alloy_x(1,i) = m1(1,i) * sd_inputs_mag_alloy(i) + av_inputs_mag_alloy(i);
    mass_val1 = opti_mass_mag_alloy_x(1,i) * beta_mass_mag_alloy_unnorm(i);
    mass_mag_alloy = mass_val1 + mass_mag_alloy + m(1);
    
    m2(1,i) = opti_mass_al_alloy_6061_norm(i) * beta_mass_al_alloy_6061(i,1);
    opti_mass_al_alloy_6061_x(1,i) = m2(1,i) * sd_inputs_al_6061(i) + av_inputs_al_6061(i);
    mass_val2 = opti_mass_al_alloy_6061_x(1,i) * beta_mass_al_alloy_6061_unnorm(i);
    mass_al_alloy_6061 = mass_val2 + mass_al_alloy_6061 + m(2);
    
    m3(1,i) = opti_mass_al_alloy_7075_norm(i) * beta_mass_al_alloy_7075(i,1);
    opti_mass_al_alloy_7075_x(1,i) = m3(1,i) * sd_inputs_al_7075(i) + av_inputs_al_7075(i);
    mass_val3 = opti_mass_al_alloy_7075_x(1,i) * beta_mass_al_alloy_7075_unnorm(i);
    mass_al_alloy_7075 = mass_val3 + mass_al_alloy_7075 + m(3);
    
    m4(1,i) = opti_mass_carbon_steel_norm(i) * beta_mass_carbon_steel(i,1);
    opti_mass_carbon_steel_x(1,i) = m4(1,i) * sd_inputs_steel(i) + av_inputs_steel(i);
    mass_val4 = opti_mass_carbon_steel_x(1,i) * beta_mass_carbon_steel_unnorm(i);
    mass_carbon_steel = mass_val4 + mass_carbon_steel + m(4);
end

disp(['Deck Magnesium Alloy (g) = ' num2str(mass_mag_alloy)])
disp(['Deck Aluminium Alloy 6061 (g) = ' num2str(mass_al_alloy_6061)])
disp(['Deck Aluminium Alloy 7075 (g) = ' num2str(mass_al_alloy_7075)])
disp(['Deck Carbon Steel (g) = ' num2str(mass_carbon_steel)])


%% Regression Function

function [rsq, beta] = regression(Input, Output)

%Shuffling the Dataset with random seed
rng(1);                                 % MATLAB random seed 1
newInd = randperm(length(Output));  % New index order for data sets

Input_New = Input(newInd,:);      % Shuffled Input
Output_New = Output(newInd,:);        % Shuffled Output  

% SPLIT DATA SETS (75%-25% RULE) & NORMALIZE WRT TRAINING SET
split_point = floor(0.75*length(Output));  % Round down using floor

%Normalisation of Training Data using mapstd
%[Input_newTrain1,PS_InpTrain1] = mapstd(Input_New(1:split_point,:)');
%[Output_newTrain1,PS_OutTrain1] = mapstd(Output_New(1:split_point,:)');

%Normalisation of Test Data using mapstd
%Input_newTest1 = mapstd('apply',Input_New(split_point+1:end,:)',PS_InpTrain1);
%Output_newTest1 = mapstd('apply',Output_New(split_point+1:end,:)',PS_OutTrain1);

beta = mvregress(Input_New,Output_New); 

rsq = 1 - norm(Input_New*beta - Output_New)^2/norm(Output_New-mean(Input_New))^2;

end

function [rsq, beta] = regression2(Input, Output)

%Shuffling the Dataset with random seed
rng(1);                                 % MATLAB random seed 1
newInd = randperm(length(Output));  % New index order for data sets

Input_New = Input(newInd,:);      % Shuffled Input
Output_New = Output(newInd,:);        % Shuffled Output  

% SPLIT DATA SETS (75%-25% RULE) & NORMALIZE WRT TRAINING SET
split_point = floor(0.75*length(Output));  % Round down using floor

%Normalisation of Training Data using mapstd
[Input_newTrain1] = (Input_New(1:split_point,:)');
[Output_newTrain1] = (Output_New(1:split_point,:)');

%Normalisation of Test Data using mapstd
Input_newTest1 = (Input_New(split_point+1:end,:)');
Output_newTest1 = (Output_New(split_point+1:end,:)');

beta = mvregress(Input_newTrain1',Output_newTrain1'); 

rsq = 1 - norm(Input_newTest1'*beta - Output_newTest1')^2/norm(Output_newTest1-mean(Input_newTest1))^2;

end

