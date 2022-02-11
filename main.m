%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main script for Gravitricity simulations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P. Majecki, 30/4/2021, Ver. 9 (model validation vs test data)

clc;clear all;close all;
restoredefaultpath
addpath('lib')
addpath('res')
addpath('huisman')
addpath('../demonstrator_data') % path to data files
warning('off','all')

model_name = 'DemonstratorSimscape';
open(model_name);

%% --- MODEL ---
% MOTOR
M.in = 1.407;	%(RM)[kg.m^2] motor plus brake plus gearbox inertias
M.effType = 0;	% Constant nominal efficiency (0) or speed-varying efficiency (1)
M.eff = 0.947;	%(RM)[-] electric motor efficiency
M.w = 1500;		%[rpm] nominal speed

% TOWER
T.ph = 12500;	%[mm] top pulley center position
T.s = 4200;		%[mm] space among tower columns
T.c = 200;		%[mm] tower column side
T.h0 = 12000;	%[mm] tower column height
T.h1 = 400;		%[mm] tower first roof short side
T.h2 = 400;		%[mm] tower second roof short side
T.hl1 = 400;	%[mm] tower first roof long side
T.hl2 = 400;	%[mm] tower second roof long side
T.s1 = 4200;	%[mm] tower first roof lenght
T.s2 = 5074;	%[mm] tower second roof length
T.b = 8000;		%[mm] tower square base long side
T.bh = 50;		%[mm] tower square base short side
T.rs = 40;		%[mm] tower second roof bars distance
T.fall = 4;		% number of falls

% CABLE
CA.d = 24;		%[mm] cable diameter
CA.mpl = 2.9;	%[kg/m] cable mass per meter
C.EA = 1*30.7e6;	%[N]cable axial stiffness
C.dr = 0.016;	%[-]damping ratio cable, average value
C.drK = 0.028;	%linear factor for time varying frequency carying damping factor

% PULLEYS / SHEAVES
P.d = 480;		%[mm] pulley diameter with cable size
P.h = 85;		%[mm] pulley height
P.up = 12500;	%[mm] upper pulley z axis position
P.m = 45;		%[kg] pulley mass
P.Bd = 720;		%[mm] B-pulley diameter with cable size
P.Bm = 50;		%[kg] B-pulley mass
P.BI = 1.6;		%[kg.m2] B-pulley inertia
P.eff = 0.988;	%[-] combined efficiency of all sheaves

% WINCH
W.d = 538;		%[mm] drum diameter
W.w = 699;		%[mm] drum width between plates
W.pd = 634;		%[mm] drum end plates diameter
W.pds = 10;		%[mm] drum end plates side
W.i = 10;		%[kg.m^2] winch drum inertia
W.nw = 2;		% number of winches
W.n0 = 1;		% initial layer
W.perc0 = 25;	% percentage of used initial layer
W.alpha_r = 0.8;%increase in winch effective radius for cable layers >=2
% Adjust winch radius to account for cable layers
W.d = 2*(W.d/2 + CA.d/2 + W.alpha_r*(W.n0-1)*CA.d);

% GEARBOX
GB.r = 39.9;	%[-] gearbox ratio
GB.eff = 0.96;	%[-] gearbox efficiency

% VARIABLE SPEED DRIVES (original ISC model)
VSD.eff = 0.97;				%[-] combined efficiency of all drives
VSD.SpdRateLimPos = 5e3;	%(RPM/s) rising slew rate limit (1e4 | 5e3)
VSD.SpdRateLimNeg = -5e3;	%(RPM/s) falling slew rate limit
VSD.PwrMax = 200;			%(kW) VSD power limit
VSD.TrqMax = 1400;			%(Nm) max motor torque limit
VSD.TrqMin = -1400;			%(Nm) min motor torque limit

% SPEED CONTROLLER (original ISC model)
C.ts = 0.01;	%[s] sample time
C.kp = 6.25;
C.ki = 6.25;
C.kd = 0.000;
C.x0i06 = 712;
C.x0d = 0;

% LOAD
L.width = 1070;		%[mm] load width
L.height = 2930;	%[mm] load height
L.depth = 2600;		%[mm] load depth
L.mass = 24920;		%[kg] half load mass (single) - original
%L.mass = 51546/2;	%[kg] half load mass (single) - to match load cell readings
%L.mass = 500000/2;	%[kg] half load mass (single) - for full-scale size
L.cpmass = 70;		%[kg] connecting plate mass (only for linked loads)
L.den = 2.86;		%[g/cm^3] average weight density
% Initial position with respect to the center of the load:
L.x0 = 0;			%[mm] load init position x
L.y0 = 0;			%[mm] load init position y
L.z0 = -2200-0e3;   %[mm] load init position z (w.r.t. the lower face of the load)
L.ap = 250;			%[mm] load anchor point with respect to the center
L.ald = P.d;		%[mm] fixed distance between load anchor point surface and center of lower anchor pulley
L.alpha0A = 0;		%[deg] initial anchor angle
L.np_up = 0.432;	%[m/s] nominal load speed [when hoisting]
L.np_dn = 0.620;	%[m/s] nominal load speed [when lowering]

% FRICTION, DAMPING and LOSSES
b_pulley = 0; %(N*m/(deg/s))
b_winch = 0; %(N*m/(deg/s))
% Individual bending cable losses on sheaves 
eta_sh = 0.015; % bending efficiency (% loss per sheave for 180deg wrap)
tau_bending = 1e-3; %(s)
eps_w_sh = 1e-1; %(rad/s) threshold for no sheave motion (1e-1)

%% Other and derived parameters
g = 9.80665; %(m/s^2)
Ip = 0.5*(P.d)*0.5*(W.d*0.5*0.001)^2; % pulley inertia
Iw = 0.5*(W.i + M.in*GB.r^2); % winch inertia incl. motor
W.mass = Iw*2/(W.d*0.5*0.001)^2;
W.dpl = W.w/(1.05*CA.d); % no. of cable loops per layer on drum
% Steady-state winch torque to hold the load (single motor and load)
W_radius_m = (W.d*0.001)*0.5;
tau_w = ((L.mass + L.cpmass/2 + P.m)*g)/2*W_radius_m;
tau_m = tau_w/GB.r;
w0_winch = 0; %(deg/s) winch speed
tau = tau_m;
%lateral horizontal fall spring single S2]
Ls2 = ((0.5*(T.s2+P.d)-2*P.d)/1000);
kl = C.EA/Ls2;
% Change GMM 7/2/2022
% Ls1 = ((T.ph-0.5*W.d)/1000);
Ls1 = ((T.ph-0.5*W.pd)/1000);   %[m] lateral vertical fall spring single S1
kh = C.EA/Ls1;
Lc0 = abs((T.ph-L.height-P.d-L.z0)*0.001);%[m]
Part.m = (L.mass + L.cpmass/2 + P.m)/4;
Part.mL = (L.mass + L.cpmass/2 + P.m)/2; % load mass per fall
L_out = (P.up-W.pd*0.5)*1e-3; %(m) outer (winch) cable length
Part.mwc = L_out*CA.mpl; %(kg) mass of winch cable
L_top = (-1.5*P.d-0.5*P.Bd+0.5*T.s2-0.5*P.d)*1e-3; %(m) top cable length
Part.mtc = L_top*CA.mpl; %(kg) mass of top cable
%in first approximation
%DeltaCable0 = (L.mass+L.cpmass/2+P.m)*1e-3*(Lc0+Ls1*2+Ls2*2)/(4*C.EA*1e-6);%[mm]
DeltaCable0 = 0;


%% GRAPHICS
% Sheave markers
r_sheave =  P.d*0.5;
d_sheave =  P.h;
r_shm = 0.15*r_sheave;
d_shm = 0.4*d_sheave;
x_shm = 0.7*r_sheave;
c_shm = [0 0 0];
phi0_shm = -90;
% Drum markers
r_drum =  W.pd*0.5;
d_drum =  2*W.pds;
r_drm = 0.15*r_drum;
d_drm = 0.5*d_drum;
x_drm = 0.7*r_drum;
c_drm = [0 0 0];
phi0_drm = -90;
% Colors
clear c
c.sheave = rgb('Cyan');
c.Rload = [0 0 1];
c.plate = [0.6 0.6 0.0];
c.opac = 0.5;


%%  --------------------------- MODEL UPDATES 2021 ------------------------

% Electric motor dynamics: 1st order lag
tau_em = 0.06; %(s) motor time constant (0.01)

% Load parameters for PLC and Drive logic
PLC_setup_Rev3

% Kp and Tn multipliers for the speed controller
%Kp_mul = 1; Tn_mul = 1; % no scaling
Kp_mul = 0.3; Tn_mul = 10; % scaling used for the report
fprintf('\nKp_mul = %3.1f, Tn_mul = %3.1f\n',Kp_mul,Tn_mul)

% Brake model (if used)
brkP = 1e5;
brkI = 1e4;
brkD = 0;
Tmax = 400e3;
tau_brk = 0.001;



%% MOTOR SPEED REFERENCE definition
% Set Interpolation off by default
hSpdRefL = [model_name '/Speed Reference/SpdRef_L_rpm'];
hSpdRefR = [model_name '/Speed Reference/SpdRef_R_rpm'];
% Default interpolation
set_param(hSpdRefL,'Interpolate','off')
set_param(hSpdRefR,'Interpolate','off')
% Default initialization
zL_init = 0; %(m)
zR_init = 0; %(m)

% Select the scenario for SpdRef (Left and Right)
SpdRefScen = 21;
switch SpdRefScen
	case 1 % Up/Down symmetric
		SpdRefL = [0 0; 5 -1070; 15 0; 20 1400; 30 0];
		SpdRefR = [0 0; 5 1070; 15 0; 20 -1400; 30 0];
		t_stop = 35;
	case 2 % Up/Down, Right motor half-speed
		SpdRefL = [0 0; 5 -1070; 15 0; 20 1400; 30 0];
		SpdRefR = [0 0; 5 1070/2; 15 0; 20 -1400/2; 30 0];
		t_stop = 35;
	case 3 % Up/Down, Right motor stopped
		SpdRefL = [0 0; 5 -1070; 15 0; 20 1400; 30 0];
		SpdRefR = [0 0; 5 1070/2; 15 0; 20 -1400/2; 30 0];
		t_stop = 35;
	case 4 % Up/Down, Right motor opposite half-speed
		SpdRefL = [0 0; 5 -1070; 15 0; 20 1400; 30 0];
		SpdRefR = [0 0];
		t_stop = 35;
	case 5 % Up/Down, Right motor opposite full-speed
		SpdRefL = [0 0; 5 -1070; 15 0; 20 1400; 30 0];
		SpdRefR = [0 0; 5 1070/2; 10 -1070; 15 0; 20 -1400/2; 25 1400; 30 0];
		t_stop = 35;
	case 6 % Experiment 2 (Hoisting)
		load Chapter_2_7_Hoisting
		Speed_lookup = [0 0 980 980 0 0];
		Speed_lookup_time=[0 3.71 3.72 8.46 8.47 9];
		Experiment_Start_lookup = [0 0 1 1];
		Experiment_Start_lookup_time = [0 3.44 3.45 5];
		Droop_scaling = 0;
		Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
		Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
		Time = 12.16;
		SpdRefL = [Speed_lookup_time' -Speed_lookup'];
		SpdRefR = [Speed_lookup_time' Speed_lookup'];
		t_stop = 15;
		zL_init = 0; %(m)
		zR_init = 0; %(m)
		tBrake_OFF = 3.44; %(s)
		tBrake_ON = 11.7; %(s)
	case 7 % Experiment 3 (Lowering)
		load Chapter_2_7_Lowering
		Speed_lookup = [0 0 -1675 -1675 0 0];
		Speed_lookup_time=[0 4.22 4.23 6.9 6.91 9];
		Experiment_Start_lookup = [0 0 1 1 ];
		Experiment_Start_lookup_time = [0 3.95 4 5];
		Droop_scaling = 0;
		Time = 10.94;
		Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
		Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
		SpdRefL = [Speed_lookup_time' -Speed_lookup'];
		SpdRefR = [Speed_lookup_time' Speed_lookup'];
		t_stop = 15;
	case 8 % Experiment 4 (???)
		load Chapter_2_10_Experiment.mat
		% 		Speed_lookup = [0 0 1000 -1000 500 -500 0 0];
		% 		Speed_lookup_time = [-2 0 2 4 7 10 12 14] + 7.14;
		Speed_lookup =		[0 0	 572.765 1000	-1000	500		-500	0	  0];
		Speed_lookup_time = [0 7.956 8.288   9.144	11.188	14.16	17.128 19.124 22];
		Experiment_Start_lookup = [0 0 1 1];
		Experiment_Start_lookup_time = [0 7.680 7.69 8];
		Droop_scaling = 0;
		Time = 27.8;
		Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
		Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
		SpdRefL = [Speed_lookup_time' -Speed_lookup'];
		SpdRefR = [Speed_lookup_time' Speed_lookup'];
		set_param(hSpdRefL,'Interpolate','on')
		set_param(hSpdRefR,'Interpolate','on')
		
		t_stop = 30;
		zL_init = 4; %(m)
		zR_init = zL_init; %(m)
		tBrake_OFF = 7.68; %7.16; %(s)
		tBrake_ON = 21.36; %(s)
		BrakeOFF = [0 0; tBrake_OFF 0; tBrake_OFF+0.1 1; tBrake_ON 1; tBrake_ON+0.1 0; 30 0];
		
		TrqMeas = [Time_ms*1e-3 Winch2_r80_Torque];
		SpdMeas = [Time_ms*1e-3 Winch2_r61_Speed];
		PwrMeas = [Time_ms*1e-3 Infeed_r82_Power_kW];
		PosMeas = [Time_ms*1e-3 Time_ms*0];
	case 9 % From Excel spreadsheet
		[num,txt,raw] = xlsread('Example Test CSV_Connected W_ JM 15Mar21.xlsx');
		XX = num(2:end,3:5);
		indx = find(XX(:,3)==-1,1,'first');
		XX = XX(1:indx-1,:);
		SpdRefL = [XX(:,3) -XX(:,1)];
		SpdRefR = [XX(:,3) XX(:,2)];
		t_stop = XX(end,3);
		zL_init = num(2,2); %(m)
		zR_init = zL_init; %(m)
	case {10,11,12} % Latest Data from Demonstrator
		switch SpdRefScen
			case 10 % test315_2021-04-19_12.45.06
				% Dynamics of system at lower position of stroke
				% Low speed lowering and hoisting
				load MAT_test315_2021_04_19_124506
				dat = test315_20210419_124506;
				t_start = 300; t_stop = 400;
				%t_start = dat.fileinfo.Time(1); t_stop = dat.fileinfo.Time(end);
			case 11 % test321_2021-04-19_14.12.27
				% Step response test (Step in joystick position, but shows 1 second
				%   acceleration and deceleration events)
				% High speed lowering and hoisting
				load MAT_test321_2021_04_19_141227
				dat = test321_20210419_141227;
				t_start = 60; t_stop = 140;
				%t_start = 60; t_stop = 100;
				%t_start = dat.fileinfo.Time(1); t_stop = dat.fileinfo.Time(end);
			case 12 % test315_2021-04-19_12.45.06
				% Mid-speed lowering
				% Dynamics at maximum position of stroke
				load MAT_test329_2021_04_19_160035
				dat = test329_20210419_160035;
				t_start = 50; t_stop = 150;
				%t_start = dat.fileinfo.Time(1); t_stop = dat.fileinfo.Time(end);
		end
		
		n0 = find(dat.fileinfo.Time>=t_start,1,'first');
		n1 = find(dat.fileinfo.Time<=t_stop,1,'last');
		t_stop = t_stop - t_start;
		
		datTime = dat.fileinfo.Time(n0:n1) - t_start; % Time (s)
		datPos = dat.Tunable_PLC.MH_IntfDbBlockPosCombActPosBlockBot.data(n0:n1); % Block bottom position from pit buffers (m)
		datSpd1 = dat.Drive_1.actual_speed.data(n0:n1); % Motor 1 Actual speed (RPM)
		datSpd2 = dat.Drive_2.actual_speed.data(n0:n1); % Motor 2 Actual speed (RPM)
		datTrq1 = dat.Drive_1.Actual_Torque.data(n0:n1); % Motor 1 Reported Torque (Nm)
		datTrq2 = dat.Drive_2.Actual_Torque.data(n0:n1); %	Motor 2 Reported Torque (Nm)
		datTrqRef1 = dat.Drive_1.Torque_Setpoint.data(n0:n1); % Motor 1 Setpoint Torque (Nm)
		datTrqRef2 = dat.Drive_2.Torque_Setpoint.data(n0:n1); % Motor 2 Setpoint Torque (Nm)
		datSpdRef1 = dat.Drive_1.Speed_controller_setpoint_r1438.data(n0:n1); % Motor 1 Speed Setpoint into speed controller (RPM) (r1438)
		datSpdRef2 = dat.Drive_2.Speed_controller_setpoint_r1438.data(n0:n1); % Motor 2 Speed Setpoint into speed controller (RPM) (r1438)
		datSpdRefPLC1 = dat.Drive_1.p20501_speed_setpoint_int16.data(n0:n1); %	Motor 1 Speed Setpoint from PLC (RPMx10) (r2050[1])
		datSpdRefPLC2 = dat.Drive_2.p20501_speed_setpoint_int16.data(n0:n1); %	Motor 2 Speed Setpoint from PLC (RPMx10) (r2050[1])
		datBrk1 = dat.Tunable_PLC.DiPvHstW1MtrBrkLiftSw.data(n0:n1); %	Winch 1 Brake Lift sensor [0/1]
		datBrk2 = dat.Tunable_PLC.DiPvHstW2MtrBrkLiftSw.data(n0:n1); %	Winch 2 Brake Lift sensor  [0/1]
		datShLoad1 = dat.Tunable_PLC.MH_IntfDbRoLoadCombLmp1.data(n0:n1); % Sheave 1 Load Measuring pin (tons)
		datShLoad2 = dat.Tunable_PLC.MH_IntfDbRoLoadCombLmp2.data(n0:n1); % Sheave 2 Load Measuring pin (tons)
		datPwrALM = dat.ALM.Active_Power.data(n0:n1); % Grid active power (measured by ALM) (kW)
		datPwrPA = dat.PAC4200.PowerTotal_101.data(n0:n1); % Grid active power (measured by Power Analyser) (kW)
		datCurALM = dat.ALM.Infeed_DC_Current.data(n0:n1); % DC current in the DC link from Active Line Module (A)
		datVol1 = dat.Drive_1.actual_DC_link_voltage.data(n0:n1); %	DC link voltage (V)
		datVol2 = dat.Drive_2.actual_DC_link_voltage.data(n0:n1); %	DC link voltage (V)
		
		
		% From Workspace blocks
		SpdRefL_inp = double([datTime datSpdRef1]);
		SpdRefR_inp = double([datTime datSpdRef2]);
		set_param(hSpdRefL,'Interpolate','on')
		set_param(hSpdRefR,'Interpolate','on')
		SpdMeasL_inp = double([datTime datSpd1]);
		SpdMeasR_inp = double([datTime datSpd2]);
		
		TrqRefL_inp = double([datTime datTrqRef1]);
		TrqRefR_inp = double([datTime datTrqRef2]);
		TrqMeasL_inp = double([datTime datTrq1]);
		TrqMeasR_inp = double([datTime datTrq2]);
		
		BrakeOffL_inp = double([datTime datBrk1]);
		PinLoadL_inp = double([datTime datShLoad1]);
		BrakeOffR_inp = double([datTime datBrk2]);
		PinLoadR_inp = double([datTime datShLoad2]);
		
		PwrMeas_inp = double([datTime datPwrALM]);
		PosMeas_inp = double([datTime datPos]);
		
		% Initial load position
		zL_init = double(datPos(1)); %(m)
		zR_init = zL_init; %(m)
	case {20,21} % Latest Data from Demonstrator with extra PLC signals
		
		switch SpdRefScen
			case 20 % test398_20210427_190936
				load MAT_test398_2021_04_27_190936
				dat = test398_20210427_190936;
				%t_start = 60; t_stop = 140;
				%t_start = 60; t_stop = 100;
				t_start = dat.fileinfo.Time(1); t_stop = dat.fileinfo.Time(end);
				
			case 21 % test398_20210427_190936
				load MAT_test402_2021_04_28_201654
				dat = test402_20210428_201654;
				t_start = dat.fileinfo.Time(1); t_stop = dat.fileinfo.Time(end);
				
				
		end
		
		n0 = find(dat.fileinfo.Time>=t_start,1,'first');
		n1 = find(dat.fileinfo.Time<=t_stop,1,'last');
		t_stop = t_stop - t_start;
		
		datTime = dat.fileinfo.Time(n0:n1) - t_start; % Time (s)
		datPos = dat.Tunable_PLC.MH_IntfDbBlockPosCombActPosBlockBot.data(n0:n1); % Block bottom position from pit buffers (m)
		datSpd1 = dat.Drive_1.actual_speed.data(n0:n1); % Motor 1 Actual speed (RPM)
		datSpd2 = dat.Drive_2.actual_speed.data(n0:n1); % Motor 2 Actual speed (RPM)
		datTrq1 = dat.Drive_1.Actual_Torque.data(n0:n1); % Motor 1 Reported Torque (Nm)
		datTrq2 = dat.Drive_2.Actual_Torque.data(n0:n1); %	Motor 2 Reported Torque (Nm)
		datTrqRef1 = dat.Drive_1.Torque_Setpoint.data(n0:n1); % Motor 1 Setpoint Torque (Nm)
		datTrqRef2 = dat.Drive_2.Torque_Setpoint.data(n0:n1); % Motor 2 Setpoint Torque (Nm)
		datSpdRef1 = dat.Drive_1.Speed_controller_setpoint_r1438.data(n0:n1); % Motor 1 Speed Setpoint into speed controller (RPM) (r1438)
		datSpdRef2 = dat.Drive_2.Speed_controller_setpoint_r1438.data(n0:n1); % Motor 2 Speed Setpoint into speed controller (RPM) (r1438)
		datSpdRefPLC1 = dat.Drive_1.p20501_speed_setpoint_int16.data(n0:n1); %	Motor 1 Speed Setpoint from PLC (RPMx10) (r2050[1])
		datSpdRefPLC2 = dat.Drive_2.p20501_speed_setpoint_int16.data(n0:n1); %	Motor 2 Speed Setpoint from PLC (RPMx10) (r2050[1])
		datBrk1 = dat.Tunable_PLC.DiPvHstW1MtrBrkLiftSw.data(n0:n1); %	Winch 1 Brake Lift sensor [0/1]
		datBrk2 = dat.Tunable_PLC.DiPvHstW2MtrBrkLiftSw.data(n0:n1); %	Winch 2 Brake Lift sensor  [0/1]
		datShLoad1 = dat.Tunable_PLC.MH_IntfDbRoLoadCombLmp1.data(n0:n1); % Sheave 1 Load Measuring pin (tons)
		datShLoad2 = dat.Tunable_PLC.MH_IntfDbRoLoadCombLmp2.data(n0:n1); % Sheave 2 Load Measuring pin (tons)
		datPwrALM = dat.ALM.Active_Power.data(n0:n1); % Grid active power (measured by ALM) (kW)
		datPwrPA = dat.PAC4200.PowerTotal_101.data(n0:n1); % Grid active power (measured by Power Analyser) (kW)
		datCurALM = dat.ALM.Infeed_DC_Current.data(n0:n1); % DC current in the DC link from Active Line Module (A)
		datVol1 = dat.Drive_1.actual_DC_link_voltage.data(n0:n1); %	DC link voltage (V)
		datVol2 = dat.Drive_2.actual_DC_link_voltage.data(n0:n1); %	DC link voltage (V)
		
		datPwrDrive1 = dat.Drive_1.Active_Power_Electric_smoothed_r822.data(n0:n1); % some power
		datPwrDrive2 = dat.Drive_2.Active_Power_Electric_smoothed_r822.data(n0:n1); % some power
		
		datPgain1 = dat.Drive_1.Speed_Controller_pgain_eff_r1468.data(n0:n1);
		datItime1 = dat.Drive_1.Speed_Controller_ITime_eff_r1469.data(n0:n1);
		datPout1 = dat.Drive_1.n_ctrl_P_M_output_r1481.data(n0:n1);
		datIout1 = dat.Drive_1.n_ctrl_I_M_output_r1482.data(n0:n1);
		datPIout1 = datPout1 + datIout1;
		datTQ01 = dat.Drive_1.M_set_before_M_suppl_r1508.data(n0:n1);
		datSupTQ1 = dat.Drive_1.M_suppl_plus_M_accel_r1516.data(n0:n1);
		datTQmax1 = dat.Drive_1.M_max_upper_effective_r1538.data(n0:n1);
		datTQmin1 = dat.Drive_1.M_max_lower_effective_r1539.data(n0:n1);
		datSpdErr1 = dat.Drive_1.n_ctrl_system_deviation_r0064.data(n0:n1);
		datSpdErrI1 = dat.Drive_1.n_ctrl_I_sys_dev_r1454.data(n0:n1);
		
		datPgain2 = dat.Drive_2.Speed_Controller_pgain_eff_r1468.data(n0:n1);
		datItime2 = dat.Drive_2.Speed_Controller_ITime_eff_r1469.data(n0:n1);
		datPout2 = dat.Drive_2.n_ctrl_P_M_output_r1481.data(n0:n1);
		datIout2 = dat.Drive_2.n_ctrl_I_M_output_r1482.data(n0:n1);
		datPIout2 = datPout2 + datIout2;
		datTQ02 = dat.Drive_2.M_set_before_M_suppl_r1508.data(n0:n1);
		datSupTQ2 = dat.Drive_2.M_suppl_plus_M_accel_r1516.data(n0:n1);
		datTQmax2 = dat.Drive_2.M_max_upper_effective_r1538.data(n0:n1);
		datTQmin2 = dat.Drive_2.M_max_lower_effective_r1539.data(n0:n1);
		datSpdErr2 = dat.Drive_2.n_ctrl_system_deviation_r0064.data(n0:n1);
		
		% From Workspace blocks
		SpdRefL_inp = double([datTime datSpdRef1]);
		SpdRefR_inp = double([datTime datSpdRef2]);
		set_param(hSpdRefL,'Interpolate','on')
		set_param(hSpdRefR,'Interpolate','on')
		SpdMeasL_inp = double([datTime datSpd1]);
		SpdMeasR_inp = double([datTime datSpd2]);
		
		TrqRefL_inp = double([datTime datTrqRef1]);
		TrqRefR_inp = double([datTime datTrqRef2]);
		TrqMeasL_inp = double([datTime datTrq1]);
		TrqMeasR_inp = double([datTime datTrq2]);
		
		BrakeOffL_inp = double([datTime datBrk1]);
		PinLoadL_inp = double([datTime datShLoad1]);
		BrakeOffR_inp = double([datTime datBrk2]);
		PinLoadR_inp = double([datTime datShLoad2]);
		
		PwrMeas_inp = double([datTime datPwrALM]);
		PosMeas_inp = double([datTime datPos]);
		
		PLCDIAG1_inp = double([datTime datPgain1 datItime1 datPIout1 datPout1 datIout1 datTQ01 datSupTQ1 datSpdErr1]);
		PLCDIAG2_inp = double([datTime datPgain2 datItime2 datPIout2 datPout2 datIout2 datTQ02 datSupTQ2 datSpdErr2]);
		
% 		SpdRefL_inp = [0 0; 5 -200; 10 0];
% 		SpdRefR_inp = [0 0; 5 200; 10 0];
		
		
		% Initial load position
		zL_init = double(datPos(1)); %(m)
		zR_init = zL_init; %(m)
end

Droop_scaling = 0;

% Demonstrator measurement data

% Initial load position (bottom of load w.r.t. ground level)
z0 = L.z0-P.d*0.5-P.d-P.d-L.height-T.h1-T.h2-32; %(mm)
zL0 = z0 + zL_init*1e3;
zR0 = z0 + zR_init;

% Low-pass moving average filter
filt_NUM = ones(1,200)*0.005;
%filt_NUM = 1;

% % Measurement Filter
% f_Yfilt = 1; %(Hz) corner frequency
% Ts_filt = 1e-3; %(s) filter sample time
% tau = 1/(2*pi*f_Yfilt); sys = tf(1,[tau 1]); sysd = c2d(sys,Ts_filt);
% [num,den] = tfdata(sysd); Yfilt_num = num{1}; Yfilt_den = den{1};


% Scope sampling
tscope1 = 0.05; % to match recorded data dampling
tscope2 = -1; % simulation scampling




%% Plot data
plot_flag = 0;
if plot_flag
	
	close all
	figure
	
	switch SpdRefScen
		
		case {10,11,12} % Demonstrator Data
			
			subplot(5,2,1)
			plot(datTime,datSpdRef1); grid on; hold on
			plot(datTime,datSpd1);
			%plot(datTime,datSpdRefPLC1/10);
			%set(gca,'YLim',[-100 1.1*max(datSpd1)])
			title('Motor 1: Speed (RPM)')
			legend('SpdRef','Spd')
			
			subplot(5,2,3)
			plot(datTime,datTrqRef1); grid on; hold on
			plot(datTime,datTrq1);
			title('Motor 1: Torque (Nm)')
			legend('TrqRef','Trq','Location','Best')
			
			subplot(5,2,5)
			plot(datTime,datBrk1); grid on
			set(gca,'YLim',[-0.2 1.2])
			title('Winch 1: Brake Lift sensor [0/1]')
			
			subplot(5,2,7)
			plot(datTime,datShLoad1); grid on; hold on
			title('Sheave 1: Load Measuring pin (tons)')
			
			subplot(5,2,2)
			plot(datTime,datSpdRef2); grid on; hold on
			plot(datTime,datSpd2);
			%plot(datTime,datSpdRefPLC2/10);
			%set(gca,'YLim',[-100 1.1*max(datSpd2)])
			title('Motor 2: Speed (RPM)')
			legend('SpdRef','Spd')
			
			subplot(5,2,4)
			plot(datTime,datTrqRef2); grid on; hold on
			plot(datTime,datTrq2);
			title('Motor 2: Torque (Nm)')
			legend('TrqRef','Trq','Location','Best')
			
			subplot(5,2,6)
			plot(datTime,datBrk2); grid on
			set(gca,'YLim',[-0.2 1.2])
			title('Winch 2: Brake Lift sensor [0/1]')
			
			subplot(5,2,8)
			plot(datTime,datShLoad2); grid on; hold on
			title('Sheave 2: Load Measuring pin (tons)')
			
			subplot(5,2,9)
			plot(datTime,datPos); grid on
			title('Block bottom position from pit buffers (m)')
			%set(gca,'YLim',[0 9])
			xlabel('time (s)')
			
			subplot(5,2,10)
			plot(datTime,datPwrALM); grid on; hold on
			plot(datTime,datPwrPA);
			title('Grid active power (kW)')
			legend('ALM','PA')
			xlabel('time (s)')
			
			set(gcf,'units','normalized','position',[0.1 0.05 0.8 0.85])
			linkaxes(findobj(gcf,'type','axes'),'x')
			
		case {20,21} % Demonstrator Data with extra signals
			
			subplot(5,2,1)
			plot(datTime,datSpdRef1); grid on; hold on
			plot(datTime,datSpd1);
			%plot(datTime,datSpdRefPLC1/10);
			%set(gca,'YLim',[-100 1.1*max(datSpd1)])
			title('Motor 1: Speed (RPM)')
			legend('SpdRef','Spd')
			
			subplot(5,2,3)
			plot(datTime,datTrqRef1); grid on; hold on
			plot(datTime,datTrq1);
			title('Motor 1: Torque (Nm)')
			legend('TrqRef','Trq','Location','Best')
			
			subplot(5,2,5)
			plot(datTime,datBrk1); grid on
			set(gca,'YLim',[-0.2 1.2])
			title('Winch 1: Brake Lift sensor [0/1]')
			
			subplot(5,2,7)
			plot(datTime,datShLoad1); grid on; hold on
			title('Sheave 1: Load Measuring pin (tons)')
			
			subplot(5,2,2)
			plot(datTime,datSpdRef2); grid on; hold on
			plot(datTime,datSpd2);
			%plot(datTime,datSpdRefPLC2/10);
			%set(gca,'YLim',[-100 1.1*max(datSpd2)])
			title('Motor 2: Speed (RPM)')
			legend('SpdRef','Spd')
			
			subplot(5,2,4)
			plot(datTime,datTrqRef2); grid on; hold on
			plot(datTime,datTrq2);
			title('Motor 2: Torque (Nm)')
			legend('TrqRef','Trq','Location','Best')
			
			subplot(5,2,6)
			plot(datTime,datBrk2); grid on
			set(gca,'YLim',[-0.2 1.2])
			title('Winch 2: Brake Lift sensor [0/1]')
			
			subplot(5,2,8)
			plot(datTime,datShLoad2); grid on; hold on
			title('Sheave 2: Load Measuring pin (tons)')
			
			subplot(5,2,9)
			plot(datTime,datPos); grid on
			title('Block bottom position from pit buffers (m)')
			%set(gca,'YLim',[0 9])
			xlabel('time (s)')
			
			subplot(5,2,10)
			plot(datTime,datPwrALM); grid on; hold on
			plot(datTime,datPwrPA);
			title('Grid active power (kW)')
			legend('ALM','PA')
			xlabel('time (s)')
			
			set(gcf,'units','normalized','position',[0.1 0.05 0.8 0.95])
			linkaxes(findobj(gcf,'type','axes'),'x')
			
		otherwise
			
			subplot(3,2,1)
			plot(Time_ms*1e-3,Infeed_r68_Infeed_Current)
			title('Infeed\_r68\_Infeed\_Current (A)')
			
			subplot(3,2,3)
			plot(Time_ms*1e-3,Infeed_r70_DC_Voltage)
			title('Infeed\_r70\_DC\_Voltage (V)')
			
			subplot(3,2,5)
			plot(Time_ms*1e-3,Infeed_r82_Power_kW)
			title('Infeed\_r82\_Power\_kW (kW)')
			
			subplot(2,2,2)
			SP_scaled = Winch1_r2060_1_SpeedSetpoint/max(Winch1_r2060_1_SpeedSetpoint)*1000;
			plot(Time_ms*1e-3,SP_scaled); hold on
			plot(Time_ms*1e-3,Winch1_r61_Speed)
			plot(Time_ms*1e-3,Winch2_r61_Speed)
			title('Winch speed (rpm)')
			legend('r2060\_SP','r61\_Winch1','r61\_Winch2')
			
			subplot(2,2,4)
			plot(Time_ms*1e-3,Winch1_r80_Torque); hold on
			plot(Time_ms*1e-3,Winch2_r80_Torque)
			title('Winch torque (Nm)')
			legend('r80\_Winch1','r80\_Winch2')
			
			set(gcf,'units','normalized','position',[0.2 0.1 0.7 0.7])
			linkaxes(findobj(gcf,'type','axes'),'x')
			
	end
	tightfig
end


%% Power check
pwrchk_flag = 0;
if pwrchk_flag
	
	m_load1 = 50e3; %(kg)
	m_load2 = 23800*2;
	DT = datTime(2) - datTime(1);
	datVel = [diff(datPos(1:2:end))/(2*DT);0];
	
	% Filtering
	f_filt = 1; %(Hz)
	tau_filt = 1/(2*pi*f_filt); Gc = tf(1,[tau_filt 1]); Gd = c2d(Gc,2*DT);
	%[num,den] = tfdata(sysd); hcon.Yfilt_num = num{1}; hcon.Yfilt_den = den{1};
	datVel = lsim(Gd,datVel);
	
	eta = 1;
	eta = P.eff*M.eff*GB.eff*VSD.eff;
	eta = 0.9;
	EtaSigned = eta.^(-sign(datVel));
	
	GravPwr1 = m_load1*g*datVel*1e-3.*EtaSigned; %(kW)
	GravPwr2 = m_load2*g*datVel*1e-3.*EtaSigned; %(kW)
	
	close all
	subplot(211)
	plot(datTime,datPos); grid on
	title('Load position (m)')
	
	subplot(212)
	%plot(datTime,[datPwrALM datPwrPA datPwrDrive1+datPwrDrive2])
	plot(datTime,[datPwrALM datPwrPA])
	hold on; grid on
	plot(datTime(1:2:end),GravPwr1)
	plot(datTime(1:2:end),GravPwr2)
	%legend('ALM','PA','r822','GravPwr1','GravPwr2')
	legend('ALM','PA','GravPwr1','GravPwr2')
	title('Power (kW)')
	xlabel('time (s)');
	
	set(gcf,'units','normalized','position',[0.1 0.05 0.8 0.85])
	linkaxes(findobj(gcf,'type','axes'),'x')
	
end


%% Load mass check
loadmass_flag = 0;
if loadmass_flag
	
	MotTrq_0 = 890; %(Nm) in steady-state conditions (no motion)
	WinchTrq_0 = MotTrq_0 * GB.r; %(Nm) winch torque
	CableTension_0 = WinchTrq_0 / (W.d/2*1e-3) *1e-3; %(kN) Cable tension
	LoadMass = CableTension_0*4/g; %(t) total load mass
	
end


%% Control
% PID control
load ltiC_PID
% IMC control
%load ltiC_IMC
sysC = minreal(ltiC{1});






