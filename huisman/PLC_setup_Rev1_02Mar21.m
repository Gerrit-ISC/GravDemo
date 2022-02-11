%clear all

%% Plant
% to be removed once control logic is moved into full model
mass = 50*10^3;
% g = 9.81;

%% Load Experiment Data
Experiment = 244; % Select data. 2= Hoisting, 3 = Lowering

if Experiment == 2
    load Chapter_2_7_Hoisting
    Speed_lookup = [0 0 980 980 0 0];
    Speed_lookup_time=[0 3.71 3.72 8.46 8.47 9];
    Experiment_Start_lookup = [0 0 1 1];
    Experiment_Start_lookup_time = [0 3.44 3.45 5];
    Droop_scaling = 0;
    Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
    Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
    Time = 12.16;
end

if Experiment == 3
    load Chapter_2_7_Lowering
    Speed_lookup = [0 0 -1675 -1675 0 0];
    Speed_lookup_time=[0 4.22 4.23 6.9 6.91 9];
    Experiment_Start_lookup = [0 0 1 1 ];
    Experiment_Start_lookup_time = [0 3.95 4 5];
    Droop_scaling = 0;
    Time = 10.94;
    Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
    Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
end

if Experiment == 4
    load Chapter_2_10_Experiment.mat
    Speed_lookup = [0 0 1000 -1000 500 -500 0 0];
    Speed_lookup_time=[-2 0 2 4 7 10 12 14]+7.14;
    Experiment_Start_lookup = [0 0 1 1];
    Experiment_Start_lookup_time = [0 7.680 7.69 8];
    Droop_scaling = 0;
    Time = 27.8;
    Winch1_r1114_1_SpeedSetpoint = zeros(1,length(Time_ms));
    Winch1_r79_TorqueSetpoint = zeros(1,length(Time_ms));
end



%% PLC parameters
Controller_sampleTime = 0.001; 
% to be moved from model into this section

%% Drive parameters

%[3040] Direction limitation and direction reversal
VFD_n_limit_setp_p1063          = 40000;    % rpm

%[3050] Skip (suppression) bandwidth and speed limiting
VFD_n_limit_pos_p1083           = 40000;    % rpm
VFD_n_max_p1082                 = 1680;     % rpm
VFD_n_limit_neg_p1086           = -40000;   % rpm

%[3070] Expanded ramp-function generator
VFD_RFG_t_RU_scal_p1138         = 1;
VFD_RFG_ramp_up_time_p1120      = 0.7;      % s
VFD_RFG_t_RD_scal_p1139         = 1;
VFD_RFG_ramp_down_time_p1121    = 0.7;      % s

%[6030] Speed setpoint filter 
VFD_n_set_filt_1_T_p1416        = 0;        % ms

% [6040] Speed controller 
VFD_n_act_T_smooth_p1442        = 4;        % ms

%[6050] Kp_n-/Tn_n adaptation 
VFD_n_ctrl_Adpt_sig_Kp_p1455    = 0;
VFD_n_ctrl_adapt_Kp_lower_p1456 = 0;
VFD_n_ctrl_adapt_Kp_upper_p1457 = 0.1;      % actually set up as 0, but lookup table needs monotonically increasing value
VFD_Adapt_factor_lower_p1458    = 1;
VFD_Adapt_factor_upper_p1459    = 1;
VFD_Kp_n_basic_p1460            = 20.174;   % Nm/RPM
VFD_n_ctr_Kp_n_up_scal_p1461    = 1;
VFD_n_ctrl_Tn_n_lower_p1462     = 96;       % ms
VFD_n_ctr_Tn_n_up_sca_p1463     = 1;
VFD_n_ctrl_n_lower_p1464        = 0;        % RPM
VFD_n_ctrl_n_upper_p1465        = 1680;     % RPM
VFD_n_ctrl_Kp_scal_p1466        = 1;

%[6060] Torque setpoint 
VFD_M_accel_T_smooth_p1517      = 4;        % ms

%[6630] Torque limit
VFD_M_max_upper_p1520           = 1345;     % Nm
VFD_M_max_up_mot_scal_p1524     = 1;
VFD_M_lim_var_fixS_src_p1551    = 1;
VFD_M_max_lower_scal_p1525      = 1;
VFD_M_max_lower_p1521           = -1345;    % Nm

%[6640] Current-power limit 
VFD_PU_I_outp_max_r0289         = 452;
VFD_Current_limit_p0640         = 430;      %
VFD_Curr_limit_var_p0641        = 1;
VFD_P_max_mot_p1530             = 125;      % kW
VFD_P_max_gen_p1531             = -375;     % kW
VFD_Id_setp_total_r1624         = 73;
VFD_M                           = inf;      %not actual value, to be updated
VFD_Iq                          = 1;        %not actual value, to be updated