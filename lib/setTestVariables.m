function TEST=setTestVariables(T)

%TEST.end simulation time
%TEST.ST torque switCG 1 or 2 motor
%TEST.PERC is % of nominal power generated
%TEST.CG 1/-1 Charging/generating power --- Charging :dz/dt>0 --- generation dz/dt<0
%TEST.NUM -  1)step test 2)round trip test 3)ramp test

% blockNames = find_system( bdroot , 'Type' , 'Block' )
% allparams('DemonstratorSimscape_v04/Pulley 1')
%allparams('DemonstratorSimscape_v04/cable mass spring/Prismatic Joint3' )


sys = [T.model '/EM System/Motor Winch Right/Revolute_Winch'];

TEST.CG = 1;
name = T.name;

switch name
	case 'SG100'
		
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=1;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SG75'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.75;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SG50'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.50;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SG25'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.25;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SG50-1'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=1;
		TEST.PERC=0.50;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'SG25-1'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=1;
		TEST.PERC=0.25;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'SC100'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=1;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SC75'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.75;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SC50'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.50;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SC25'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.25;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SC50-1'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=1;
		TEST.PERC=0.50;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'SC25-1'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=1;
		TEST.PERC=0.25;
		TEST.CG=1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'SG100-Assist'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=1;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'SG50-Assist'
		%TEST.end=430;
		TEST.end=10;
		
		TEST.ST=2;
		TEST.PERC=0.50;
		TEST.CG=-1;
		TEST.NUM=1;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'E100'
		TEST.end=35;
		TEST.ST=2;
		TEST.PERC=1;
		TEST.NUM=2;
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'E75'
		TEST.end=35;
		TEST.ST=2;
		TEST.PERC=0.75;
		TEST.NUM=2;
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'E50'
		TEST.end=35;
		TEST.ST=2;
		TEST.PERC=0.50;
		TEST.NUM=2;
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'E25'
		TEST.end=35;
		TEST.ST=2;
		TEST.PERC=0.25;
		TEST.NUM=2;
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'E50-1'
		TEST.end=35;
		TEST.ST=1;
		TEST.PERC=0.50;
		TEST.NUM=2;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'E25-1'
		TEST.end=35;
		TEST.ST=1;
		TEST.PERC=0.25;
		TEST.NUM=2;
		
		set_param(sys,'TorqueActuationMode','ComputedTorque')
		set_param(sys,'MotionActuationMode','InputMotion')
	case 'RG25'
		TEST.end=20;
		TEST.ST=2;
		TEST.PERC=1;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RG50'
		TEST.end=15;
		TEST.ST=2;
		TEST.PERC=2;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RG100'
		TEST.end=12.5;
		TEST.ST=2;
		TEST.PERC=4;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RG150'
		TEST.end=11.6667;
		TEST.ST=2;
		TEST.PERC=6;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RG250'
		TEST.end=11;
		TEST.ST=2;
		TEST.PERC=10;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RGXXX'
		TEST.end=11;
		TEST.ST=2;
		TEST.PERC=15;
		TEST.CG=-1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RC25'
		TEST.end=20;
		TEST.ST=2;
		TEST.PERC=1;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RC50'
		TEST.end=15;
		TEST.ST=2;
		TEST.PERC=2;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RC100'
		TEST.end=12.5;
		TEST.ST=2;
		TEST.PERC=4;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RC150'
		TEST.end=11.6667;
		TEST.ST=2;
		TEST.PERC=6;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RC250'
		TEST.end=11;
		TEST.ST=2;
		TEST.PERC=10;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
	case 'RCXXX'
		TEST.end=11;
		TEST.ST=2;
		TEST.PERC=15;
		TEST.CG=1;
		TEST.NUM=3;
		
		set_param(sys,'TorqueActuationMode','InputTorque')
		set_param(sys,'MotionActuationMode','ComputedMotion')
end


switch T.SimType
	case 0
		TEST.mod = 'normal';
	case 1
		TEST.mod='accelerator';
	case 2
		TEST.mod='rapid';
end

