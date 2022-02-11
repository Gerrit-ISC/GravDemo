function [sys,x0,str,ts] = sfun_wcable(t,x,u,flag,n0,wperc0,wdpl,alpha,campl,DD,CD)
%%SFUN_WCABLE  -- compute winch cable mass and #layer
%
% Input:
%			1. Winch angle [deg]
% Outputs:
%			1. Cable mass [kg]
%			2. Cable layer [-]
%			3. Winch radius [mm]
% Parameters:
%			1. wn = initial layer [-]
%			2. wperc0 = percentage of used initial layer [%]
%			3. wdpl = number of cable loops per layer [-]
%			4. alpha = increase in winch effective radius for cable layers >=2
%			5. campl = cable mass per unit legth [kg/m]
%			6. DD = drum diameter [mm]
%			7. CD = cable diameter [mm]

switch flag
	
	case 3 %outputs
		
		nrot = u/360; %[-] No. of winch rotations (cable loops) from t_start
		n1 = nrot/wdpl; % No. of full cable layers w.r.t. initial
		I = round(n1); % as above but rounded
		D = n1 - I; % remainder after rounding
		%n0 = initial layer
		%wperc0 = percentage of incomplete initial layer
		n = n0 + I; % No. of full layers
		DeltaD = wperc0 + D*100; % Percent fraction of current layer
		if DeltaD>100
			n = n + 1;
			DeltaD = DeltaD - 100;
		elseif DeltaD<0
			n = n - 1;
			DeltaD = 100 + DeltaD;
		end
		
		% Winch radius
		rw = DD/2 + CD/2 + alpha*(n-1)*CD;
				
		% Cable mass
		mass = 0;		
		% Full layers
		for i=1:n-1
			rw = DD/2 + CD/2 + alpha*(i-1)*CD; % winch radius in layer i
			L = 2*pi*(rw*1e-3)*wdpl; % total cable length in layer i
			mass = mass + L*campl;
		end
		% Top layer
		rw = DD/2 + CD/2 + alpha*(n-1)*CD; % winch radius in layer n
		L = 2*pi*(rw*1e-3)*(DeltaD/100)*wdpl; % total cable length in layer n
		mass = mass + L*campl;
		
		% Return outputs
		sys = [mass; n; rw];
		
	case 0 %initialization
		
		sys = [0 0 3 1 0 1 1];%xc xd y u 0 1 1
		ts = 0.01;
		str = [];
		x0 = [];
		
	otherwise
		sys=[];
		
end