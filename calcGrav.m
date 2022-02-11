%% (00) Save simulation data
%RES001 = E100, Single, Basic, Filter ON
%RES002 = E100, Double, Basic, Filter ON
%RES003 = E100, Double, Basic, Filter ON, with cable mass & winch radius

%RES010 = E100, Double, ALL ON, z0=0, Flexi1
%RES011 = E100, Double, ALL ON, z0=-50, Flexi1
%RES012 = E100, Double, ALL ON, z0=-100, Flexi1
%RES013 = E100, Double, ALL ON, z0=-150, Flexi1

%%RES020 = E100, Double, ALL ON, z0=0, Flexi2
%RES021 = E100, Double, ALL ON, z0=-50, Flexi2

%RES030 = E100, Double, ALL ON, z0=0, Flexi1, no filter, TrqLim=1400
%RES031 = E100, Double, ALL ON, z0=0, Flexi1, faster filter, TrqLim=1400
%RES032 = E50, Double, ALL ON, z0=0, Flexi1, faster filter, TrqLim=1400

save res/RES021 ScopeMain ScopeLoad

%% (01) Load and plot simulation data
close all
plot_layout = 41; % 41, 21

% Select design cases to plot
DCS = {'010'}; LGNDS = {'1'};
%DCS = {'011','021'}; LGNDS = {'Flexi1','Flexi2'};
%DCS = {'010','030','031'}; LGNDS = {'Original filter','No filter','Faster filter'};
%DCS = {'010','011','012','013'}; LGNDS = {'z = 0m','z = -50m','z = -100m','z = -150m'};

cc = 'brgky'; % plots colors
cref = 'c'; sref = '--'; sref2 = '-.';


figure; clear h bfig
for i=1:length(DCS)
	
	eval(['load RES' DCS{i}])
	
	switch plot_layout
		
		case 41
			
			subplot(4,1,1)
			plot(ScopeMain.time,ScopeMain.signals(1).values(:,1),'Color',cref,'LineStyle',sref)
			grid on; hold on
			h(i)=plot(ScopeMain.time,ScopeMain.signals(1).values(:,2),'Color',cc(i));
			title('Electrical Power (kW)')
			
			subplot(4,1,2)
			plot(ScopeMain.time,ScopeMain.signals(2).values(:,2),'Color',cc(i));
			grid on; hold on
			title('Motor speed (rpm)')
			
			subplot(4,1,3)
			plot(ScopeMain.time,ScopeMain.signals(3).values,'Color',cc(i));
			grid on; hold on
			title('Motor torque (Nm)')
			
			subplot(4,1,4)
			plot(ScopeMain.time,ScopeMain.signals(4).values,'Color',cc(i));
			grid on; hold on
			title('Load position (m)')
			xlabel('time (sec)')
			
			subplot(4,1,4);
			fpos = [0.1 0.1 0.6 0.8];

		case 21
			
			subplot(2,1,1)
			plot(ScopeMain.time,ScopeMain.signals(1).values(:,1),'Color',cref,'LineStyle',sref)
			grid on; hold on
			h(i)=plot(ScopeMain.time,ScopeMain.signals(1).values(:,2),'Color',cc(i));
			%grid on; hold on
			
% 			plot(ScopeMain.time([1 end]),[1.02*62.5]*[1;1],'Color',cref,'LineStyle',sref)
% 			plot(ScopeMain.time([1 end]),[0.98*62.5]*[1;1],'Color',cref,'LineStyle',sref)
% 			plot(ScopeMain.time([1 end]),[-1.02*62.5]*[1;1],'Color',cref,'LineStyle',sref)
% 			plot(ScopeMain.time([1 end]),[-0.98*62.5]*[1;1],'Color',cref,'LineStyle',sref)
			
			title('Electrical Power (kW)')
			
			subplot(2,1,2)
			plot(ScopeLoad.time,ScopeLoad.signals(2).values,'Color',cc(i));
			grid on; hold on
			title('Load velocity (m/s)')
			
			subplot(2,1,2);
			fpos = [0.3 0.3 0.4 0.6];
			
			
	end
	
	
	%%% Compute benchmark figures %%%
	% Distance (km)
	bfig(i,1) = ScopeMain.signals(1).values(end);
	
end

legend(h,LGNDS{:},'Location','Best')
set(gcf,'units','normalized','position',fpos)
linkaxes(findobj(gcf,'type','axes'),'x')
%tightfig

%%%% Display benchmark figures %%%
%clc
fprintf('\n----- Benchmark Figures -----\n')
fprintf(' PwrEnd = Last Pwr (kW)\n')
fprintf('\t\t\t\tDist\n')
for i=1:length(DCS)
	fprintf('%12s\t%4.3f\n',LGNDS{i},bfig(i,:))
end

%% (02) Calculation of natural frequency
EA = 30.7e6;	%[N]cable axial stiffness

z_v = [-10, -10, -750]'; %(m) load position
m_v = []'; %(kg) load mass
nc_v = []'; %(-) No of lifting cables


L_out = (P.up-W.pd*0.5)*1e-3; %(m) outer (winch) cable length
L_in = (P.up-L.height-P.d-z_v*1e3)*1e-3; %(m) inner (tower) cable length

m_out = L_out*CA.mpl; %(kg) outer cable mass
m_in = 2*(L_in*CA.mpl); %(kg) inner cable mass x2
m_cab = m_out + m_in; %(kg) total cable mass
m_load = L.mass + P.m + L.cpmass/2; %(kg) load mass incl. sheave and half-plate
m_tot = m_cab + m_load; %(kg) total load mass

k_out = C.EA/L_out; %(N/m) outer cable stiffness
k_in = 2*(C.EA./L_in); %(N/m) inner cables stiffness (parallel connection)
k_cab = 1*k_out + k_in; %(N/m) total cable stiffness (parallel connection)
%k_cab = k_out*k_in./(k_out+k_in);

fn_v = 1/(2*pi)*sqrt(k_cab./m_tot) %(Hz) system natural frequency

%% (02A) Calculation of natural frequency
EA = 30.7e6;	%[N]cable axial stiffness

z_v = [10, 750, 750]'; %(m) load position
m_v = [50e3, 50e3, 500e3]'; %(kg) load mass
nc_v = [4, 8, 8]'; %(-) No of lifting cables
nn = length(z_v);

mc_v = nc_v.*(z_v.*CA.mpl); %(kg) cable mass
mL_v = m_v + mc_v;
k_v = nc_v.*(EA./z_v); %(N/m) effective cable stiffness (parallel connection)
fn_v = 1/(2*pi).*sqrt(k_v./mL_v) %(Hz) system natural frequency

%% (03) Holding torque calculation
z0 = 0;
% Winch radius (m)
rw0 = 0.001*(W.d/2); % no cable
rw = 0.001*(W.d/2 + CA.d/2 + W.alpha_r*(W.n0-1)*CA.d);
% Load mass (kg)
mL = L.mass + L.cpmass/2 + P.m;
% 
Tauw0 = (mL*g)/2*rw0; Taum0 = Tauw0/GB.r
Tauw = (mL*g)/2*rw; Taum = Tauw/GB.r

%% ----------------- NEW WORK (September 2021) ----------------------------

%% (04) Deep Shaft System Analysis - Model #1 with TFs
EA = 30.7e6;	%[N] cable axial stiffness
mc_pm = 2.9;	%[kg/m] cable mass per meter
rw0 = 0.562/2;	%[m] winch radius (no cable)
tau_em = 0.06;  %(s) motor time constant
GBr = 39.9;		%(-) gearbox ratio
b = 1e4;			%(N/(m/s)) damping

% Parameters
z_v = [150, 450, 750]'; %(m) rope length / load position
m_v = [50e3, 50e3, 50e3]'; %(kg) load mass
nc_v = [4, 4, 4]'; %(-) No of lifting cables
nn = length(z_v);

% Natural frequency (Hz)
mc_v = nc_v.*(z_v.*mc_pm); %(kg) cable mass
mL_v = m_v + mc_v; %(kg) total load mass
k_v = nc_v.*(EA./z_v); %(N/m) effective cable stiffness (parallel connection)
fn_v = 1/(2*pi).*sqrt(k_v./mL_v) %(Hz) system natural frequency

% Balancing (static) motor torque per winch (Nm)
tauM0_v = (m_v + mc_v)*g/2./(nc_v/2)*rw0/GBr; 

% LTI models: stiffness
G0 = tf(1,1)*ones(nn,1);
s = tf('s');
opt = stepDataOptions('StepAmplitude',100*GBr/rw0*2);
for i=1:nn
	G0(i) = (1/mL_v(i))/(s^2 + b/mL_v(i)*s + k_v(i)/mL_v(i));
end

step(G0,20,opt)

% Closed-loop
Kp = 0.1;
Ti = 1000;
Td = 2;
C0 = Kp*(1 + 1/Ti/s + Td*s);

G1 = tf(1,1)*ones(nn,1);
H = tf(1,1)*ones(nn,1);
for i=1:nn
	G1(i) = G0(i)/s;
	%G1(i) = dcgain(G0(i))/s;
	C1 = C0 / dcgain(G0(i));
	H(i) = G1(i)*C1/(1 + G1(i)*C1);
	H(i) = minreal(H(i));
end

close all
step(H,50)

%% (05) Deep Shaft System Analysis - Model #2, M1-K-M2 
EA = 30.7e6;	%[N] cable axial stiffness
mc_pm = 2.9;	%[kg/m] cable mass per meter
rw0 = 0.562/2;	%[m] winch radius (no cable)
tau_em = 0.06;  %(s) motor time constant
GBr = 39.9;		%(-) gearbox ratio
b = 1e4;		%(N/(m/s)) damping

% Parameters
z_v = [150, 450, 750]'; %(m) rope length / load position
m_v = [50e3, 50e3, 50e3]'; %(kg) load mass
nc_v = [4, 4, 4]'; %(-) No of lifting cables
nn = length(z_v);

% Natural frequency (Hz)
mc_v = nc_v.*(z_v.*mc_pm); %(kg) cable mass
mL_v = m_v + mc_v; %(kg) total load mass
k_v = nc_v.*(EA./z_v); %(N/m) effective cable stiffness (parallel connection)
fn_v = 1/(2*pi).*sqrt(k_v./mL_v) %(Hz) system natural frequency

% Balancing (static) motor torque per winch (Nm)
tauM0_v = (m_v + mc_v)*g/2./(nc_v/2)*rw0/GBr; 

% LTI models: stiffness
G0 = tf(1,1)*ones(nn,1);
s = tf('s');
opt = stepDataOptions('StepAmplitude',100*GBr/rw0*2);
for i=1:nn
	G0(i) = (1/mL_v(i))/(s^2 + b/mL_v(i)*s + k_v(i)/mL_v(i));
end

step(G0,20,opt)

% Closed-loop
Kp = 0.1;
Ti = 1000;
Td = 2;
C0 = Kp*(1 + 1/Ti/s + Td*s);

G1 = tf(1,1)*ones(nn,1);
H = tf(1,1)*ones(nn,1);
for i=1:nn
	G1(i) = G0(i)/s;
	%G1(i) = dcgain(G0(i))/s;
	C1 = C0 / dcgain(G0(i));
	H(i) = G1(i)*C1/(1 + G1(i)*C1);
	H(i) = minreal(H(i));
end

close all
step(H,50)


%% (06) GravSimple model definition and analysis - single depth
EA = 30.7e6;	%[N] cable axial stiffness
mc_pm = 2.9;	%[kg/m] cable mass per meter
rw0 = 0.562/2;	%[m] winch radius (no cable) = sheave radius
tau_em = 60e-3; %(s) motor time constant
GBr = 39.9;		%(-) gearbox ratio
b = 5e3;		%(N/(m/s)) damping
Jm = 1.407;		%[kg.m^2] motor plus brake plus gearbox inertias
Jw = 10;		%[kg.m^2] winch drum inertia
Js = 1.6;		%[kg.m2] sheaves inertia

D = 50; %(m) rope length / load position
M = 50e3; %(kg) load mass
nc = 4; %(-) No of lifting cables
nm = 2; % number of motors

mc = nc*(D*mc_pm); %(kg) cable mass
alfa = 0.5; % cable mass split coef
M1 = alfa*mc; % upper mass
M2 = (1-alfa)*mc + M; % lower mass

% Cable stiffness and Natural frequency
mL = M + mc; %(kg) total load mass
K = nc*(EA/D); %(N/m) effective cable stiffness (parallel connection)
fn = 1/(2*pi).*sqrt(K/mL) %(Hz) system natural frequency
fn2 = 1/(2*pi).*sqrt(K.*(M1+M2)./(M1.*M2)) %(Hz) system natural frequency

% Balancing (static) motor torque per motor (Nm)
tauM0 = (M + mc)*g/nm/(nc/2)*rw0/GBr; 
tau_max = tauM0;

% Controller tuning
Kp = 0.1;
Ti = 1000;
Td = 0;

PID.P = Kp;
PID.I = Kp/Ti;
PID.D = Kp*Td;

s = tf('s');
sysC = 0.30306*(s+26.45)*(s+1.154)/s/(s+156.9);
%sysC = C_imc;

% Linear Analysis
% x = [Tm, z1, z1dot, z2, z2dot]
% u = Tref
% y = (30/pi)*(4*GBr/rw0)*z1dot
% Params: tau_em, K, M1, M2, b, k0=4*GBr/rw0
k0 = 4*GBr/rw0;
ssA = [-1/tau_em, 0, 0, 0, 0;
	0, 0, 1, 0, 0;
	k0/M1, -K/M1, -b/M1, K/M1, 0;
	0, 0, 0, 0, 1;
	0, K/M2, 0, -K/M2, -b/M2];
ssB = [1/tau_em; 0; 0; 0; 0];
ssC = [0, 0, 30/pi*k0, 0, 0];
ssD = [0];
ssGrav01 = ss(ssA,ssB,ssC,ssD);

% close all
% figure; step(ssGrav01,20)
% figure; step(feedback(ssGrav01*sysC,1),20)

%% (07) GravSimple model definition and analysis - multiple depths
format shortg; iC=1;
%Ds = [50]; iC=1; %(m) operating depths / load positions
%Ds = [400]; iC=2; %(m) operating depths / load positions
%Ds = [750]; iC=3; %(m) operating depths / load positions
%Ds = [50; 400; 750]; %(m) operating depths / load positions
Ds = [7; 50; 400; 750]; %(m) operating depths / load positions
%Ds = [50:50:750]'; %(m) operating depths / load positions
M = 50e3; %(kg) load mass
nc = 4; %(-) No of lifting cables
nm = 2; % number of motors

EA = 30.7e6;	%[N] cable axial stiffness
mc_pm = 2.9;	%[kg/m] cable mass per meter
rw0 = 0.562/2;	%[m] winch radius (no cable) = sheave radius
rs0 = 0.480/2;	%[m] sheave radius
tau_em = 60e-3; %(s) motor time constant
GBr = 39.9;		%(-) gearbox ratio
b = 5e3;		%(N/(m/s)) viscous friction coef ~ overall effective damping

Jm = 1.407;		%[kg.m^2] motor plus brake plus gearbox inertias
Jw = 10;		%[kg.m^2] winch drum inertia
Js = 1.6;		%[kg.m2] sheaves inertia

% Calculations
Mmot = 2*Jm*GBr^2/rw0^2;
Mwin = 2*Jw/rw0^2;
Msh = 7*Js/rs0^2;

mc = nc*(Ds*mc_pm);	%(kg) cable mass
alfa = 0.5;			% cable mass split coef
M1 = alfa*mc + Mmot + Mwin + Msh;	% upper mass
M2 = (1-alfa)*mc + M; % lower mass


% Equivalent cable stiffness and Natural frequency
K = nc*(EA./Ds); %(N/m) effective cable stiffness (parallel connection)
fn0 = 1/(2*pi)*sqrt(K./M2); %(Hz) system natural frequency (fixed top)
fn = 1/(2*pi).*sqrt(K.*(M1+M2)./(M1.*M2)); %(Hz) system natural frequency (free top)
zeta = b./(2*sqrt(K.*M2)); %(-) damping

% Balancing (static) motor torque per motor (Nm)
k0 = 4*GBr/rw0;
tauM0 = (M + mc)*g/k0; 
tau_max = tauM0;

% Controller tuning
Kp = 0.1;
Ti = 1000;
Td = 0;

PID.P = Kp;
PID.I = Kp/Ti;
PID.D = Kp*Td;

s = tf('s');
%ltiC = 0.30306*(s+26.45)*(s+1.154)/s/(s+156.9);
%ltiC = C_imc;

% Linear Analysis
% x = [Tm, z1, z1dot, z2, z2dot]
% u = Tref
% y = (30/pi)*(4*GBr/rw0)*z1dot
% Params: tau_em, K, M1, M2, b, k0=4*GBr/rw0
clear ltiG
for i=1:length(Ds)
	ssA = [-1/tau_em, 0, 0, 0, 0;
		0, 0, 1, 0, 0;
		k0/M1(i), -K(i)/M1(i), -b/M1(i), K(i)/M1(i), 0;
		0, 0, 0, 0, 1;
		0, K(i)/M2(i), 0, -K(i)/M2(i), -b/M2(i)];
	ssB = [1/tau_em; 0; 0; 0; 0];
	ssC = [0, 0, 30/pi*k0, 0, 0];
	ssD = [0];
	ltiG{i} = ss(ssA,ssB,ssC,ssD);
end

% PLOTS & CALC
clc; close all
% (1) Relative Inertias
fprintf('\nEquivalent inertias (kg):\n')
fprintf('  2 Motors + gearboxes:\t%3.1f\n',2*Jm*GBr^2/rw0^2)
fprintf('  2 Winch drums:\t\t%3.1f\n',2*Jw/rw0^2)
fprintf('  7 sheaves:\t%3.1f\n',7*Js/rs0^2)
for i=1:length(Ds)
	fprintf('  Cable @%dm:\t\t%3.1f\n',Ds(i),mc(i))
end
fprintf('  Load:\t\t\t\t%3.1f\n',M)

% (2) Natural frequencies etc.
n = length(Ds);
fprintf('\nSystem properties:\n')
fprintf(['  Depth (m):\t\t\t\t' repmat('\t\t%d',1,n) '\n'],Ds)
fprintf(['  Cable mass (t):\t\t\t' repmat('\t\t%3.1f',1,n) '\n'],mc*1e-3)
fprintf(['  Holding motor torque (Nm):\t' repmat('\t%3.1f',1,n) '\n'],tauM0)
fprintf(['  Equiv. cable stiffness (kN/m):' repmat('\t%3.1f',1,n) '\n'],K*1e-3)
fprintf(['  Natural frequency fn (Hz):\t' repmat('\t%3.2f',1,n) '\n'],fn)
fprintf(['  Natural frequency fn0 (Hz):\t' repmat('\t%3.2f',1,n) '\n'],fn0)

if 0
	figure(1)
	plot(Ds,fn0,'bo-'); grid on
	title('System natural frequency (fixed top)')
	ylabel('fn [Hz]'); xlabel('Depth [m]')
	
	figure(101)
	plot(Ds,fn,'bo-'); grid on
	title('System natural frequency (free top)')
	ylabel('fn [Hz]'); xlabel('Depth [m]')

end

% (3) Open-loop responses
if 0
	clear strLeg
	for i=1:n
		figure(2);bode(ltiG{i},{1e-3,1e3}); hold on; grid on
		figure(3); impulse(ltiG{i},5); hold on; grid on
		strLeg{i} = ['D = ' num2str(Ds(i)) 'm'];
	end
	figure(2);legend(strLeg{:},'Location','SouthWest')
	figure(3);legend(strLeg{:})
end

% Eigenvalues
if 0
	eigG = zeros(5,n);
	for i=1:n
		eigG(:,i) = eig(ltiG{i});
	end
	eigG
end


% Controller
clear ltiC

% PID control
load ltiC_PID
% IMC control
%load ltiC_IMC

if 0

	clear strLeg
	for i=1:n
		if (n==1)
			CC = ltiC{iC};
		else
			CC = ltiC{i};
			%CC = ltiC{2};			
		end
		
		figure(10); step(feedback(ltiG{i}*CC,1),10); hold on; grid on
		figure(11); bodemag(feedback(ltiG{i}*CC,1),{1e-3,1e3}); hold on; grid on		
		figure(12); nichols(ltiG{i}*CC); hold on; grid on				
		figure(13); bodemag(feedback(CC,ltiG{i}),{1e-3,1e3}); hold on; grid on	
		strLeg{i} = ['D = ' num2str(Ds(i)) 'm'];
	end
	figure(10);legend(strLeg{:})
	figure(11);legend(strLeg{:},'Location','SouthWest')
	figure(12);legend(strLeg{:})
	figure(13);legend(strLeg{:},'Location','SouthWest')
	
end

% Simulations (GravSimple model)
sysC = ltiC{iC};
if 0
	
	
end



%% (10) Deep Shaft analysis and control

clear ltiC
%load ltiC_PID % PID control
load ltiC_IMC % IMC control

format shortg; iC=1;
%Ds = [7]; iC=1; %(m) operating depths / load positions
%Ds = [50]; iC=2; %(m) operating depths / load positions
%Ds = [400]; iC=3; %(m) operating depths / load positions
Ds = [750]; iC=4; %(m) operating depths / load positions
%Ds = [50; 400; 750]; %(m) operating depths / load positions
%Ds = [7; 50; 400; 750]; %(m) operating depths / load positions
%Ds = [50:50:750]'; %(m) operating depths / load positions
M = 50e3; %(kg) load mass
nc = 4; %(-) No of lifting cables
nm = 2; % number of motors

EA = 30.7e6;	%[N] cable axial stiffness
mc_pm = 2.9;	%[kg/m] cable mass per meter
rw0 = 0.562/2;	%[m] winch radius (no cable) = sheave radius
rs0 = 0.480/2;	%[m] sheave radius
tau_em = 60e-3; %(s) motor time constant
GBr = 39.9;		%(-) gearbox ratio
b = 1e4;		%(N/(m/s)) viscous friction coef ~ overall effective damping

Jm = 1.407;		%[kg.m^2] motor plus brake plus gearbox inertias
Jw = 10;		%[kg.m^2] winch drum inertia
Js = 1.6;		%[kg.m2] sheaves inertia

% Calculations
Mmot = 2*Jm*GBr^2/rw0^2;
Mwin = 2*Jw/rw0^2;
Msh = 6*Js/rs0^2;

mc = nc*(Ds*mc_pm);	%(kg) cable mass
alfa = 0.5;			% cable mass split coef
M1 = alfa*mc + Mmot + Mwin + Msh;	% upper mass
M2 = (1-alfa)*mc + M; % lower mass

% Equivalent cable stiffness and Natural frequency
K = nc*(EA./Ds); %(N/m) effective cable stiffness (parallel connection)
fn0 = 1/(2*pi)*sqrt(K./M2); %(Hz) system natural frequency (fixed top)
fn = 1/(2*pi).*sqrt(K.*(M1+M2)./(M1.*M2)); %(Hz) system natural frequency (free top)
zeta = b./(2*sqrt(K.*M2)); %(-) damping

% Balancing (static) motor torque per motor (Nm)
k0 = GBr/rw0;
tauM0 = (M + mc)*g/4/k0; 
tau_max = tauM0;
% Static pin load (t)
Fpin0 = 4*k0*tauM0/2/(g*1000);


% Controller tuning
Kp = 0.1;
Ti = 1000;
Td = 0;

PID.P = Kp;
PID.I = Kp/Ti;
PID.D = Kp*Td;

s = tf('s');
%ltiC = 0.30306*(s+26.45)*(s+1.154)/s/(s+156.9);
%ltiC = C_imc;

% Linear Analysis
% x = [Tm, z1, z1dot, z2, z2dot]
% u = Tref
% y = (30/pi)*(4*GBr/rw0)*z1dot
% Params: tau_em, K, M1, M2, b, k0=4*GBr/rw0
clear ltiG ltiG2
for i=1:length(Ds)
	ssA = [-1/tau_em, 0, 0, 0, 0;
		0, 0, 1, 0, 0;
		4*k0/M1(i), -K(i)/M1(i), -b/M1(i), K(i)/M1(i), 0;
		0, 0, 0, 0, 1;
		0, K(i)/M2(i), 0, -K(i)/M2(i), -b/M2(i)];
	ssB = [1/tau_em; 0; 0; 0; 0];
	ssC = [0, 0, 30/pi*2*k0, 0, 0];
	ssC2 = [0, 0, 0, 0, 30/pi*2*k0];
	ssD = [0];
	ltiG{i} = ss(ssA,ssB,ssC,ssD);
	ltiG2{i} = ss(ssA,ssB,ssC2,ssD);	
end

% PLOTS & CALC
clc; close all
% (1) Relative Inertias
fprintf('\nEquivalent inertias (kg):\n')
fprintf('  2 Motors + gearboxes:\t%3.1f\n',2*Jm*GBr^2/rw0^2)
fprintf('  2 Winch drums:\t\t%3.1f\n',2*Jw/rw0^2)
fprintf('  6 sheaves:\t\t%3.1f\n',6*Js/rs0^2)
for i=1:length(Ds)
	fprintf('  Cable @%dm:\t\t%3.1f\n',Ds(i),mc(i))
end
fprintf('  Load:\t\t\t\t%3.1f\n',M)

% (2) Natural frequencies etc.
n = length(Ds);
fprintf('\nSystem properties:\n')
fprintf(['  Depth (m):\t\t\t\t' repmat('\t\t%d',1,n) '\n'],Ds)
fprintf(['  Cable mass (t):\t\t\t' repmat('\t\t%3.1f',1,n) '\n'],mc*1e-3)
fprintf(['  Holding motor torque (Nm):\t\t' repmat('%3.0f\t\t',1,n) '\n'],tauM0)
fprintf(['  Static Pin load (t):\t\t\t' repmat('\t%3.1f',1,n) '\n'],Fpin0)
fprintf(['  Equiv. cable stiffness (kN/m)\t\t' repmat('%3.0f\t',1,n) '\n'],K*1e-3)
fprintf(['  Natural frequency fn (Hz):\t' repmat('\t%3.2f',1,n) '\n'],fn)
fprintf(['  Natural frequency fn0 (Hz):\t' repmat('\t%3.2f',1,n) '\n'],fn0)
fprintf(['  Damping ratio zeta (-):\t\t' repmat('\t%3.3f',1,n) '\n'],zeta)


%%% Natural frequency plots
if 0
	figure(1)
	plot(Ds,fn0,'bo-'); grid on
	title('System natural frequency (fixed top)')
	ylabel('fn [Hz]'); xlabel('Depth [m]')
	
	figure(101)
	plot(Ds,fn,'bo-'); grid on
	title('System natural frequency (free top)')
	ylabel('fn [Hz]'); xlabel('Depth [m]')

end


%%% Eigenvalues
if 0
	eigG = zeros(5,n);
	for i=1:n
		eigG(:,i) = eig(ltiG{i});
	end
	eigG
end


%%% Open-loop responses
if 0
	close all
	clear strLeg
	for i=1:n
		figure(2);h2=bodeplot(ltiG{i},{1e-3,1e4}); hold on; grid on
		figure(3);h3=stepplot(ltiG{i},10); hold on; grid on
		figure(4);h4=impulseplot(ltiG{i},10); hold on; grid on
		strLeg{i} = ['D = ' num2str(Ds(i)) 'm'];
	end
	
	figure(2); hax = findobj(gcf,'type','axes');
	axes(hax(2)); legend(strLeg{:},'Location','SouthWest')
	options = getoptions(h2);
	options.FreqUnits = 'Hz';
	options.XLim = [1e-3 1e2];
	options.Title.String = 'Bode Diagram: From Tref (Nm) to n (rpm)';
	setoptions(h2,options);
	
	figure(3);legend(strLeg{:},'Location','SouthEast')
	options = getoptions(h3);
	options.Title.String = 'Step Response: From Tref (Nm) to n (rpm)';
	options.YLabel.String = 'Amplitude (rpm)';
	setoptions(h3,options);
	
	figure(4);legend(strLeg{:},'Location','NorthEast')
	options = getoptions(h4);
	options.Title.String = 'Impulse Response: From Tref (Nm) to n (rpm)';
	options.YLabel.String = 'Amplitude (rpm)';	
	setoptions(h4,options);
	
	
	if 0
		figure(22)
		bode(ltiG{3},ltiG2{3},{1e-3,1e3});  hold on; grid on
		legend('z1','z2')
		figure(23)
		impulse(ltiG{3},ltiG2{3},5);  hold on; grid on
		legend('z1','z2')
		
	end
end


%%% Closed-loop control
if 0

	MotTF = 1/(tau_em*s+1);
	clear strLeg
	for i=1:n
		if (n==1)
			%CC = ltiC{iC};
			CC = minreal(2*ltiC{iC});
		else
			CC = minreal(2*ltiC{i})
			%CC = ltiC{3};			
		end
		figure(11); h1=stepplot(feedback(ltiG{i}*CC,1),10); hold on; grid on
		figure(12); h2=stepplot(feedback(CC,ltiG{i}),10); hold on; grid on		
		figure(13); h3=bodeplot(feedback(ltiG{i}*CC,1),{1e-3,1e4}); hold on; grid on		
		figure(14); h4=bodeplot(feedback(CC,ltiG{i}),{1e-3,1e4}); hold on; grid on	
		figure(15); h5=nicholsplot(ltiG{i}*CC); hold on; grid on				
		figure(16); h6=stepplot(feedback(CC,ltiG{i})*MotTF,10); hold on; grid on				
		strLeg{i} = ['D = ' num2str(Ds(i)) 'm'];
	end
	
	figure(11);legend(strLeg{:},'Location','SouthEast')
	options = getoptions(h1);
	options.Title.String = 'Step Response: From nref (rpm) to n (rpm)';
	options.YLabel.String = 'Amplitude (rpm)';
	setoptions(h1,options);

	figure(12);legend(strLeg{:},'Location','NorthEast')
	options = getoptions(h2);
	options.Title.String = 'Step Response: From nref (rpm) to Tref (Nm)';
	options.YLabel.String = 'Amplitude (Nm)';
	setoptions(h2,options);

	figure(16);legend(strLeg{:},'Location','NorthEast')
	options = getoptions(h6);
	options.Title.String = 'Step Response: From nref (rpm) to Tm (Nm)';
	options.YLabel.String = 'Amplitude (Nm)';
	setoptions(h6,options);
	
	figure(13);	hax = findobj(gcf,'type','axes');
	axes(hax(2)); legend(strLeg{:},'Location','SouthWest')
	setoptions(h3,'FreqUnits','Hz','PhaseVisible','off'); set(gca,'XLim',[1e-3 1e2])	
	options = getoptions(h3);
	options.Title.String = 'Bode Diagram: From nref (rpm) to n (rpm)';
	setoptions(h3,options);

	figure(14);hax = findobj(gcf,'type','axes');
	axes(hax(2)); legend(strLeg{:},'Location','NorthWest')
	setoptions(h4,'FreqUnits','Hz','PhaseVisible','off');set(gca,'XLim',[1e-3 1e2])	
	options = getoptions(h4);
	options.Title.String = 'Bode Diagram: From nref (rpm) to Tref (Nm)';
	setoptions(h4,options);

	figure(15);legend(strLeg{:})	
	options = getoptions(h5);
	options.Title.String = 'Nichols Chart of open-loop transfer function G(s)*IMC(s)';
	setoptions(h5,options);
	
	
end


%%% Simulations (GravSimple model)
%sysC = ltiC{iC};
sysC = minreal(2*ltiC{iC});
if 0
	
	
end


%% Validation LTI step vs GravSimple
figure
subplot(121)
plot(h1.Responses.Data.Time,h1.Responses.Data.Amplitude*1000)
hold on
plot(out.ScopeData.time,out.ScopeData.signals(1).values(:,2))
set(gca,'XLim',[0 10])
legend('LTI step','sim','Location','South')
title('Motor speed (rpm)'); xlabel('time (sec)')
%figure
subplot(122)
plot(h2.Responses.Data.Time,h2.Responses.Data.Amplitude*1000)
hold on
plot(out.ScopeData.time,(out.ScopeData.signals(2).values(:,2)-tauM0))
set(gca,'XLim',[0 10])
title('Motor torque deviations from nominal (Nm)')
xlabel('time (sec)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.2, 0.5, 0.4]);








