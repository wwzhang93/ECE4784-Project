%ECE 4784
%Project Phase 1
%William Zhang
%Due September, 29th 2014

%Constants Provided:
simTtot = 100; %100 ms total simulation time
t = 0:.01:simTtot; %create time vector using .01 steps between data points

gK = 36; %36 mS/cm^2
gNa = 120; %120 mS/cm^2
gL = 0.3; %0.3 mS/cm^2
EK = -12; %-12 mV
ENa = 115; %115 mV
EL = 10.6; %10.6 mV
VRest = -70; %-70 mV resting potential of membrane.

%Equations
%Gating Variables:
V = 0; %Initial assumption that V = 0;
alpham = 0.1.*((25-V)/(exp((25-V)/10)-1));
betam = 4*exp(-V/18);
alphan = 0.01.*((10-V)/(exp((10-V)/10)-1));
betan = 0.125*exp(-V/80);
alphah = 0.07*exp(-V/20);
betah = 1/((exp(30-V)/10)+1);

%Initial Values:
m(1) = alpham/(alpham-betam);
n(1) = alphan/(alphan-betan);
h(1) = alphah/(alphah-betah);

%Currents:
% I = zeros(length(t));
%I(1:50) = 5;
%curStim = 10;
%I = ones(length(t)) .* curStim;

%iNa = m.^3.*gNa.*h.*(V-ENa);
%iK = n.^4*gK*(V-EK);
%iL = gL*(V-EL);
%iIon = I - iK - iNa - iL;
