%% ECE 4784 %%
%Project Phase 1
%William Zhang
%Due September, 29th 2014

%% Constants Provided:
simTtot = 100; %100 ms total simulation time
step = input('Enter desired time step resolution for simulation [must be below .05ms]: '); %degree of freedom for resolution
t = 0 : step : simTtot; %create time vector using inputted step between data points
I = zeros(1, length(t)); % Current vector must be same size as time vector

%Maximum Conductances
gKBAR = 36; %36 mS/cm^2
gNaBAR = 120; %120 mS/cm^2
gLBAR = 0.3; %0.3 mS/cm^2
EK = -12; %-12 mV
ENa = 115; %115 mV
EL = 10.6; %10.6 mV
VRest = -70; %-70 mV resting potential of membrane.
Cm = 1.0; %uF/cm^2

%%Equations
%Gating Variables:
V(1) = 0; %Initial assumption that V = 0;
alpham(1) = 0.1.*((25-V(1))/(exp((25-V(1))/10)-1));
betam(1) = 4*exp(-V(1)/18);
alphan(1) = 0.01.*((10-V(1))/(exp((10-V(1))/10)-1));
betan(1) = 0.125*exp(-V(1)/80);
alphah(1) = 0.07*exp(-V(1)/20);
betah(1) = 1/(exp((30-V(1))/10)+1);

%Initial Values:
m(1) = alpham(1)/(alpham(1)+betam(1));
n(1) = alphan(1)/(alphan(1)+betan(1));
h(1) = alphah(1)/(alphah(1)+betah(1));

%% Current - Input based for stimulation purposes
a = input('Desired amplitude of current? (in uA/cm^2): '); %inputs asking for desired characteristics
d = input('Desired duration of current? must be between 0 and 100 ms): ');
for i = 1:d/step % filling in vector where current is desired.
I(i) = a;
end


%% Euler's Method
j = 1; %initialize loop index
while j < length(t)
    alpham = 0.1*((25-V(j))/(exp((25-V(j))/10)-1));
    betam = 4*exp(-V(j)/18);
    alphan = 0.01*((10-V(j))/(exp((10-V(j))/10)-1));
    betan = 0.125*exp(-V(j)/80);
    alphah = 0.07*exp(-V(j)/20);
    betah(1) = 1/(exp((30-V(j))/10)+1);
    
    gK(j) = n(j)^4*gKBAR; %Membrane conductances - finding based n, m, and h
    gNa(j) = m(j)^3*gNaBAR *h(j);
    gL(j) = gLBAR;
    
    iK(j) = gK(j)*(V(j)-EK); %Current equation definitions
    iNa(j) = gNa(j)*(V(j)-ENa);
    iL(j) = gL(j)*(V(j)-EL);
    iIon(j) = I(j)-iNa(j)-iK(j)-iL(j);
    
    derivm = alpham*(1-m(j))-betam*m(j); %Euler's method - determining derivative
    derivn = alphan*(1-n(j))-betan*n(j); %based on curent alpha and beta values
    derivh = alphah*(1-h(j))-betah*h(j);
    
    m(j+1) = m(j)+step*derivm; %Stepping m, n, and h values based on the derivation
    n(j+1) = n(j)+step*derivn;
    h(j+1) = h(j)+step*derivh;
    
    derivV(j) = iIon(j)/Cm;
    V(j+1) = V(j)+step*derivV(j); %final voltage calculation based on equation given to us
    
    j=j+1;
end

%% Plots of Data
%Voltage
V = V + VRest; %sets action potential starting voltage to -70mV 
subplot(2,1,1) %subplot for both graphs on one figure
plot(t,V) 
title('Membrane Potential'); %labels
xlabel('Time [ms]');
ylabel('Membrane Voltage [mV]');
legend('Vm');
axis([0, 100, -100, 40]);

%Conductances
t2 = t(1:length(t)-1); % Time vector must match length of conductances
subplot(2,1,2);
plot(t2, gNa, 'g', t2, gK, 'b', t2, gL, 'r')
title('gK, gNa, gL') %lables
xlabel('Time [ms])')
ylabel('Conductance [mS/cm^2]')
legend('gNa','gK', 'gL')
