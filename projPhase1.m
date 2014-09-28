%ECE 4784
%Project Phase 1
%William Zhang
%Due September, 29th 2014

%%Constants Provided:
simTtot = 100; %100 ms total simulation time
step = .02;
t = 0 : step : simTtot; %create time vector using .01 steps between data points

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
betah(1) = 1/((exp(30-V(1))/10)+1);

%Initial Values:
m(1) = alpham(1)/(alpham(1)+betam(1));
n(1) = alphan(1)/(alphan(1)+betan(1));
h(1) = alphah(1)/(alphah(1)+betah(1));

%%Injected Current
amp = input('Enter the amplitude of the injected current (in uA/cm^2): ');
dur = input('Enter the duration of the injected current [between 0 and 100 ms]: ');
I = zeros(1, length(t)); % Setting up I to be the same size as the time array
for y = 1:dur/step % injected current goes from 0 to (duration entered / time step), elsewhere is 0
I(y) = amp;
end


for j = 1 : length(t)-1
    alpham = 0.1.*((25-V(j))/(exp((25-V(j))/10)-1));
    betam = 4*exp(-V(j)/18);
    alphan = 0.01.*((10-V(j))/(exp((10-V(j))/10)-1));
    betan = 0.125*exp(-V(j)/80);
    alphah = 0.07*exp(-V(j)/20);
    betah = 1/((exp(30-V(j))/10)+1);
    
    gK(j) = n(j)^4*gKBAR;
    gNa(j) = m(j)^3*gNaBAR *h(j);
    gL(j) = gLBAR;
    
    iK(j) = gK(j)*(V(j)-EK);
    iNa(j) = gNa(j)*(V(j)-ENa);
    iL(j) = gL(j)*(V(j)-EL);
    iIon(j) = I(j)-iNa(j)-iK(j)-iL(j);
    
    derivm = alpham*(1-m(j))-betam*m(j);
    derivn = alphan*(1-n(j))-betan*n(j);
    derivh = alphah*(1-h(j))-betah*h(j);
    
    m(j+1) = m(j)+step.*derivm;
    n(j+1) = n(j)+step.*derivn;
    h(j+1) = h(j)+step.*derivh;
    
    derivV(j) = iIon(j)/Cm;
    V(j+1) = V(j)+step*derivV(j);
    
end
