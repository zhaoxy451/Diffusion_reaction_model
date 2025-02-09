% units: time ~ 1 ps = 1, length ~ 1 Å = 1;

dt = 0.1;
dx = 3;
% diffusion constant air in water 1.2E-5 cm^2·s^-1
D1 = 0.12;
D2 = 0.12;
a1 = D1*dt/(dx*dx);
a2 = D2*dt/(dx*dx);

% equilibrium concentration
c0 = 0.5;
kevp = 0.03;

%% setting the external terms
% incoming flux as boundary condition, effusive beam pressure 0.5 Torr, 
% 2 inches away from an Ø1 mm orifice, 7E19 #·m^-2·s^-1
F = 0.007;
ext1 = zeros(33,1);
ext2 = zeros(33,1);

% incoming flux added to the first box
ext1(1) = F*dt/dx;


%% initial conentration table: all zero
conctable1 = zeros(33,1);

%% setting the initial propagator
ppg1 = zeros(33,33);

for i = 1:32
    ppg1(i,i+1) = a1;
    ppg1(i+1,i) = a1;
    ppg1(i+1,i+1) = 1-2*a1;

end

ppg1(1,1) = 1-a1;

% secondary approximation on conctable(34)
ppg1(33,31) = 0.5*a1;
ppg1(33,32) = -a1;
ppg1(33,33) = 1+0.5*a1;

%% propagation

for j = 1:1E4
    conctable1 = ppg1*conctable1 + ext1;
    plot(conctable1);
 
    t = dt*j/1000;
    tstring = num2str(t);
    tstring = strcat(tstring,' ns');
    lgd = legend;
    lgd.Title.String = tstring;
    
    pause(0.00005);
end

