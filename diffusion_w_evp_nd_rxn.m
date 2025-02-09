%% physical constants and algorithm constants
% units: time ~ 1 ps = 1, length ~ 1 Å = 1;

dt = 1;
dx = 3;
N = 33;

% propagating step 1 ps, diffusion layer thickness 3 Å, number of diffusion
% layers in the container
% diffusion constant air in water 1.2E-5 cm^2·s^-1
De = 0.5;
Dh = 0.8;
Dh2 = 0.45;
% propagation term between the adjacent layers 
ae = De*dt/(dx*dx);
ah = Dh*dt/(dx*dx);
ah2 = Dh2*dt/(dx*dx);

% equilibrium concentration-evaporation rate constant:
% bulk concentration to evaporative flux, as a reference CO2~1E-11
kevph2 = 1E-8;
% k(e+h3o+) = 0.55E10 M^-1·s^-1, k(e+h2o) = 1E3 s^-1
% pH = 7, ke = 3E-9; pH = 6, ke = 2E-8;
% pH = 8, ke = 1E-9; pH = 9, ke = 1E-9;
% pH = 1
ke = 2E-3;
kee = 9.4;
khh = 8.5;

%% setting the initial external terms
% incoming flux as boundary condition, effusive beam pressure 0.5 Torr, 
% 2 inches away from an Ø1 mm orifice, 7E19 #·m^-2·s^-1, convert unite to
% Å^-2·ps^-1, equivalent pressure ~0.5 Torr * 5E-5
FNa = 7E-13;
exte = zeros(N,1);
exth = zeros(N,1);
exth2 = zeros(N,1);


%% initial conentration table: all zero
econc = zeros(N,1);
hconc = zeros(N,1);
h2conc = zeros(N,1);

%% setting the initial propagator for electron
ppge = zeros(N,N);

% the triple diagonal propagation matrix
for i = 1:N-1
    ppge(i,i+1) = ae;
    ppge(i+1,i) = ae;
    ppge(i+1,i+1) = 1-2*ae;

end

% the outer surface of the container
ppge(1,1) = 1-ae;

% the container stops at the N-th box, no diffusion backward
ppge(N,N) = 1-ae;

%% setting the initial propagator for hydrogen radical
ppgh = zeros(N,N);

% the triple diagonal propagation matrix
for i = 1:N-1
    ppgh(i,i+1) = ah;
    ppgh(i+1,i) = ah;
    ppgh(i+1,i+1) = 1-2*ah;

end

% the outer surface of the container
ppgh(1,1) = 1-ah;

% the container stops at the N-th box, no diffusion backward
ppgh(N,N) = 1-ah;

%% setting the initial propagator for hydrogen radical
ppgh2 = zeros(N,N);

% the triple diagonal propagation matrix
for i = 1:N-1
    ppgh2(i,i+1) = ah2;
    ppgh2(i+1,i) = ah2;
    ppgh2(i+1,i+1) = 1-2*ah2;

end

% the outer surface of the container
ppgh2(1,1) = 1-ah2;

% the container stops at the N-th box, no diffusion backward
ppgh2(N,N) = 1-ah2;

%% propagation

for j = 1:1E7
%   reset external terms every time
    exte = zeros(N,1);
    exth = zeros(N,1);
    exth2 = zeros(N,1);
    
%   incoming flux added to the first box
    exte(1) = FNa*dt/dx;
%   modify electron concentration with incident flux & rxn
    econc = econc + exte - dt*ke*econc - dt*kee*econc.*econc;
%   apply fick's second law to electron concentration table
    econc = ppge*econc;
    
%   rxn: modify hydrogen radical concentration and dihydrogen concentration
%   add chemically generated concentration to external terms
    exth = dt*ke*econc;
    hconc = hconc + exth;
    exth2 = dt*kee*econc.*econc + dt*khh*hconc.*hconc;
%   dihydrogen evaporation from its first layer
    exth2(1) = exth2(1) - kevph2*dt*h2conc(1)/dx;
    h2conc = h2conc + exth2;
    
%   apply diffusion to hydrogen radical and dihydrogen concentration tables
    hconc = ppgh*hconc;
    h2conc = ppgh2*h2conc;
    

%     t = dt*j/1000;
%     tstring = num2str(t);
%     tstring = strcat(tstring,' ns');
%     lgd = legend;
%     lgd.Title.String = tstring;
    
%     pause(0.00005);
end
%     clf;
    conctable1mol = h2conc;
%     * 1/(6E-4);
    depth = dx*(1:N);
    line(depth, conctable1mol);