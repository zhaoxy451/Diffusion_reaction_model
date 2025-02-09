% units: time ~ 1 ps = 1, length ~ 1 Å = 1;

dt = 1;
dx = 3;
N = 33;
% propagating step 1 ps, diffusion layer thickness 3 Å, number of diffusion
% layers in the container
% diffusion constant air in water 1.2E-5 cm^2·s^-1
D = 0.12;

% propagation term between the adjacent layers 
ae = D*dt/(dx*dx);


% equilibrium concentration-evaporation rate constant:
% bulk concentration to evaporative flux, as a reference CO2~1E-11
kevp = 1E-11;


%% setting the external terms
% incoming flux as boundary condition, effusive beam pressure 0.5 Torr, 
% 2 inches away from an Ø1 mm orifice, 7E19 #·m^-2·s^-1, convert unite to
% Å^-2·ps^-1, equivalent pressure ~0.5 Torr * 5E-5
F = 7E-13;
exte = zeros(N,1);
exth = zeros(N,1);
exth2 = zeros(N,1);

% incoming flux added to the first box
exte(1) = F*dt/dx;

%% initial conentration table: all zero
econc = zeros(N,1);
hconc = zeros(N,1);
h2conc = zeros(N,1);

%% setting the initial propagator
ppgh = zeros(N,N);

% the triple diagonal propagation matrix for each species
for i = 1:N-1
    ppgh(i,i+1) = ae;
    ppgh(i+1,i) = ae;
    ppgh(i+1,i+1) = 1-2*ae;

end

% the outer surface of the container
ppgh(1,1) = 1-ae;

% the container stops at the N-th box, no diffusion backward
ppgh(N,N) = 1-ae;



%% propagation

for j = 1:1E4
%   apply fick's second law to the concentration table
    econc = ppgh*econc + exte;
%   evaporation from the first layer
    econc(1) = econc(1) - kevp*econc(1)*dt/dx;

%     t = dt*j/1000;
%     tstring = num2str(t);
%     tstring = strcat(tstring,' ns');
%     lgd = legend;
%     lgd.Title.String = tstring;
    
%     pause(0.00005);
end
%     clf;
    conctable1mol = econc * 1/(6E-4);
    depth = dx*(1:N);
    line(depth, conctable1mol);