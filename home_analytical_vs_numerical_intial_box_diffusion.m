% units: time ~ 1 ps = 1, length ~ 1 Å = 1;
% diffusion base unite 1E-4 cm^2·s^-1 = 1;
% concentration base unit 1667 mol/L or 1 per Å^3 = 1;

dt = 10;
dx = 2;
K = 100;
% propagating step 1 ps, diffusion layer thickness 2 Å, number of diffusion
% layers in the container
% diffusion constant: air in water 1.2E-5 cm^2·s^-1
D = 0.12;

% propagation term between the adjacent layers 
a = D*dt/(dx*dx);

% initial concentration 1M/1667 = 6E-4
% length of the box 20 Å = number of box hbx * box width 2 Å
c0 = 6E-4;
hbx = 20;
h = hbx*dx;

%% initial conentration table, box range C0
conc = zeros(K,1);

for i = 1:hbx
    conc(i) = c0;
end


%% setting the initial propagator
ppg = zeros(K,K);

% the triple diagonal propagation matrix for each species
for i = 1:K-1
    ppg(i,i+1) = a;
    ppg(i+1,i) = a;
    ppg(i+1,i+1) = 1-2*a;

end

% the outer surface of the container
ppg(1,1) = 1-a;

% the container stops at the K-th box, no diffusion backward
ppg(K,K) = 1-a;


%% analytical evaluation parameters
t = 0;
concana = zeros(K,1);
x = dx*(0:K-1)';


conctable1mol = conc * 1/(6E-4);
line(x, conctable1mol, 'Color', 'red');

% analytical evaluation

sigma = sqrt(D*t);
concana = c0/2*(erf(0.5*(h-x)/sigma) + erf(0.5*(h+x)/sigma));

% xd = x - dx;
conctable1mola = concana * 1/(6E-4);
line(x, conctable1mola, 'Color', 'black');

%% numerical propagation & analytical evaluation

for j = 1:4E2
%   apply fick's second law to the concentration table
    conc = ppg*conc;
%   timer for later analytical use
    t = t + dt;
end

%     clf;
conctable1mol = conc * 1/(6E-4);
line(x, conctable1mol, 'Color', 'red');

% analytical evaluation

sigma = sqrt(D*t);
concana = c0/2*(erf(0.5*(h-x)/sigma) + erf(0.5*(h+x)/sigma));

% xd = x - dx;
conctable1mola = concana * 1/(6E-4);
line(x, conctable1mola, 'Color', 'black');