% units: time ~ 1 ps = 1, length ~ 1 Å = 1;

dt = 4;
dx = 3;
N = 45;
% propagating step 1 ps, diffusion layer thickness 3 Å, number of diffusion
% layers in the container
% diffusion constant air in water 1.2E-5 cm^2·s^-1
D = 0.12;

% propagation term between the adjacent layers 
a = D*dt/(dx*dx);

% initial concentration and length of the box, 10 M, 15 Å
c0 = 6E-3;
hbx = 5;
h = hbx*dx;

%% initial conentration table, box range C0
conc = zeros(N,1);

for i = 1:hbx
    conc(i) = c0;
end


%% setting the initial propagator
ppgh = zeros(N,N);

% the triple diagonal propagation matrix for each species
for i = 1:N-1
    ppgh(i,i+1) = a;
    ppgh(i+1,i) = a;
    ppgh(i+1,i+1) = 1-2*a;

end

% the outer surface of the container
ppgh(1,1) = 1-a;

% the container stops at the N-th box, no diffusion backward
ppgh(N,N) = 1-a;


%% analytical evaluation parameters
t = 0;
conca = zeros(N,1);
x = dx*(0:N-1)';

%% numerical propagation & analytical evaluation


for j = 1:1E2
%   apply fick's second law to the concentration table
    conc = ppgh*conc;
%   timer for later analytical use
    t = t + dt;
end

%     clf;
conctable1mol = conc * 1/(6E-4);
line(x, conctable1mol, 'Color', 'red');

% analytical evaluation

sigma = sqrt(D*t);
conca = c0*(erfc(0.5*x/sigma) - erfc(0.5*(h-x)/sigma));

conctable1mola = conca * 1/(6E-4);
line(x, conctable1mola, 'Color', 'black');