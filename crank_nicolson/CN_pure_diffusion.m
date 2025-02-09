% units: time ~ 1 ps = 1, length ~ 1 Å = 1;
% diffusion base unite 1E-4 cm^2·s^-1 = 1;
% concentration base unit 1667 mol/L or 1 per Å^3 = 1;

dt = 1000;
dx = 2;
K = 1000;
% propagating step 1 ps, diffusion layer thickness 2 Å, number of diffusion
% layers in the container
% diffusion constant: air in water 1.2E-5 cm^2·s^-1
D = 0.12;

% propagation term between the adjacent layers 
a = D*dt/(dx*dx);

% initial concentration 1M/1667 = 6E-4
% length of the box 40 Å = number of box hbx 20 * box width 2 Å
c0 = 6E-4;
hbx = 20;
h = hbx*dx;

%% initial conentration table, box range C0
conc = zeros(K,1);

for i = 1:hbx
    conc(i) = c0;
end


%% setting the propagators for cn and cn+1
ppgn = zeros(K,K);
ppgn1 = zeros(K,K);

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgn(i,i+1) = a/2;
    ppgn(i+1,i) = a/2;
    ppgn(i+1,i+1) = 1-a;

end

% start diagonal term of cn
ppgn(1,1) = 1-a/2;

% end diagonal term of cn
ppgn(K,K) = 1-a/2;

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgn1(i,i+1) = -a/2;
    ppgn1(i+1,i) = -a/2;
    ppgn1(i+1,i+1) = 1+a;

end

% start diagonal term of cn
ppgn1(1,1) = 1+a/2;

% end diagonal term of cn
ppgn1(K,K) = 1+a/2;

ppgcn = ppgn1\ppgn;

%% analytical evaluation parameters
t = 0;
concana = zeros(K,1);
x = dx*(0:K-1)';

%     clf;
conctable1mol = conc * 1/(6E-4);
line(x, conctable1mol, 'Color', 'red');

% analytical evaluation

sigma = sqrt(D*t);
concana = c0/2*(erf(0.5*(h-x)/sigma) + erf(0.5*(h+x)/sigma));

conctable1mola = concana * 1/(6E-4);
line(x, conctable1mola, 'Color', 'black');

%% numerical propagation & analytical evaluation

for j = 1:4E2
%   apply fick's second law to the concentration table
    conc = ppgcn*conc;
%   timer for later analytical use
    t = t + dt;
end

%     clf;
conctable1mol = conc * 1/(6E-4);
line(x, conctable1mol, 'Color', 'red');

% analytical evaluation

sigma = sqrt(D*t);
concana = c0/2*(erf(0.5*(h-x)/sigma) + erf(0.5*(h+x)/sigma));

conctable1mola = concana * 1/(6E-4);
line(x, conctable1mola, 'Color', 'black');