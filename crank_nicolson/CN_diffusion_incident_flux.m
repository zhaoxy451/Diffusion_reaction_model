% units: time ~ 1 ps = 1, length ~ 1 Å = 1;
% diffusion base unite 1E-4 cm^2·s^-1 = 1;
% concentration base unit 1667 mol/L or 1 per Å^3 = 1;

dt = 100;
dx = 2;
K = 400;
% propagating step 1 ps, diffusion layer thickness 2 Å, number of diffusion
% layers in the container
% diffusion constant: hydrated electron in water 5E-5 cm^2·s^-1
De = 0.5;
Dh = 0.8;
Dhh = 0.45;

% propagation term between the adjacent layers 
ae = De*dt/(dx*dx);
ah = Dh*dt/(dx*dx);
ahh = Dhh*dt/(dx*dx);

%% initial conentration table, box range C0
% initial concentration 1M/1667 = 6E-4

conce = zeros(K,1);
conch = zeros(K,1);
conchh = zeros(K,1);

pH = 1;
c0p = 10^(-pH);
concp = ones(K,1) * c0p;

% for i = 1:K
%     conc(i) = c0;
% end

%%  incoming flux 
%   650 K, Na incident flux on microjet 3E19 #·m^-2·s^-1
FNa = 3E-13;


%%  define chemistry

% kee = 6E9 M^-1·s^-1
% khh = 5E9 M^-1·s^-1
% kep = 2E10 M^-1·s^-1
% keh = 2E10 M^-1·s^-1
% ke = 1E3 s^-1
% unit conversion 1 M^-1·s^-1 = 6E-16 Å^-3·ps^-1
% unit conversion 1 s^-1 = 1E-12 ps^-1

kee = 3.6E-6;
khh = 3E-6;
kep = 1.2E-5;
keh = 1.2E-5;
ke = 1E-9;

%% setting the propagators for cn and cn+1 (es, hydrated electron)
ppgne = zeros(K,K);
ppgn1e = zeros(K,K);

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgne(i,i+1) = ae/2;
    ppgne(i+1,i) = ae/2;
    ppgne(i+1,i+1) = 1-ae;

end

% start diagonal term of cn
ppgne(1,1) = 1-ae/2;

% end diagonal term of cn
ppgne(K,K) = 1-ae/2;

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgn1e(i,i+1) = -ae/2;
    ppgn1e(i+1,i) = -ae/2;
    ppgn1e(i+1,i+1) = 1+ae;

end

% start diagonal term of cn
ppgn1e(1,1) = 1+ae/2;

% end diagonal term of cn
ppgn1e(K,K) = 1+ae/2;

ppgcne = ppgn1e\ppgne;

%% setting the propagators for cn and cn+1 (H·, hydrogen radical)
ppgnh = zeros(K,K);
ppgn1h = zeros(K,K);

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgnh(i,i+1) = ah/2;
    ppgnh(i+1,i) = ah/2;
    ppgnh(i+1,i+1) = 1-ah;

end

% start diagonal term of cn
ppgnh(1,1) = 1-ah/2;

% end diagonal term of cn
ppgnh(K,K) = 1-ah/2;

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgn1h(i,i+1) = -ah/2;
    ppgn1h(i+1,i) = -ah/2;
    ppgn1h(i+1,i+1) = 1+ah;

end

% start diagonal term of cn
ppgn1h(1,1) = 1+ah/2;

% end diagonal term of cn
ppgn1h(K,K) = 1+ah/2;

ppgcnh = ppgn1h\ppgnh;

%% setting the propagators for cn and cn+1 (H2, dihydrogen gas)
ppgnhh = zeros(K,K);
ppgn1hh = zeros(K,K);

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgnhh(i,i+1) = ahh/2;
    ppgnhh(i+1,i) = ahh/2;
    ppgnhh(i+1,i+1) = 1-ahh;

end

% start diagonal term of cn
ppgnhh(1,1) = 1-ahh/2;

% end diagonal term of cn
ppgnhh(K,K) = 1-ahh/2;

% the triple diagonal propagation matrix for cn
for i = 1:K-1
    ppgn1hh(i,i+1) = -ahh/2;
    ppgn1hh(i+1,i) = -ahh/2;
    ppgn1hh(i+1,i+1) = 1+ahh;

end

% start diagonal term of cn
ppgn1hh(1,1) = 1+ahh/2;

% end diagonal term of cn
ppgn1hh(K,K) = 1+ahh/2;

ppgcnhh = ppgn1hh\ppgnhh;

%% define operating physical variables
t = 0;
x = dx*(0:K-1)' + dx/2;

tiledlayout(3,1)

%% numerical propagation

for j = 1:1E4
%   add incoming flux into the top layer in concentration table
    conce(1) = conce(1) + FNa*dt/dx;
    
%   carry out chemistry
    ree = kee* conce .* conce;
    rhh = khh* conch .* conch;
    rep = kep* conce .* concp;
    reh = keh* conce .* conch;
    re = ke* conce;

    Re = -2*ree - rep - reh - re;
    Rh = -2*rhh - reh + re + rep;
    Rhh = ree + rhh + reh;
    
    conce = conce + Re*dt;
    conch = conch + Rh*dt;
    conchh = conchh + Rhh*dt;
    
%   apply fick's second law to the concentration table
    conce = ppgcne*conce;
    conch = ppgcnh*conch;
    conchh = ppgcnhh*conchh;
    
%   helium evaporation


%   timer for later analytical use
    t = t + dt;

%   
end

conctable1molh = conch * 1/(6E-4);
conctable1molhh = conchh * 1/(6E-4);
conctable1mole = conce * 1/(6E-4);



%     clf;
nexttile;
line(x, conctable1mole, 'Color', 'red');
nexttile;
line(x, conctable1molh, 'Color', 'blue');
nexttile;
line(x, conctable1molhh, 'Color', 'black');

