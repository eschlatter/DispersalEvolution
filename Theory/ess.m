clear

syms v1 v2 v1_ v2_ %proportion of larvae allocated to each distance
v0 = 1-v1-v2;
v0_ = 1-v1_-v2_;
assume(v0>0 & v1>0 & v2>0 & v0<1 & v1<1 & v2<1 & ...
    v0_>0 & v1_>0 & v2_>0 & v0_<1 & v1_<1 & v2_<1)


n2 = 8; % total number of spots reachable at distance 2
n1 = 4;
n0 = 1;

g2 = 4; % habitable number of spots reachable at distance 2
g1 = 3;
g0 = 1;

% chance of a larva from a (v1_,v2_) parent winning the competition at
% each distance, in a population of (v1,v2) parents

w2 = g2*((1/n2)*v2_)/(v0+(g1/n1)*v1+((g2-1)/n2)*v2+(1/n2)*v2_);
w1 = g1*((1/n1)*v1_)/(v0+((g1-1)/n1)*v1+(g2/n2)*v2+(1/n1)*v1_);
w0 = v0_/((g1/n1)*v1+(g2/n2)*v2+v0_);

w = w0+w1+w2;

%% Plot fitness surface

xy = [];

for i = 0:0.01:1
    for ii = 0:0.01:(1-i)
        xy = [xy; i ii];
    end
end

% xyn = find(xy(:,1)~=0);
% xy = xy(xyn,:);
% xyn = find(xy(:,2)~=0);
% xy = xy(xyn,:);

w_surf_gen = subs(w,{v1 v2},{.5 .5}); %arbitrary choice of existing strategy
w_surf = subs(w_surf_gen,{v1_ v2_},{xy(:,1) xy(:,2)});

plot3(xy(:,1),xy(:,2),w_surf,'*')
xlabel('v1')
ylabel('v2')
zlabel('w')

% computationally find optimum value
[M,I]=max(w_surf);
opt = xy(I,:); % v1* v2*

%% %% Look for ESS

% find both partial derivatives wrt v1_ and v2_
dwdv1_ = diff(w,'v1_');
dwdv2_ = diff(w,'v2_');

% set d1_=d1=d1star, d2_=d2=d2star
dwdv1_new = subs(dwdv1_,[v1_ v2_],[v1 v2]);
dwdv2_new = subs(dwdv2_,[v1_ v2_],[v1 v2]);

%solve for v1star and v2star
[v1star, v2star] = solve(dwdv1_new==0,dwdv2_new==0);

%% 

v = 0:0.1:5;
y = (9/8)*v.^3-10*v.^2+26*v-16;
plot(v,y)

%% nbins=2, nmax=0 fitness surface


xy = [];
for i = 0:0.01:1
    for ii = 0:0.01:1
        xy = [xy; i ii];
    end
end

v = xy(:,1);
vprime = xy(:,2);

w = (1-vprime)./(1-vprime+.75*v)+.75*vprime./(vprime+4-2*v);

plot3(v,vprime,w,'*')
xlabel('v')
ylabel('vprime')
zlabel('w')
