clear

syms d1 d2 d1_ d2_ %proportion of larvae allocated to each distance
d0 = 1-d1-d2;
d0_ = 1-d1_-d2_;

n2 = 8; % total number of spots reachable at distance 2
g2 = 4; % habitable number of spots reachable
n1 = 4;
g1 = 3;
n0 = 1;
g0 = 1;

% chance of a larva from a (d1_,d2_) parent winning the competition at
% each distance, in a population of (d1,d2) parents
w2 = g2*((d2_/n2)/((d2_/n2)+(d0+(g1/n1)*d1+(g2/n2)*d2)-(d2/n2)));
w1 = g1*((d1_/n1)/((d1_/n1)+(d0+(g1/n1)*d1+(g2/n2)*d2)-(d1/n1)));
w0 = d0_/(d0_+(d0+(g1/n1)*d1+(g2/n2)*d2)-d0);

% total fitness (expected number of larvae successfully dispersing and
% winning their competition) for a parent playing (d1_,d2_) in a population
% of parents playing (d1,d2)
w = w2+w1+w0;

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

w_surf_gen = subs(w,{d1 d2},{.7 0}); %arbitrary choice of existing strategy
w_surf = subs(w_surf_gen,{d1_ d2_},{xy(:,1) xy(:,2)});

plot3(xy(:,1),xy(:,2),w_surf,'*')
xlabel('d1')
ylabel('d2')
zlabel('w')

% computationally find optimum value
[M,I]=max(w_surf);
opt = xy(I,:); % d1* d2*

%% Look for ESS

% find both partial derivatives wrt d1_ and d2_
dwdd1_ = diff(w,'d1_');
dwdd2_ = diff(w,'d2_');

% set d1_=d1=d1star, d2_=d2=d2star
dwdd1_new = subs(dwdd1_,[d1_ d2_],[d1 d2]);
dwdd2_new = subs(dwdd2_,[d1_ d2_],[d1 d2]);

%solve for d1star and d2star
[d1star, d2star] = solve(dwdd1_new==0,dwdd2_new==0);
%only one option with sum=<1
d1star = double(d1star(1));
d2star = double(d2star(1));

w_opt = double(subs(w, [d1 d1_ d2 d2_], [d1star d1star d2star d2star]));
w_test = double(subs(w, [d1 d1_ d2 d2_], [d1star d1star+0.1 d2star d2star+0.1]));
%w_opt<w_test, so [d1star d2star] is a minimum, and not an ESS.