%% set up fitness equation

clear
nmax = 2;
v = sym('v',[1 nmax]); %population displacement strategy
v_ = sym('v_',[1 nmax]); %focal displacement strategy
assume(sum(v)==1 & sum(v_)==1); %probabilities sum to 1

n = [1 4]; %total number of sites at each distance
g = [1 3]; %number of habitable sites at each distance

w = sym('w',[1 nmax]);

for k = 1:nmax
    w(k) = g(k)*(v_(k)/n(k))/((v_(k)-v(k))/n(k)+sum(g.*v./n));
end

%% choose population displacement strategy and plot fitness surface
% xy = [];
% for i = 0:0.01:1
%     for ii = 0:0.01:(1-i)
%         xy = [xy; i ii];
%     end
% end

x = 0:0.01:1;
y = 1-x;

w_surf_gen = subs(w,v,[.5 .5]); %arbitrary choice of existing strategy
w_surf = subs(w_surf_gen,{v_(1) v_(2)},{x y}); %plug in x and y values -- not working

plot(x,w_surf,'*')
