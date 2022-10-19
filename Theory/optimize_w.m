clear

Dprime = sym('Dprime',[1 3]);
D = sym('D',[1 3]);
assume(Dprime(1)+Dprime(2)+Dprime(3)==1);
assume(D(1)+D(2)+D(3)==1);

w(Dprime,D) = Dprime(3)^2/(16*(7*D(3)/8+Dprime(3)/8+D(2)+D(1)))+...
                3*Dprime(2)^2/(16*(D(3)+3*D(2)/4+Dprime(2)/4+D(1)))+...
                Dprime(1)^2/(D(3)+D(2)+Dprime(1));

Aprime = [.25, .25, .5];
A = [.5, .25, .25];

w(.25, .25, .5,.5, .25, .25)

%diff(w,'Dprime')