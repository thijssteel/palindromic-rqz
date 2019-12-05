function a = randlog()
%RANDLOG generate a logarithmically distributed rand number
t = rand(1);
s = 2*randi([0,1],1)-1;
a = s*10.^(-15 + 15*t);
end

