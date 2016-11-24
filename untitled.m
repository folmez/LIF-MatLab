clear all;
tic;
x=rand(1e3,1);
y=rand(1e3,1);
z=rand(1e3,1);
n=1e5;
evk = @(a,v,b) (-1).*a.*v + b;
for i=1:n
    evaluate_k(x,y,z);
%     evk(x,y,z);
end
toc;