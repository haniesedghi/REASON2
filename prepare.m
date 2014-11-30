function prepare

clc;
clear all;

sigma=0.5;
p=2000;%200;
m=p;
n=p;
size=p;
r=100;%10; 
lmb=0.05;
%SampleSize=10000;
%Data1=zeros(p^2,SampleSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r = 100 % rank
incoh = 1.6; % incoherence (atleast 1)
%tune_lambda = 0;
randn_mat = zeros(m, n); % random iid gaussian submatrix for incoherence
z_u = floor(m/incoh^2);
z_v = floor(n/incoh^2);
randn_mat(1:z_u, 1:z_v) = randn(z_u, z_v);
[U, Sig, V] = svds(randn_mat, r);

for i = 1:r-1
    Sig(i, i) = 1;%mean*(1 + vari*mean*randn);
end
Sig(r, r) = 1;
L_star = U*Sig*V';
%rank(L_star)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_star=zeros(p,p);

 s=500; %50 sparsity level
for m=1:s
    if rand <0.5
        S_star(randi(p),randi(p))=1;
    else
       S_star(randi(p),randi(p))=-1; 
    end
end



S_Star=S_star*r/p;



X=S_star+L_star;

% 
% for i=1:SampleSize
%     D=S_star+L_star+normrnd(0,sigma*eye(p),[p p]);
%     Data1(:,i)=D(:);
% end
    
s=nnz(S_star)
savefile='linearnoise.mat';
%save(savefile,'Data1','L_star','S_star');
save(savefile,'X','L_star','S_star');
