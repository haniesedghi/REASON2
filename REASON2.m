function REASON2 
clc
clear all

addpath PROPACK;


p=2000; %matrix dimension is p by p
m=p;
n=p;
size=p;
R=20; %initial radius
temp = open('linearnoise.mat'); %input, noisy version of M*=S*+L*

X=temp.X;

lambda=1/sqrt(m);
tol=1e-8;
maxIter=10;

if nargin < 2 %nargin Number of function input arguments
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 4
    maxIter = 100;
elseif maxIter == -1
    maxIter = 100;
end



M = X;
EpochCenter=temp.S_star;

S = zeros( m, n);
L = zeros( m, n);
W = zeros( m, n);

Z = zeros( m, n);

G = zeros( m, n);

Sold=zeros(p,p);


TIME=50; %time in seconds for algorith
total_svd = 0;
alg_iter=12; %number of iterations of the algorithm. can be changed by the user



multiplier=1.5;
dnorm = norm(X, 'fro');
tolProj = 1e-6 * dnorm; 
total_svd = 0;

% initialize
Z = X;
norm_two = lansvd(Z, 1, 'L');
norm_inf = norm( Z(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Z = Z / dual_norm;
rho = .5/norm_two ;
rho_bar = rho * 1e7;

tau=rho; %dont change

tauK=1; %dont change


KK=100;

iter = 0;
converged = false;
stopCriterion = 1;
sv = 10;
svp = sv;
mu=1; 
sigma=1;

tic
  tstart=tic; 
%while ~converged 
for k=1:alg_iter
        telapsed=toc(tstart);
        if telapsed > TIME
      stopCriterion = norm(X-S-L, 'fro') / norm(X, 'fro')
    
      test_err = norm(L-temp.L_star, 'fro') / norm(temp.L_star, 'fro')
        test_err2 = norm(S-temp.S_star, 'fro') / norm(temp.S_star, 'fro')
            break;
        end   
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;


   int_iter=0;
   while primal_converged == false
        int_iter=int_iter+1;

     G=M-S-L+(1/rho)*Z;   %G update for decomposing L,S update
        softi=S+tauK*G;
    thresh=(tauK*lambda)/rho;
    W=max( softi - thresh,0)+ min( softi + thresh,0); %L1 update
    if norm(W-EpochCenter,1) > R
    v(:,1)=W(:);%v=vector(W);
    e_cent(:,1)=EpochCenter(:);
                
                
    L1norm=norm(v,1);
                    

                % Find kappa=rho
                
                % if already in the ball, do not proceed with the projection, it is done
                if L1norm <= R
                    theta=v+e_cent;
                else %else do projection
                    
                    v_sort=sort(v,'descend');
                    i=1;
                    sum=0;
                    f=0;
                    kappa=0;
                    zeta=0;
                    while (f==0) && (i<= size)
                        sum=sum+v_sort(i,1);
                        if v_sort(i,1)-(sum-R)/(i+1) <= 0
                            f=1;
                            kappa=max(i,0);
                            sum=sum-v_sort(i,1);
                        else
                        i=i+1;
                        end
                        
                    end
                    
                    if (f==0)
                        kappa=size;
                    end
                    
                    if (kappa==0)
                        zeta=0;
                    else
                     zeta=(sum-R)/kappa;
                    end
                    for i=1:size
                        if (v(i,1)-zeta >0)
                            v(i,1)=sign(v(i,1)-zeta)*v(i,1);
                        else v(i,1)=0;
                        end
                    end
                end
                    
                    
                    % Recentralize here done!
                
                    theta=v+e_cent;
                    %map back to matrix
                
    
        if int_iter==4
            primal_converge=1;
        end
    
                S=reshape(theta(:,1),[p p]);
    else
                    S=W;
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    Gb=M-S+1/rho*Z;%L+tauK*G;
        if choosvd(n, sv) == 1
            [Ub Sb V] = lansvd(Gb, sv, 'L');
       else
            [Ub Sb V] = svd(Gb, 'econ');
        end
        diagS = diag(Sb);
        svp = length(find(diagS > mu*tauK/rho));
        if svp < sv
            sv = min(svp + 1, n);
        else
            sv = min(svp + round(0.05*n), n);
        end
        
        
            
    % thresh=(mu)/rho;
    % Sb = max( Sb - thresh,0) + min( Sb + thresh,0);

    % diagS = diag(Sb);
    
    
    
        temp_L = Ub(:,1:svp)*diag(diagS(1:svp)-mu*tauK/rho)*V(:,1:svp)';   %honey double check soft threshold e ha! diag(diagS(1:svp)-mu/rho) 
        
        if norm(L - temp_L, 'fro') < tolProj  && norm(S - Sold, 'fro') < tolProj
            primal_converged = true;
        end
        L = temp_L;
        Sold=S;
%        W2=(diagS(1:svp)-mu*tauK/rho);
%        if norm(W2,1) > R %if norm(W2-EpochCenter2,1) > R
%               v2=W2;%v=vector(W);
%               %e_cent2(:,1)=EpochCenter2(:);
%               %c1=zeros(1,2);
%               %c1=size(v2);
%                e_cent2=zeros(svp,1);
%            
%                 
%                 L1norm2=norm(v2,1);
%                     
%                 size2=svp;
%                 % Find kappa=rho
%                 
%                 % if already in the ball, do not proceed with the projection, it is done
%                 if L1norm2 <= R
%                     theta2=v2+e_cent2;
%                 else
%                     
%                     v_sort2=sort(v2,'descend');
%                     i=1;
%                     sum=0;
%                     f=0;
%                     kappa=0;
%                     zeta=0;
%                     while (f==0) && (i<= size2)
%                         sum=sum+v_sort2(i,1);
%                         if v_sort2(i,1)-(sum-R)/(i+1) <= 0
%                             f=1;
%                             kappa=max(i,0);
%                             sum=sum-v_sort2(i,1);
%                         else
%                         i=i+1;
%                         end
%                         
%                     end
%                     
%                     if (f==0)
%                         kappa=size2;
%                     end
%                     
%                     if (kappa==0)
%                         zeta=0;
%                     else
%                      zeta=(sum-R)/kappa;
%                     end
%                     for i=1:size2
%                         if (v2(i,1)-zeta >0)
%                             v2(i,1)=sign(v2(i,1)-zeta)*v2(i,1);
%                         else v2(i,1)=0;
%                         end
%                     end
%                 end
%                     
%                     
%                     % Recentralize here done!
%                 
%                     theta2=v2+e_cent2;
%                     %map back to matrix
%                 
%     
% 
%     
%                 S2=eye(svp)*theta2;
%                 L = Ub(:,1:svp)*diag(S2)*V(:,1:svp)';
%                  
%     end
%            
           
           
           
           
           
           
           
           
           
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
   
     Z=Z+tau*(M-S-L);
     rho = min(multiplier*rho, rho_bar);
     rho=rho*multiplier;
       tau=rho;


    


    
 

    %% stop Criterion    
    stopCriterion = norm(X-S-L, 'fro') / norm(X, 'fro');%dnorm;
    if stopCriterion < tol
        converged = true;
    end    
    
    disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' r(A) ' num2str(svp)...
        ' |E|_0 ' num2str(length(find(abs(M-S-L)>0)))...
        ' stopCriterion ' num2str(stopCriterion)]);
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
     test_err = norm(L-temp.L_star, 'fro') / norm(temp.L_star, 'fro')
        test_err2 = norm(S-temp.S_star, 'fro') / norm(temp.S_star, 'fro')
end
toc
%S_avg=S_sum/(KK-100);

         test_err = norm(L-temp.L_star, 'fro') / norm(temp.L_star, 'fro')  %Frobenius error for low rank part
        test_err2 = norm(S-temp.S_star, 'fro') / norm(temp.S_star, 'fro')   %Frobenius error for sparse part
   %     rank(L)

% savefile='linearnoiseResults.mat';
% save(savefile,'S','L');



