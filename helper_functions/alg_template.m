function Yhat=alg_template(Y,predictionMethod,test_ind,left_out)
%alg_template predicts DTIs based on the prediction method selected in
%start.m or sensitivity_analysis.m
%
% INPUT:
%  Y:                     interaction matrix
%  predictionMethod:      method to use for prediction
%  test_indices:          indices of the test set instances
%  left_out:              in case of S1: left_out = test_indices
%                         in case of S2: left_out = left out drugs
%                         in case of S3: left_out = left out targets
%
% OUTPUT:
%  Yhat:                  prediction scores matrix

    % Parameters
    global Sd St
    predFn = str2func(['alg_'  predictionMethod]);
    Yhat = predFn(Y,Sd,St,test_ind,left_out);

end

function Yhat=alg_rls_wnn(Y,ka,kb,~,~)
%alg_rls_wnn predicts DTIs based on the algorithm described in the following paper: 
% Twan van Laarhoven, Elena Marchiori,
% (2013) Predicting drug–target interactions for new drug compounds using a
%           weighted nearest neighbor profile 
% 
% Code below is adapted from the code available at this website:
% http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/

    %%%%%%%%%%%
    %%% WNN %%%
    %%%%%%%%%%%

    global eta
    %eta = 0.7;     %default
    Y = preprocess_WNN(Y,ka,kb,eta);


    %%%%%%%%%%%
    %%% GIP %%%
    %%%%%%%%%%%

    global alpha
    %alpha = 0.5;   %default
    ka = alpha*ka + (1-alpha)*getGipKernel(Y);
    kb = alpha*kb + (1-alpha)*getGipKernel(Y');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Regularized Least Squares (RLS-kron) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global sigma
	%sigma = 1;     %default
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	m1 = va' * Y * vb;
	m2 = m1 .* l;
	Yhat = va * m2 * vb';
end

function y3=alg_grmf(Y,Sd,St,test_ind,~)
%alg_grmf predicts DTIs based on the algorithm described in the following paper: 
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2016) Drug-target interaction prediction with graph-regularized matrix factorization
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  y3:          prediction matrix
%

    % parameters
    global num_iter p k lambda_l lambda_d lambda_t

    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    Sd = preprocess_PNN(Sd,p);
     St = preprocess_PNN(St,p);

    % Laplacian Matrices
    Dd = diag(sum(Sd));
    Ld = Dd - Sd;
     Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));
    Dt = diag(sum(St));
    Lt = Dt - St;
     Lt = (Dt^(-0.5))*Lt*(Dt^(-0.5));

    % (W)GRMF
    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W);    % update A & B

    % compute prediction matrix
    y3 = A*B';
end

function [A,B]=alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_grmf_predict performs alternating least squares for GRMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Ld:          drug graph Laplacian
%  Lt:          target graph Laplacian
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Ld = lambda_d*Ld;          % to avoid 
    lambda_t_Lt = lambda_t*Lt;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))

        %%%%%%%%%%%%
        %%% GRMF %%%
        %%%%%%%%%%%%

        for z=1:num_iter
            A = (Y*B  - lambda_d_Ld*A) / (B'*B + lambda_l_eye_K);
            B = (Y'*A - lambda_t_Lt*B) / (A'*A + lambda_l_eye_K);
        end
        
        
    else

        %%%%%%%%%%%%%
        %%% WGRMF %%%
        %%%%%%%%%%%%%

        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             for i=1:size(A,1)
%                 A(i,:) = (H(i,:)*B - lambda_d*Ld(i,:)*A_old) / (B'*diag(W(i,:))*B + lambda*eye(k));
%             end
%             B_old = B;
%             for j=1:size(B,1)
%                 B(j,:) = (H(:,j)'*A - lambda_t*Lt(j,:)*B_old) / (A'*diag(W(:,j))*A + lambda*eye(k));
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_minus_alpha_Ld_A_old = H*B - lambda_d_Ld*A_old;
            for a=1:size(A,1)
                A(a,:) = HB_minus_alpha_Ld_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_K);
            end
            
            B_old = B;
            HtA_minus_beta_Lt_B_old = H'*A - lambda_t_Lt*B_old;
            for b=1:size(B,1)
                B(b,:) = HtA_minus_beta_Lt_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_K);
            end
        end
    end
    
end

function [S,p]=preprocess_PNN(S,p)
%preprocess_PNN sparsifies S by keeping, for each drug/target, the "p"
% nearest neighbors (NNs) and discarding the rest. 

    NN_mat = zeros(size(S));
    for j=1:length(NN_mat)
        [~,indx] = sort(S(j,:),'descend');
        indx = indx(1:p+1);     % keep drug/target j and its "p" NNs
        NN_mat(j,indx) = 1;
    end
    NN_mat = (NN_mat+NN_mat')/2;
    S = NN_mat .* S;

end

function [A,B]=initializer(Y,k)
%initializer initializes the A and B latent feature matrices for either
% of the CMF or GRMF algorithms.
%
% INPUT:
%  Y:   interaction matrix
%  k:   number of latent features
%
% OUTPUT:
%  A:   latent feature matrix for drugs
%  B:   latent feature matrix for targets
%

    [u,s,v] = svds(Y,k);
    A = u*(s^0.5);
    B = v*(s^0.5);

%     % Alternative: Use non-negative matrix factorization
%     k = min(k, min(size(Y)));
%     [A,B] = nnmf(Y,k);
%     B = B';

end

function y3=alg_cmf(Y,Sd,St,test_ind,~)
%alg_cmf predicts DTIs based on the algorithm described in the following paper:
% Xiaodong Zheng, Hao Ding, Hiroshi Mamitsuka and Shanfeng Zhu
% (2013) Collaborative Matrix Factorization with Multiple Similarities for Predicting Drug-Target Interactions
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  y3:          prediction matrix
%

    % parameters
    global num_iter k lambda_l lambda_d lambda_t

    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W);     % update A & B

    % compute prediction matrix
    y3 = A*B';

end

function [A,B]=alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_cmf_predict performs alternating least squares for CMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Sd = lambda_d*Sd;          % to avoid 
    lambda_t_St = lambda_t*St;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d_Sd*A)  / (BtB + lambda_l_eye_K + lambda_d*(AtA));
            AtA = A'*A;
            B = (Y'*A + lambda_t_St*B) / (AtA + lambda_l_eye_K + lambda_t*(BtB));
            BtB = B'*B;
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             lambda_d_A_oldt_A_old = lambda_d*(A_old'*A_old);
%             for a=1:size(A,1)
%                 A(a,:) = (H(a,:)*B + lambda_d_Sd(a,:)*A_old) / (B'*B + lambda_l_eye_k + lambda_d_A_oldt_A_old);
%             end
%             B_old = B;
%             lambda_t_B_oldt_B_old = lambda_t*(B_old'*B_old);
%             for b=1:size(B,1)
%                 B(b,:) = (H(:,b)'*A + lambda_t_St(b,:)*B_old) / (A'*A + lambda_l_eye_k + lambda_t_B_oldt_B_old);
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d_Sd*A_old;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = lambda_l_eye_K + lambda_d*(A_old'*A_old);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_k_plus_lambda_d_A_oldt_A_old);
            end
            B_old = B;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t_St*B_old;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = lambda_l_eye_K + lambda_t*(B_old'*B_old);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_k_plus_lambda_t_B_oldt_B_old);
            end

        end
    end
    
end

function Yhat=alg_mc(Y,Sd,St,test_ind,~)

global rank

  W = ones(size(Y));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
    
IDX = find(W);
M = opRestriction(numel(W),IDX);
y = M(Y(:),1);
r=rank;
[Yhat] = IST_MC(y,M,size(W),r); %mask changed and lansvd changed, with NN constraint

end

function Yhat=alg_mgrnnm(matrix,Sd,St,test_ind,~)

global lamda nu1 nu2 mu1 mu2 pp method;

Z=matrix';
Y=matrix;

    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    
    Sd = preprocess_PNN(Sd,pp);
    St = preprocess_PNN(St,pp);

    
% Laplacian Matrices    
Dd = diag(sum(Sd)); Lr = Dd - Sd;  
if(det(Dd)==0)
    Dd=0.1*eye(size(Dd))+Dd;
end
Lr = (Dd^(-0.5))*Lr*(Dd^(-0.5));

Dt = diag(sum(St)); Lc = Dt - St;
if(det(Dt)==0)
    Dt=0.1*eye(size(Dt))+Dt;
end
Lc = (Dt^(-0.5))*Lc*(Dt^(-0.5));
    %}
  W = ones(size(matrix));          % weight matrix W
  W(test_ind) = 0;            % set W=0 for test instances
  IDX = find(W);
   
   A = opRestriction(numel(W),IDX); %ii1=sqrt(nu1)*ones(size(matrix));
   A2=opDiag(numel(W),sqrt(nu1));%opMask(ii1(:));%opDirac(numel(W));%opDiag(size(matrix,1), sqrt(nu1));%ii2=sqrt(nu2)*ones(size(matrix));
   A3=opDiag(numel(W),sqrt(nu2));%opMask(ii2(:));%A3= opDiag(numel(matrix), sqrt(nu2));
   AA=opStack(A,A2,  A3  );

   AA_1BMC=[W;W;W];%[W; sqrt(nu1)*eye([size(W,2) size(W,2)]); sqrt(nu2)*eye([size(W,2) size(W,2)]) ];
   
   ys=[]; ys2=[]; ys3=[];
  
for i=1:20
    %YY=AA(matrix(:),1);  
    
    r2=sqrt(nu1).*Z'; r3=sqrt(nu2).*Y;
    YY=[A(matrix(:),1); r2(:); r3(:) ];
    YY_1BMC=[W.*matrix; W.*matrix;  W.*matrix];%[W.*matrix; sqrt(nu1)*Z';  sqrt(nu2)*Y];
    
     if strcmp(method,'NNM')
         [X] = SVS( YY , AA,[ 1*size(matrix,1) size(matrix,2)],0,lamda); 
    
     elseif strcmp(method,'1BMC')
        %temp=ones(size(matrix));t2=AA(temp(:),1);AAnew=AA(t2,2);
        
        [X] = OBMC( matrix ,W, [ 1*size(matrix,1) size(matrix,2)]); 
     end
   Z=sylvester(nu1*eye(size(Z,1)),mu1*Lr,nu1*X');
   Y=sylvester(nu2*eye(size(Y,1)),mu2*Lc,nu2*X);
   
   ys=[ys; norm(YY-AA(X(:),1),'fro')];
   ys2=[ys2; lamda*norm(X(:),1)];
   ys3=[ys3; mu1*trace(X'*Lr*X)+mu2*trace(X*Lc*X')];
end

  %figure;  plot(ys+ys2+ys3);
  
Yhat=X; 
   
end


