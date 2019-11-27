%%%%%%
%%%%% File name: main.m 
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% Purpose: This is the file where researchers need to formulate the
%%%%% problem which is in our case: calculating the frontier of the
%%%%% \cite{latexTarantola} paper's geophysics ill-posed nonlinear inverse 
%%%%% problem. This is the main background of the coding. 
%%%%%%%%%%%%%%%%%%%%%%%
%
%Regarding the computing artifact we have considered \cite{latexTarantola} paper. Considering \\
%discrete data with continuous parameter case, authors of \cite{latexTarantola} use the following \\
%back propagation where they have used $z(w)$ with two iterations $k$-th and $(k+1)$-th but without \\
%discretizing the $w$. This back propagation is estimating the frontier
%\hat{z}_{k+1}(w)= z_0 + \int dw' \sum_i \sum_j C_{p_0}(w,w') G_{k}^{i}(w') (S^{-1})^{ij} \\
%[ d_0^j - g^j(\hat{\bfz}_k) + \int dw'' G_k^j(w'') . [\hat{z}_k(w'') - z_0(w'')] \\
%
function z_hat=main(w,x,d,sigma_d,sigma,theta,K,Gfun,ffun)
m=length(w); % Number of grid points
W=max(w); % Maximum of the grid
n= length(x); % Number of data
z= zeros(1,m) ; % Initial values for the frontier
H = 10; % Deepth between the surface and the subsurface
%%
%%%%%%%%%%%%
cov = covfun(sigma,theta,w); % Covariance matrix,C_p(w,w')=\sigma^2 exp[-\frac{1}{2} \frac{(w-w')^2}{\Delta^2}]
for k=1:K ;  %K= Number of iteration and it is changable
    G = Gfun(x, H, z, w); %Jacobian
    for i=1:n ;
        for j=1:n ;
            G1 = cov.*G(j,:); % Matrix multiplication of the covariance  and the jacobian
            I3 = trapz(w,G1,2); % Trapzoidal integration of G1
            G2 = G(i,:)'.* I3; % Matrix multiplication the jacobian and the Trapzoidal integration of G1                  
            I4= trapz(w,G2);  % Trapzoidal integration of G2    
            % S matrix calculation
            S(i,j)=I4; %S matrix 
            if i==j
                S(i,j)=sigma_d+ I4;
            end
        end
    end
    %Inverse of S
    Sinv= inv(S); % Inverse of S matrix
    for j=1:n ;
        G3 = G(j,:).*z ; 
        I5(j) = trapz(w,G3); % Trapzoidal integration of G3
        f = ffun(x(j),w,H,z);
        % little g=integration of 'f'
        g(j)= trapz(w,f); % g= Trapzoidal integration of f w.r.t. w
        % Calculating part2(j)=[ d_0^j - g^j(\hat{\bfz}_k) + \int dw'' G_k^j(w'')\\
        %.[\hat{z}_k(w'')-z_0(w'')]
        part2(j)= d(j)-g(j)+I5(j);
    end
    % Calculating \hat{z}_{k+1}(w)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Let's consider all = \hat{z}_{k+1}(w)
    all = zeros(m,m);
    for i=1:n ;
        for j=1:n ;
            all = all + cov.*G(i,:).*Sinv(i,j)*part2(j) ;
        end
    end
    % Calculating the integration of \hat{z}_{k+1}(w)
    I6= trapz(w,all,2) ;
    z=I6' ;
    z_hat(:,k)= I6 ;
end
%%%%%%%%%%%%%
%%
% Plotting z_hat i.e. the results of the integration of \hat{z}_{k+1}(w)
figure(8); clf; hold on;
for i=1:K
plot(w,z_hat(:,i),'LineWidth', 5) %setting the line width of a plot
set(gca, 'Fontsize', 14) %setting a title of the plot
set(gca,'FontWeight','bold') %%setting a font weight of the plot
title('Line plot of frontier') %setting a title of the plot
xlabel('w') %labeling x-axis
ylabel('z(w)') %labeling x-axis
end
hold off
end % end of the z_hat function
%%%%%%%%%%%%%%%%
%covariance function, cov
function cov = covfun(sigma,theta,w)
[W1, W2]= meshgrid(w);
cov = sigma.^2 *exp((-0.5 * ((W1-W2).^2 )/theta)); %calculating the covariance
end
