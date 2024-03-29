%%%%%%
%%%%% File name: Gfun.m 
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% Purpose: This is the file where researchers need to formulate the
%%%%% 'G', the jacobian. This also helps to find out the frontier of 
%%%%% the geophysics inverse problem.
%%%%%%%%%%%%%%%%%%%%%%%
function G= Gfun(x, H, z, w)
% x are in n data point & it should be in a column
 %%%%%%%
 m=length(w);
 n= length(x);
 w= repmat(w,n,1);
 x=repmat(x,1,m);
 z=repmat(z,n,1);
 G= (2 *(H-z))./((x-w).^2 + (H-z).^2);
end

