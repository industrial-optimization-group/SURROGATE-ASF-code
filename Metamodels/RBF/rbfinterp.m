function [f] = rbfinterp(x, options)
tic;
phi       = options.('rbfphi');
rbfconst  = options.('RBFConstant');
nodes     = options.('x');%by author
rbfcoeff  = (options.('rbfcoeff'))';


[dim              n] = size(nodes);
[dimPoints  nPoints] = size(x);%By author


if (dim~=dimPoints)
   error(sprintf('x should have the same number of rows as an array used to create RBF interpolation'));% by author
    
end;

f = zeros(1, nPoints);
r = zeros(1, n);

for i=1:1:nPoints
	s=0;
    r =  (x(:,i)*ones(1,n)) - nodes;%by author  
    r = sqrt(sum(r.*r, 1));
%     for j=1:n
%          r(j) =  norm(x(:,i) - nodes(:,j));
%     end
%     s =  sum(rbfcoeff(1:n).*feval(phi, r, rbfconst));%by me
%     temp1 = feval(phi, r, rbfconst);
    
    s_radial = sum(rbfcoeff(1:n).*feval(phi, r, rbfconst));
%     s_aug = sum(rbfcoeff(n+1:end).*constant(x(:,i)));
    s_aug = sum(rbfcoeff(n+1:end).*linear(x(:,i)));
%    s_aug = sum(rbfcoeff(n+1:end).*quad(x(:,i)));% see line 46
%     s_aug = sum(rbfcoeff(n+1:end).*cubic(x(:,i)));
%     s_aug = sum(rbfcoeff(n+1:end).*quatric(x(:,i)));
%     s_aug = sum(rbfcoeff(n+1:end).*degreefive(x(:,i)));
%     s = s_radial;
    s = s_radial + s_aug;
%     s = rbfcoeff(n+1) + sum(rbfcoeff(1:n).*feval(phi, r, rbfconst));% by authors
%  
% 	for k=1:dim
%        s=s+rbfcoeff(k+n+1)*x(k,i);     % linear part
% 	end
	f(i) = s;
end;

if (strcmp(options.('Stats'),'on'))
    fprintf('Interpolation at %d points was computed in %e sec\n', length(f), toc);    
end;

function aug = constant(x)
aug = ones(1,1);

function aug = linear(x)%num of additional term:3
aug = [ones(1,1) x'];

function aug = quad(x)

aug = [1 x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2];

function aug = cubic(x)

aug = [1 x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3];
   
function aug = quatric(x)
aug = [1 x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3 ...
       x(1,:)' .* x(2,:)'.^3 x(1,:)'.^2 .* x(2,:)'.^2 x(1,:)'.^3 .* x(2,:)' ...
       x(1,:)'.^4 x(1,:)'.^4];
   
function aug = degreefive(x)
aug = [1 x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3 ...
       x(1,:)' .* x(2,:)'.^3 x(1,:)'.^2 .* x(2,:)'.^2 x(1,:)'.^3 .* x(2,:)' x(1,:)'.^4 ...
       x(1,:)'.^4 ...
       x(1,:)' .* x(2,:)'.^4 x(1,:)'^2 .* x(2,:)'.^3 x(1,:)'^3 .* x(2,:)'.^2 ...
       x(1,:)'^4 .* x(2,:)' x(1,:)'^5 x(2,:)'.^5];