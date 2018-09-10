function [options xBasis yBasis selectedBasisInd] = rbfcreate(x, y, varargin)
%RBFCREATE Creates an RBF interpolation
%   OPTIONS = RBFSET(X, Y, 'NAME1',VALUE1,'NAME2',VALUE2,...) creates an   
%   radial base function interpolation 
%   
%   RBFCREATE with no input arguments displays all property names and their
%   possible values.
%   
%RBFCREATE PROPERTIES
% 

%
% Alex Chirokov, alex.chirokov@gmail.com
% 16 Feb 2006
tic;
% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('               x: [ dim by n matrix of coordinates for the nodes ]\n');
  fprintf('               y: [   1 by n vector of values at nodes ]\n');
  fprintf('     RBFFunction: [ gaussian  | thinplate | cubic | multiquadrics | {linear} ]\n');
  fprintf('     RBFConstant: [ positive scalar     ]\n');
  fprintf('       RBFSmooth: [ positive scalar {0} ]\n');
  fprintf('           Stats: [ on | {off} ]\n');
  fprintf('\n');
  return;
end
Names = [
    'RBFFunction      '
    'RBFConstant      '
    'RBFSmooth        '
    'Stats            '
    'xBasis           '%by Mohammad
    'yBasis           '%By Mohammad
    'selectedBasisInd '%By Mohammad
];
[m,n] = size(Names);
names = lower(Names);

options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end

%**************************************************************************
%Check input arrays 
%**************************************************************************
[nXDim nXCount]=size(x);
[nYDim nYCount]=size(y);

if (nXCount~=nYCount)
  error(sprintf('x and y should have the same number of rows'));
end;

if (nYDim~=1)
  error(sprintf('y should be n by 1 vector'));
end;

options.('x')           = x;
options.('y')           = y;
%**************************************************************************
%Default values 
%**************************************************************************
options.('RBFFunction') = 'linear';
%**********The radius of the basises***************************************
d_max = max(pdist(x','euclidean'));
options.('RBFConstant') = d_max / ((nXDim * nYCount)^(1/nXDim));% Kitayama
% options.('RBFConstant') = (prod(max(x')-min(x'))/nXCount)^(1/nXDim); %approx. average distance between the nodes 
%**************************************************************************
options.('RBFSmooth')   = 0;
options.('Stats')       = 'off';

%**************************************************************************
% Argument parsing code: similar to ODESET.m
%**************************************************************************

i = 1;
% A finite state machine to parse name-value pairs.
if rem(nargin-2,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin-2
  arg = varargin{i};
    
  if ~expectval
    if ~isstr(arg)
      error(sprintf('Expected argument %d to be a string property name.', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;      
  end
  i = i + 1;
end

if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end

    
%**************************************************************************
% Creating RBF Interpolatin
%**************************************************************************

switch lower(options.('RBFFunction'))
      case 'linear'          
        options.('rbfphi')   = @rbfphi_linear;
      case 'cubic'
        options.('rbfphi')   = @rbfphi_cubic;
      case 'multiquadric'
        options.('rbfphi')   = @rbfphi_multiquadrics;
      case 'thinplate'
        options.('rbfphi')   = @rbfphi_thinplate;
      case 'gaussian'
        options.('rbfphi')   = @rbfphi_gaussian;
    otherwise
        options.('rbfphi')   = @rbfphi_linear;
end

phi       = options.('rbfphi');

A=rbfAssemble(x, phi, options.('RBFConstant'), options.('RBFSmooth'));

%**************************************************************************
%Select basis for RBF by Mohammad
%**************************************************************************
% [xBasis, yBasis, selectedBasisInd, rbfcoeff] = selectBasis(x,y, A, options);
% options.('xBasis') = xBasis;
% options.('yBasis') = yBasis;
% options.('selectedBasisInd') = selectedBasisInd;

%**************************************************************************

% b = y'; %by me
% save('A','A');
b=[y'; zeros(nXDim+1,1)]; % by authors zeros(numTerm,1)
% b = y';
rbfcoeff = A\b;
% rbfcoeff = A^(-1) * b; % by me

% cap_landa = 1e-3 * eye(size(A)); % by me inspired by Kitayama
% rbfcoeff = (A' * A + cap_landa) ^(-1) * A' * b; % by me inspired by Kitayama
%inverse
%%%%%By me%%%%%%%%%%%%%
% landa = (1e-3) * eye(nYDim);
% rbfcoeff = (A' * A +landa)^(-1) * (A'*b);
% rbfcoeff = (A' * A)^(-1) * (A'*b);
%%%%%By the authors%%%%
% rbfcoeff=A\b;
%%%%%%%%%%%%%%%%%%%%%%%
% coeff_old = rbfcoeff;% by me
% save('coeff_old','coeff_old');%by me
%SVD
% [U,S,V] = svd(A);
% 
% for i=1:1:nXCount+1
%     if (S(i,i)>0) S(i,i)=1/S(i,i); end;   
% end;    
% rbfcoeff = V*S'*U*b;


options.('rbfcoeff') = rbfcoeff;


if (strcmp(options.('Stats'),'on'))
    fprintf('%d point RBF interpolation was created in %e sec\n', length(y), toc);  
    fprintf('\n');
end;

end

function [A]=rbfAssemble(x, phi, const, smooth)
% x = x';
[dim, n]=size(x);
%%%%%%%%%By author%%%%%%%%%
% A=zeros(n,n);
% for i=1:n   
%     A(i,i) = A(i,i) - smooth;
%     for j=1:i
%         r=norm(x(:,i)-x(:,j));        
%         temp=feval(phi,r, const);
%         A(i,j)=temp;
%         A(j,i)=temp;
%     end      
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%From
%%%%%%%%%%%%%%%http://www.princeton.edu/~kung/ele571/571-MatLab/571svm/svmkernel.m%%%%%%%
x = x';
dist2 = repmat(sum((x.^2)', 1), [n 1])' + ...
repmat(sum((x.^2)',1), [n 1]) - ...
            2*x*(x');
x = x';
dist = sqrt(abs(dist2));
A = feval(phi,dist, const);
% dif = sum(sum(A - ANew))
%%%%%%By me%%%%%%%%%%%%
% phi_old = A;% by me
% save('phi_old','phi_old');
% P = constant(x,n);%add term: 1
P = linear(x,n);%add term: 3
% P = quadratic(x,n);%add term: 6
% P = cubic(x,n);%add term: 10
% P = quatric(x,n);%add term: 15
% P = degreefive(x,n); %add term: 21
A = [ A      P % by author
       P' zeros(dim+1,dim+1)];%by author
% A = [P A];
%%%%%%By authors%%%%%%%
% % Polynomial part
% P=[ones(n,1) x'];It means that the augmented part is p(x) = 1 + x(1) + x(2)
% A = [ A      P
%       P' zeros(dim+1,dim+1)];
end%rbfAssemble
%**************************************************************************
% Radial Base Functions
%************************************************************************** 
function u=rbfphi_linear(r, const)
u=r;
end
function u=rbfphi_cubic(r, const)
u=r.*r.*r;
end
function u=rbfphi_gaussian(r, const)
u=exp(-0.5*r.*r/(const*const));
end
function u=rbfphi_multiquadrics(r, const)
u=sqrt(1+r.*r/(const*const));
end
function u=rbfphi_thinplate(r, const) 
u=r.*r.*log(r); % by authors
u(isnan(u)) = 0;
% u = zeros(size(r));
% nonvanish = find (r > 1e-200)
% for i = 1 : length(nonvanish)% This for is time-consuming
%     u(nonvanish(i)) = r(nonvanish(i))*r(nonvanish(i))*log(r(nonvanish(i)));
% end% i

% if r < 1e-200
%     u = 0;
% else
%     u=r.*r.*log(r); % by me
end
%*****************************************************************************
% Augmenting functions
%*****************************************************************************
function aug = constant(x,n)%num of additional term:1
aug = ones(n,1);
end
function aug = linear(x,n)%num of additional term:3
aug = [ones(n,1) x'];
end
function aug = quadratic(x,n)%num of additional term:6
aug = [ones(n,1) x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2];
end
function aug = cubic(x,n)%num of additional term:10
aug = [ones(n,1) x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3];
 end  
function aug = quatric(x,n)%num of additional term:15
aug = [ones(n,1) x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3 ...
       x(1,:)' .* x(2,:)'.^3 x(1,:)'.^2 .* x(2,:)'.^2 x(1,:)'.^3 .* x(2,:)' x(1,:)'.^4 ...
       x(1,:)'.^4];
end   
function aug = degreefive(x,n)%num of additional term:21
aug = [ones(n,1) x(1,:)' x(2,:)' x(1,:)' .* x(2,:)' x(1,:)'.^2 x(2,:)'.^2 ...
       x(1,:)'.^2 .* x(2,:)' x(1,:)' .* x(2,:)'.^2 x(1,:)'.^3 x(2,:)'.^3 ...
       x(1,:)' .* x(2,:)'.^3 x(1,:)'.^2 .* x(2,:)'.^2 x(1,:)'.^3 .* x(2,:)' x(1,:)'.^4 ...
       x(1,:)'.^4 ...
       x(1,:)' .* x(2,:)'.^4 x(1,:)'.^2 .* x(2,:)'.^3 x(1,:)'.^3 .* x(2,:)'.^2 ...
       x(1,:)'.^4 .* x(2,:)' x(1,:)'.^5 x(2,:)'.^5];
end