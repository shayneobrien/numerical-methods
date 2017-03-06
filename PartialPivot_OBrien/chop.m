%% Algorithm 6.2 and 6.3: Chop function

% Name: Shayne O'Brien
% Course: MATH 345 (Dr. Haddad)
% Due Date: Wednesday, 11/2/16 by 11:30 am
% Content: Chop function

%% Taken from: http://matlab-coding.blogspot.com/2008/04/chop.html

function Xout=chop(Xin,n)

% Usage: y = chop(x,n)

% rounds elements of x to n decimal digits. 

% For example, chop(24.567,4) returns 24.57

% chop(24.564,4) returns 24.56

% chop(24.564,2) returns 25.00

% chop(24.564,1) returns 20.00

%

% Users of more recent versions of Matlab 4.1c+ should

% use the system provided command "chop" instead.



% Abort if only one input argumentis is given.

if nargin<2, help chop, return, end



% abs(Xin)=A*10^E, A=0.ddddd... & E=floor(log10(|Xin|))+1

% round the decimals in abs(Xin)/10^K & K=E-n.

zero = (Xin==0);

Xout=abs(Xin)+zero; size_x = size(Xout);

E_X= (10*ones(size_x)).^(floor(log10(Xout))-n+1);

Xout=round(Xout./E_X).*E_X; Xout=sign(Xin).*Xout.*(1-zero);