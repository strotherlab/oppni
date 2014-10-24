function DL = get_legendre( ord, N )
%
x  = linspace(-1,1,N)';
%
DL = [];
%
if (ord>=0)
    DL = [];
    for i = 0:ord
        X = generate_legendre(i,N);
        DL = [DL X];
    end
end


% if( ord>=0 ) d0 = ones(N,1);                                       DL=[DL d0]; end
% if( ord>=1 ) d1 = x;                                               DL=[DL d1]; end
% if( ord>=2 ) d2 = 0.5    *(3*x.^2   - 1);                          DL=[DL d2]; end
% if( ord>=3 ) d3 = 0.5    *(5*x.^3   - 3*x);                        DL=[DL d3]; end
% if( ord>=4 ) d4 = 0.125  *(35*x.^4  - 30*x.^2  + 3);               DL=[DL d4]; end
% if( ord>=5 ) d5 = 0.125  *(63*x.^5  - 70*x.^3  + 15*x);            DL=[DL d5]; end
% if( ord>=6 ) d6 = 0.0625 *(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);    DL=[DL d6]; end
% if( ord>=7 ) d7 = 0.0625 *(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x); DL=[DL d7]; end
%

function X = generate_legendre(ord,N)

x  = linspace(-1,1,N)';
pk = LegendrePoly(ord);
xp = bsxfun(@power,x,[ord:-1:0]);    
X = bsxfun(@times,xp,pk');    
X = sum(X,2);


