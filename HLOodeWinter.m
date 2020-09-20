function dxdt=HLOodeWinter(t,x)

global b a d f m u g
dxdt=[-x(1).*x(2)./(b+x(1))-a.*x(1).*x(3)./(d+x(1)); f.*x(1).*x(2)./(b+x(1))-m.*x(2); g.*x(1).*x(3)./(d+x(1))-u.*x(3)];
end