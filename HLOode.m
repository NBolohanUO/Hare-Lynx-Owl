function dxdt=HLOode(t,x)

global b B a d f m h v u g Ts s
dxdt=[Ts.*(x(1).*(1-x(1))-x(1).*x(2)./(b+x(1))-x(1).^2.*x(3)./(B+x(1).^2))+(1-Ts).*(-x(1).*x(2)./(b+x(1))-a.*x(1).*x(3)./(d+x(1))); f.*x(1).*x(2)./(b+x(1))-m.*x(2); Ts.*(h.*x(1).^2.*x(3)./(B+x(1).^2)+s.*x(3)./(v+x(3))-u.*x(3))+(1-Ts).*(g.*x(1).*x(3)./(d+x(1))-u.*x(3))];
end