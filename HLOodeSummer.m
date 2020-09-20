function dxdt=HLOodeSummer(t,x)

global b B f m h v u s
dxdt=[x(1).*(1-x(1))-x(1).*x(2)./(b+x(1))-x(1).^2.*x(3)./(B+x(1).^2); f.*x(1).*x(2)./(b+x(1))-m.*x(2); h.*x(1).^2.*x(3)./(B+x(1).^2)+s.*x(3)./(v+x(3))-u.*x(3)];
end