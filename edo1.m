function dx = edo1 (t,y,Par,F,ke,kw)

    Par(6) = max(y(3), Par(6));
    mi = (Par(3)*(y(2)/(Par(4)+y(2)+((y(2)^2)/Par(5)))))*((1-(y(3)/Par(6)))^Par(7));
    
    dx = zeros(4,1);
    
    aux = (F - (((ke*y(3))+(kw*(994-y(3))))*(y(4)/994)));
    dx(1) = (mi-(aux/y(4)))*y(1);    
    dx(2) = (371.4*(F/y(4)))- ((mi*y(1)/Par(1))+(aux*y(2)/y(4)));
    dx(3) = (mi*(Par(2)/Par(1))*y(1))- (((ke)+(aux/y(4)))*y(3));
    dx(4) = aux;
    
end
