function F = PID(SP,t,x_Error,F)
       

%     disp(error);
    T = 1e-3;  
%    SP = 46.2780;
    Kc = 0.9;
    tal_i = (50e-3)/2;
    tal_d = (2e-3)/3;
    error = SP-x_Error ;
    delta_c = Kc*(1+(T/tal_i)+(tal_d/T))*error(end) - Kc*error(end-1)*(1+(2*tal_d/T))+ Kc*error(end-2)*(tal_d/T);
%     disp(delta_c);
    

%     
    F = F + delta_c;
    
end
