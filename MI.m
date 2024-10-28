function mi = MI (y,Par)

    Par(6) = max(y(3), Par(6));
    mi = (Par(3)*(y(2)/(Par(4)+y(2)+((y(2)^2)/Par(5)))))*((1-(y(3)/Par(6)))^Par(7));

end