i=1;
j=7;

A = -(xMidpoint(i)-xVertex(j))*cos(phi(j))-(yMidpoint(i)-yVertex(j))*sin(phi(j));
B = (xMidpoint(i)-xVertex(j))^2 + (yMidpoint(i)-yVertex(j))^2;
C = sin(phi(i)-phi(j));
D = (yMidpoint(i)-yVertex(j))*cos(phi(i))-(xMidpoint(i)-xVertex(j))*sin(phi(i));
Sj = sqrt((xVertex(j+1)-xVertex(j))^2 + (yVertex(j+1)-yVertex(j))^2);
E = sqrt(B-A^2);

disp([A B C D Sj E])

disp((C/2)*log((Sj^2+2*A*Sj+B)/B)+((D-A*C)/E)*(atan2((Sj+A),E)-atan2(A,E)))