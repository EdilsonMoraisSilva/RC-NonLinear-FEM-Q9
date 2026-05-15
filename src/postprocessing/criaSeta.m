

  
function [ux, uy] = criaSeta(Dados,q, coordCentro, escalaVetor,z);
 

switch Dados.InicioSeta
    case 'Botton-Up'
           normq = norm(q);
           % create the arrow horizontally
           vx(1) = 0.0;
           vy(1) = 0.0;

           vx(2) = normq;
           vy(2) = 0;

           vx(3) = normq;
           vy(3) = -0.1*normq;

           vx(4) = 1.2*normq;
           vy(4) = 0;

           vx(5) = normq;
           vy(5) = 0.1*normq;

           vx(6) = normq;
           vy(6) = 0;
           
    case 'Top-down'
        
           normq = norm(q);
           % create the arrow horizontally
           vx(1) = 0.0;
           vy(1) = 0.0;

           vx(2) = -(1.2*normq-normq);
           vy(2) = -0.1*normq;

           vx(3) = -(1.2*normq-normq);
           vy(3) = 0;

           vx(4) = -(1.2*normq-normq);
           vy(4) = 0.1*normq;
            
           vx(5) = 0;
           vy(5) = 0;
           
           vx(6) = -(1.2*normq-normq);
           vy(6) = -0.1*normq;
           
           vx(7) = -(1.2*normq-normq);
           vy(7) = 0;
                     
           vx(8) = -(1.2*normq-normq);
           vy(8) = 0;

           vx(9) = -1.2*normq;
           vy(9) = 0;
end
        
        
        

   % rotate the arrow to the correct inclination
   cost = q(1)/normq;
   sent = q(2)/normq;
   
   r = [cost -sent;
        sent cost];
    
   u = r * [vx;
            vy];
   
   ux = u(1,:);
   uy = u(2,:);
   
   ux = ux*escalaVetor + coordCentro(1);
   uy = uy*escalaVetor + coordCentro(2);
   %uz(1:6) = z;
   

return
      
   
   