function [RAoutTMP,RAinTMP,RTMP,RARTMP,BPTMP,RABPTMP,HoxTMP,KroxTMP] = reactions(state,reaction,RAout,RAin,R,RAR,BP,RABP,Hox,Krox,production,X,f0,beta,kp,gamma,kdeg,mon,moff,jalpha,jbeta,bpdeg1, ...
   Vr,Vbp,rdeg1,rdeg2,DiffRA,d,e,lambda,cH,cK,kappaH,kappaK,dH,dK)

        %Solve Reactions
        if reaction == 1
            RAoutTMP  = production - beta*RAout + kp*RAin + DiffRA;

            RAinTMP   = beta*RAout - kp*RAin -((kdeg*RAR) ./(gamma+RAR)).*RAin ...
                    - mon*RAin.*BP + moff*RABP - rdeg1*RAin;
            RTMP      = Vr - rdeg2*R- jalpha*RABP.*R + jbeta*BP.*RAR;

            RARTMP    = jalpha*RABP.*R - jbeta*BP.*RAR;
            BPTMP     = Vbp - bpdeg1*BP - mon*RAin.*BP + moff*RABP ...
                    + jalpha*RABP.*R - jbeta*BP.*RAR + ((d.*RAR)./(e+RAR));
            RABPTMP   = mon*RAin.*BP - moff*RABP -jalpha*RABP.*R ...
                    + jbeta*BP.*RAR;
        elseif reaction == 2
            RASignal  = (lambda*RAR).*(lambda*RAR); 
            Cyp       = kdeg*RASignal./(1+RASignal + f0*exp(-.00001*(400-X)));
            
            RAoutTMP  = production - beta*RAout + kp*RAin + DiffRA;

            RAinTMP   = beta*RAout - kp*RAin -Cyp.*RAin ...
                    - mon*RAin.*BP + moff*RABP - rdeg1*RAin;
            RTMP      = Vr - rdeg2*R- jalpha*RABP.*R + jbeta*BP.*RAR;

            RARTMP    = jalpha*RABP.*R - jbeta*BP.*RAR;
            BPTMP     = Vbp - bpdeg1*BP - mon*RAin.*BP + moff*RABP ...
                    + jalpha*RABP.*R - jbeta*BP.*RAR + ((d.*RAR)./(e+RAR));
            RABPTMP   = mon*RAin.*BP - moff*RABP -jalpha*RABP.*R ...
                    + jbeta*BP.*RAR;
        end
            
        if state == 3
            HoxTMP    = (cH.*Hox.*Hox   + kappaH.*RAin.*RAin)./(1+cH.*Hox.*Hox + cK.*Krox.*Krox + kappaH.*RAin.*RAin) - dH.*Hox;
            KroxTMP   = (cK.*Krox.*Krox + kappaK.*RAin.*RAin)./(1+cH.*Hox.*Hox + cK.*Krox.*Krox + kappaH.*RAin.*RAin) - dK.*Krox;
        else
            HoxTMP = Hox;
            KroxTMP= Krox;
        end
end