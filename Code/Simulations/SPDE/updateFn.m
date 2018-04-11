function [RAout,RAin,R,RAR,BP,RABP,Hox,Krox] = updateFn(state,dt,RAout,RAoutN1,RAoutN2,RAoutN3,RAoutN4,RAin,RAinN1,RAinN2,RAinN3,RAinN4,R,RN1,RN2,RN3,RN4,...
    RAR,RARN1,RARN2,RARN3,RARN4,BP,BPN1,BPN2,BPN3,BPN4,RABP,RABPN1,RABPN2,RABPN3,RABPN4,Hox,Hox1,Hox2,Hox3,Hox4,Krox,Krox1,Krox2,Krox3,Krox4,RAoutNoise,RAinNoise,RARNoise,HoxNoise,KroxNoise)

                RAout = RAout + (dt/6)*(RAoutN1 + 2*RAoutN2 + 2*RAoutN3 + RAoutN4) + RAoutNoise;
                RAin  = RAin  + (dt/6)*(RAinN1 + 2*RAinN2 + 2*RAinN3 + RAinN4) + RAinNoise ;
                R     = R     + (dt/6)*(RN1 + 2*RN2 + 2*RN3 + RN4);
                RAR   = RAR   + (dt/6)*(RARN1 + 2*RARN2 + 2*RARN3 + RARN4) + RARNoise;
                BP    = BP    + (dt/6)*(BPN1 + 2*BPN2 + 2*BPN3 + BPN4);
                RABP  = RABP  + (dt/6)*(RABPN1 + 2*RABPN2 + 2*RABPN3 + RABPN4);
                if state == 3
                    Hox   = Hox   + (dt/6)*(Hox1 + 2*Hox2 + 2*Hox3 + Hox4)  + HoxNoise;
                    Krox  = Krox  + (dt/6)*(Krox1+ 2*Krox2+ 2*Krox3+ Krox4) + KroxNoise;
                end
end