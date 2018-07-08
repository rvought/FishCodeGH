            SET    1,1,0           ; Get rate & scaling OK

            VAR    V45,LoopC=0     ; Define variable for section loops
            VAR    V46,RampC=0     ; Define variable for ramp loops
            VAR    V47,DelayC=0    ; Define variable for delay loops
            VAR    V48,Delay2=0    ;  and another one
            VAR    V49,Delay3=0    ;  and another one
            VAR    V50,Delay4=0    ;  and another one
            VAR    V51,Delay5=0    ;  and another one

E0:         DIGOUT [......00]
            DAC    0,0
            DAC    1,0
            DELAY  s(0.492)-1
            SZ     0,0.2           ; Set cosine amplitude
            OFFSET 0,0             ;  cosine centre
            ANGLE  0,90            ;  cosine phase
            RATE   0,Hz(0.5)       ;  set rate and start cosine off
            DELAY  s(9.999)-1
            RATE   0,0             ; Stop cosine output
            DAC    0,0
            DELAY  s(0.499)-1
            HALT                   ; End of this sequence section

