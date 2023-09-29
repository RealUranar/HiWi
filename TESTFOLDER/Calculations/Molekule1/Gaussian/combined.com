%Chk=F1a_ortho.chk
#P RHF/6-31G* Opt

COMMENT

0  1
C          -6.63120        -1.10050        -0.53780
C          -6.78870         0.27180        -0.56920
C          -7.56700         0.21510         1.65190
C          -7.40940        -1.16270         1.68500
C          -6.93630        -1.84130         0.58270
H          -6.25820        -1.66700        -1.39610
H          -6.53550         0.81060        -1.47270
H          -7.94150         0.73150         2.53180
H          -7.64450        -1.76010         2.55780
H          -6.81830        -2.92240         0.62460
N          -7.43960         2.37180         0.54190
N          -6.49990         3.19100         0.30440
C          -5.20710         2.66600         0.01260
C          -4.56120         1.86590         0.93330
C          -4.55820         2.92710        -1.17410
C          -3.31050         1.34680         0.67330
H          -5.05530         1.64740         1.87580
C          -3.30930         2.39930        -1.41480
H          -5.05260         3.55770        -1.91470
C          -2.64600         1.59810        -0.50990
H          -2.81270         0.71410         1.41480
H          -2.79120         2.59820        -2.34340
C          -7.26170         0.95790         0.53020
C          -0.78490         0.26200         0.22080
H          -1.46920        -0.61210         0.37440
H          -0.70620         0.79460         1.18930
C           0.54660        -0.18430        -0.26690
C           1.47300        -0.98670         0.36880
H           1.36000        -1.42230         1.34860
N           1.12150         0.11810        -1.45860
N           2.28940        -0.45550        -1.52540
N           2.53990        -1.14230        -0.42230
C           3.75690        -1.91670        -0.14160
H           4.15970        -2.21910        -1.14200
H           3.53270        -2.83310         0.41530
C           4.78850        -1.06700         0.53160
H           5.01570        -0.18760        -0.08910
H           4.38620        -0.79080         1.53370
C           6.01300        -1.91320         0.79770
O           7.02540        -1.13980         1.44210
H           5.72790        -2.78300         1.40510
C           6.66950        -2.37250        -0.50440
H           7.61110        -2.88800        -0.17830
H           6.02780        -3.01580        -1.10610
O           7.00890        -1.23100        -1.22880
C           7.41090         0.00000        -0.62630
C           7.41800         0.04200         0.70670
C           7.81090         1.24410         1.31380
H           7.89720         1.51970         2.37580
S           8.17690         2.33460        -0.09020
C           7.77190         1.08550        -1.34090
H           7.83900         1.27360        -2.39880
O          -1.38250         1.07030        -0.76590


--Link1--
%Chk=F1a_ortho.chk
#P HF/6-31G* SCF=Tight Geom=AllCheck Guess=Read
Pop=MK IOp(6/33=2, 6/41=10, 6/42=17)