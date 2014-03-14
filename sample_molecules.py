charge=["X",
"H","He",
"Li","Be","B","C","N","O","F","Ne",
"Na","Mg","Al","Si","P","S","Cl","Ar"
]


def atom(*args, **kwargs):
    A, = args
    assert A in charge
    return """%s"""%(A)

def diatom(*args, **kwargs ):
    A, B, R = args
    try:
        assert A in charge
        assert B in charge
        _R = float(R)
    except():
       raise SystemExit(1)
    return """%s
%s 1 R
 Variables:
R=%s
"""%(A, B, R)

def triatom_AB2(*args, **kwargs ):
    """Symmetric three-atom AB2 molecule"""
    A, B, r, a = args
    try:
        assert A in charge
        assert B in charge
        _r = float(r)
        _a = float(a)
    except():
       raise SystemExit(1)
    basis = kwargs.get("basis", "6-31+G**")
    return """
%s
%s 1 R
%s 1 R 2 A
 Variables:
R=%s
A=%s
"""%(A, B, B, r, a)


molinp={
   "h":atom("H"),
   "o":atom("O"),
   "f":atom("F"),
   "h2":diatom("H","H",0.7361),
   "oh":diatom("O","H",0.9665),
   "hcl":diatom("Cl","H",1.2744),
   "hf":diatom("F","H",0.9151),
   "hs":diatom("S","H",1.3394),
   "nh":diatom("N","H",1.0317),
   "o2":diatom("O","H",1.3394),
   "sh":diatom("S","H",1.3405),
   "nh2":triatom_AB2("N","H",1.019,104.270),
   "h2o":triatom_AB2("O","H",0.95,105),
   "sh2":triatom_AB2("S","H",1.3361,93.19)
   }


molinp["nh2"]="""N
H 1 R
H 1 R 2 A
 Variables:
R=1.0191
A=104.27
"""

molinp["n2h"]="""N
H 1 RNH
N 1 RNN 2 A 
 Variables:
RNH=1.041 
RNN=1.169
A=116.566
"""

molinp["ph2"]="""P
H 1 R
H 1 R 2 A
 Variables:
R=1.4149
A=92.4849
"""

molinp["nh3"]="""N
X 1 1.0
H 1 R 2 A
H 1 R 2 A 3 120
H 1 R 2 A 4 120
 Variables:
R=1.0058
A=110.2596
"""
#X help atom, A calculated from
#sin A = 2/sqrt(3)*sin(HCH/2), HCH=108.6715

molinp["ph3"]="""P
X 1 1.0
H 1 R 2 A
H 1 R 2 A 3 120
H 1 R 2 A 4 120
 Variables:
R=1.410 
A=122.0756
"""

molinp["h3"]="""H
H 1 R
H 2 R 1 180
 Variables:
R=0.925
"""

molinp["h2h"]="""H
H 1 R1
H 2 R2 1 180
 Variables:
R1=0.9246
R2=0.9246
"""

molinp["ch3"]="""C
H 1  R
H 1  R  2  120
H 1  R  2  120 3 180
 Variables:
R=1.074
"""

molinp["ch4"]="""C 
H 1 R
H 1 R 2 109.4712206
H 1 R 2 109.4712206 3 120
H 1 R 2 109.4712206 3 240
 Variables:
R=1.085
"""

molinp["ch4_par"]="""C 
H 1 R
H 1 R 2 A
H 1 R 2 A 3 D1
H 1 R 2 A 3 D2
 Variables:
R=1.085
A=109.4712206
D1=120
D2=240
"""

molinp["ch3oh"]="""C
O 1 RCO
H 1 RCH1  2 OCH1
H 1 RCH2  2 OCH2  3 D1
H 2 ROH   1 COH   3 D2
H 1 RCH1  2 OCH1  3 D3
 Variables:
RCO=1.408
RCH1=1.090
OCH1=112.020
RCH2=1.084
OCH2=106.854
D1=241.4
ROH=0.953
COH=109.807
D2=298.6
D3=122.8
"""

molinp["ch2oh"]="""C
H 1 RCH1
H 1 RCH2   2 AHCH
O 1 RCO    2 AHCO  3 DHHCO
H 4 ROH    1 ACOH  2 DHCOH
 Variables:
RCH1=1.074
RCH2=1.078
AHCH=120.906
RCO=1.356
AHCO=113.472  
DHHCO=150.2
ROH=0.953
ACOH=110.744
DHCOH=182.6
"""

molinp["c2h6"]="""C 
C 1 RCC  
H 1 RCH  2 A
H 1 RCH  2 A  3 120.0
H 1 RCH  2 A  3 240.0
H 2 RCH  1 A  3  60.0
H 2 RCH  1 A  3 180.0
H 2 RCH  1 A  3 300.0
 Variables:
RCC=1.522
RCH=1.087
A=111.295
"""

molinp["c2h5"]="""C
C 1 RCC
H 1 R1 2 A1
H 2 R2 1 A2 3 -D2
H 2 R2  1 A2 3 D2
H 1 R3  2 A3 4 D3
H 1 R3  2 A3 5 -D3
 Variables:
RCC=1.485
R1=1.094
A1=111.744
R2=1.077
A2=120.892
D2=85.534
R3=1.088
A3=111.661
D3=33.871
"""

molinp["n2h2"]="""N
H 1 RNH
N 1 RNN  2 ANNH
H 3 RNH  1 ANNH  2 180.0
 Variables:
RNH=1.023
RNN=1.227
ANNH=107.343
"""

molinp["TS01"]="""Cl
H   1 R1
H   2 R2  1 180
 Variables:
R1=1.371
R2=1.114
"""

molinp["TS02"]="""O   
H   1 R1
H   1 R2  2 A2
H   3 R3  1 A3  2   0.0
 Variables:
R1=0.963
R2=1.264
A2=100.451
R3=0.848
A3=166.210
"""

molinp["TS03"]="""C
H   1 R1   
H   1 R1  2 A1
H   1 R1  2 A1  3 D1
H   1 R2  2 A2  3 D2
H   2 R3  1 A3  3 D3
 Variables:
R1=1.080
A1=114.677
D1=135.8
R2=1.381  
A2=103.565
D2=247.9
R3=2.742
A3=53.912
D3=112.1
"""

molinp["TS04"]="""C   1
H   1 R1
H   1 R1  2 A1
H   1 R2  2 A2  3 D1
O   4 R3  1 A3  2 D2
H   5 R4  4 A4  1   0.0
H   1 R1  2 A1  3 D3
 Variables:
R1=1.082
R2=1.244
R3=1.251
R4=0.952
A1=113.227
A2=104.276
A3=174.761
A4=100.548
D1=116.7
D2=239.5
D3=229.4
"""

molinp["TS05"]="""C   1
O   1 R1
H   1 R2  2 A1
H   1 R3  2 A2  3 D1
H   2 R4  1 A3  3 D2
H   1 R5  2 A4  3 D3
H   2 R6  1 A5  3 D4
 Variables:
R1=1.374
R2=1.087
R3=1.081
R4=0.958
R5=1.330
R6=3.032
A1=115.089
A2=109.674
A3=110.524
A4=110.196
A5=45.235
D1=230.6
D2=314.3
D3=116.8
D4=243.3
"""

molinp["TS06"]="""H
H 1 R
H 2 R 1 180
 Variables:
R=0.925
"""

molinp["TS07"]="""H
O 1 R1
H 2 R2  1 A1
N 3 R3  2 A2  1 D1
H 4 R4  3 A3  2 D2
H 4 R5  3 A4  2 D3
 Variables:
R1=0.959
R2=1.248
R3=1.151
R4=1.011
R5=1.011
A1=108.803
A2=150.779
A3=110.384
A4=108.054
D1=305.5
D2=67.3
D3=309.5
"""

molinp["TS08"] = """Cl
H 1 R1
C 2 R2 1 180
H 3 R3 2 A 1 0
H 3 R3 2 A 1 120
H 3 R3 2 A 1 -120
 Variables:
R1 = 1.4265
R2 = 1.3981
R3 = 1.0796
A = 101.66
"""

molinp["TS09"] = """C
C 1 RCC
H 1 RCH1 2 ALPHA
H 1 RCH2 2 BETA 3 TAU
H 1 RCH2 2 BETA 3 -TAU
H 2 RCH4 1 GAMMA 3 180
H 2 RCH5 1 DELTA 3 CHI
H 2 RCH5 1 DELTA 3 -CHI
O 6 ROH4 2 EPS 1 180
H 9 ROH7 6 ZETA 2 0
 Variables:
RCC = 1.508
RCH1 = 1.0866
RCH2 = 1.084
RCH4 = 1.1670
RCH5 = 1.0834
ROH4 = 1.3601
ROH7 = 0.9675
ALPHA = 110.65
BETA = 111.01
GAMMA = 107.01
DELTA = 113.45
EPS = 175.42
ZETA = 95.90
TAU = 119.87
CHI = 116.51
"""


molinp["TS10"]="""F
H   1 R1
H   2 R2  1 180
 Variables:
R1=1.440 
R2=0.772
"""

molinp["TS11"] = """O
H 1 R1
C 2 R2 1 A1
H 3 R3 2 A2 1 180
H 3 R4 2 A3 4 D
H 3 R4 2 A3 4 -D
Variables:
R1 = 1.179
R2 = 1.289
R3 = 1.082
R4 = 1.081
A1 = 181.05
A2 = 103.9
A3 = 104.5
D = 120
"""

molinp["TS12"] = """H
H 1 1.2587
P 2 1.4851 1 170.61
H 3 R 2 A  1 D
H 3 R 2 A  1 -D
Variables:
R = 1.4152
A = 92.92
D = 45
"""

molinp["TS13"] = """Cl
H 1 R
H 1 R 2 180
 Variables:
R = 1.4850
"""

molinp["TS14"]="""O
H   1 R1
H   2 R2  1 180
 Variables:
R1=1.2145
R2=0.8941
"""

molinp["TS15"] = """H
N 1 1.06
N 2 1.25 1 105.8
H 3 1.17 2 106.2 1 0
H 4 1.20 3 173.5 2 0
"""


molinp["TS16"] = """H 
H 1 RH1H2
S 2 RSH2  1 H1H2S
H 3 RSH3  2 H2SH3 1 D
 Variables: 
RH1H2 = 1.2385
RSH2 = 1.4046
RSH3 = 1.3372
H1H2S = 175.26
H2SH3 = 95.27
D = 0
"""

molinp["TS17"] = """O
H 1 ROH
Cl 2 RClH 1 OHCl
 Variables:
ROH = 2.44
RClH = 2.66
OHCl = 131.4
"""

molinp["TS18"] = """C
H 1 CH2
H 1 CH3 2 H3CH2
H 1 CH4 2 H4CH2 3 H4CH2H3
H 1 CH4 3 H4CH2 2 -H4CH2H3
N 2 NH2 1 NH2C  3 0
H 6 NH7 2 H7NH2 1 0
 Variables:
CH2 = 1.4447 
CH3 = 1.0769
CH4 = 1.0773
NH2 = 1.2091
NH7 = 1.0244
H3CH2 = 100.77
H4CH2 = 104.51
NH2C = 175.14
H7NH2 = 103.66
H4CH2H3 = 119.42
"""

molinp["TS19"] = """C
H 1 H2C1
C 1 C3C1 2 C3C1H2
H 1 H4C1 3 H4C1C3 2 D1
H 1 H4C1 3 H4C1C3 2 -D1
N 2 N6H2 1 N6H2C1 3 0
H 6 H7N6 2 H7N6H2 1 0
H 3 H8C3 1 H8C3C1 2 0
H 3 H9C3 1 H9C3C1 2 D2
H 3 H9C3 1 H9C3C1 2 -D2
 Variables:
H2C1 = 1.3960
C3C1 = 1.5118
H4C1 = 1.0820
N6H2 = 1.2326
H7N6 = 1.0181
H8C3 = 1.0887
H9C3 = 1.0850
C3C1H2 = 105.85
H4C1C3 = 103.28
N6H2C1 = 177.05
H7N6H2 = 99.33
H8C3C1 = 111.76
H9C3C1 = 107.77
D1 = 113.50
D2 = 60.19
"""

molinp["TS20"] = """C
C 1 CC
H 1 H3C 2 H3CC
N 3 NH3 1 NH3C 2 0
H 2 H5C 3 H5CC 1 H5CCH3
H 2 H5C 3 H5CC 1 -H5CCH3
H 2 H7C 1 H7CC 3 180
H 1 H8C 2 H8CC 7 H8CCH7 
H 1 H8C 2 H8CC 7 -H8CCH7 
H 4 H10N 3 H10NH3 1 H10NH3C
H 4 H10N 3 H10NH3 1 -H10NH3C
Variables:
CC = 1.517
H3C = 1.254
H3CC = 111.088
NH3 = 1.313
NH3C = 175.404
H5C = 1.093
H5CC = 111.088
H5CCH3 = 60
H7C = 1.096 
H7CC = 111.088
H8C = 1.092
H8CC = 114.157 
H8CCH7 = 60
H10N = 1.024
H10NH3 = 100.007
H10NH3C = 60
"""

