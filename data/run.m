b[nF_, nH_] := {
    20 nF/9 + nH/6, -22/3 + 4 nF/3 + nH/6, -11 + 4 nF/3};

b2HSM = b[3, 2];

BetagY = M D[GGy[M], M] - b2HSM[[1]]/(16 Pi^2)*(GGy[M])^3;
Betag2 = M D[GG2[M], M] - b2HSM[[2]]/(16 Pi^2)*(GG2[M])^3;
Betag3 = M D[GG3[M], M] - b2HSM[[3]]/(16 Pi^2)*(GG3[M])^3;
Betayt = M D[YT[M],
      M] - (9/2 YT[M]^2 - 17 GGy[M]^2/12 - 9 GG2[M]^2/4 - 8 GG3[M]^2)*
     YT[M]/(16 Pi^2);
Betal1 = M D[La1[M],
      M] - (12 La1[M]^2 + 4 La3[M]^2 + 4 La3[M]*La4[M] + 2 La4[M]^2 +
       2 La5[M]^2 + 24 La6[M]^2 +
       3 (3 GG2[M]^4 + GGy[M]^4 + 2 GG2[M]^2*GGy[M]^2)/4 -
       3 La1[M] (3 GG2[M]^2 + GGy[M]^2))/(16 Pi^2);
Betal2 = M D[La2[M],
      M] - (12 La2[M]^2 + 4 La3[M]^2 + 4 La3[M]*La4[M] + 2 La4[M]^2 +
       2 La5[M]^2 + 24 La7[M]^2 +
       3 (3 GG2[M]^4 + GGy[M]^4 + 2 GG2[M]^2*GGy[M]^2)/4 -
       3 La2[M] (3 GG2[M]^2 + GGy[M]^2 - 4 (YT[M]*Sin[beta])^2) -
       12 (YT[M]*Sin[beta])^4)/(16 Pi^2);
Betal3 = M D[La3[M],
      M] - ((La1[M] + La2[M]) (6 La3[M] + 2 La4[M]) + 4 La3[M]^2 +
       2 La4[M]^2 + 2 La5[M]^2 + 4 La6[M]^2 + 4 La7[M]^2 +
       16 La6[M]*La7[M] +
       3 (3 GG2[M]^4 + GGy[M]^4 - 2 GG2[M]^2*GGy[M]^2)/4 -
       3 La3[
         M] (3 GG2[M]^2 + GGy[M]^2 -
          2 (YT[M]*Sin[beta])^2))/(16 Pi^2);
Betal4 = M D[La4[M],
      M] - (2 (La1[M] + La2[M]) La4[M] + 8 La3[M]*La4[M] +
       4 La4[M]^2 + 8 La5[M]^2 + 10 La6[M]^2 + 10 La7[M]^2 +
       4 La6[M]*La7[M] + 3 GG2[M]^2*GGy[M]^2 -
       3 La4[M] (3 GG2[M]^2 + GGy[M]^2 -
          2 (YT[M]*Sin[beta])^2))/(16 Pi^2);
Betal5 = M D[La5[M],
      M] - (La5[M] (2 La1[M] + 2 La2[M] + 8 La3[M] + 12 La4[M]) +
       10 (La6[M]^2 + La7[M]^2) + 4 La6[M]*La7[M] -
       3 La5[M] (3 GG2[M]^2 + GGy[M]^2 -
          2 (YT[M]*Sin[beta])^2))/(16 Pi^2);
Betal6 = M D[La6[M],
      M] - (La6[M] (12 La1[M] + 6 La3[M] + 8 La4[M]) +
       La7[M] (6 La3[M] + 4 La4[M]) + 10 La5[M]*La6[M] +
       2 La5[M]*La7[M] -
       3 La6[M] (3 GG2[M]^2 + GGy[M]^2 -
          2 (YT[M]*Sin[beta])^2))/(16 Pi^2);
Betal7 = M D[La7[M],
      M] - (La7[M] (12 La2[M] + 6 La3[M] + 8 La4[M]) +
       La6[M] (6 La3[M] + 4 La4[M]) + 10 La5[M]*La7[M] +
       2 La5[M]*La6[M] -
       3 La7[M] (3 GG2[M]^2 + GGy[M]^2 -
          2 (YT[M]*Sin[beta])^2))/(16 Pi^2);

l1bound = 0.9031500 (*2*0.448538*);
l2bound = 0.2849600 (*2*0.138807*);
l3bound = 0.2233890 (*2.22612*);
l4bound = -0.7657490 (*-0.765326*);
l5bound = -0.7657490 (*-0.765326*);
l6bound = 0;
l7bound = 0;
Mbc = 10^7;
gybound = 0.384231;
gwbound = 0.595574;
gsbound = 0.765717;
ytbound = 0.7617750 (*0.761686*);
beta = ArcTan[2.246540 (*1.88701*)];

RUN = NDSolve[{BetagY == 0, Betag2 == 0, Betag3 == 0, Betayt == 0,
    Betal1 == 0, Betal2 == 0, Betal3 == 0, Betal4 == 0, Betal5 == 0,
    Betal6 == 0, Betal7 == 0, GGy[Mbc] == gybound, GG2[Mbc] == gwbound,
     GG3[Mbc] == gsbound, YT[Mbc] == ytbound, La1[Mbc] == l1bound,
    La2[Mbc] == l2bound, La3[Mbc] == l3bound, La4[Mbc] == l4bound,
    La5[Mbc] == l5bound, La6[Mbc] == l6bound,
    La7[Mbc] == l7bound}, {GGy[M], GG2[M], GG3[M], YT[M], La1[M],
    La2[M], La3[M], La4[M], La5[M], La6[M], La7[M]}, {M, 1, Mbc}]

Print[Evaluate[La1[M] /. RUN /. M->125]];
Print[Evaluate[La2[M] /. RUN /. M->125]];
