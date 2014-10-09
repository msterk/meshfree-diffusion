clear domain;
global domain;
domain.nodes = [
0.000000000000000 0.000000000000000 1 0.000000000000000;
0.250000000000000 0.000000000000000 1 0.000000000000000;
0.500000000000000 0.000000000000000 1 0.000000000000000;
0.750000000000000 0.000000000000000 1 0.000000000000000;
1.000000000000000 0.000000000000000 1 0.000000000000000;
0.000000000000000 0.250000000000000 1 0.000000000000000;
0.250000000000000 0.250000000000000 0 1.000000000000000;
0.500000000000000 0.250000000000000 0 1.000000000000000;
0.750000000000000 0.250000000000000 0 1.000000000000000;
1.000000000000000 0.250000000000000 1 0.000000000000000;
0.000000000000000 0.500000000000000 1 0.000000000000000;
0.250000000000000 0.500000000000000 0 1.000000000000000;
0.500000000000000 0.500000000000000 0 1.000000000000000;
0.750000000000000 0.500000000000000 0 1.000000000000000;
1.000000000000000 0.500000000000000 1 0.000000000000000;
0.000000000000000 0.750000000000000 1 0.000000000000000;
0.250000000000000 0.750000000000000 0 1.000000000000000;
0.500000000000000 0.750000000000000 0 1.000000000000000;
0.750000000000000 0.750000000000000 0 1.000000000000000;
1.000000000000000 0.750000000000000 1 0.000000000000000;
0.000000000000000 1.000000000000000 1 0.000000000000000;
0.250000000000000 1.000000000000000 1 0.000000000000000;
0.500000000000000 1.000000000000000 1 0.000000000000000;
0.750000000000000 1.000000000000000 1 0.000000000000000;
1.000000000000000 1.000000000000000 1 0.000000000000000;
];
clear mfreeParams; global mfreeParams; 
mfreeParams.numNodes = 25;
mfreeParams.randSeed = 17;
domain.rHole = 0.000000000000000;
domain.xHole = 0.600000000000000;
domain.yHole = 0.350000000000000;
mfreeParams.C = 0.001000000000000;
mfreeParams.irregularity = 0.300000000000000;
mfreeParams.baseFtype = 'mls';
mfreeParams.mlsDegree = 2;
mfreeParams.quadDegree = 5;
mfreeParams.timeStep = 0.001000000000000;
mfreeParams.nI = 13;
mfreeParams.betaQuad = 0.700000000000000;
mfreeParams.preFindQuadSupport = 0;
mfreeParams.numRepetitions = 1;
mfreeParams.dataDistrType = 0;
reg25.A=spalloc(25, 25, 1);
reg25.A(1,1)=-9.324016336347762;
reg25.A(1,2)=4.111671958034440;
reg25.A(1,3)=-0.490434456937199;
reg25.A(1,4)=-0.149045458549265;
reg25.A(1,6)=4.111671958034458;
reg25.A(1,7)=4.254935183754265;
reg25.A(1,8)=-0.812004012738423;
reg25.A(1,9)=0.000869752349042;
reg25.A(1,11)=-0.490434456937202;
reg25.A(1,12)=-0.812004012738426;
reg25.A(1,13)=-0.253034411723704;
reg25.A(1,16)=-0.149045458549265;
reg25.A(1,17)=0.000869752349042;
reg25.A(2,1)=-1.583701417813938;
reg25.A(2,2)=-2.889499261014200;
reg25.A(2,3)=-1.403533904428681;
reg25.A(2,4)=-0.123265416743180;
reg25.A(2,6)=2.136172586194484;
reg25.A(2,7)=3.729843576056023;
reg25.A(2,8)=1.999068449102314;
reg25.A(2,9)=0.134915388647185;
reg25.A(2,11)=-0.536407081406246;
reg25.A(2,12)=-0.888536575964703;
reg25.A(2,13)=-0.547342283750757;
reg25.A(2,14)=-0.027714058878296;
reg25.A(2,17)=0.000000000000000;
reg25.A(3,1)=-0.035572430038485;
reg25.A(3,2)=-1.326650743558669;
reg25.A(3,3)=-3.275553652805677;
reg25.A(3,4)=-1.326650743558675;
reg25.A(3,5)=-0.035572430038486;
reg25.A(3,6)=-0.022469566128146;
reg25.A(3,7)=1.859433842320772;
reg25.A(3,8)=4.326071447614749;
reg25.A(3,9)=1.859433842320764;
reg25.A(3,10)=-0.022469566128146;
reg25.A(3,12)=-0.300615114095545;
reg25.A(3,13)=-1.398769771808908;
reg25.A(3,14)=-0.300615114095544;
reg25.A(4,1)=0.000000000000000;
reg25.A(4,2)=-0.123265416743178;
reg25.A(4,3)=-1.403533904428673;
reg25.A(4,4)=-2.889499261014199;
reg25.A(4,5)=-1.583701417813949;
reg25.A(4,7)=0.134915388647191;
reg25.A(4,8)=1.999068449102307;
reg25.A(4,9)=3.729843576056007;
reg25.A(4,10)=2.136172586194498;
reg25.A(4,12)=-0.027714058878295;
reg25.A(4,13)=-0.547342283750758;
reg25.A(4,14)=-0.888536575964711;
reg25.A(4,15)=-0.536407081406240;
reg25.A(5,2)=0.000921175499896;
reg25.A(5,3)=-0.025132646958045;
reg25.A(5,4)=0.081558557218736;
reg25.A(5,5)=0.940802574133600;
reg25.A(5,7)=0.000929164605909;
reg25.A(5,8)=0.006964792690617;
reg25.A(5,9)=-0.083901494197756;
reg25.A(5,10)=0.081558557218660;
reg25.A(5,13)=0.012616833949940;
reg25.A(5,14)=0.006964792690647;
reg25.A(5,15)=-0.025132646958017;
reg25.A(5,19)=0.000929164605910;
reg25.A(5,20)=0.000921175499900;
reg25.A(6,1)=-1.583701417813950;
reg25.A(6,2)=2.136172586194482;
reg25.A(6,3)=-0.536407081406244;
reg25.A(6,6)=-2.889499261014196;
reg25.A(6,7)=3.729843576056019;
reg25.A(6,8)=-0.888536575964701;
reg25.A(6,11)=-1.403533904428677;
reg25.A(6,12)=1.999068449102309;
reg25.A(6,13)=-0.547342283750755;
reg25.A(6,16)=-0.123265416743181;
reg25.A(6,17)=0.134915388647184;
reg25.A(6,18)=-0.027714058878296;
reg25.A(6,21)=0.000000000000000;
reg25.A(7,1)=0.014136310198338;
reg25.A(7,2)=0.035460566451191;
reg25.A(7,3)=0.015877392263175;
reg25.A(7,4)=-0.000849569707348;
reg25.A(7,5)=-0.000644326564222;
reg25.A(7,6)=0.035460566451191;
reg25.A(7,7)=0.161518496483850;
reg25.A(7,8)=0.075869417154706;
reg25.A(7,9)=0.006457201243340;
reg25.A(7,10)=-0.000703921646040;
reg25.A(7,11)=0.015877392263175;
reg25.A(7,12)=0.075869417154706;
reg25.A(7,13)=0.015750379019093;
reg25.A(7,14)=0.000425818412815;
reg25.A(7,15)=-0.000155399981270;
reg25.A(7,16)=-0.000849569707348;
reg25.A(7,17)=0.006457201243341;
reg25.A(7,18)=0.000425818412814;
reg25.A(7,19)=-0.000053415881286;
reg25.A(7,20)=-0.000016656909212;
reg25.A(7,21)=-0.000644326564222;
reg25.A(7,22)=-0.000703921646040;
reg25.A(7,23)=-0.000155399981270;
reg25.A(7,24)=-0.000016656909212;
reg25.A(8,1)=-0.000552598197074;
reg25.A(8,2)=0.011961280597280;
reg25.A(8,3)=0.044918913130209;
reg25.A(8,4)=0.011961280597277;
reg25.A(8,5)=-0.000552598197074;
reg25.A(8,6)=0.009059211150560;
reg25.A(8,7)=0.053530709715060;
reg25.A(8,8)=0.186237998243230;
reg25.A(8,9)=0.053530709715052;
reg25.A(8,10)=0.009059211150560;
reg25.A(8,11)=0.002035650920755;
reg25.A(8,12)=0.016328780522327;
reg25.A(8,13)=0.075011978626682;
reg25.A(8,14)=0.016328780522322;
reg25.A(8,15)=0.002035650920757;
reg25.A(8,16)=-0.000517566771717;
reg25.A(8,17)=0.000872247812292;
reg25.A(8,18)=0.007453105459234;
reg25.A(8,19)=0.000872247812292;
reg25.A(8,20)=-0.000517566771715;
reg25.A(8,21)=-0.000157926658175;
reg25.A(8,22)=-0.000251180364950;
reg25.A(8,23)=-0.000635535814857;
reg25.A(8,24)=-0.000251180364950;
reg25.A(8,25)=-0.000157926658175;
reg25.A(9,1)=-0.000644326564222;
reg25.A(9,2)=-0.000849569707347;
reg25.A(9,3)=0.015877392263178;
reg25.A(9,4)=0.035460566451189;
reg25.A(9,5)=0.014136310198339;
reg25.A(9,6)=-0.000703921646040;
reg25.A(9,7)=0.006457201243342;
reg25.A(9,8)=0.075869417154711;
reg25.A(9,9)=0.161518496483844;
reg25.A(9,10)=0.035460566451193;
reg25.A(9,11)=-0.000155399981270;
reg25.A(9,12)=0.000425818412815;
reg25.A(9,13)=0.015750379019089;
reg25.A(9,14)=0.075869417154706;
reg25.A(9,15)=0.015877392263175;
reg25.A(9,16)=-0.000016656909212;
reg25.A(9,17)=-0.000053415881286;
reg25.A(9,18)=0.000425818412815;
reg25.A(9,19)=0.006457201243342;
reg25.A(9,20)=-0.000849569707348;
reg25.A(9,22)=-0.000016656909212;
reg25.A(9,23)=-0.000155399981270;
reg25.A(9,24)=-0.000703921646039;
reg25.A(9,25)=-0.000644326564222;
reg25.A(10,3)=-0.013029845984023;
reg25.A(10,4)=-0.066167344054254;
reg25.A(10,5)=0.160539267612127;
reg25.A(10,7)=0.000000000000000;
reg25.A(10,8)=0.015748825270937;
reg25.A(10,9)=0.110576487688360;
reg25.A(10,10)=0.629648454319150;
reg25.A(10,13)=-0.002153251676373;
reg25.A(10,14)=-0.003160665040734;
reg25.A(10,15)=0.249340149438668;
reg25.A(10,18)=-0.000565727610419;
reg25.A(10,19)=-0.041248478593576;
reg25.A(10,20)=-0.039527871369858;
reg25.A(11,1)=-0.035572430038484;
reg25.A(11,2)=-0.022469566128146;
reg25.A(11,6)=-1.326650743558666;
reg25.A(11,7)=1.859433842320771;
reg25.A(11,8)=-0.300615114095545;
reg25.A(11,11)=-3.275553652805682;
reg25.A(11,12)=4.326071447614746;
reg25.A(11,13)=-1.398769771808908;
reg25.A(11,16)=-1.326650743558674;
reg25.A(11,17)=1.859433842320766;
reg25.A(11,18)=-0.300615114095544;
reg25.A(11,21)=-0.035572430038485;
reg25.A(11,22)=-0.022469566128146;
reg25.A(12,1)=-0.000552598197074;
reg25.A(12,2)=0.009059211150560;
reg25.A(12,3)=0.002035650920755;
reg25.A(12,4)=-0.000517566771717;
reg25.A(12,5)=-0.000157926658175;
reg25.A(12,6)=0.011961280597280;
reg25.A(12,7)=0.053530709715060;
reg25.A(12,8)=0.016328780522327;
reg25.A(12,9)=0.000872247812292;
reg25.A(12,10)=-0.000251180364950;
reg25.A(12,11)=0.044918913130209;
reg25.A(12,12)=0.186237998243230;
reg25.A(12,13)=0.075011978626682;
reg25.A(12,14)=0.007453105459234;
reg25.A(12,15)=-0.000635535814857;
reg25.A(12,16)=0.011961280597277;
reg25.A(12,17)=0.053530709715051;
reg25.A(12,18)=0.016328780522322;
reg25.A(12,19)=0.000872247812292;
reg25.A(12,20)=-0.000251180364950;
reg25.A(12,21)=-0.000552598197074;
reg25.A(12,22)=0.009059211150560;
reg25.A(12,23)=0.002035650920757;
reg25.A(12,24)=-0.000517566771715;
reg25.A(12,25)=-0.000157926658175;
reg25.A(13,1)=-0.000483234263835;
reg25.A(13,2)=0.000081184287111;
reg25.A(13,3)=0.007371805085845;
reg25.A(13,4)=0.000081184287112;
reg25.A(13,5)=-0.000483234263835;
reg25.A(13,6)=0.000081184287110;
reg25.A(13,7)=0.011308749291760;
reg25.A(13,8)=0.058322424366573;
reg25.A(13,9)=0.011308749291755;
reg25.A(13,10)=0.000081184287111;
reg25.A(13,11)=0.007371805085845;
reg25.A(13,12)=0.058322424366570;
reg25.A(13,13)=0.181754505331748;
reg25.A(13,14)=0.058322424366565;
reg25.A(13,15)=0.007371805085846;
reg25.A(13,16)=0.000081184287112;
reg25.A(13,17)=0.011308749291755;
reg25.A(13,18)=0.058322424366565;
reg25.A(13,19)=0.011308749291759;
reg25.A(13,20)=0.000081184287111;
reg25.A(13,21)=-0.000483234263835;
reg25.A(13,22)=0.000081184287111;
reg25.A(13,23)=0.007371805085846;
reg25.A(13,24)=0.000081184287113;
reg25.A(13,25)=-0.000483234263836;
reg25.A(14,1)=-0.000157926658175;
reg25.A(14,2)=-0.000517566771717;
reg25.A(14,3)=0.002035650920757;
reg25.A(14,4)=0.009059211150557;
reg25.A(14,5)=-0.000552598197073;
reg25.A(14,6)=-0.000251180364950;
reg25.A(14,7)=0.000872247812294;
reg25.A(14,8)=0.016328780522359;
reg25.A(14,9)=0.053530709715037;
reg25.A(14,10)=0.011961280597280;
reg25.A(14,11)=-0.000635535814857;
reg25.A(14,12)=0.007453105459239;
reg25.A(14,13)=0.075011978626696;
reg25.A(14,14)=0.186237998243169;
reg25.A(14,15)=0.044918913130228;
reg25.A(14,16)=-0.000251180364950;
reg25.A(14,17)=0.000872247812291;
reg25.A(14,18)=0.016328780522321;
reg25.A(14,19)=0.053530709715053;
reg25.A(14,20)=0.011961280597292;
reg25.A(14,21)=-0.000157926658175;
reg25.A(14,22)=-0.000517566771717;
reg25.A(14,23)=0.002035650920758;
reg25.A(14,24)=0.009059211150552;
reg25.A(14,25)=-0.000552598197067;
reg25.A(15,4)=-0.007479499976368;
reg25.A(15,5)=-0.023867882316688;
reg25.A(15,8)=-0.001569843560556;
reg25.A(15,9)=-0.039918604317484;
reg25.A(15,10)=0.166877977050278;
reg25.A(15,13)=0.003139687121019;
reg25.A(15,14)=0.094796208587868;
reg25.A(15,15)=0.713979810532753;
reg25.A(15,18)=-0.001569843560556;
reg25.A(15,19)=-0.039918604317479;
reg25.A(15,20)=0.166877977050270;
reg25.A(15,24)=-0.007479499976368;
reg25.A(15,25)=-0.023867882316689;
reg25.A(16,1)=0.000000000000000;
reg25.A(16,6)=-0.123265416743179;
reg25.A(16,7)=0.134915388647194;
reg25.A(16,8)=-0.027714058878295;
reg25.A(16,11)=-1.403533904428681;
reg25.A(16,12)=1.999068449102316;
reg25.A(16,13)=-0.547342283750755;
reg25.A(16,16)=-2.889499261014204;
reg25.A(16,17)=3.729843576056006;
reg25.A(16,18)=-0.888536575964711;
reg25.A(16,21)=-1.583701417813937;
reg25.A(16,22)=2.136172586194495;
reg25.A(16,23)=-0.536407081406243;
reg25.A(17,1)=-0.000644326564222;
reg25.A(17,2)=-0.000703921646040;
reg25.A(17,3)=-0.000155399981270;
reg25.A(17,4)=-0.000016656909212;
reg25.A(17,6)=-0.000849569707347;
reg25.A(17,7)=0.006457201243342;
reg25.A(17,8)=0.000425818412815;
reg25.A(17,9)=-0.000053415881286;
reg25.A(17,10)=-0.000016656909212;
reg25.A(17,11)=0.015877392263178;
reg25.A(17,12)=0.075869417154710;
reg25.A(17,13)=0.015750379019089;
reg25.A(17,14)=0.000425818412815;
reg25.A(17,15)=-0.000155399981270;
reg25.A(17,16)=0.035460566451189;
reg25.A(17,17)=0.161518496483844;
reg25.A(17,18)=0.075869417154706;
reg25.A(17,19)=0.006457201243342;
reg25.A(17,20)=-0.000703921646039;
reg25.A(17,21)=0.014136310198339;
reg25.A(17,22)=0.035460566451192;
reg25.A(17,23)=0.015877392263176;
reg25.A(17,24)=-0.000849569707348;
reg25.A(17,25)=-0.000644326564222;
reg25.A(18,1)=-0.000157926658175;
reg25.A(18,2)=-0.000251180364950;
reg25.A(18,3)=-0.000635535814857;
reg25.A(18,4)=-0.000251180364950;
reg25.A(18,5)=-0.000157926658175;
reg25.A(18,6)=-0.000517566771717;
reg25.A(18,7)=0.000872247812294;
reg25.A(18,8)=0.007453105459239;
reg25.A(18,9)=0.000872247812291;
reg25.A(18,10)=-0.000517566771717;
reg25.A(18,11)=0.002035650920757;
reg25.A(18,12)=0.016328780522358;
reg25.A(18,13)=0.075011978626696;
reg25.A(18,14)=0.016328780522321;
reg25.A(18,15)=0.002035650920758;
reg25.A(18,16)=0.009059211150557;
reg25.A(18,17)=0.053530709715037;
reg25.A(18,18)=0.186237998243170;
reg25.A(18,19)=0.053530709715052;
reg25.A(18,20)=0.009059211150552;
reg25.A(18,21)=-0.000552598197073;
reg25.A(18,22)=0.011961280597280;
reg25.A(18,23)=0.044918913130228;
reg25.A(18,24)=0.011961280597291;
reg25.A(18,25)=-0.000552598197067;
reg25.A(19,2)=-0.000016656909212;
reg25.A(19,3)=-0.000155399981270;
reg25.A(19,4)=-0.000703921646039;
reg25.A(19,5)=-0.000644326564222;
reg25.A(19,6)=-0.000016656909212;
reg25.A(19,7)=-0.000053415881286;
reg25.A(19,8)=0.000425818412817;
reg25.A(19,9)=0.006457201243343;
reg25.A(19,10)=-0.000849569707349;
reg25.A(19,11)=-0.000155399981270;
reg25.A(19,12)=0.000425818412817;
reg25.A(19,13)=0.015750379019096;
reg25.A(19,14)=0.075869417154704;
reg25.A(19,15)=0.015877392263170;
reg25.A(19,16)=-0.000703921646039;
reg25.A(19,17)=0.006457201243344;
reg25.A(19,18)=0.075869417154712;
reg25.A(19,19)=0.161518496483849;
reg25.A(19,20)=0.035460566451184;
reg25.A(19,21)=-0.000644326564222;
reg25.A(19,22)=-0.000849569707349;
reg25.A(19,23)=0.015877392263169;
reg25.A(19,24)=0.035460566451180;
reg25.A(19,25)=0.014136310198355;
reg25.A(20,8)=-0.000565727610419;
reg25.A(20,9)=-0.041248478593548;
reg25.A(20,10)=-0.039527871369860;
reg25.A(20,13)=-0.002153251676374;
reg25.A(20,14)=-0.003160665040653;
reg25.A(20,15)=0.249340149438604;
reg25.A(20,17)=0.000000000000000;
reg25.A(20,18)=0.015748825270894;
reg25.A(20,19)=0.110576487688397;
reg25.A(20,20)=0.629648454319082;
reg25.A(20,23)=-0.013029845984057;
reg25.A(20,24)=-0.066167344054244;
reg25.A(20,25)=0.160539267612188;
reg25.A(21,6)=0.000921175499896;
reg25.A(21,7)=0.000929164605909;
reg25.A(21,11)=-0.025132646958045;
reg25.A(21,12)=0.006964792690618;
reg25.A(21,13)=0.012616833949940;
reg25.A(21,16)=0.081558557218738;
reg25.A(21,17)=-0.083901494197756;
reg25.A(21,18)=0.006964792690647;
reg25.A(21,19)=0.000929164605910;
reg25.A(21,21)=0.940802574133602;
reg25.A(21,22)=0.081558557218659;
reg25.A(21,23)=-0.025132646958018;
reg25.A(21,24)=0.000921175499900;
reg25.A(22,7)=0.000000000000000;
reg25.A(22,11)=-0.013029845984023;
reg25.A(22,12)=0.015748825270938;
reg25.A(22,13)=-0.002153251676373;
reg25.A(22,14)=-0.000565727610419;
reg25.A(22,16)=-0.066167344054254;
reg25.A(22,17)=0.110576487688357;
reg25.A(22,18)=-0.003160665040733;
reg25.A(22,19)=-0.041248478593576;
reg25.A(22,21)=0.160539267612127;
reg25.A(22,22)=0.629648454319149;
reg25.A(22,23)=0.249340149438669;
reg25.A(22,24)=-0.039527871369858;
reg25.A(23,12)=-0.001569843560556;
reg25.A(23,13)=0.003139687121018;
reg25.A(23,14)=-0.001569843560556;
reg25.A(23,16)=-0.007479499976368;
reg25.A(23,17)=-0.039918604317485;
reg25.A(23,18)=0.094796208587866;
reg25.A(23,19)=-0.039918604317479;
reg25.A(23,20)=-0.007479499976368;
reg25.A(23,21)=-0.023867882316688;
reg25.A(23,22)=0.166877977050279;
reg25.A(23,23)=0.713979810532753;
reg25.A(23,24)=0.166877977050271;
reg25.A(23,25)=-0.023867882316689;
reg25.A(24,12)=-0.000565727610419;
reg25.A(24,13)=-0.002153251676373;
reg25.A(24,14)=0.015748825270894;
reg25.A(24,15)=-0.013029845984057;
reg25.A(24,17)=-0.041248478593548;
reg25.A(24,18)=-0.003160665040653;
reg25.A(24,19)=0.110576487688394;
reg25.A(24,20)=-0.066167344054247;
reg25.A(24,21)=0.000000000000000;
reg25.A(24,22)=-0.039527871369860;
reg25.A(24,23)=0.249340149438604;
reg25.A(24,24)=0.629648454319082;
reg25.A(24,25)=0.160539267612190;
reg25.A(25,9)=0.000929164605902;
reg25.A(25,10)=0.000921175499883;
reg25.A(25,13)=0.012616833949832;
reg25.A(25,14)=0.006964792690417;
reg25.A(25,15)=-0.025132646958013;
reg25.A(25,17)=0.000929164605900;
reg25.A(25,18)=0.006964792690348;
reg25.A(25,19)=-0.083901494197416;
reg25.A(25,20)=0.081558557219214;
reg25.A(25,22)=0.000921175499873;
reg25.A(25,23)=-0.025132646958073;
reg25.A(25,24)=0.081558557219385;
reg25.A(25,25)=0.940802574132732;
reg25.B=spalloc(25, 25, 1);
reg25.B(7,1)=0.014139030952195;
reg25.B(7,2)=0.035460658116755;
reg25.B(7,3)=0.015878834231769;
reg25.B(7,4)=-0.000848759571398;
reg25.B(7,5)=-0.000644303435595;
reg25.B(7,6)=0.035460658116756;
reg25.B(7,7)=0.161511913632444;
reg25.B(7,8)=0.075867271838627;
reg25.B(7,9)=0.006458076149392;
reg25.B(7,10)=-0.000704016059195;
reg25.B(7,11)=0.015878834231769;
reg25.B(7,12)=0.075867271838628;
reg25.B(7,13)=0.015750717135235;
reg25.B(7,14)=0.000426550985590;
reg25.B(7,15)=-0.000155448526212;
reg25.B(7,16)=-0.000848759571398;
reg25.B(7,17)=0.006458076149392;
reg25.B(7,18)=0.000426550985590;
reg25.B(7,19)=-0.000053243238040;
reg25.B(7,20)=-0.000016667343519;
reg25.B(7,21)=-0.000644303435595;
reg25.B(7,22)=-0.000704016059195;
reg25.B(7,23)=-0.000155448526212;
reg25.B(7,24)=-0.000016667343519;
reg25.B(8,1)=-0.000551826681343;
reg25.B(8,2)=0.011962441663827;
reg25.B(8,3)=0.044920681581363;
reg25.B(8,4)=0.011962441663824;
reg25.B(8,5)=-0.000551826681343;
reg25.B(8,6)=0.009059867029122;
reg25.B(8,7)=0.053528285400930;
reg25.B(8,8)=0.186232759394393;
reg25.B(8,9)=0.053528285400922;
reg25.B(8,10)=0.009059867029123;
reg25.B(8,11)=0.002036332002843;
reg25.B(8,12)=0.016328443378969;
reg25.B(8,13)=0.075011777873706;
reg25.B(8,14)=0.016328443378964;
reg25.B(8,15)=0.002036332002846;
reg25.B(8,16)=-0.000517380391111;
reg25.B(8,17)=0.000872731042006;
reg25.B(8,18)=0.007454584686697;
reg25.B(8,19)=0.000872731042005;
reg25.B(8,20)=-0.000517380391109;
reg25.B(8,21)=-0.000157952222901;
reg25.B(8,22)=-0.000251178713936;
reg25.B(8,23)=-0.000635651455758;
reg25.B(8,24)=-0.000251178713936;
reg25.B(8,25)=-0.000157952222901;
reg25.B(9,1)=-0.000644303435595;
reg25.B(9,2)=-0.000848759571397;
reg25.B(9,3)=0.015878834231772;
reg25.B(9,4)=0.035460658116754;
reg25.B(9,5)=0.014139030952197;
reg25.B(9,6)=-0.000704016059195;
reg25.B(9,7)=0.006458076149393;
reg25.B(9,8)=0.075867271838632;
reg25.B(9,9)=0.161511913632437;
reg25.B(9,10)=0.035460658116757;
reg25.B(9,11)=-0.000155448526212;
reg25.B(9,12)=0.000426550985591;
reg25.B(9,13)=0.015750717135231;
reg25.B(9,14)=0.075867271838627;
reg25.B(9,15)=0.015878834231769;
reg25.B(9,16)=-0.000016667343519;
reg25.B(9,17)=-0.000053243238040;
reg25.B(9,18)=0.000426550985591;
reg25.B(9,19)=0.006458076149393;
reg25.B(9,20)=-0.000848759571398;
reg25.B(9,22)=-0.000016667343519;
reg25.B(9,23)=-0.000155448526212;
reg25.B(9,24)=-0.000704016059195;
reg25.B(9,25)=-0.000644303435595;
reg25.B(12,1)=-0.000551826681343;
reg25.B(12,2)=0.009059867029122;
reg25.B(12,3)=0.002036332002843;
reg25.B(12,4)=-0.000517380391111;
reg25.B(12,5)=-0.000157952222901;
reg25.B(12,6)=0.011962441663827;
reg25.B(12,7)=0.053528285400930;
reg25.B(12,8)=0.016328443378969;
reg25.B(12,9)=0.000872731042006;
reg25.B(12,10)=-0.000251178713936;
reg25.B(12,11)=0.044920681581363;
reg25.B(12,12)=0.186232759394393;
reg25.B(12,13)=0.075011777873706;
reg25.B(12,14)=0.007454584686697;
reg25.B(12,15)=-0.000635651455758;
reg25.B(12,16)=0.011962441663824;
reg25.B(12,17)=0.053528285400922;
reg25.B(12,18)=0.016328443378964;
reg25.B(12,19)=0.000872731042005;
reg25.B(12,20)=-0.000251178713936;
reg25.B(12,21)=-0.000551826681343;
reg25.B(12,22)=0.009059867029123;
reg25.B(12,23)=0.002036332002846;
reg25.B(12,24)=-0.000517380391109;
reg25.B(12,25)=-0.000157952222901;
reg25.B(13,1)=-0.000483089914478;
reg25.B(13,2)=0.000081580711140;
reg25.B(13,3)=0.007372895067013;
reg25.B(13,4)=0.000081580711141;
reg25.B(13,5)=-0.000483089914479;
reg25.B(13,6)=0.000081580711140;
reg25.B(13,7)=0.011308388978283;
reg25.B(13,8)=0.058321481761030;
reg25.B(13,9)=0.011308388978278;
reg25.B(13,10)=0.000081580711140;
reg25.B(13,11)=0.007372895067013;
reg25.B(13,12)=0.058321481761027;
reg25.B(13,13)=0.181751608293496;
reg25.B(13,14)=0.058321481761022;
reg25.B(13,15)=0.007372895067015;
reg25.B(13,16)=0.000081580711141;
reg25.B(13,17)=0.011308388978278;
reg25.B(13,18)=0.058321481761022;
reg25.B(13,19)=0.011308388978282;
reg25.B(13,20)=0.000081580711140;
reg25.B(13,21)=-0.000483089914479;
reg25.B(13,22)=0.000081580711140;
reg25.B(13,23)=0.007372895067014;
reg25.B(13,24)=0.000081580711142;
reg25.B(13,25)=-0.000483089914479;
reg25.B(14,1)=-0.000157952222901;
reg25.B(14,2)=-0.000517380391111;
reg25.B(14,3)=0.002036332002845;
reg25.B(14,4)=0.009059867029120;
reg25.B(14,5)=-0.000551826681342;
reg25.B(14,6)=-0.000251178713936;
reg25.B(14,7)=0.000872731042008;
reg25.B(14,8)=0.016328443379001;
reg25.B(14,9)=0.053528285400908;
reg25.B(14,10)=0.011962441663827;
reg25.B(14,11)=-0.000635651455758;
reg25.B(14,12)=0.007454584686702;
reg25.B(14,13)=0.075011777873720;
reg25.B(14,14)=0.186232759394332;
reg25.B(14,15)=0.044920681581382;
reg25.B(14,16)=-0.000251178713936;
reg25.B(14,17)=0.000872731042005;
reg25.B(14,18)=0.016328443378964;
reg25.B(14,19)=0.053528285400923;
reg25.B(14,20)=0.011962441663838;
reg25.B(14,21)=-0.000157952222901;
reg25.B(14,22)=-0.000517380391111;
reg25.B(14,23)=0.002036332002846;
reg25.B(14,24)=0.009059867029115;
reg25.B(14,25)=-0.000551826681336;
reg25.B(17,1)=-0.000644303435595;
reg25.B(17,2)=-0.000704016059195;
reg25.B(17,3)=-0.000155448526212;
reg25.B(17,4)=-0.000016667343519;
reg25.B(17,6)=-0.000848759571397;
reg25.B(17,7)=0.006458076149393;
reg25.B(17,8)=0.000426550985591;
reg25.B(17,9)=-0.000053243238040;
reg25.B(17,10)=-0.000016667343519;
reg25.B(17,11)=0.015878834231772;
reg25.B(17,12)=0.075867271838632;
reg25.B(17,13)=0.015750717135231;
reg25.B(17,14)=0.000426550985591;
reg25.B(17,15)=-0.000155448526212;
reg25.B(17,16)=0.035460658116754;
reg25.B(17,17)=0.161511913632438;
reg25.B(17,18)=0.075867271838627;
reg25.B(17,19)=0.006458076149393;
reg25.B(17,20)=-0.000704016059195;
reg25.B(17,21)=0.014139030952196;
reg25.B(17,22)=0.035460658116757;
reg25.B(17,23)=0.015878834231769;
reg25.B(17,24)=-0.000848759571398;
reg25.B(17,25)=-0.000644303435595;
reg25.B(18,1)=-0.000157952222901;
reg25.B(18,2)=-0.000251178713936;
reg25.B(18,3)=-0.000635651455758;
reg25.B(18,4)=-0.000251178713936;
reg25.B(18,5)=-0.000157952222901;
reg25.B(18,6)=-0.000517380391111;
reg25.B(18,7)=0.000872731042008;
reg25.B(18,8)=0.007454584686702;
reg25.B(18,9)=0.000872731042005;
reg25.B(18,10)=-0.000517380391111;
reg25.B(18,11)=0.002036332002845;
reg25.B(18,12)=0.016328443379001;
reg25.B(18,13)=0.075011777873720;
reg25.B(18,14)=0.016328443378964;
reg25.B(18,15)=0.002036332002846;
reg25.B(18,16)=0.009059867029120;
reg25.B(18,17)=0.053528285400908;
reg25.B(18,18)=0.186232759394332;
reg25.B(18,19)=0.053528285400923;
reg25.B(18,20)=0.009059867029115;
reg25.B(18,21)=-0.000551826681342;
reg25.B(18,22)=0.011962441663827;
reg25.B(18,23)=0.044920681581382;
reg25.B(18,24)=0.011962441663838;
reg25.B(18,25)=-0.000551826681336;
reg25.B(19,2)=-0.000016667343519;
reg25.B(19,3)=-0.000155448526212;
reg25.B(19,4)=-0.000704016059195;
reg25.B(19,5)=-0.000644303435595;
reg25.B(19,6)=-0.000016667343519;
reg25.B(19,7)=-0.000053243238040;
reg25.B(19,8)=0.000426550985592;
reg25.B(19,9)=0.006458076149395;
reg25.B(19,10)=-0.000848759571399;
reg25.B(19,11)=-0.000155448526212;
reg25.B(19,12)=0.000426550985593;
reg25.B(19,13)=0.015750717135237;
reg25.B(19,14)=0.075867271838625;
reg25.B(19,15)=0.015878834231764;
reg25.B(19,16)=-0.000704016059195;
reg25.B(19,17)=0.006458076149395;
reg25.B(19,18)=0.075867271838634;
reg25.B(19,19)=0.161511913632443;
reg25.B(19,20)=0.035460658116749;
reg25.B(19,21)=-0.000644303435595;
reg25.B(19,22)=-0.000848759571399;
reg25.B(19,23)=0.015878834231763;
reg25.B(19,24)=0.035460658116744;
reg25.B(19,25)=0.014139030952213;
reg25.fstar=[];
reg25.fstar(1)=0.000000000000000;
reg25.fstar(2)=0.000000000000000;
reg25.fstar(3)=0.000000000000000;
reg25.fstar(4)=0.000000000000000;
reg25.fstar(5)=0.000000000000000;
reg25.fstar(6)=0.000000000000000;
reg25.fstar(7)=0.000000000000000;
reg25.fstar(8)=0.000000000000000;
reg25.fstar(9)=0.000000000000000;
reg25.fstar(10)=0.000000000000000;
reg25.fstar(11)=0.000000000000000;
reg25.fstar(12)=0.000000000000000;
reg25.fstar(13)=0.000000000000000;
reg25.fstar(14)=0.000000000000000;
reg25.fstar(15)=0.000000000000000;
reg25.fstar(16)=0.000000000000000;
reg25.fstar(17)=0.000000000000000;
reg25.fstar(18)=0.000000000000000;
reg25.fstar(19)=0.000000000000000;
reg25.fstar(20)=0.000000000000000;
reg25.fstar(21)=0.000000000000000;
reg25.fstar(22)=0.000000000000000;
reg25.fstar(23)=0.000000000000000;
reg25.fstar(24)=0.000000000000000;
reg25.fstar(25)=0.000000000000000;