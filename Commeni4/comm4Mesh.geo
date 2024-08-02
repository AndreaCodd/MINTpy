//COMMENI 4

xL = -30000.;
x1 =  -6000.;
x2 =   2000.;
x3 =   5000.;
xR =  30000.;

yB = -50000.;
y1 = -25000.;
y2 =  -4000.;
y3 =  -3000.;
y4 =  -1000.;
y5 =   -500.;
yG =      0.;
yA =  50000.;

elB = 5000.;
el1 = 2000.;
el2 =   50.;
el3 =  100.;
el4 =  100.;
elG =  100.; 
elA = 10000.; 
elS =   20.;


Point(1)  = {xL, yB, 0., elB};
Point(2)  = {xR, yB, 0., elB};
Point(3)  = {xL, y1, 0., el1};
Point(4)  = {xR, y1, 0., el1};
Point(5)  = {xL, y3, 0., el3};
Point(6)  = {x1, y3, 0., el2};
Point(7)  = {x1, y2, 0., el2};
Point(8)  = {x2, y2, 0., el2};
Point(9)  = {x3, y4, 0., el2};
Point(10) = {xR, y4, 0., el3};
Point(11) = {xL, y5, 0., el4};
Point(12) = {x1, y5, 0., el2};
Point(13) = {xR, y5, 0., el4};
Point(14) = {xL, yG, 0., elG};
Point(15) = {xR, yG, 0., elG};
Point(16) = {xL, yA, 0., elA};
Point(17) = {xR, yA, 0., elA};

// sensors
Point(20) = {-20000., 0., 0., elS};
Point(21) = {-10000., 0., 0., elS};
Point(22) = { -7000., 0., 0., elS};
Point(23) = { -6000., 0., 0., elS};
Point(24) = { -5000., 0., 0., elS};
Point(25) = {     0., 0., 0., elS};
Point(26) = {  2000., 0., 0., elS};
Point(27) = {  5000., 0., 0., elS};
Point(28) = {  8000., 0., 0., elS};
Point(29) = { 16000., 0., 0., elS};

// across first
Line(1)  = { 1,  2};
Line(2)  = { 3,  4};
Line(3)  = { 5,  6};
Line(4)  = { 7,  8};
Line(5)  = { 8,  9};
Line(6)  = { 9, 10};
Line(7)  = {11, 12};
Line(8)  = {12, 13};

// ground level and sensors
Line(9)  = {14, 20};
Line(10) = {16, 17};
Line(11) = { 1,  3};
Line(12) = { 3,  5};
Line(13) = { 5, 11};
Line(14) = {11, 14};
Line(15) = {14, 16};
Line(16) = { 7,  6};
Line(17) = { 6, 12};
Line(18) = { 2,  4};
Line(19) = { 4, 10};
Line(20) = {10, 13};
Line(21) = {13, 15};
Line(22) = {15, 17};
Line(25) = {20, 21};
Line(26) = {21, 22};
Line(27) = {22, 23};
Line(28) = {23, 24};
Line(29) = {24, 25};
Line(30) = {25, 26};
Line(31) = {26, 27};
Line(32) = {27, 28};
Line(33) = {28, 29};
Line(34) = {29, 15};

Line Loop(1) = {1, 18, -2, -11};
Plane Surface(1) = {1};
Physical Surface("Base")= {1};
Line Loop(2) = {2, 19, -6, -5, -4, 16, -3, -12};
Plane Surface(2) = {2};
Physical Surface("Middle") = {2};
Line Loop(3) = {3, 17, -7, -13};
Plane Surface(3) = {3};
Physical Surface("Blob1")= {3};
Line Loop(4) = {4, 5, 6, 20, -8, -17, -16};
Plane Surface(4) = {4};
Physical Surface("Blob2")= {4};
Line Loop(5) = {7, 8, 21, -34:-25, -9, -14};
Plane Surface (5) = {5};
Physical Surface("Top")= {5};
Line Loop(6) = {9, 25:34, 22,-10, -15};
Plane Surface (6) = {6};
Physical Surface("air")= {6};































 











