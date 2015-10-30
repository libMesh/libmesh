lc = 0.65;
lco = 0.3*lc;
lc1 = 0.75*lc;

d = 5; 
l = 5; 
b = 3; 
a = 0.5; 
dgap = 0.001; 
R0 = 12; R1 = Sqrt(R0*R0+(a+dgap)*(a+dgap));

Point(10000) = {R0,0.,0.,lc};
Point(1) = {0.,a,0.,lco};
Point(2) = {d,a,0.,lco};
Point(3) = {d+l,b,0.,lco};
Point(4) = {d,a+dgap,0.,lc};
Point(5) = {0,a+dgap,0.,lc};
Point(6) = {R0,R1,0.,lc};
Point(7) = {R0+R1,0,0.,lc1};
Point(8) = {R0,-R1,0.,lc};
Point(9) = {0,-(a+dgap),0.,lc};
Point(10) = {d,-(a+dgap),0.,lc};
Point(11) = {d+l,-b,0.,lco};
Point(12) = {d,-a,0.,lco};
Point(13) = {0,-a,0.,lco};
Point(14) = {d,0,0.,lco};
Point(15) = {d+l,0,0.,lco};
Point(16) = {d+l,b+a,0.,lc};
Point(17) = {d,b+a,0.,lc};
Point(18) = {d+l,-(b+a),0.,lc};
Point(19) = {d,-(b+a),0.,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Circle(5) = {5,10000,6};
Circle(6) = {6,10000,7};
Circle(7) = {7,10000,8};
Circle(8) = {8,10000,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,1};

Line(14) = {2,14};
Line(15) = {14,12};

Line(16) = {10,19};
Line(17) = {19,18};
Line(18) = {18,11};
Line(19) = {4,17};
Line(20) = {17,16};
Line(21) = {16,3};
Line(22) = {14,15};
Line(23) = {15,3};
Line(24) = {15,11};

Line Loop(1) = {1,14,15,12,13};
Line Loop(2) = {4,5,6,7,8,9,16,17,18,-24,23,-21,-20,-19};
Line Loop(3) = {-14,-22,-23,2};
Line Loop(4) = {22,24,11,-15};
Line Loop(5) = {19,20,21,3};
Line Loop(6) = {-16,10,-18,-17};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Physical Line(1) = {13};
Physical Line(2) = {5,6,7,8};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
