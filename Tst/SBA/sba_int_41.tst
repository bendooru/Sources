LIB "tst.lib";
tst_init();

ring r = (integer),(x,y,z),dp;
ideal i = 22xy+94xz,74x2,83yz;
ideal stdi = std(i);
ideal sbai = sba(i,1,0);
std(reduce(stdi,sbai));
std(reduce(sbai,stdi));


tst_status(1);
$
