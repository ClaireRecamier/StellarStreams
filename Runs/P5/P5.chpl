//import streakline;
use H5Helper;
var Dom = {1..0,1..0,1..0};
var arr : [Dom] real;
var ff = openH5File("phi_000000.h5");
readRealArray(ff,'phi',Dom,arr);
writeln(Dom);
