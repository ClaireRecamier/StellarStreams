use IO;
use Math;
use Random;
use HDF5;
use H5Helper;
use /Code/streakline;
//DIMENSIONLESS VERSION
//math constants
const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kmkpc = 3.24078e-17;
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const yrsec: real = 60*60*24*365.24; //seconds to years
const secyr: real = 1 / yrsec;
const kpckm = 1.0 / kmkpc;
const G = 6.67430e-11; //m^3 / (kg sec^2)
const G_code = 6.67430e-11 * ((1e-3)**3) * (yrsec)**2 * kmkpc**3 * (1/toCodeTime(1))**2  * toCodeLength(1)**3 * (1/toCodeMass(1)); //G in km^3 / (kg * sec^2)
const hbar = 1.0545718e-34; //((m^2)*(kg/s))
const parsec = 3.08571e16 ; // m
const H0 =(67.7e3)/(parsec*(1e6)) ; // 1/sec
const omegaM0 = 0.31 ;
const maxion = 1e-22 * 1.783 * 1e-36;

//setup constants
config const r0 = toCodeLength(15.0); // radius of ring, in km
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const Ne = 1; //number of particles in ring
config const nbins = 4;
var period: real = 0.0;
/* config const po = 27492.260803351877;
config const rc = 0.05359134269304645; */
config const po = (0.5**4) * 27492.260803351877;
config const rc = 2.0*0.05359134269304645;
var pot = 6; //set galactic potential to that of imported soliton
config const integrator = 0;   //set Integrator to LF
config const N = 100000;//6000 set number of timesteps to equal total of 6 Gyr
var dt = toCodeTime(1000000.0); //set timestep to seconds per Myr


var Dom = {1..0,1..0,1..0};
//var Dom = {{1..0,1..0,1..0},{1..0,1..0,1..0},{1..0,1..0,1..0}};
var arr : [Dom] real;


proc main(){

  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real; //stores angular momentum of each particle, updated at each timestep
  var SD: [0..N] real; //stores std dev at each timestep
  var PS: [1..nbins] real; //whole array gets updated at each timestep
  //hardcode galactic potential parameters
  var calcpar: [0..5] real;

  var chfile = open("bp1.csv",iomode.cw); //create test.csv and open
  var SDfile = open("SD19.csv",iomode.cw); //create test.csv and open
  var AMfile = open("AM19.csv",iomode.cw); //create test.csv and open
  var PSfile = open("PS19.csv",iomode.cw); //create test.csv and open
  var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
  var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
  var WritingChannel = chfile.writer(); //open writing channel to test.csv
  var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv

  load_pot(calcpar,pot);
  //set initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar, WritingChannel,AMWritingChannel,PSWritingChannel);
  //dt = dt / period;
  writeln("initialized ring");
  writeln("pos 1 ",pos[1]);
  writeln("vel 1", vel[1]);
  //
  /* writeln("calling search force ",search_force((0.5,0.5,0.5))); */
  //writeln("calling search force ",search_force((0.4,0.4,0.4)));

  fwd_orbit_ring(pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel,SDWritingChannel,AMWritingChannel,PSWritingChannel);
  WritingChannel.close();
  SDWritingChannel.close();
  AMWritingChannel.close();
  PSWritingChannel.close();
  PSfile.close();
  AMfile.close();
  SDfile.close();
  chfile.close();
}
