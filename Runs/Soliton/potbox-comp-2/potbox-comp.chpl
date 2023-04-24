use streakline;
use IO;
use Math;
use Random;
use HDF5;
use H5Helper;

/* This module compares the csvbox created from an analytic potential
in the soliton-imp-csvbox-1 run, to the hdf5 potential. To do so,
it converts the hdf5 box values of potential to acceleration, then
outputs the difference of hdf5 acceleration and analytic
acceleration to a file. */


//setup constants
/* config const r0 = toCodeLength(50.0); // radius of ring - kpc converted to codeunits */
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const Ne = 1; //number of particles in ring
config const nbins = 4;
var period: real = 0.0;
config const po = (0.5**4) * 27492.260803351877;
config const rc = 2.0*0.05359134269304645;
config const integrator = 0;   //set Integrator to LF
config const N = 100000;//6000 set number of timesteps to equal total of 6 Gyr
var dt = toCodeTime(1000000.0); //set timestep to seconds per Myr
config const grid_size: real = 127.0;
config const Lbox: real = 2.0;
const G_code = 6.67430e-11 * ((1e-3)**3) * (yrsec)**2 * kmkpc**3 * (1/toCodeTime(1))**2  * toCodeLength(1)**3 * (1/toCodeMass(1)); //G in km^3 / (kg * sec^2)


proc pot_to_acc (x,y,z, pot_array){
  var acc: 3*real;
  var a: real;
  var h: real = 2.0 * (Lbox / grid_size); //codeunits per gridpoint
  var x_pos: int = x;
  var y_pos: int = y;
  var z_pos: int = z;


  if x > 126 {
    x_pos = 126;
  }
  else if x < 2 {
    x_pos = 2;
  }
  if y > 126 {
    y_pos = 126;
  }
  else if y < 2{
    y_pos = 2;
  }
  if z > 126 {
    z_pos = 126;
  }
  else if z < 2{
    z_pos = 2;
  }
  //acceleration is gradient divided by mass, total mass of which equals 25
  /* acc(0) = -1*(pot_array[x_pos + 1, y_pos,z_pos] - pot_array[x_pos-1,y_pos,z_pos])/ (h*25.0);
  acc(1) = -1*(pot_array[x_pos, y_pos + 1,z_pos] - pot_array[x_pos,y_pos-1,z_pos])/ (h*25.0);
  acc(2) = -1*(pot_array[x_pos, y_pos,z_pos + 1] - pot_array[x_pos,y_pos,z_pos-1])/ (h*25.0); */

  acc(0) = -1*(pot_array[x_pos + 1, y_pos,z_pos] - pot_array[x_pos-1,y_pos,z_pos])/ (h);
  acc(1) = -1*(pot_array[x_pos, y_pos + 1,z_pos] - pot_array[x_pos,y_pos-1,z_pos])/ (h);
  acc(2) = -1*(pot_array[x_pos, y_pos,z_pos + 1] - pot_array[x_pos,y_pos,z_pos-1])/ (h);
  a = len(acc);
  /* writeln(pot_array[39,63,64],"\n"); */
  return a;
}

proc plot_accelerations(pot_array){
  //setup writing to csv file - z=0 plane
  var plot_acc_file_1 = open("plot-acc1.csv",iomode.cw);
  var plot_acc_1 = plot_acc_file_1.writer();
  //open second file for x = 0 plane
  var plot_acc_file_2 = open("plot-acc2.csv",iomode.cw);
  var plot_acc_2 = plot_acc_file_2.writer();
  //setup reading from imported analytic soliton
  var nfwbox = open("../soliton-imp-csvbox-1/box-pot.csv",iomode.r);
  var boxchannel = nfwbox.reader();
  //iterate through csv box of accelerations
  for i in 0..(grid_size: int){
    for j in 0..(grid_size: int){
      for k in 0..(grid_size: int){
        //read in acceleration at this point
        var tmp: string;
        boxchannel.read(tmp);
        var acceleration = tmp: real;
        //convert location in box to cartesian coordinates
        var dr: real = Lbox/grid_size; //codeunits per gridpoint
        var x: real = (i - (grid_size / 2)) * dr;
        var y: real = (j - (grid_size / 2)) * dr;
        var z: real = (k - (grid_size / 2)) * dr;
        //convert cartesian coordinates
        var r: real = sqrt(x**2 + y**2 + z**2);
        /* var theta: real = acos(z/r);
        var phi: real = atan(y/x);
        //only collect data on x axis
        if (theta == Pi/2 && phi == 0){
          plot_acc.write(r, ',', acceleration,'\n');
        } */
        var theta: real = atan(y/x);
        if (k == 63){
          plot_acc_1.write(r, ',', acceleration,',');
          plot_acc_1.write(r, ',', pot_to_acc (i,j,k, pot_array),'\n');
        }
        if (i == 63){
          plot_acc_2.write(r, ',', acceleration,',');
          plot_acc_2.write(r, ',', pot_to_acc (i,j,k, pot_array),'\n');
        }

      }
    }
  }
  plot_acc_1.close();
  plot_acc_file_1.close();
  plot_acc_2.close();
  plot_acc_file_2.close();

  boxchannel.close();
  nfwbox.close();

}
proc compare_accelerations(pot_array,impanly1,comp,hdf5_box_pot){
  var impanly1_pot: real;
  var acc: real;
  var anly_acc: real;
  var tmp: string;
  var tmp1: 3*real;

  //iterate through box
  for i in 0..(grid_size: int){
    for j in 0..(grid_size: int){
      for k in 0..(grid_size: int){
        //convert potential at this point from hdf5 file to acc
        acc = pot_to_acc(i,j,k,pot_array);
        //obtain acc at this point from csv file
        impanly1.read(tmp);
        anly_acc = tmp: real;
        //write difference of accelerations to comparison file
        comp.write(acc - anly_acc, ",");
        //write hdf5 potential converted to acceleration to file
        hdf5_box_pot.write(acc, ",");
        //print out this specific range of points to terminal
        if (i < 56) && (i > 50){
          if j == 60 {
            if k == 60 {
              writeln("from hdf5 box ", acc);
              writeln("from analytic box ",anly_acc);
            }
          }
        }
      }
      comp.write("\n");
      hdf5_box_pot.write("\n");
    }
  }
}

proc plot_masses(calcpar,pot_array){
  //setup writing of analytic mass to csv file
  var plot_mass_file = open("plot-mass.csv",iomode.cw);
  var plot_mass = plot_mass_file.writer();
  //setup writing of hdf5 mass to csv file
  var plot_mass_file_2 = open("plot-mass2.csv",iomode.cw);
  var plot_mass_2 = plot_mass_file_2.writer();
  //plot pot calculated from analytic soliton
  var plot_pot_file = open("plot-pot.csv",iomode.cw);
  var plot_pot = plot_pot_file.writer();

  var plot_pot_file_2 = open("plot-pot2.csv",iomode.cw);
  var plot_pot_2 = plot_pot_file_2.writer();
  //calc analytic mass along radius
  var po: real = calcpar[0];
  var rc: real = calcpar[1];
  var r: real = rc * 10.0;
  var height: real = 0.0;
  var dr: real = 0.001; //changed from 0.1 which wasn't working
  var mass: real = 0.0;
  var bound: int = (r/dr): int;
  var a: real = 0.0;
  for i in 0..bound {
    const r1 = i * dr;
    height = po * 4 * Pi * (r1**2) / (1 + 0.091 * (i * dr/rc)**2)**8; //
    //height = po / (1 + 0.091 * (i * dr/rc)**2)**8; //
    mass = mass + (height * dr);
    //a += -G * height * dr / (i*dr)**2;
    plot_mass.write(r1, ',', mass,'\n');
    plot_pot.write(r1, ',', -1.0 * G_code * mass / r1 ,'\n');
  }
  //calc hdf5 mass along radius
  writeln(pot_array[64,64,64]);
  writeln(pot_array[63,63,63]);
  mass = 0.0;
  for i in 0..(grid_size: int) {
    for j in 0..(grid_size: int) {
      for k in 0..(grid_size: int){
        //convert ijk to xyz and get r
        var dr: real = Lbox / grid_size; //codeunits per gridpoint
        var x: real = (i - (grid_size / 2)) * dr;
        var y: real = (j - (grid_size / 2)) * dr;
        var z: real = (k - (grid_size / 2)) * dr;
        var r: real = sqrt(x**2 + y**2 + z**2);
        //convert potential to mass
        var pot: real = pot_array[i,j,k];
        var mass: real = (-1.0 * pot * r) / G_code;
        //fix mass calculation
        //mass = (-1.0 * pot * r);
        plot_mass_2.write(r, ',', mass,'\n');
        plot_pot_2.write(r, ',', pot,'\n');
      }
    }


  }



  plot_mass.close();
  plot_mass_file.close();
  plot_mass_2.close();
  plot_mass_file_2.close();
  plot_pot.close();
  plot_pot_file.close();
  plot_pot_2.close();
  plot_pot_file_2.close();

}

proc main(){
  //open files with analytic soliton acceleration box
  var impanly1file = open("../soliton-imp-csvbox-1/box-pot.csv",iomode.r); //create test.csv and open
  var impanly1 = impanly1file.reader(); //open writing channel to test.csv
  //open file to write comparison of accelerations to
  var compfile = open("comparison.csv",iomode.cw); //create test.csv and open
  var comp = compfile.writer(); //open writing channel to test.csv
  //open file to write in hdf5 potential converted to acceleration
  var impfile = open("hdf5-box-pot.csv",iomode.cw);
  var hdf5_box_pot = impfile.writer(); //open writing channel to test.csv

  //load hdf5 potential into pot_array
  var Dom = {1..0,1..0,1..0};
  var pot_array : [Dom] real;
  var imp1_file = openH5File("../../../Code/phi_000000.h5");
  readRealArray(imp1_file,'phi',Dom,pot_array);

  var calcpar: [0..1] real = [po, rc];

  compare_accelerations(pot_array,impanly1,comp,hdf5_box_pot);
  plot_accelerations(pot_array);
  plot_masses(calcpar,pot_array);

  impanly1.close();
  comp.close();
  hdf5_box_pot.close();
  impanly1file.close();
  compfile.close();
  impfile.close();
}
