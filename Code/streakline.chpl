//uses chplultra units
use IO;
use Math;
use Random;
use HDF5;
use H5Helper;

//math constants
const Pi = 3.141592653589793;
const Msun = 4 * Pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kmkpc = 3.24078e-17;
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const yrsec: real = 60*60*24*365.24; //seconds to years
const secyr: real = 1 / yrsec;
const kpckm = 1.0 / kmkpc;

//physics constants
const G = 6.67430e-11; //m^3 / (kg sec^2)
const G_code = 6.67430e-11 * ((1e-3)**3) * (yrsec)**2 * kmkpc**3 * (1/toCodeTime(1))**2  * toCodeLength(1)**3 * (1/toCodeMass(1)); //G in km^3 / (kg * sec^2)
const hbar = 1.0545718e-34; //((m^2)*(kg/s))
const parsec = 3.08571e16 ; // m
const H0 =(67.7e3)/(parsec*(1e6)) ; // 1/sec
const omegaM0 = 0.31;
const maxion = 1e-22 * 1.783 * 1e-36;

//setup constants
/* config const r0 = toCodeLength(15.0); // radius of ring, in km
config const integrator = 0;   //set Integrator to LF
config const N = 100000;//number of timesteps ~ total of 6 Gyr
var dt = toCodeTime(1000000.0); //set timestep to seconds per Myr
config const mcli: real = 20000 * Msun; //initial mass of cluster
config const mclf: real = 20000 * Msun; //final mass of cluster
config const M = 1; //particles are released every Mth timestep
config const Rcl = 20 * 0.001 * kpcau; //radius of plummer core
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const Ne = 1; //number of particles in ring THIS IS DIFF in DIFF VERSION
config const nbins = 4; //power spectra setup
var pot = 6; //set galactic potential: 5 = analytic soliton,6 = imported potential
var period: real = 0.0;
var k: int = 0; //record how many particles have been released
var dm = (mcli - mclf)/N; //amount of mass released per timesteps
var mcl: real = mcli; //current mass of cluster */
/* var Dom = {1..0,1..0,1..0}; //domain of box potential array
var pot_array : [Dom] real; //box potential array */
/* var offset: [0..1] real = [0.2,0.2]; */
var randStream = new RandomStream(real);


proc pot_box(Lbox, dt, calcpar,pot, grid_size,pot_array,box_param) { //creates acceleration box in km/seconds squared
  //var nfwbox = open("../TestingHDF5/mathematica/nfwbox1.csv",iomode.cw); //create test.csv and open
  var box = open("box-pot.csv",iomode.cw); //create test.csv and open
  var boxchannel = box.writer(); //open writing channel to test.csv

  var center_offset: 3*real = (box_param[0],box_param[1],box_param[2]);
  var grid_size: real = box_param[3];
  var Lbox: real = box_param[4];
  //var pot_array = box_param[5];

  /* var par: [0..5] real = [417.0, 36.54 * kpckm, 90.0 * Pi / 180, 1.0, 1.0 ,1.0];
  if pot == 3 { //if using triaxial NFW potential
    //assuming par = [V, rhalo, phi, q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM in km cubed / seconds squared
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1];
    var acc: 3*real;
    //acc = force(r0 * pos[j],pot,calcpar); //convert pos to dimensions and get force
    writeln("GM in km cubed / seconds squared ", calcpar[0]);
    acc = force((r0,0.0,0.0),pot,calcpar);
    writeln("centripetal acceleration in km/sec^2 ",acc[0]," ",acc[1]," ",acc[2]);
    calcpar[5] = par[1] / r0; //Rhalo becomes dimensionless
    period = 2 * Pi * r0 / sqrt(len(acc) * r0); //in seconds; 2Pi * r/v

    //period = 2 * Pi * r0 / magVel; //in seconds; 2Pi * r/v
    calcpar[0] = calcpar[0] * (period ** 2)/ (r0 ** 3); //GM becomes dimensionless
    dt = dt / period; //divide dt by period to make it dimensionless
    writeln("period ", period);
    writeln("dt ", dt);

    }
  else if pot == 0 {
    calcpar[0] = G * 1000000000.0 * (1.989e30); //GM in km^3 / sec^2
  } */
  var dr: real = Lbox / grid_size; //codeUnits per gridpoint
  var x, y, z: real = 0.0;
  /* var x, y, z: real = 0.0; */
  var acc: 3*real;
  /* for i in -63..63 {//x */
  for i in -63..63 {//x
    x = i * dr; //position in codeunits
    /* for j in -64..63 {//y */
    for j in -63..63 {//y
      y = j * dr;
      /* for k in -64..63 {//z */
      for k in -63..63 {//z
        z = k * dr; //in codeunits
        /* acc = force((x,y,z),pot,calcpar,box_param); */
        acc = force((x,y,z),pot,calcpar,pot_array,box_param);
        boxchannel.write(len(acc)," ");
      }
      boxchannel.write("\n");
    }
  }
  boxchannel.close();
  box.close();
}

proc load_anly_pot(r0,par,calcpar,pot,ref dt) { //creates acceleration box in km/seconds squared
  if pot == 3 {
    //triaxial NFW potential
    //assuming par = [V (km/sec), rhalo (km), phi (radians), q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    //converts calcpar to dimensionless units
    //(only deals with GM and rhalo - wb c1 - c4???)
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM in km cubed/seconds squared
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1]; //rhalo in km

    /* this is no longer needed if using chpl ultra units?
    //convert calcpar to dimensionless units
    var acc: 3*real;
    //acc = force(r0 * pos[j],pot,calcpar); //convert pos to dimensions and get force
    acc = force((r0,0.0,0.0),pot,calcpar); //obtain force at r0 in units
    calcpar[5] = par[1] / r0; //Rhalo becomes dimensionless
    period = 2 * Pi * r0 / sqrt(len(acc) * r0); //in seconds; 2Pi * r/v
    calcpar[0] = calcpar[0] * (period ** 2)/ (r0 ** 3); //GM becomes dimensionless
    dt = dt / period; //divide dt by period to make it dimensionless */

    }
  else if pot == 0 {
    //using point potential (Newton)
    calcpar[0] = G_code * 1000000000.0 * (1.989e30); //GM in km^3 / sec^2
  }
  /* else if pot == 5 {
    calcpar[0] = par[0];
    calcpar[1] = par[1];
  } */
}
proc load_imp_pot(pot,ref Dom, ref pot_array, pot_file) {
  //loads the specified file into the specified array
  if pot == 6 {
    //importing potential in pot_array
    readRealArray(pot_file,'phi',Dom,pot_array);
  }
}

proc force(pos,pot,calcpar,pot_array,box_param){
//proc force(pos,pot,calcpar,box_param){
  var acc: 3*real;
  var r: real;
  if pot == 0 { //if using point mass potential
    r = len(pos);
    //acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
    //acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r**3;
    //writeln("acc ",acc);
    acc = (-1,-1,-1)*(calcpar[0] * pos)/r**3;

  }
  else if pot == 3 { //NFW triaxial potential, C version
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    r = sqrt(calcpar[1]*pos[0]*pos[0] + calcpar[2]*pos[1]*pos[1] + calcpar[3]*pos[0]*pos[1] + calcpar[4]*pos[2]*pos[2]);
    var aux: real = 0.5 * calcpar[0] / (r**3) * (1.0/(1.0 + calcpar[5]/r)-log(1.0+r/calcpar[5]));

    acc[0]=aux*(2*calcpar[1]*pos[0] + calcpar[3]*pos[1]);
    acc[1]=aux*(2*calcpar[2]*pos[1] + calcpar[3]*pos[0]);
    acc[2]=aux*(2*calcpar[4]*pos[2]);
    //writeln("acc magnitude ", len(acc));
  }
  else if pot == 5 { //soliton approximation
    var po: real = calcpar[0];
    var rc: real = calcpar[1];

    var center_offset: 3*real = (box_param[0],box_param[1],box_param[2]);

    var shifted_pos: 3*real;
    shifted_pos[0] = pos[0] + center_offset[0];
    shifted_pos[1] = pos[1] + center_offset[1];
    shifted_pos[2] = pos[2] + center_offset[2];
    r = len(shifted_pos);
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
    }
    //writeln(mass);
    acc = (- G_code * mass / r**2) * (pos/r);
    //acc = a * (pos/r);
  }
  else if pot == 6 { //search hdf5 box
    var center_offset: 3*real = (box_param[0],box_param[1],box_param[2]);
    var grid_size: real = box_param[3];
    var Lbox: real = box_param[4];
    var dr: real = grid_size / Lbox;
    //HDF5 box goes from 0 to 127
    var x_pos: int = floor(abs(-Lbox/2.0 + center_offset[0] - pos[0]) * dr): int;
    var y_pos: int = floor(abs(-Lbox/2.0 + center_offset[1] - pos[1]) * dr): int;
    var z_pos: int = floor(abs(-Lbox/2.0 + center_offset[2] - pos[2]) * dr): int;
    //writeln("x pos ",abs(-1.0 - pos[0]) * dr-0.25," y pos ",abs(-1.0 - pos[1]) * dr-0.25," z pos ",abs(-1.0 - pos[2]) * dr-0.25);
    /* writeln("x pos ",x_pos," y pos ",y_pos," z pos ",z_pos); */
    //var a: [0..2] real;
    var a: real;
    a = pot_array[x_pos,y_pos,z_pos];
    /* writeln("a ",a); */
    //take the derivative of a
    var h: real = 2.0 * (Lbox / grid_size);

    //derivative blows up at edges, check for edge
    if x_pos > 126 {
      x_pos = 126;
    }
    if y_pos > 126 {
      y_pos = 126;
    }
    if z_pos > 126 {
      z_pos = 126;
    }

    acc(0) = -1*(pot_array[x_pos + 1, y_pos,z_pos] - pot_array[x_pos-1,y_pos,z_pos])/ h;
    acc(1) = -1*(pot_array[x_pos, y_pos + 1,z_pos] - pot_array[x_pos,y_pos-1,z_pos])/ h;
    acc(2) = -1*(pot_array[x_pos, y_pos,z_pos + 1] - pot_array[x_pos,y_pos,z_pos-1])/ h;
    /* writeln("f(x+1) ",pot_array[x_pos+1, y_pos,z_pos]);
    writeln(x_pos," ",y_pos," ",z_pos,"\n"); */

    /* writeln("f(x-1) ",pot_array[x_pos-1, y_pos,z_pos]); */
    /* writeln("f(x+1) - f(x-1) ",(pot_array[x_pos+1, y_pos,z_pos] - pot_array[x_pos-1, y_pos,z_pos])/(2*0.015625)); */


    //acc = acc * pos /len_pos;

    /* acc(0) = a[0];
    acc(1) = a[1];
    acc(2) = a[2]; */
  }
  else { //search csv box
    var center_offset: 3*real = (box_param[0],box_param[1],box_param[2]);
    var grid_size: real = box_param[3];
    var Lbox: real = box_param[4];
    //searches up potential from box and converts to dimensionless units
    var box = open("box-pot.csv",iomode.r); // open box
    var boxchannel = box.reader(); //open reading channel
    var dr: real = grid_size / Lbox; //gridpts per codeunit
    //writeln(dr);
    var tmp: string;
    var tmp2: string;
    //determine column number from z position
    var col: int = floor(abs(-Lbox/2.0 + center_offset[2] - pos[2]) * dr): int;
    var row: real = floor(abs(-Lbox/2.0 + center_offset[0] - pos[0]) * dr);
    //determine row number from x and y positions
    row = 128 * row; //row = 10 * row;
    row += (floor(abs(-Lbox/2.0 + center_offset[1] - pos[1]) * dr));
    var frow: int = row: int; //final row number
    //writeln(frow);
    //writeln(col);
    //writeln(pos);
    for i in 1..frow{
      boxchannel.readln(tmp);
    }
    for i in 0..col{
      /* boxchannel.read(tmp2); */
      boxchannel.read(tmp);
    }
    /* writeln(tmp2); */
    /* var a: real = tmp2: real; */
    var a: real = tmp: real;
    //acc = -(pos/len(pos)) * a * (period**2)/r0; //acc made dimensionless
    acc = -(pos/len(pos)) * a;
    //writeln(acc);
  }
  return acc;
}

//initialize particles in a ring
proc init_ring (pos, vel, AM, SD, PS, calcpar,WritingChannel,AMWritingChannel,PSWritingChannel, pot, Ne, r0,nbins,pot_array,box_param,scale_magVel) {
  //assign counterclockwise from rightmost point on x axis.
  var posAngle, velAngle, magVel: real;
  var acc: 3*real;
  acc = force((r0,0.0,0.0),pot,calcpar,pot_array,box_param);
  //acc = force((1.0,0.0,0.0),pot,calcpar);
  /* acc = search_force((r0,0.0,0.0)); */
  //magVel = sqrt(len(acc) * r0) * period / r0;
  magVel = scale_magVel*sqrt(len(acc) * r0);
  //magVel = sqrt(len(acc));

  //var velAngle: real;
  for j in 1..Ne { //Ne=number of desired particles
    posAngle = j * 2.0 * Pi / Ne; //cycle thru angles in a circle

    pos[j][0] = r0 * cos(posAngle);
    pos[j][1] = r0 * sin(posAngle);
    pos[j][2] = 0.0;
    /* pos[j][0] = cos(posAngle); //dimensionless
    pos[j][1] = sin(posAngle);
    pos[j][2] = 0.0; */

    velAngle = (Pi / 2.0) - posAngle; //angles are complementary

    vel[j][0] = magVel * cos(velAngle); //dimensionless
    vel[j][1] = -1.0 * magVel * sin(velAngle);
    vel[j][2] = 0.0;

    AM[j] = angMom(pos[j],vel[j]);
    //writeln("initial pos ", pos[j]);
    //writeln("initial vel ", vel[j]);
    WritingChannel.write(pos[j][0]*vel[j][0] + pos[j][1]*vel[j][1] + pos[j][2]*vel[j][2],",0.0,0.0,",len(vel[j]),",",len(pos[j]),",");
    WritingChannel.write(pos[j][0],",",pos[j][1],",",pos[j][2],",");
    WritingChannel.write(vel[j][0],",",vel[j][1],",",vel[j][2],"\n");
    AMWritingChannel.write(AM[j],",");
  }
  pow_spec(pos,vel,PS,PSWritingChannel,nbins,Ne);
  AMWritingChannel.write("\n");

  //iterate through 2Pi radians divided by 10 (Ne) particles.
  //for each, multiply r0 by sin to get y, r0 by cosine to get x
  //calculate magnitude of centripetal velocity vector for each particles
  //angle velocity vector makes with axis is 90 - theta.
  //multiply magnitude of velocity vector with sin and cosine of angle.
}

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc fwd_orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail, calcpar,k, M, myReadingChannelv,myReadingChannelp, offsetOption, myWritingChannel) {
    if integrator == 0  { //if leapfrog
      //move velocity forward half a timestep
      halfstep(pos[0],vel[0],pot,dt,1.0, calcpar);
      //myWritingChannel.write(aukpc * pos[0][0],",",aukpc * pos[0][1],",",aukpc * vel[0][0],",",aukpc * vel[0][1],"\n");
      writeln("pos cluster after first halfstep: ", aukpc * pos[0]);
      writeln("vel cluster after first halfstep: ",aukpc * vel[0]);

      for i in 1..N-1 {//make N full steps in pos and vel forwards
        //mcl -= dm;//decrease mass

        leapfrog(pos,vel,i,dt,1.0, calcpar);
        //myWritingChannel.write(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1],",");

        for j in 0..k-1 {
          stream_step(pos_lead[j], vel_lead[j], pos[i], dt, calcpar);
          stream_step(pos_trail[j], vel_trail[j], pos[i], dt, calcpar);
          //myWritingChannel.write(",",pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1]);
          //writeln("stars ",j," at ",i," step ",aukpc * pos_lead[j][0],",",aukpc * pos_lead[j][1],",",aukpc * vel_lead[j][0],",",aukpc * vel_lead[j][1],",",aukpc * pos_trail[j][0],",",aukpc * pos_trail[j][1],",",aukpc * vel_trail[j][0],",",aukpc * vel_trail[j][1]);
        }

        if i % M == 0 { //EJECT

          //load dvl, dvt, r1, and r2 based on offsetOption
          var dvl, dvt: real;
          var r1, r2: 3*real;
          loadOffsets(dvl,dvt,r1,r2,myReadingChannelv,myReadingChannelp,randStream);

          eject(pos[i],vel[i],pos_lead[k], vel_lead[k], pos_trail[k], vel_trail[k], dvl, dvt, r1, r2, calcpar);
          //myWritingChannel.write(dvl,",",dvt,",",aukpc * pos_lead[k][0],",",aukpc * pos_lead[k][1],",",aukpc * vel_lead[k][0],",",aukpc * vel_lead[k][1]);
          //myWritingChannel.write(",",aukpc * pos_trail[k][0],",",aukpc * pos_trail[k][1],",",aukpc * vel_trail[k][0],",",aukpc * vel_trail[k][1],"\n");


          //writeln("after ejecting ", vel_lead[k]);
          //writeln(pos_trail[k]);
          k+=1;
        }
        /* //print snapshot at 5000th timestep
        if i == N-1000 {
          for j in 1..k {
            myWritingChannel.write(pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1],"\n");
          }
        }
        */
        //myWritingChannel.write("\n");

      }
      //move velocity backward half a timestep
      halfstep(pos[N-1],vel[N-1],pot,dt,-1.0, calcpar);
      //position of cluster at last timestep
      writeln("pos cluster after last timestep ", aukpc * pos[N-1]);
      writeln("vel cluster after last timestep ", aukpc * vel[N-1]);

      //myWritingChannel.write(aukpc * pos[N-1][0],",",aukpc * pos[N-1][1],",",aukpc * vel[N-1][0],",",aukpc * vel[N-1][1],"\n");
      if offsetOption == 1 {
        myReadingChannelv.close();
        myReadingChannelp.close();
        /* cfile1.close();
        cfile2.close(); */
      }

      //myWritingChannel.write("\n");

      //positions of streams at last timestep
      for j in 0..k-1 {
        myWritingChannel.write("\n",aukpc * pos[j][0],",",aukpc * pos[j][1],",",aukpc * pos[j][2],",",aukpc * vel[j][0],",",aukpc * vel[j][1],",",aukpc * vel[j][2],",");
        myWritingChannel.write(aukpc * pos_lead[j][0],",",aukpc * pos_lead[j][1],",",aukpc * vel_lead[j][0],",",aukpc * vel_lead[j][1],",",aukpc * pos_trail[j][0],",",aukpc * pos_trail[j][1],",",aukpc * vel_trail[j][0],",",aukpc * vel_trail[j][1]);
      }


      //myWritingChannel.write(aukpc * pos[N-1][0],",",aukpc * pos[N-1][1],",",aukpc * vel[N-1][0],",",aukpc * vel[N-1][1],",");
  }

}


//moves each particle in ring forward
proc fwd_orbit_ring (pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel, SDWritingChannel,AMWritingChannel,PSWritingChannel, Ne, r0, nbins,pot_array, box_param) {
    if integrator == 0  { //if leapfrog
      //move velocity of all particles a halfstep forward

      for j in 1..Ne {
        halfstep(pos[j],vel[j],pot,dt,1.0, calcpar,pot_array,box_param);
        AM[j] = angMom(pos[j],vel[j]);
        //writeln("vel after halfstep ",vel[j]);
        //writeln("pos after halfstep ",pos[j]);

        //WritingChannel.write(pos[j][0] * aukpc,",",pos[j][1] * aukpc,",",pos[j][2] * aukpc,",");
        //WritingChannel.write(vel[j][0] * aukpc,",",vel[j][1] * aukpc,",",vel[j][2] * aukpc,"\n");
        //AMWritingChannel.write(AM[j],",");
      }
      writeln("first halfstep forward complete");
      //AMWritingChannel.write("\n");
      for i in 0..N { //for each subsequent timestep
        //move velocity of all particles a fullstep forwards
        for j in 1..Ne {
          step(pos[j], vel[j], dt, calcpar,WritingChannel, pot,pot_array,box_param);
          //writeln("streamstep complete");
          AM[j] = angMom(pos[j],vel[j]);
          //writing position and velocity of to file (one particle per line)
          WritingChannel.write(pos[j][0],",",pos[j][1],",",pos[j][2],",");
          WritingChannel.write(vel[j][0],",",vel[j][1],",",vel[j][2],"\n");
          //writing angular momentum to a diff file
          AMWritingChannel.write(AM[j],",");
          }
        AMWritingChannel.write("\n");
        SD[i] = stdev(AM,Ne);
        pow_spec(pos,vel,PS,PSWritingChannel,nbins,Ne);
        //writeln("fullstep complete");
        /* //print snapshot at 5000th timestep
        if i == N-1000 {
          for j in 1..k {
            myWritingChannel.write(pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1],"\n");
          }
        }
        */
      }
      //for final timestep, move everything back
      for j in 1..Ne {
        halfstep(pos[j],vel[j],pot,dt,-1.0, calcpar,pot_array,box_param);
      }
      for i in 0..N {
        //initially Sd[i] is multiplied by r0/period?
        //SDWritingChannel.write(i,",",SD[i] * (r0/period),"\n");
        SDWritingChannel.write(i,",",SD[i],"\n");
      }
  }
}

//back orbit procedure for cluster and stream setup
proc back_orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail, calcpar)
{
  if integrator == 0  { //if leapfrog
    //move velocity forward half a timestep
    //writeln("pos cluster before first halfstep: ",aukpc * pos[0]);
    //writeln("vel cluster before first halfstep: ",aukpc * vel[0]);
    halfstep(pos[0],vel[0],pot,dt,-1.0, calcpar);
    //myWritingChannel.write(aukpc * pos[0][0],",",aukpc * pos[0][1],",",aukpc * vel[0][0],",",aukpc * vel[0][1],"\n");
    //writeln("pos cluster after first halfstep: ",aukpc * pos[0]);
    //writeln("vel cluster after first halfstep: ",aukpc * vel[0]);
  for i in 1..N-1 {//make N-1 full steps in pos and vel forwards
    //mcl -= dm;//decrease mass

    leapfrog(pos,vel,i,dt,-1.0, calcpar);
    //myWritingChannel.write(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1],"\n");
    //writeln(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1]);
  }
  //move velocity backward half a timestep
  halfstep(pos[N-1],vel[N-1],pot,dt,1.0, calcpar);
  //myWritingChannel.write(aukpc * pos[N-1][0],",",aukpc * pos[N-1][1],",",aukpc * vel[N-1][0],",",aukpc * vel[N-1][1],"\n");
  //writeln("ending pos ",pos[N-1]," ending vel ",vel[N-1]);
  }
}


//procedure to advance particles by a timestep. only subject to force from central pot
proc step(ref pos, ref vel, dt, calcpar,WritingChannel,pot,pot_array,box_param) {

  //update position of jth particle
  pos += dt * vel;
  //writeln(vel);
  //search acceleration
  var a: 3*real;
  //writeln("searching for force");
  a = force(pos,pot,calcpar,pot_array,box_param); //returns dimensionless acc
  //writeln("found force");

  //writeln(a);
  //update velocity of jth particle
  vel += dt * a;
  //writeln(vel);

  //writeln("position ",pos);
  //writeln("velocity ", vel);
  WritingChannel.write(pos[0]*vel[0] + pos[1]*vel[1] + pos[2]*vel[2],",");
  WritingChannel.write(a[0]*vel[0] + a[1]*vel[1] + a[2]*vel[2],",");
  WritingChannel.write(len(a),",");
  WritingChannel.write(len(vel),",");
  WritingChannel.write(len(pos),",");
  //writeln("magnitude velocity: ",len(vel));
  //writeln("magnitude position: ",len(pos));
  //writeln(a[0]*vel[0] + a[1]*vel[1] + a[2]*vel[2]);
}

//procedure to advance ejected particles by a timestep
proc stream_step(ref pos, ref vel, pos_cl, dt, calcpar, pot,Rcl) {
  //update position of jth particle
  pos += dt * vel;
  //calculate acceleration from Msun
  var a: 3*real;
  a = force(pos,pot, calcpar);
  //writeln("pos ",pos);
  //writeln("a ",a);
  //writeln("force plummer ",force_plummer(pos_cl - pos, Rcl));
  //calculate plummer acceleration
  a += force_plummer(pos_cl - pos, Rcl);
  //update velocity of jth particle
  vel += dt * a;
}

//SHOULD NO LONGER BE IN USE; CONSOLIDATED WITH FORCE
proc search_force(pos, pot, calcpar,pot_array){

  //search acceleration
  var acc: 3*real;
  if pot == 6 {
    var dr: real = 128 / 2.0;
    var x_pos: int = floor(abs(-0.9921875 - pos[0]) * dr): int;
    var y_pos: int = floor(abs(-1.00585938 - pos[1]) * dr): int;
    var z_pos: int = floor(abs(-1.0 - pos[2]) * dr): int;
    //writeln("x pos ",abs(-1.0 - pos[0]) * dr-0.25," y pos ",abs(-1.0 - pos[1]) * dr-0.25," z pos ",abs(-1.0 - pos[2]) * dr-0.25);
    /* writeln("x pos ",x_pos," y pos ",y_pos," z pos ",z_pos); */
    //var a: [0..2] real;
    var a: real;
    a = calcpar[x_pos,y_pos,z_pos];
    /* writeln("a ",a); */
    //take the derivative of a
    var h: real = 2.0 * (2.0 / 128);

    //derivative blows up at edges, check for edge
    if x_pos > 126 {
      x_pos = 126;
    }
    if y_pos > 126 {
      y_pos = 126;
    }
    if z_pos > 126 {
      z_pos = 126;
    }

    acc(0) = -1*(pot_array[x_pos + 1, y_pos,z_pos] - pot_array[x_pos-1,y_pos,z_pos])/ h;
    acc(1) = -1*(pot_array[x_pos, y_pos + 1,z_pos] - pot_array[x_pos,y_pos-1,z_pos])/ h;
    acc(2) = -1*(pot_array[x_pos, y_pos,z_pos + 1] - pot_array[x_pos,y_pos,z_pos-1])/ h;
    /* writeln("f(x+1) ",pot_array[x_pos+1, y_pos,z_pos]);
    writeln("f(x-1) ",pot_array[x_pos-1, y_pos,z_pos]);
    writeln("f(x+1) - f(x-1) ",(pot_array[x_pos+1, y_pos,z_pos] - pot_array[x_pos-1, y_pos,z_pos])/(2*0.015625)); */


    //acc = acc * pos /len_pos;

    /* acc(0) = a[0];
    acc(1) = a[1];
    acc(2) = a[2]; */
  }
  else {
    //searches up potential from box and converts to dimensionless units
    var nfwbox = open("nfwbox.csv",iomode.r); // open box
    var boxchannel = nfwbox.reader(); //open reading channel
    var dr: real = 128.0 / toCodeLength(40); //nsteps / range(made dimensionless)

    //writeln(dr);
    var tmp: string;
    //determine column number from z position
    var col: int = floor(abs(toCodeLength(-20.0) - pos[2]) * dr): int;
    var row: real = floor(abs(toCodeLength(-20.0) - pos[0]) * dr);
    //determine row number from x and y positions
    row = 128 * row; //row = 10 * row;
    row += (floor(abs(toCodeLength(-20.0) - pos[1]) * dr));
    var frow: int = row: int; //final row number
    //writeln(frow);
    //writeln(col);
    //writeln(pos);
    for i in 1..frow{
      boxchannel.readln(tmp);
    }
    for i in 0..col{
      boxchannel.read(tmp);
    }
    //writeln(tmp);
    var a: real = tmp: real;
    //acc = -(pos/len(pos)) * a * (period**2)/r0; //acc made dimensionless
    acc = -(pos/len(pos)) * a;
    //writeln(acc);
  }
  return acc;

}
//shifts velocity by a halfstep in direction of sign
proc halfstep(p, ref v, pot, dt, sign, calcpar,pot_array,box_param) {
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = force(p,pot,calcpar,pot_array,box_param);
  //writeln("dot force and vel ",a[0]*v[0] + a[1]*v[1] + a[2]*v[2]);

  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
}

proc angMom(pos,vel){
  var omega: 3*real;
  omega = cross(pos, vel); //r x v
  //omega = omega / (len(pos)**2); //(r x v)/r^2
  return len(omega);
}
proc cross((x1,y1,z1),(x2,y2,z2)) {
  var cross_product: 3*real;
  cross_product[0] = (y1 * z2) - (z1 * y2);
  cross_product[1] = (z1 * x2) - (x1 * z2);
  cross_product[2] = (x1 * y2) - (y1 * x2);
  return cross_product;
}
proc stdev(AM, Ne){
    var mean, SD: real = 0.0;
  for j in 1..Ne { //iterate over all stream particles
    mean += AM[j];
  }
  mean = mean / Ne;
  for j in 1..Ne {
    SD += (AM[j] - mean)**2;
  }
  SD = SD / Ne;
  SD = sqrt(SD);
  //return SD/r0;
  return SD;
}
proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}
proc pow_spec (pos,vel, PS,PSWritingChannel , nbins, Ne) {//gets computed at each timstep
  //keep track of index of particle
  var j: int = 1;
  var avg: real = Ne / nbins;
  //loop through number of bins, aka sections of angle
  var low_bound, up_bound, posAngle: real = 0.0;
  posAngle = acos( pos[j][0] / len(pos[j]) );
  for i in 1..nbins{
    low_bound = (i - 1.0) * (2.0 * Pi / nbins);
    up_bound = i * (2.0 * Pi / nbins);
    //writeln("low: ",low_bound);
    //writeln("high ", up_bound);
    //calculate angle of position of current particle
    //while angle of position is in current bin, iterate over particles
    PS[i] = 0;
    while posAngle > low_bound && posAngle <= up_bound{
      //writeln("j ", j);
      //writeln(pos[j]);
      //writeln(posAngle);
      PS[i] += 1; //add one to current bin
      j += 1; //move on to next particle
      posAngle = acos( pos[j][0] / len(pos[j]) ); //update angle
      if pos[j][1] < 0{  //if y is neg, aka in lower quad,
        posAngle = posAngle + 2 * (Pi - posAngle);
        }
    }
    //PSWritingChannel.write(i,",",PS[i] / avg,",");
    PSWritingChannel.write(PS[i] / avg,",");

  }
  PSWritingChannel.write("\n");
}
//update positions at every timestep based on potential box

proc toreal (str) {
  /*helper function to load in offsets from c code to replicate results of random
  radial offsets from c in chapel. converts offset string into a real*/
  var numb: [0..1] int;
  var n = str.split(".");
  numb[0] = abs(n[0]:int);
  numb[1] = n[1]:int;
  var n1: real = 1.0 * numb[1] / 10**n[1].size + numb[0];
  if n[0][0] == '-' then n1 = n1 * -1.0;
  //writeln(n1);
  return n1;
}

proc norm_rand (ref r1, ref r2, randStream) {
  /*function to recreate random radial offsets in chapel.
  fills two tuples with normally distributed random numbers */
  var arr1: [0..2] real;
  var arr2: [0..2] real;
  randStream.fillRandom(arr1);
  randStream.fillRandom(arr2);

  r1[0] = sqrt(-2 * log(arr1[0])) * cos(2*Pi*arr1[1]);
  r1[1] = sqrt(-2 * log(arr1[0])) * sin(2*Pi*arr1[1]);

  r1[2] = sqrt(-2 * log(arr1[2])) * cos(2*Pi*arr2[0]);
  r2[0] = sqrt(-2 * log(arr1[2])) * sin(2*Pi*arr2[0]);

  r2[1] = sqrt(-2 * log(arr2[1])) * cos(2*Pi*arr2[2]);
  r2[2] = sqrt(-2 * log(arr2[1])) * sin(2*Pi*arr2[2]);
}

proc readTuple(ReadingChannel, ref r){
  /* helper function to load in offsets from c code to replicate results of random
  radial offsets from c in chapel.fills tuple r with numbers from reading
  channel and calls toreal to convert from string to real type*/
  var tmp: string;
  for i in 0..2 {
    ReadingChannel.read(tmp);
    r(i) = toreal(tmp);
  }
}

proc loadOffsets(ref dvl, ref dvt, ref r1, ref r2, ReadingChannelv, ReadingChannelp, randStream, offsetOption, offset,WritingChannel) {
  /*central function controlling offsets.
  offsetOption 0 = no radial offset
  offsetOption 1 = use C radial offsets to replicate results in chpl
  offsetOption 2 = create new random radial offsets entirely in chpl
  calls all helper offset loading functions */

  if offsetOption == 0 { //if adding no random radial offset
    dvl = 0;
    dvt = 0;
    r1 = (0.0,0.0,0.0);
    r2 = (0.0,0.0,0.0);
    offset = [0.0,0.0];
  }
  else if offsetOption == 1 { //if using C's random radial offsets
    readTuple(ReadingChannelv, r1);
    dvl = len(r1) * offset[1]/3;
    readTuple(ReadingChannelv,r2);
    dvt = len(r2) * offset[1]/3;
    readTuple(ReadingChannelp,r1);
    readTuple(ReadingChannelp,r2);
  }
  else if offsetOption == 2 { //if adding chapel generated random radial offsets
    norm_rand(r1,r2,randStream); //get new tuple of normalized random numbers
    dvl = len(r1)*offset[1]/3; //convert to maxwell distribution
    dvt = len(r2)*offset[1]/3; //convert to maxwell distribution
    norm_rand(r1, r2, randStream);//get new tuple of normalized random numbers
  }
  WritingChannel.write(dvl,",",dvt,",");
}


proc toSeconds (l) { //not in use
  var a: real = ((8 * Pi)/(3* H0**2 *omegaM0))**0.5;
  return a / l**2;
}

proc toMeters (l) { //not in use
  var a: real = (8 * Pi * (hbar**2))/(3 * (maxion**2) * (H0**2) * omegaM0);
  return (a**0.25) / l;
}

proc toCodeTime(yrs){
  return yrs / (75.4212 * 10**9);
}
proc toCodeLength(kpc){
  return kpc / 38.3609;
}
proc toCodeMass(kg) {
  return kg/4.42837e36;
}

proc force_plummer(r, Rcl, mcl) {
  // dist from cluster to particle, cluster radius
  var dist = len(r);
  var raux: real = sqrt(dist**2 + Rcl**2);
  var acc: 3*real;
  acc = (mcl / raux**3 ) * r;
  //write("{",acc[0],",",acc[1],"},");
  return acc;
}

proc tidal_radius(pos_cl,vel_cl,omega,calcpar, pot, mcl) {
  //writeln("pos cl ",pos_cl," vel cl ",vel_cl," omega ",omega);
  var un_vec: 3*real = pos_cl/len(pos_cl);
  var delta = 0.02 * kpcau * un_vec;
  var x1: 3*real = pos_cl - delta;
  //writeln("x1 " , x1);
  var x2: 3*real = pos_cl + delta;
  //writeln("x2 " , x2);
  var a1: 3*real = force(x1,pot,calcpar);
  var a2: 3*real = force(x2,pot,calcpar);
  //writeln("a1 " , a1);
  //writeln("a2 " , a2);

  var dpot: real = (len(a1) - len(a2))/(len(x1-x2));
  //var dpot: 3*real = (a1 - a2)/(len(x1-x2));
  //writeln("len x1-x2 ",len(x1-x2));
  //writeln("a1-a2 ",a1-a2);
  //writeln("mcl ",mcl);
  return un_vec*((mcl/abs(omega*omega + dpot))**(1.0/3.0));
}

//procedure to eject particles
proc eject(pos_cl, vel_cl, ref pos_lead, ref vel_lead, ref pos_trail, ref vel_trail, dvl, dvt, r1, r2, calcpar)
{

  //calculate angular velocity of Cluster
  //writeln("pos cl at ejection ",pos_cl, " vel cl at ejection", vel_cl);
  var omega: 3*real;
  omega = cross(pos_cl, vel_cl); //r x v
  //writeln("omega of cluster: ",omega);
  var r = len(pos_cl);
  omega = omega / (r**2); //(r x v)/r^2
  var om: real = len(omega);

  var Rj: 3*real = tidal_radius(pos_cl, vel_cl, om, calcpar);//calculate tidal radius
  var dRRj: real = len(Rj) * offset[0];
  //initial position of particle is position of cluster plus or minus tidal radius
  //writeln("Rj: ",Rj);
  //pos_lead = pos_cl - (Rj * (pos_cl/r)); //multiplied by unit vector
  //pos_trail = pos_cl + (Rj * (pos_cl/r));

  pos_lead = pos_cl - Rj + (dRRj * r1);
  pos_trail = pos_cl + Rj + (dRRj * r2);

  //pos_lead = pos_cl - Rj;
  //pos_trail = pos_cl + Rj;

  var mag_vcl: real = len(vel_cl);
  var uv_vel: 3*real = vel_cl / mag_vcl; // unit vector of vel vector
  var uv_pos: 3*real = pos_cl / r; //unit vector of pos vector

  vel_lead = uv_vel * (mag_vcl - len(Rj) * om) - (dvl * uv_pos);
  vel_trail = uv_vel * (mag_vcl + len(Rj) * om) + (dvt * uv_pos);
  /*
  var mag_vcl: real = len(vel_cl);
  var uv_vel: 3*real = vel_cl / mag_vcl; // unit vector of vel vector
  vel_lead = uv_vel * (mag_vcl - len(Rj) * om);
  vel_trail = uv_vel * (mag_vcl + len(Rj) * om);
  */

  //writeln("cluster position at ejection ",pos_cl);
  //writeln("cluster velocity at ejection ",vel_cl);
  //writeln("pos_lead at ejection ",pos_lead);
  //writeln("vel_lead at ejection", vel_lead);
  //writeln("pos_trail at ejection ",pos_trail);
  //writeln("vel_trail at ejection", vel_trail);
}
/*
proc main_ring(){

  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real; //stores angular momentum of each particle, updated at each timestep
  var SD: [0..N] real; //stores std dev at each timestep
  var PS: [1..nbins] real; //whole array gets updated at each timestep
  //hardcode galactic potential parameters
  var par: [0..5] real = [417.0, 36.54 * kpckm, 90.0 * Pi / 180, 1.0, 1.0 ,1.0];
  var calcpar: [0..5] real;

  var chfile = open("mathematica/bp19.csv",iomode.cw); //create test.csv and open
  var SDfile = open("mathematica/SD19.csv",iomode.cw); //create test.csv and open
  var AMfile = open("mathematica/AM19.csv",iomode.cw); //create test.csv and open
  var PSfile = open("mathematica/PS19.csv",iomode.cw); //create test.csv and open
  var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
  var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
  var WritingChannel = chfile.writer(); //open writing channel to test.csv
  var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv

  load_pot(par,calcpar,pot,Dom,arr);
  //set initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar, WritingChannel,AMWritingChannel,PSWritingChannel);
  //dt = dt / period;
  writeln("initialized ring");
  writeln("pos 1 ",pos[1]);
  writeln("vel 1", vel[1]);
  //
   writeln("calling search force ",search_force((0.5,0.5,0.5)));
  //writeln("calling search force ",search_force((0.4,0.4,0.4)));

  fwd_orbit(pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel,SDWritingChannel,AMWritingChannel,PSWritingChannel);
  WritingChannel.close();
  SDWritingChannel.close();
  AMWritingChannel.close();
  PSWritingChannel.close();
  PSfile.close();
  AMfile.close();
  SDfile.close();
  chfile.close();
}

proc main_stream () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;

  var pos_lead: [0..Ne] 3*real;
  var pos_trail: [0..Ne] 3*real;
  var vel_lead: [0..Ne] 3*real;
  var vel_trail: [0..Ne] 3*real;

  //hardcode galactic potential parameters
  var par: [0..5] real = [417.0, 36.54 * kpckm, 90.0 * Pi / 180, 1.0, 1.0 ,1.0];
  var calcpar: [0..5] real;
  load_pot(par,calcpar,pot,Dom,arr);


  //hardcode initial position and velocity of cluster
  pos[0]=(50.0*kpcau,0,0);
  vel[0]=(0,0.5* sqrt(calcpar[0]/len(pos[0])),0.0); //(if pot 3) set velocity equal to half of centripetal velocity
  //vel[0]=(0,0.5* sqrt(1000000000*Msun/len(pos[0])),0); // (if pot 0)set velocity equal to half of centripetal velocity

  myWritingChannel.write("x cluster,y cluster,x cluster vel,y cluster vel,x lead trail,y lead tail,x vel lead tail, y vel lead tail,x trail tail,y trail tail,x vel trail tail, y vel trail tail\n");

  //writeln("initial energy ", energy(pos,vel,0));
  back_orbit(pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail, calcpar);
  pos[0] = pos[N-1];
  vel[0] = vel[N-1];
  fwd_orbit(pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail, calcpar);
  //writeln("final energy ", energy(pos,vel,N));
  myWritingChannel.close();
  chfile.close();
}
*/
