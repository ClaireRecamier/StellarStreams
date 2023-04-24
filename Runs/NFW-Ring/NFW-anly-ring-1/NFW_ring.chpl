use streakline;
use IO;
use Math;
use Random;
use HDF5;
use H5Helper;
//DIMENSIONLESS VERSION

//setup constants
config const r0 = toCodeLength(50.0); // radius of ring - kpc converted to codeunits
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const Ne = 1; //number of particles in ring
config const nbins = 4;
var period: real = 0.0;
config const po = (0.5**4) * 27492.260803351877;
config const rc = 2.0*0.05359134269304645;
var pot = 3; //set galactic potential to analytic NFW
config const integrator = 0;   //set Integrator to LF
config const N = 100000;//6000 set number of timesteps to equal total of 6 Gyr
var dt = toCodeTime(1000000.0); //set timestep to seconds per Myr
config const grid_size: int = 128;

proc main(){
  //create array of tuples to hold position and velocity and velocity at each timestep

  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real; //stores angular momentum of each particle, updated at each timestep
  var SD: [0..N] real; //stores std dev at each timestep
  var PS: [1..nbins] real; //whole array gets updated at each timestep

  var Dom = {1..0,1..0,1..0};
  var pot_array : [Dom] real; //placeholder
  var center_offset: 3*real = (0.0,0.0,0.0);

  var chfile = open("pos-vel.csv",iomode.cw); //create test.csv and open
  var SDfile = open("standard-dev.csv",iomode.cw); //create test.csv and open
  var AMfile = open("ang-mom.csv",iomode.cw); //create test.csv and open
  var PSfile = open("pow-spec.csv",iomode.cw); //create test.csv and open
  var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
  var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
  var WritingChannel = chfile.writer(); //open writing channel to test.csv
  var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv

  //hardcode galactic potential parameters UNITS
  var par: [0..5] real = [toCodeLength(417.0 * kmkpc)/ toCodeTime(secyr), toCodeLength(36.54), 90.0 * pi / 180, 1.0, 1.0 ,0.94];
  var calcpar: [0..5] real;
  var scale_magVel: real = 0.75;

  //pot_box parameter pot set to 5, so soliton approximation is used in creation of box
  load_anly_pot(r0,par,calcpar,pot,dt);
  writeln("initialized potential");
  var acc: 3*real = force((r0,0.0,0.0),pot,calcpar,pot_array,center_offset,grid_size);
  var period: real = 2 * Pi * r0 / sqrt(len(acc) * r0); //in seconds; 2Pi * r/v
  writeln("period /dt ",period/dt,"\n");

  //initialize initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar,WritingChannel,AMWritingChannel,PSWritingChannel, pot, Ne, r0,nbins,pot_array,center_offset,grid_size,scale_magVel);
  //dt = dt / period;
  writeln("initialized ring");
  /* writeln(pos[1]); */
  /* writeln(vel[1]); */
  /* force((toCodeLength(-20),toCodeLength(-20),toCodeLength(-20)),pot,calcpar,pot_array,center_offset,grid_size);    */
  fwd_orbit_ring(pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel, SDWritingChannel,AMWritingChannel,PSWritingChannel, Ne, r0, nbins,pot_array, center_offset, grid_size);
  WritingChannel.close();
  SDWritingChannel.close();
  AMWritingChannel.close();
  PSWritingChannel.close();
  PSfile.close();
  AMfile.close();
  SDfile.close();
  chfile.close();
}


/*
proc main () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;

  var pos_lead: [0..Ne] 3*real;
  var pos_trail: [0..Ne] 3*real;
  var vel_lead: [0..Ne] 3*real;
  var vel_trail: [0..Ne] 3*real;

  //hardcode galactic potential parameters
  var calcpar: [0..5] real;
  var par: [0..5] real = [toCodeLength(417.0 * (10**3) * mau)/ toCodeTime(secyr), toCodeLength(36.54), 90.0 * pi / 180, 1.0, 1.0 ,0.94];
  if pot == 3 { //if using triaxial NFW potential
    //assuming par = [V, rhalo, phi, q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1];
    //writeln("initializing param ",calcpar);
  }

  //hardcode initial position and velocity
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

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc fwd_orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail, calcpar) {
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
        //print snapshot at 5000th timestep
        if i == N-1000 {
          for j in 1..k {
            myWritingChannel.write(pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1],"\n");
          }
        }

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
        cfile1.close();
        cfile2.close();
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

//procedure to advance ejected particles by a timestep
proc stream_step(ref pos, ref vel, pos_cl, dt, calcpar) {
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

  //writeln("cluster position at ejection ",pos_cl);
  //writeln("cluster velocity at ejection ",vel_cl);
  //writeln("pos_lead at ejection ",pos_lead);
  //writeln("vel_lead at ejection", vel_lead);
  //writeln("pos_trail at ejection ",pos_trail);
  //writeln("vel_trail at ejection", vel_trail);
}

proc cross((x1,y1,z1),(x2,y2,z2)) {
  var cross_product: 3*real;
  cross_product[0] = (y1 * z2) - (z1 * y2);
  cross_product[1] = (z1 * x2) - (x1 * z2);
  cross_product[2] = (x1 * y2) - (y1 * x2);
  return cross_product;
}

//shifts velocity by a halfstep in direction of sign
proc halfstep(p, ref v, pot, dt, sign, calcpar) {
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = force(p,pot,calcpar);
  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
}

//shifts velocity and position by a full step forward
proc leapfrog(pos, vel, i, dt, sign, calcpar) {
  var a: 3*real;
  pos[i] = pos[i-1] + (vel[i-1] * dt * sign);
  a = force(pos[i],pot,calcpar);
  vel[i] = vel[i-1] + (dt * sign * a);

}


proc force(pos,pot,calcpar){
  var acc: 3*real;
  var r: real;
  if pot == 0 { //if using point mass potential
    r = len(pos);
    //acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
    acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r**3;
    //writeln("acc ",acc);
  }
  else if pot == 3 { //NFW triaxial potential, C version
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    r = sqrt(calcpar[1]*pos[0]*pos[0] + calcpar[2]*pos[1]*pos[1] + calcpar[3]*pos[0]*pos[1] + calcpar[4]*pos[2]*pos[2]);
    var aux: real = 0.5 * calcpar[0] / (r**3) * (1.0/(1.0 + calcpar[5]/r)-log(1.0+r/calcpar[5]));

    acc[0]=aux*(2*calcpar[1]*pos[0] + calcpar[3]*pos[1]);
    acc[1]=aux*(2*calcpar[2]*pos[1] + calcpar[3]*pos[0]);
    acc[2]=aux*(2*calcpar[4]*pos[2]);
    //writeln("acc ",acc);
  }
  return acc;
}

proc force_plummer(r, Rcl) {
  // dist from cluster to particle, cluster radius
  var dist = len(r);
  var raux: real = sqrt(dist**2 + Rcl**2);
  var acc: 3*real;
  acc = (mcl / raux**3 ) * r;
  //write("{",acc[0],",",acc[1],"},");
  return acc;
}

proc tidal_radius(pos_cl,vel_cl,omega,calcpar) {
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

proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}


proc energy(pos,vel,i){
  var energy: real = 0.0;

  energy += 0.5 * mcl * len(vel[i])**2;
  energy -= Msun * mcl / len(pos[i]);

  return energy;
}

proc toreal (str) {
  //writeln(str);
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
  //fills two tuples with normally distributed random numbers
  var arr1: [0..2] real;
  var arr2: [0..2] real;
  randStream.fillRandom(arr1);
  randStream.fillRandom(arr2);

  r1[0] = sqrt(-2 * log(arr1[0])) * cos(2*pi*arr1[1]);
  r1[1] = sqrt(-2 * log(arr1[0])) * sin(2*pi*arr1[1]);

  r1[2] = sqrt(-2 * log(arr1[2])) * cos(2*pi*arr2[0]);
  r2[0] = sqrt(-2 * log(arr1[2])) * sin(2*pi*arr2[0]);

  r2[1] = sqrt(-2 * log(arr2[1])) * cos(2*pi*arr2[2]);
  r2[2] = sqrt(-2 * log(arr2[1])) * sin(2*pi*arr2[2]);
}

proc readTuple(myReadingChannel, ref r){
  //fills tuple r with numbers from reading channel
  var tmp: string;
  for i in 0..2 {
    myReadingChannel.read(tmp);
    r(i) = toreal(tmp);
  }
}

proc loadOffsets(ref dvl, ref dvt, ref r1, ref r2, myReadingChannelv, myReadingChannelp, randStream) {
  if offsetOption == 0 { //if adding no random radial offset
    dvl = 0;
    dvt = 0;
    r1 = (0.0,0.0,0.0);
    r2 = (0.0,0.0,0.0);
    offset = [0.0,0.0];
  }
  else if offsetOption == 1 { //if using C's random radial offsets
    readTuple(myReadingChannelv, r1);
    dvl = len(r1) * offset[1]/3;
    readTuple(myReadingChannelv,r2);
    dvt = len(r2) * offset[1]/3;
    readTuple(myReadingChannelp,r1);
    readTuple(myReadingChannelp,r2);
  }
  else if offsetOption == 2 { //if adding chapel generated random radial offsets
    norm_rand(r1,r2,randStream); //get new tuple of normalized random numbers
    dvl = len(r1)*offset[1]/3; //convert to maxwell distribution
    dvt = len(r2)*offset[1]/3; //convert to maxwell distribution
    norm_rand(r1, r2, randStream);//get new tuple of normalized random numbers
  }
  myWritingChannel.write(dvl,",",dvt,",");
}
*/
