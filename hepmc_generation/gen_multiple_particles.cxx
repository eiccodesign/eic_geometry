#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <random>


using namespace HepMC3;

// Generate single particle as input to the Insert simulation.
// --
// We generate events with  polar angle with respect to
// the proton direction and then rotate so that the events are given
// in normal lab coordinate system
// --
void gen_multiple_particles(
                    int n_events = 1000, 
                    const char* out_fname = "gen_particles.hepmc", 
                    TString particle_name = "e-",
                    double th_min = 3., // Minimum polar angle, in degrees
		                double th_max = 3., // Maximum polar angle, in degrees
		                double phi_min = 0., // Minimum azimuthal angle, in degrees
                    double phi_max = 360., // Maximum azimuthal angle, in degrees
                    double p_low = 10.,  // Momentum in GeV/c,
                    double p_high = 100.,
		                TString dist = "log10continuous",  // Momentum distribution: fixed, uniform, 
                                              // Gaussian, log10continuous, discrete
                    int min_num_particles = 1,
                    int max_num_particles = 10
                  )
{ 
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);
  if(dist == "fixed" && p_low != p_high)
  {
    std::cerr << "Energies need to match for fixed distribution!" << std::endl;
    exit(1);
  }

  // Random number generator
  TRandom3 *r1 = new TRandom3(0); //Use time as random seed
  
  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();
  TParticlePDG *particle = pdg->GetParticle(particle_name);
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {

    //Set the event number
    evt.set_event_number(events_parsed);

    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 111 - pi0
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    //If wanted, set non-zero vertex
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    double vt = 0.;

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    evt.shift_position_by(FourVector(vx, vy, vz, vt));
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    // For each event, have a random number of particles between min_num_particles and max_num_particles
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(min_num_particles, max_num_particles); // define the range
    const int num_event_particles = distr(gen);
    
    // Generating num_event_particles particles
    // Each is distributed uniformly in theta and phi, with energy distribution according to argument
    for(int i_particle = 0; i_particle < num_event_particles; i_particle++)
    {
      // Define momentum with respect to proton direction
      double phi   = r1->Uniform(phi_min*TMath::DegToRad(),phi_max*TMath::DegToRad());
      double th    = r1->Uniform(th_min*TMath::DegToRad(),th_max*TMath::DegToRad());

      //Total momentum distribution
      double pevent = -1;
      if(dist=="fixed"){ //fixed
        pevent = p_low;
      }
      else if(dist == "uniform"){ //Uniform: random between p_low and p_high
        pevent = r1->Uniform(p_low, p_high);
      }
      else if(dist == "gaussian"){  //Gaussian: Sigma = 0.1*mean
        while(pevent<0) //Avoid negative values
        {
          pevent = r1->Gaus(p_low, 5);
        } 
      }
      else if(dist == "log10continuous")
      {
        // For continuous in log10
        // Set to between 1 and p_high GeV by default
        double num_log_uniform_energies = r1->Uniform(log10(p_low), log10(p_high));      
        pevent = pow(10, num_log_uniform_energies);
        
      }
      else if(dist == "discrete")
      {
        // For discrete in log10
        const int num_loguniform_energies = 36;
        const int random_power = (int) r1->Uniform(0, num_loguniform_energies);
        double random_pow = (random_power*0.0423)+1;
        pevent = (int) pow(10, random_pow);
      }

      double px    = pevent * std::cos(phi) * std::sin(th);
      double py    = pevent * std::sin(phi) * std::sin(th);
      double pz    = pevent * std::cos(th);
      TVector3 pvec(px,py,pz); 

      //Rotate to lab coordinate system
      double cross_angle = -25./1000.; //in Rad
      TVector3 pbeam_dir(sin(cross_angle),0,cos(cross_angle)); //proton beam direction
      // pvec.RotateY(-pbeam_dir.Theta()); // Theta is returned positive, beam in negative X
      // type 1 is final state
      // pdgid 11 - electron 0.510 MeV/c^2
      GenParticlePtr p_ith_particle = std::make_shared<GenParticle>(
          FourVector(
              pvec.X(), pvec.Y(), pvec.Z(),
              sqrt(pevent*pevent + (mass * mass))),
          pdgID, 1);

      v1->add_particle_out(p_ith_particle);
    }
    
    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 100 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}