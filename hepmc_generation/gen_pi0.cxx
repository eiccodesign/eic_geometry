#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <iostream>
#include <random>
#include <TMath.h>

using namespace HepMC3;

std::tuple<double, int, double> GetParticleInfo(TDatabasePDG* pdg, TString particle_name)
{
  TParticlePDG *particle = pdg->GetParticle(particle_name);
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();
  const double lifetime = particle->Lifetime();
  return std::make_tuple(mass, pdgID, lifetime);
}
// Calculates the decay length of a particle. Samples from an exponential decay.
double GetDecayLength(TRandom3* r1, double lifetime, double mass, double momentum_magnitude)
{ 
  double c_speed = TMath::C() * 1000.; // speed of light im mm/sec
  double average_decay_length = (momentum_magnitude/mass) * lifetime * c_speed;
  return r1->Exp(average_decay_length);
}

// Generate single sigma baryons and decay them to a neutron + 2 photons
void gen_pi0(int n_events = 10000, 
                      const char* out_fname = "gen_sigma_decay.hepmc", 
                      double th_min = 0, // Minimum polar angle, in deg
                      double th_max = 0.3, // Maximum polar angle, in deg
                      double phi_min = 0., // Minimum azimuthal angle, in degrees
                      double phi_max = 360., // Maximum azimuthal angle, in degrees
                      double p_low = 10.,  // Momentum in GeV/c,
                      double p_high = 300.,
                      TString dist = "log10continuous"  // Momentum distribution: fixed, uniform, 
                                              // Gaussian, log10continuous, discrete
                  )
{
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(0);
  cout<<"Random number seed is "<<r1->GetSeed()<<"!"<<endl;

  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();

  auto pi0_info = GetParticleInfo(pdg, "pi0");
  double pi0_mass = std::get<0>(pi0_info);
  int pi0_pdgID = std::get<1>(pi0_info);

  auto photon_info = GetParticleInfo(pdg, "gamma");
  double photon_mass = std::get<0>(photon_info);
  int photon_pdgID = std::get<1>(photon_info);

  for (events_parsed = 0; events_parsed < n_events; events_parsed++)
  {
    //Set the event number
    evt.set_event_number(events_parsed);

    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum with respect to proton direction
    double phi   = r1->Uniform(phi_min*TMath::DegToRad(), phi_max*TMath::DegToRad());
    double theta = r1->Uniform(th_min*TMath::DegToRad(), th_max*TMath::DegToRad());

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
          pevent = r1->Gaus(p_low,0.1*p_low);
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
      double log_min = log10(p_low);
      double log_max = log10(p_high);
      double random_pow = log_min + random_power * (log_max - log_min) / (num_loguniform_energies - 1);
      pevent = (int) pow(10, random_pow);
    }

    // Define momentum with respect to EIC proton beam direction
    Double_t pi0_p     = pevent;
    Double_t pi0_phi   = phi;
    Double_t pi0_th    = theta; // Divide by 1000 for radians
    Double_t pi0_px    = pi0_p * TMath::Cos(pi0_phi) * TMath::Sin(pi0_th);
    Double_t pi0_py    = pi0_p * TMath::Sin(pi0_phi) * TMath::Sin(pi0_th);
    Double_t pi0_pz    = pi0_p * TMath::Cos(pi0_th);
    Double_t pi0_E     = TMath::Sqrt(pi0_p*pi0_p + pi0_mass*pi0_mass);

    TVector3 pi0_pvec(pi0_px, pi0_py, pi0_pz);

    double cross_angle = 0; // in Rad
    TVector3 pbeam_dir(TMath::Sin(cross_angle), 0, TMath::Cos(cross_angle)); //proton beam direction
    pi0_pvec.RotateY(cross_angle); // Theta is returned positive, beam in negative X

    // type 2 is state that will decay
    GenParticlePtr p_pi0 = std::make_shared<GenParticle>(
    FourVector(pi0_pvec.X(), pi0_pvec.Y(), pi0_pvec.Z(), pi0_E), pi0_pdgID, 2 );
    // Generating pi0 particle, will be generated at origin
    // Must have input electron + proton for vertex
    GenVertexPtr pi0_initial_vertex = std::make_shared<GenVertex>();
    pi0_initial_vertex->add_particle_in(p1);
    pi0_initial_vertex->add_particle_in(p2);
    pi0_initial_vertex->add_particle_out(p_pi0);
    evt.add_vertex(pi0_initial_vertex);

    TLorentzVector pi0_lab(pi0_px, pi0_py, pi0_pz, pi0_E);
    TVector3 boost_vec = pi0_lab.BoostVector();

    // Generate gammas in pi0 rest frame
    TLorentzVector gamma1_rest, gamma2_rest;

    double cost_gamma1rest = r1->Uniform(-1,1);
    double th_gamma1rest = std::acos(cost_gamma1rest);
    double sint_gamma1rest = std::sin(th_gamma1rest);

    double phi_gamma1rest = r1->Uniform(-1.*TMath::Pi(),1.*TMath::Pi());
    double cosp_gamma1rest = std::cos(phi_gamma1rest);
    double sinp_gamma1rest = TMath::Sin(phi_gamma1rest);

    gamma1_rest.SetE(pi0_mass/2.);
    gamma1_rest.SetPx( (pi0_mass/2.)*sint_gamma1rest*cosp_gamma1rest );
    gamma1_rest.SetPy( (pi0_mass/2.)*sint_gamma1rest*sinp_gamma1rest );
    gamma1_rest.SetPz( (pi0_mass/2.)*cost_gamma1rest );

    gamma2_rest.SetE(pi0_mass/2.);
    gamma2_rest.SetPx( -gamma1_rest.Px() );
    gamma2_rest.SetPy( -gamma1_rest.Py() );
    gamma2_rest.SetPz( -gamma1_rest.Pz() );

    //Boost gammas to lab frame
    TLorentzVector gamma1_lab = gamma1_rest;
    gamma1_lab.Boost(boost_vec);
    TLorentzVector gamma2_lab = gamma2_rest;
    gamma2_lab.Boost(boost_vec);

    // type 1 is final state
    // pdgid 22 - gamma
    GenParticlePtr p_gamma1 = std::make_shared<GenParticle>(
        FourVector(gamma1_lab.Px(), gamma1_lab.Py(), gamma1_lab.Pz(), gamma1_lab.E()),
        photon_pdgID, 1 );

    GenParticlePtr p_gamma2 = std::make_shared<GenParticle>(
        FourVector(gamma2_lab.Px(), gamma2_lab.Py(), gamma2_lab.Pz(), gamma2_lab.E()),
        photon_pdgID, 1 );

    GenVertexPtr v_pi0_decay= std::make_shared<GenVertex>();
    v_pi0_decay->add_particle_in(p_pi0);

    v_pi0_decay->add_particle_out(p_gamma1);
    v_pi0_decay->add_particle_out(p_gamma2);
    evt.add_vertex(v_pi0_decay);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
