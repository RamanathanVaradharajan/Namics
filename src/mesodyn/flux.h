#ifndef FLUX_H
#define FLUX_H

#include <thread>
#include "factory.h"
#include "../lattice.h"
#include "lattice_object.h"
#include "neighborlist.h"
#include "component.h"
#include "gaussian_noise.h"
#include "lattice_accessor.h"

class IFlux;

namespace Flux {
  typedef Factory_template<IFlux, Dimensionality, Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>> Factory;
}

class IFlux {
  public:

    IFlux(Lattice* lat_, shared_ptr<IComponent> A_, shared_ptr<IComponent> B_);
    ~IFlux() {}

    // Calculate flux at site depending on property differenceof IComponent A and B and store in J.
    // Assumes the flux generated by difference between B and A is opposite and equal
    virtual void flux() = 0;
    Lattice_object<Real> J;
    Lattice* m_lat;

    shared_ptr<IComponent> component_a;
    shared_ptr<IComponent> component_b;
};

class ILangevin_flux : public IFlux {
  public:
    ILangevin_flux(Lattice*, Real, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>);

  protected:
    int onsager_coefficient(Lattice_object<Real>&, Lattice_object<Real>&);
    int potential_difference(Lattice_object<Real>&, Lattice_object<Real>&);

    Lattice_object<Real> L;
    Lattice_object<Real> mu;

    const Real D;
    shared_ptr<Gaussian_noise> gaussian;
};

class Flux1D : public ILangevin_flux {
public:
  Flux1D(Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>);
  virtual ~Flux1D();

  virtual void flux() override;

  enum error {
    ERROR_SIZE_INCOMPATIBLE,
    ERROR_NOT_IMPLEMENTED,
  };

  Lattice_object<Real> J_plus;

protected:
  void attach_neighborlists(shared_ptr<Neighborlist>, const Offset_map&);
  int langevin_flux(const Offset_map&);

  Lattice_object<Real> t_L;
  Lattice_object<Real> t_mu;

private:
  Offset_map offset;
};

class Flux2D : public Flux1D {
public:
  Flux2D(Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>);
  virtual ~Flux2D();

  virtual void flux() override;

private:
  Offset_map offset;
};

class Flux3D : public Flux2D {
public:
  Flux3D(Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>);
  ~Flux3D();

  virtual void flux() override;

private:
  Offset_map offset;
};

class Flux3D_extended_stencil : public ILangevin_flux {
public:
  Flux3D_extended_stencil(Lattice*, Real, const Lattice_object<size_t>&, shared_ptr<IComponent>, shared_ptr<IComponent>, shared_ptr<Gaussian_noise>);
  ~Flux3D_extended_stencil();

  virtual void flux() override;
  int langevin_flux_forward(const Offset_map&);
  int langevin_flux_backward(const Offset_map&);

private:
  std::vector<Offset_map> forward_offsets;
  std::vector<Offset_map> backward_offsets;
  Lattice_object<Real> t_L;
  Lattice_object<Real> t_mu;
  Lattice_object<Real> t_J;
};

#endif