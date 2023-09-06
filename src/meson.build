# Source files
_srcs = [
  'constants.f90',
  'dirac.f90',
  'energies.f90',
  'fe.f90',
  'feutils.f90',
  'graphs.f90',
  'graphs_potential.f90',
  'hartree_screening.f90',
  'lapack.f90',
  'linalg.f90',
  'mesh.f90',
  'mixings.f90',
  'quadrature.f90',
  'schroed_dirac_solver.f90',
  'schroed_glob.f90',
  'solvers.f90',
  'states.f90',
  'string_utils.f90',
  'types.f90',
  'xc.f90',
]

# Dependencies
gauss_jacobilib = library('gjp_gw',
                 ['gjp_gw_single.f90'],
                 dependencies: _deps,
                 cpp_args: _args,
                 link_with: _linkto,
                 include_directories: _incdirs,
                 install: true)

_linkto += [gauss_jacobilib]

# Create library
featomlib = library('featom',
                 _srcs,
                 dependencies: _deps,
                 cpp_args: _args,
                 link_with: _linkto,
                 include_directories: _incdirs,
                 install: true)

featom_dep = declare_dependency(
  include_directories: _incdirs,
  link_with: _linkto + [featomlib],
  dependencies: _deps)
