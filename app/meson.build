# --------------------- Executables

_deps += featom_dep
conv = executable('conv', 'conv.f90',
                  dependencies: _deps,
                  include_directories: _incdirs,
                  cpp_args: _args,
                  link_with: _linkto,
                  install: false)

conv_potential = executable('conv_potential',
                            'conv_potential.f90',
                            dependencies: _deps,
                            include_directories: _incdirs,
                            cpp_args: _args,
                            link_with: _linkto,
                            install: false)

gpd_coulomb_schroed_nelements = executable('gpd_coulomb_schroed_nelements',
                                           'gpd_coulomb_schroed_nelements.f90',
                                           dependencies: _deps,
                                           include_directories: _incdirs,
                                           cpp_args: _args,
                                           link_with: _linkto,
                                           install: false)
