project('featom', 'fortran',
        version:'0.1.0',
       default_options: ['wrap_mode=default',
                         'warning_level=3'])

# Add compiler options
_args = [] # Extra arguments
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries
_incdirs = [] # All the includes

fc = meson.get_compiler('fortran')

# Platform detection
host_system = host_machine.system()
is_windows = host_system == 'windows'
is_mingw = is_windows and fc.get_id() == 'gcc'

cpu_family = host_machine.cpu_family()

if is_mingw
  # For mingw-w64, don't use LTO
  add_project_arguments('-fno-use-linker-plugin', language: ['fortran'])
endif

# --------------------- Dependencies
lapack = dependency('lapack', required: true)
_deps += [lapack]

# --------------------- Library
subdir('src') # defines featom_dep
_deps += [featom_dep]

# --------------------- Executable apps
if get_option('with_app')
  subdir('app')
endif

# ------------------------ Tests

if get_option('with_tests') and not is_windows
test_array = [#
  ['CoulombDirac', 'testCoulombDirac', 'test_coulomb_dirac.f90'],
  ['CoulombSchroed', 'testCoulombSchroed', 'test_coulomb_schroed.f90'],
  ['DftDirac', 'testDftDirac', 'test_dft_dirac.f90'],
  ['DftSchroed', 'testDftSchroed', 'test_dft_schroed.f90'],
  ['HarmonicSchroed', 'testHarmonicSchroed', 'test_harmonic_schroed.f90'],
  ['HarmonicDirac', 'testHarmonicDirac', 'test_harmonic_dirac.f90'],
  ['DftDiracFast', 'testDftDiracFast', 'test_dft_dirac_fast.f90'],
  ['DftSchroedFast', 'testDftSchroedFast', 'test_dft_schroed_fast.f90'],
             ]
foreach test : test_array
  test(test.get(0),
       executable(test.get(1),
          sources : ['test/'+test.get(2)],
          dependencies : _deps,
          link_with : _linkto,
          cpp_args : _args,
          include_directories: _incdirs,
                 ),
      )
endforeach
endif
