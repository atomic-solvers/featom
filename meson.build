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

# --------------------- Executable
subdir('app')
