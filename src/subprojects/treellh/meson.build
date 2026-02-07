project('treellh', ['cpp'],
  default_options: [
    'cpp_std=c++20', 
    'warning_level=3', 
    'werror=true'])


tskit_proj = subproject('tskit')
tskit_dep = tskit_proj.get_variable('tskit_dep')
tskit_inc = tskit_proj.get_variable('tskit_inc')

lib_deps = [tskit_dep]

# Subprojects should compile in the static library for simplicity.
treellh_inc = [tskit_inc, include_directories('.')]
treellh = static_library('treellh', 'treellh.cpp', dependencies: lib_deps)
treellh_dep = declare_dependency(link_with : treellh, include_directories: treellh_inc)

