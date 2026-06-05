project('ARGLy', 'cpp',
  default_options: [
    'cpp_std=c++23']
    )

treellh_proj = subproject('treellh')
treellh_dep = treellh_proj.get_variable('treellh_dep')

openmp_dep = dependency('openmp')

executable('main',
  'main.cpp',
   dependencies : [treellh_dep, openmp_dep],
   install : true)

executable('argly',
  'argly_cli.cpp',
   dependencies : [treellh_dep],
   install : true)

