project('ARGLy', 'cpp',
  default_options: [
    'cpp_std=c++20',
    'optimization=3',
    'debug=false',
    'b_ndebug=true']
    )

treellh_proj = subproject('treellh')
treellh_dep = treellh_proj.get_variable('treellh_dep')

executable('main',
  'main.cpp',
   dependencies : [treellh_dep],
   install : true)

nlopt_dep = dependency('nlopt', required: false)
nlopt_inc = []
if not nlopt_dep.found()
  brew = find_program('brew', required: false)
  if brew.found()
    nlopt_prefix = run_command(brew, '--prefix', 'nlopt', check: true).stdout().strip()
    nlopt_inc = include_directories(nlopt_prefix / 'include')
    nlopt_dep = meson.get_compiler('cpp').find_library('nlopt', dirs: nlopt_prefix / 'lib', required: true)
  else
    nlopt_dep = meson.get_compiler('cpp').find_library('nlopt', required: true)
  endif
endif

executable('argly',
  'argly_cli.cpp', 'fit_search.cpp',
   dependencies : [treellh_dep, nlopt_dep],
   include_directories : nlopt_inc,
   install : true)

cache_test = executable('cache_test', 'cache_test.cpp', dependencies: [treellh_dep])
test('likelihood cache', cache_test)

executable('argly-epoch-fit', 'epoch_cli.cpp', 'epoch_model.cpp', 'epoch_fit.cpp',
  dependencies: [treellh_dep, nlopt_dep], include_directories: nlopt_inc, install: true)
epoch_test = executable('epoch_test', 'epoch_test.cpp', 'epoch_model.cpp', 'epoch_fit.cpp',
  dependencies: [treellh_dep, nlopt_dep], include_directories: nlopt_inc)
test('epoch likelihood and ECME', epoch_test)
