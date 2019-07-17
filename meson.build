project('gallium-oxide', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])

gallium_oxide = executable('gallium_oxide',
                           'gallium_oxide.cc',
                           dependencies : [subproject('probable').get_variable('probable_dep')])
