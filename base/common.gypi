{
  'target_defaults': {
    'default_configuration': 'Release',
    'cflags': ['-std=c++11', '-Wall'],
    'include_dirs': ['<(DEPTH)'],

    'configurations': {
      'Debug': {
        'cflags+': ['-g', '-O0'],
        'msvs_settings': {
          'VCCLCompilerTool': {
            'Optimization': '0', # /Od
            'RuntimeLibrary': '1', # /MTd
          },
          'VCLinkerTool': {
            'LinkIncremental': '2',
          },
        },
        'xcode_settings': {
          'GCC_OPTIMIZATION_LEVEL': '0', # -O0
        },
      }, # Debug
      'Release': {
        'cflags+': ['-O3'],
        'msvs_settings':{
          'VCCLCompilerTool': {
            'Optimization': '2', # /O2
            'InlineFunctionExpansion': '2',
            'RuntimeLibrary': '0', # /MT
          },
          'VCLinkerTool': {
            'LinkIncremental': '1',
            'OptimizeReferences': '2',
          },
        },
        'xcode_settings': {
          'GCC_OPTIMIZATION_LEVEL': '3', # -O3
        },
      }, # Release
    },
    'variables': {
      'component%': 'static_library',
    },
  }, # target_defaults
}
