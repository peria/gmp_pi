{
  'includes': [
    'base/common.gypi'
  ],
  'targets': [
    {
      'target_name': 'chudnovsky',
      'type': 'none',
      'dependencies': [
        'chudnovsky/chudnovsky.gyp:chudnovsky_all',
      ]
    }
  ]
}
