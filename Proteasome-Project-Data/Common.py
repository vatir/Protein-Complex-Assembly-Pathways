#!/data/koanb/vpy-2.7.3/bin/python2.7
# -*- coding: utf-8 -*-

def Chunks(l, k):
  """ Yield k successive chunks from l. With largest chunks first."""
  if k < 1:
    yield []
    raise StopIteration
  n = len(l)
  avg = n/k
  remainders = n % k
  start, end = 0, avg
  while start < n:
    if remainders > 0:
      end = end + 1
      remainders = remainders - 1
    yield l[start:end]
    start, end = end, end+avg