#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import re,os,sys
f = open(sys.argv[1], 'r')

seq_length = 0
seq_count = 0

for line in f:
    if line.startswith(">"):
        seq_id = line.strip()[1:]
        seq_count += 1
        seq_length = 0
    else:
        seq_length += len(line.strip())
        print(seq_id, seq_length)


