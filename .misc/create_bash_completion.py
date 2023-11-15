#!/usr/bin/env python
import os
import re
import sys
import wxgen.driver
import wxgen.metric
import wxgen.output
from cStringIO import StringIO


"""
This script creates a bash completion script by parsing command-line options from wxgen's
description as well as reading all classes from wxgen.metric and wxgen.output
"""

def get_options(message):
   message = message.split('\n')
   for start in range(len(message)):
      if message[start] == "":
         break
   options = list()
   for i in range(start+1, len(message)):
      value = [q for q in message[i].replace(',','').split(' ') if q != ""]
      if len(value) > 0 and value[0][0] == '-':
         options += [value[0].strip()]
   return options

parser,sp = wxgen.driver.get_parsers()

print("# START wxgen completion")
print("_wxgen()")
print("{")
print("local cur prev opts mode")

print('COMPREPLY=()')
print('cur="${COMP_WORDS[COMP_CWORD]}"')
print('prev="${COMP_WORDS[COMP_CWORD-1]}"')
print('mode="${COMP_WORDS[1]}"')

#print 'if [ "$cur" = "" ] || [[ "$cur" =~ -* ]]; then'
#print "   COMPREPLY=( $( compgen -f -W '",
#for s in sp.keys():
#   print s,
#print "' -- $cur ) )"
#print 'fi'
old_stdout = sys.stdout
sys.stdout = message = StringIO()

sp["sim"].print_help()
sim_flags = get_options(message.getvalue())
sys.stdout = message = StringIO()
sp["truth"].print_help()
truth_flags = get_options(message.getvalue())
sys.stdout = message = StringIO()
sp["verif"].print_help()
verif_flags = get_options(message.getvalue())
sys.stdout = old_stdout

aggregators = [mod for mod in wxgen.driver.get_module_names(wxgen.aggregator) if mod != "quantile"]
transforms = wxgen.driver.get_module_names(wxgen.transform)
plots = wxgen.driver.get_module_names(wxgen.plot)
metrics = wxgen.driver.get_module_names(wxgen.metric)

def get_option_string(flag, options, use_elif=True):
   s = '   '
   if use_elif:
      s += 'elif '
   else:
      s += 'if '
   s += '[ "$prev" = "' + flag + '" ]; then\n'
   s += '      COMPREPLY=( $( compgen -W "'
   for option in options:
      s += '%s ' % option
   s += '" -- $cur ) )'
   return s

def get_flag_string(flags, use_elif=True):
   s = '   '
   if use_elif:
      s += 'elif '
   else:
      s += 'if '
   s += '[ "$cur" = "" ] || [[ "$cur" =~ -* ]]; then\n'
   s += '      COMPREPLY=( $( compgen -f -W "'
   for flag in flags:
      s += '%s ' % flag
   s += '" -- $cur ) )'
   return s

print('if [ "$mode" = "" ] ; then')
print('   COMPREPLY=( $( compgen -W " sim truth verif" -- $cur ) )')

print('elif [ "$mode" = "sim" ]; then')
print('   if [ "$prev" = "-db" ]; then')
print('      COMPREPLY=( $( compgen -f -W -- $cur ) )')
print(get_option_string("-m", metrics))
print(get_flag_string(sim_flags))
print('   fi')

print('elif [ "$mode" = "truth" ]; then')
print('   if [ "$prev" = "-db" ]; then')
print('      COMPREPLY=( $( compgen -f -W -- $cur ) )')
print(get_flag_string(truth_flags))
print('   fi')

print('elif [ "$mode" = "verif" ]; then')
print(get_option_string('-a', aggregators, False))
print(get_option_string('-tr', transforms))
print(get_option_string('-m', plots))
print(get_flag_string(verif_flags))
print('   fi')
print('fi')

print('return 0')
print('}')
print('complete -F _wxgen wxgen')
print('# END wxgen completion')
