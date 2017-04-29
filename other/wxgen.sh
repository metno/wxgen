# START wxgen completion
_wxgen()
{
local cur prev opts mode
COMPREPLY=()
cur="${COMP_WORDS[COMP_CWORD]}"
prev="${COMP_WORDS[COMP_CWORD-1]}"
mode="${COMP_WORDS[1]}"
if [ "$mode" = "" ] ; then
   COMPREPLY=( $( compgen -W " sim truth verif" -- $cur ) )
elif [ "$mode" = "sim" ]; then
   if [ "$prev" = "-db" ]; then
      COMPREPLY=( $( compgen -f -W -- $cur ) )
   elif [ "$prev" = "-m" ]; then
      COMPREPLY=( $( compgen -W "exp rmsd " -- $cur ) )
   elif [ "$cur" = "" ] || [[ "$cur" =~ -* ]]; then
      COMPREPLY=( $( compgen -W "-h -n -t -m -rs -w -i -j -b -p -db -dbtype -o -wl -v --debug -s -lat -lon " -- $cur ) )
   fi
elif [ "$mode" = "truth" ]; then
   if [ "$prev" = "-db" ]; then
      COMPREPLY=( $( compgen -f -W -- $cur ) )
   elif [ "$cur" = "" ] || [[ "$cur" =~ -* ]]; then
      COMPREPLY=( $( compgen -W "-h -sd -ed -n -t -db -dbtype --db -o -wl -v --debug -s -lat -lon " -- $cur ) )
   fi
elif [ "$mode" = "verif" ]; then
   if [ "$prev" = "-m" ]; then
      COMPREPLY=( $( compgen -W "consecutive iqr max mean median min range std sum variance " -- $cur ) )
   elif [ "$prev" = "-m" ]; then
      COMPREPLY=( $( compgen -W "dryday frostday nothing summerday wetday " -- $cur ) )
   elif [ "$prev" = "-m" ]; then
      COMPREPLY=( $( compgen -W "covarmap histogram jump map sortstat timestat timeseries variance " -- $cur ) )
   elif [ "$cur" = "" ] || [[ "$cur" =~ -* ]]; then
      COMPREPLY=( $( compgen -W "-h -fs -m -o -truth -xlim -xlog -ylim -ylog -r -tr -a -clim -cmap -tm -ts -v --debug -s -lat -lon " -- $cur ) )
   fi
fi
return 0
}
complete -F _wxgen wxgen
# END wxgen completion
