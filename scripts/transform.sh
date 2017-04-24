infile=$1
collapse=$2
if [ "$2" == "collapse" ]
    then python vwize.py -i $infile -C 1 | python collapse_subtypes.py
else
    python vwize.py -n -i $infile -C 1
fi

