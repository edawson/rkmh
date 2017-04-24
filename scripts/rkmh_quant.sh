cat $1 | grep -v "FAIL" | cut -f 2 | cut -f 2 -d " " | cut -f 3 -d "|" | sort | uniq -c
