cat $1 | shuf | vw --passes=25 --cache_file cache.f --binary --ect 5 --interactions vvv -f multiway.model --ignore s
