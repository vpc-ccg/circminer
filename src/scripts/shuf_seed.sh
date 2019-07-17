#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "usage: ./shuf_seed.sh FILE"
    exit 1
fi

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

shuf --random-source=<(get_seeded_random 7) $1
