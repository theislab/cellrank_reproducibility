#!/usr/bin/env bash

set -e

SIZE=$2
SPLIT=$3
N_JOBS=$4

case $1 in
"gpcca")
  python3 scripts/utils.py gpcca $SIZE $SPLIT --n-jobs=$N_JOBS
  ;;
"velocyto")
  rm -rf velocyto.loom
  python3 scripts/utils.py velocyto $SIZE $SPLIT --n-jobs=$N_JOBS
  ;;
"palantir")
  python3 scripts/utils.py palantir $SIZE $SPLIT --n-jobs=$N_JOBS
  ;;
"stemnet")
  rm -rf *.rds
  python3 scripts/utils.py stemnet $SIZE $SPLIT
  Rscript --vanilla scripts/benchmark_stemnet.r $SIZE $SPLIT
  ;;
"fateid")
  rm -rf *.rds
  python3 scripts/utils.py fateid $SIZE $SPLIT
  Rscript --vanilla scripts/benchmark_fateid.r $SIZE $SPLIT
  ;;
*)
  echo "Invalid method \`$1\`."
  exit 1
esac
