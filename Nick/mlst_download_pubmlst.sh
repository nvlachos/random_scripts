#!/bin/bash

set -e

OUTDIR=/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts
if [[ ! -d "$OUTDIR"]]; then
  mkdir -p "$OUTDIR"
fi
wget -P "$OUTDIR" http://pubmlst.org/data/dbases.xml
#wget --no-clobber -P "$OUTDIR" http://pubmlst.org/data/dbases.xml

for URL in $(grep '<url>' $OUTDIR/dbases.xml); do
#  echo $URL
  URL=${URL//<url>}
  URL=${URL//<\/url>}
#  echo ${URL: -4}
  if [ ${URL:(-4)} = ".txt" ]; then
    PROFILE=$(basename $URL .txt)
    echo "# $PROFILE "
    PROFILEDIR="$OUTDIR/$PROFILE"
    echo "mkdir -p '$PROFILEDIR'"
    echo "(cd '$PROFILEDIR' && echo "$URL" && wget -q '$URL')"
  elif [ ${URL:(-4)} = ".tfa" ]; then
    echo "(cd '$PROFILEDIR' && echo "$URL" && wget -q '$URL')"
  fi
done

# delete fungi schemes
echo rm -frv "$OUTDIR"/{tvaginalis,ckrusei,ctropicalis,calbicans,sparasitica}
