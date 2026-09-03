#!/usr/bin/env bash

# git-archive-submodules - A script to produce an archive tar.gz file of the a git module including all git submodules
# based on https://ttboj.wordpress.com/2015/07/23/git-archive-with-submodules-and-tar-magic/


usage()
{
  echo >&2 "git-archive-submodules - A script to produce an archive tar.gz file of the
                       a git module including all git submodules"
  echo >&2 "requires: git, sed, gzip, and tar (or gtar on macos)"
  echo >&2 "usage: $0 [destination]"
}

# requires gnu-tar on mac
export TARCOMMAND=tar
case "$OSTYPE" in
  darwin*)    export TARCOMMAND=gtar;;
  linux-gnu*) export TARCOMMAND=tar;;
  *)        echo "unknown: $OSTYPE" && exit 1;;
esac
command -v ${TARCOMMAND} >/dev/null 2>&1 || { usage; echo >&2 "ERROR: I require ${TARCOMMAND} but it's not installed.  Aborting."; exit 1; }

# reqiures git
command -v git >/dev/null 2>&1 || { usage; echo >&2 "ERROR:I require git but it's not installed.  Aborting."; usage; exit 1; }

# requires sed
command -v sed >/dev/null 2>&1 || { usage; echo >&2 "ERROR:I require sed but it's not installed.  Aborting."; usage; exit 1; }

# requires gzip
command -v gzip >/dev/null 2>&1 || { usage; echo >&2 "ERROR:I require gzip but it's not installed.  Aborting."; usage; exit 1; }

export TARMODULE=`basename \`git rev-parse --show-toplevel\``
export TARVERSION=`git describe --tags --abbrev=0 | sed 's/v//g'`
export TARPREFIX="${TARMODULE}-${TARVERSION}"

# create module archive
git archive --prefix=${TARPREFIX}/ -o ${TMPDIR}${TARPREFIX}.tar v${TARVERSION}
if [[ ! -f "${TMPDIR}${TARPREFIX}.tar" ]]; then
  echo "ERROR: base sourcecode archive was not created. check git output in log above."
  usage
  exit 1
fi

# force init submodules
git submodule update --init --recursive

# tar each submodule recursively
git submodule foreach --recursive 'git archive --prefix=${TARPREFIX}/${displaypath}/ HEAD > ${TMPDIR}tmp.tar && ${TARCOMMAND} --concatenate --file=${TMPDIR}${TARPREFIX}.tar ${TMPDIR}tmp.tar'

# compress tar file
gzip -9 ${TMPDIR}${TARPREFIX}.tar
if [[ ! -f "${TMPDIR}${TARPREFIX}.tar.gz" ]]; then
  echo "ERROR: gzipped archive was not created. check git output in log above."
  usage
  exit 1
fi

# copy file to final name and location if specified
if [[ -z "$1" ]]; then
  destination=${TARPREFIX}-including-submodules.tar.gz
else
  destination=$1
fi
cp ${TMPDIR}${TARPREFIX}.tar.gz ${destination}
if [[ -f "${TMPDIR}${TARPREFIX}.tar.gz" ]]; then
  rm ${TMPDIR}${TARPREFIX}.tar.gz
  echo "created ${destination}"
else
  echo "ERROR copying ${TMPDIR}${TARPREFIX}.tar.gz to ${destination}"
  usage
  exit 1
fi
