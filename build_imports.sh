#!/usr/bin/env bash

# wdl/tools/picard/collect_rna_seq_metrics.wdl

set -e

function cleanup() {
    cd ${prev_dn}
    rm -rf "${working_dn}"
}
trap cleanup EXIT

if [ -z "${1}" ]; then
    echo USAGE: ${0} WDL && exit 2
fi

wdl="${1}"
if [ ! -s "${wdl}" ]; then
    echo Given WDL "${wdl}" does not exist! && exit 2
fi
wdl_dn=$(readlink -f "${wdl}" | xargs -I% dirname %)
echo  "${wdl_dn}"

working_dn=$(mktemp -d 2>/dev/null || mktemp -d -t 'build-imports-XXX')
echo WDL: "${wdl}"
for import_wdl in $(grep import "${wdl}" | awk -F\" '{print $2}'); do
    dest_bn=$(dirname ${import_wdl})
    dest_dn="${working_dn}/${dest_bn}"
    mkdir -p "${dest_dn}"
    #cp "../${import_wdl}" ${dest_dn}
    #import_wdl_imports="../${import_wdl/.wdl/.imports.zip}"
    cp "${import_wdl}" ${dest_dn}
    import_wdl_imports="${import_wdl/.wdl/.imports.zip}"
    if [ -e "${import_wdl_imports}" ]; then
        unzip -qo "${import_wdl_imports}" -d "${dest_dn}"
    fi
done

prev_dn=$(pwd -P)
cd ${working_dn}
zip_fn=$(basename "${wdl}" | sed 's/\.wdl/.imports.zip/')
zip -qr "${zip_fn}" .
wdl_dn=$(dirname "${wdl}")
cd "${prev_dn}"
cd "${wdl_dn}"
cp -f "${working_dn}/${zip_fn}" .
unzip -l "${zip_fn}"
echo Created imports file: "${wdl_dn}/${zip_fn}"

exit
