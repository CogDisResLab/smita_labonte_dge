#!/usr/bin/env bash

# This script does find and replace for the files given in the arguments

# Save the arguments into a list
files=("$@")

# fix pval_column
for file in "${files[@]}"; do
    echo "Fixing the errors in $file"
    sed -i bak 's/pval_column/pvalColumn/g' $file
    sed -i bak 's/gene_column/geneColumn/g' $file
    sed -i bak 's/logfc_column/logfcColumn/g' $file
    sed -i bak 's/get_concordants/getConcordants/g' $file
    sed -i bak 's/filter_signature/filterSignature/g' $file
    sed -i bak 's/prepare_signature/prepareSignature/g' $file
    sed -i bak 's/5%/25%/g' $file
    sed -i bak 's/10%/50%/g' $file
    sed -Ei bak 's/(.*_(up|dn)_sig),.*$/\1)/' $file
done