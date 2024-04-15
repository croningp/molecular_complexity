#!/usr/bin/env fish

echo Predicting Spectra of 5-Aminoisophthalic acid

curl 'http://127.0.0.1:8091/v1/ir?smiles=C1=C(C=C(C=C1C(=O)O)N)C(=O)O&method=GFNFF' >> ir_example_GFNFF.json
curl 'http://127.0.0.1:8091/v1/ir?smiles=C1=C(C=C(C=C1C(=O)O)N)C(=O)O&method=GFN1xTB' >> ir_example_GFN1xTB.json
curl 'http://127.0.0.1:8091/v1/ir?smiles=C1=C(C=C(C=C1C(=O)O)N)C(=O)O&method=GFN2xTB' >> ir_example_GFN2xTB.json
