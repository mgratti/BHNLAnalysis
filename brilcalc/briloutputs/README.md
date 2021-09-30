## Outputs of BRILCALC

Documentation of the BRILCALC tool: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html

Example command:

```
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json  -i recommended_json_2018.txt --hltpath "HLT_Mu7_IP4*part*" -u /fb -o output_byls_ HLT_Mu7_IP4.csv
```

The json file `recommended_json_2018.txt` is [this one](https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt)
