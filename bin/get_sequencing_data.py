#!/usr/bin/env python

import requests
import subprocess as sp
import sys
import glob

run_id = sys.argv[1]

sra_prefetch = sp.run(['prefetch', f"{run_id}"], stdout=sp.PIPE)
if sra_prefetch.returncode == 0:
    fastq_dump = sp.run(["fasterq-dump", "--split-files", f"{run_id}"], stdout=sp.PIPE)
    fastq_files = glob.glob(f"{run_id}*.fastq")
    if fastq_files:
        sp.run(["gzip"] + fastq_files, check=True)
    else:
        print("[ERROR] I could not find any .fastq to compress")
        sys.exit(81)
elif (sra_prefetch.returncode == 3 ) and (run_id.startswith("ERR")):
    response = requests.get(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_id}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes")
    print("status_code_request_ena: {}".format(response.status_code))
    if response.status_code == 200:
        payload = response.content.decode('utf-8').split("\n")
        payload = list(filter(None, payload))
        if len(payload) == 2:
            read_set = payload[1].split("\t")[1].split(";")
            for lnk in read_set:
                d = sp.run(["curl", "-O", "--retry", "5", f"http://{lnk}"])
        else:
            print("[ERROR] ENA response: payload of unexpected length. Possible explanation: no public data is associated with this run")
            print(payload)
            sys.exit(80)
    else:
        print(response.status_code)
        print(response.content)      
else:
    print("[ERROR] Something went wrong! Check the prefetch exit code")
